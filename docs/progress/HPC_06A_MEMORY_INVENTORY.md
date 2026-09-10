# HPC-06A 記憶體生命週期盤點

狀態：首次 source inventory 與一項輸出複製移除；尚未完成大型記憶體驗收。
日期：2026-09-10。以下由 `c920d9e` 的實作盤點，輸出陣列修正另列於下方。

令 P 為 communicator ranks、N 為控制點數、F 為 transport fields、E 為元素數、
R 為資料庫所有 rank required-element 索引的總長、Q 為可視化 Bezier 點數。
大小公式只計有效陣列 payload，假定 real64／int32；不含 capacity、allocator、
PETSc scatter／matrix、MPI buffers 或容器 metadata，不能當作 RSS 量測。

| 路徑／來源 | 分布與大小 | 生命週期與必要性 |
|---|---|---|
| [GatherAllPetscReal](../../solvers/cpu/include/PetscGather.hpp) | 每 rank 同時有完整 PETSc Vec 與回傳 vector；D rows 至少 16D bytes | 每次呼叫的短期 gather，加上 caller 保留回傳值的時間。是顯式全域診斷介面，並非 owned/ghost assembly；不能把它誤判成每次 transport assembly 都執行。 |
| [TransientTransportRuntime::GatherState](../../solvers/cpu/include/TransientTransportRuntime.hpp) | D=N×F；以共用 gather 複製 | 呼叫者要求完整 state 時才發生。assembly 使用 `GatherRequiredState`，integral 使用 fields 大小的 Allreduce。後續須依實際 call graph／profile 分辨測試與 production 使用頻率。 |
| [OneDGetVectorAll](../../solvers/one_d/include/OneDImplicit.hpp) | 每 rank 完整解向量；gather Vec 與 candidate 至少 16D bytes，舊 `values` 可能仍存在 | 線性求解結果、非線性 residual 與 Jacobian callback 都會呼叫。Residual 另建立 D 長度 values。小型 replicated network 可保留；大型 network 必須量測 callback 頻率與 owner／分散式替代成本。 |
| [WriteFlowOutput](../../solvers/cpu/src/iga_navier_stokes.cpp) | root gather 4N doubles，再分離 3N velocity 與 N pressure，至少 64N bytes | 每次場輸出；scatter 與 gather Vec 持續到 writer 返回。分離陣列為目前 VTU／VTKHDF API 輸入。下方移除 initializer-list 額外複製；root 完整場仍存在。 |
| [iga_solve 輸出](../../solvers/cpu/src/iga_solve.cpp)、[iga_transport 輸出](../../solvers/cpu/src/iga_transport.cpp) | `VecScatterCreateToZero`，root 完整場 | 輸出期間存在；不能由 owned-row assembly 推論輸出已分散。HPC-06B 須提供另一條不收集完整場的路徑。 |
| [Database](../../solvers/cpu/include/IgaDatabase.hpp) | 每 rank offsets 8(E+1)、owners 4E、rank_offsets 8(P+1)、rank_elements 8R bytes | Database 生命週期持有完整索引。`LoadRequired` 只載入該 rank 所需元素，非每 rank 完整元素 payload；另短期複製該 rank indices。R 可大於 E，不能漏掉多 rank required 重複。 |
| [TransientFlowRuntime](../../solvers/cpu/include/TransientFlowRuntime.hpp) | required/owned 元素、局部 PETSc state／ghost，加 replicated boundary/label arrays | trial 期間還保存 committed boundaries、outlets 等，以支持 rollback。checkpoint capture 另複製 boundary 狀態；不可直接刪去 transaction 必需快照。 |
| [MovingImmersedTransientFlowRuntime](../../solvers/cpu/include/MovingImmersedTransientFlowRuntime.hpp) | committed 與 candidate/trial 各有 geometry、catalog 與 inner runtime；BeginTrial 另建 pressure、extension、seed、map | BeginTrial 建構成功前保留 committed；FinalizeCommit 交換後刪除舊 epoch；AbortTrial 丟棄 trial。雙 epoch 是目前強例外安全與 rollback 契約的一部分。需量測各 phase，而非把雙份狀態直接視為 leak。 |
| [TemporalVtkHdfWriter](../../include/TemporalVtkHdf.hpp)、[BezierVisualization](../../include/BezierVisualization.hpp) | writer 擁有 mesh 副本，含 coordinates、connectivity、signatures；Append 另抽取 Q×components doubles | mesh 存活至 writer 解構；抽取陣列存活至 Append 返回，與控制點陣列同時存在。HDF5 compression 並不移除這些 host allocations。 |
| [MemoryReport](../../solvers/cpu/include/MemoryReport.hpp) | root Gather 6P doubles | 每個 memory stage 的小型 metadata；相較全場 gather 屬可保留的簡單實作，輸出包含每 rank metrics。 |

## 已移除的輸出複製

`WriteFlowOutput` 原先以 initializer-list 建立 `vector<VtkPointArray>`。
即使先將 velocity／pressure move 到 initializer-list 元素，這些元素是 const，
vector 建構仍複製其中的 payload。現在 reserve 兩個 entries，再分別 push_back
含 moved vector 的暫存物件，讓 payload 所有權移入最終 arrays。

這移除陣列組裝期間額外的 4N doubles（32N bytes）；例如 N=1,000,000 時是
32,000,000 bytes 的一次 payload 複製。這是由容器語義推得的局部配置差異，
不是已量測的整個程序峰值 RSS 降幅。後續 Bezier extraction 可能才是總峰值。
文字場格式、順序、VTU／VTKHDF arrays schema、collective cleanup 不變。

既有 `flow_output_failure_test` 直接包含 production writer；重建後在 world 3 ranks
及 split 1+2 groups，各 11 種故障／重試，共 33 種條件通過。核對文字速度／壓力、
VTU、來源 Vec 不變與 collective 錯誤傳遞。這批不包含成功的 VTKHDF dataset 讀回，
不宣稱完成其格式驗收。三份 rank reports 均 exit 0、無 timeout。
來源與 binary／logs hashes：`outputs/hpc06/memory/flow-writer-audit.json`。

```bash
make -C solvers/cpu flow_output_failure_test \
  PETSC_DIR=/usr/lib/petscdir/petsc3.15/x86_64-linux-gnu-real
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 timeout --kill-after=5s 120s \
  mpiexec -np 3 solvers/cpu/flow_output_failure_test NEW_OUTPUT
```

## 下一步與驗收缺口

先以既有 per-rank memory reports 與 rank wrapper RSS 做固定案例 phase 比較，
區分 database indices、required-element payload、assembly workspace、accepted/trial
與 output 的重疊峰值。現有 concurrent 大方管作業可提供功能驗收時 RSS 觀察，
不能作為無干擾 scaling 證據。沒有實測瓶頸證據前，不改小型 1D／metadata 路徑。

CPU flow 已提供 PVTU 分片輸出與時間索引，並通過實際 PDE／ParaView 比較，
詳見 [HPC-06B 進度](HPC_06B_PARTITIONED_VTK_PROGRESS.md)。大型 RSS 與其他求解
路徑仍待驗證，HPC-06A／整份清單保持未勾選。

## CPU flow 階段量測接線

CPU `iga_navier_stokes` 新增可選 `--memory-report PATH`，沿用
`DistributedMemoryRecorder` 的 JSONL schema。各階段記錄全部 ranks 的目前 RSS、
歷史峰值 RSS、PETSc allocator／process 計數，不只 root。未啟用時不執行量測。
開檔與收尾共同協調；關閉成功後才印出 solver 成功摘要。各 rank 的 report path
亦須一致。

量測點包括 database_open、initialized_state、visualization_geometry、
visualization_ready、flow_trial_solved、output_begin、flow_closed。
PVTU 另有 parallel_piece_extracted／parallel_piece_published；序列輸出則在
writer 返回後記錄 serial_output_released。Geometry 量測在 root mesh 釋放前，
ready 則在 PVTU 模式釋放後；HDF 模式維持原 mesh 生命週期。

這是階段邊界採樣。RSS current 不保證涵蓋階段內瞬時峰值；VmHWM 是到該時點
為止的程序歷史峰值，不能當作該階段獨立 peak。各 rank 峰值總和也不是同步
aggregate peak。PVTU published 時 piece 仍存活，writer 返回後會釋放；序列
released 則已釋放 gather 暫存，兩者生命週期不同，不能直接相減當成收益。

1／2-rank HDF 與 PVTU 共四次實際小案例求解，所有階段有各 rank 有效 RSS，
每一步 output records 符合既有頻率；與前批無量測求解的 checkpoint.state
SHA256 完全相等。另以 report path 目錄阻擋驗證兩 rank 共同拒絕、無 timeout
且沒有成功摘要。這批仍是 one-cell 接線測試，不是大型記憶體效益證據。

```bash
python3 scripts/hpc_flow_memory_regression.py \
  --reference-root outputs/hpc06/flow-pvtu-v1/cli \
  --output-dir NEW_OUTPUT
```

證據：`outputs/hpc06/flow-memory-v1/regression/acceptance.json` 與 `audit.json`。
