# HPC-02：浸入式暫態體積積分 OpenMP 整合

日期：2026-09-08。狀態：實作、實際元素錯誤回復及 1／2／4／8-thread 完整 FSI
參考場比較、固定幾何 Newton、隔離效能矩陣與修正後無 OpenMP 回歸已通過；
貼體同 rank 單 thread 記憶體基準亦已完成。
HPC-02A／B／C 已完成；整份清單的其他階段仍保持原範圍。

## 實作範圍

[ImmersedTransientFlowRuntime.hpp](../../solvers/cpu/include/ImmersedTransientFlowRuntime.hpp)
的 `Assemble` 已接入有界元素批次執行器。呼叫執行緒 materialize 元素並
透過 PETSc 取得 nodal fields；worker 只執行 `BuildVolume`，使用凍結的
body force、唯讀 history／layout／quadrature 與私有局部矩陣。
所有 worker 結束且該批次無錯誤後，依原 cell 順序插入 PETSc。
壁面、保守 trace、port、ghost penalty 與 gauge 工作仍在呼叫執行緒。
Moving runtime 與 FSI 使用相同體積積分路徑；靜態浸入式 runtime 尚未接入。
貼體 flow 後續已另行接入，見 [MPI／OpenMP 進度](HPC_02_HYBRID_PROGRESS.md)。

[ElementAssemblyExecution.hpp](../../solvers/cpu/include/ElementAssemblyExecution.hpp)
在 runtime 建立時凍結執行設定，不加入物理 input hash：

- OpenMP build 預設使用 `omp_get_max_threads()`，可由 `OMP_NUM_THREADS` 配置。
- `IGA_ASSEMBLY_THREADS` 可明確覆寫；無 OpenMP build 只接受 1。
- 批次預設為 `min(threads, 8)` 個元素，`IGA_ASSEMBLY_BATCH_SIZE` 可覆寫。
- 多 threads 要求 MPI 提供至少 `MPI_THREAD_FUNNELED`；constructor／assembly
  必須由 MPI initialization thread 呼叫。assembly 在任何 PETSc 操作前拒絕外來執行緒。
- `IGA_PROFILE=1` 記錄每次組裝實際 team size、cells、batches 與峰值 resident items。
  元素數上限不等同 byte 預算，最終仍須檢查 RSS。

## 已發現並修正的失敗回復問題

實際 runtime 測試在第二批的 cell 10／12 完成積分後拋出例外，確認最小
cell ID 的錯誤由呼叫執行緒接收；第一批此時已插入 PETSc。測試暴露出
原有 rollback 不會完成 pending MatSetValues，重試的 `MatZeroEntries`
因此回報 PETSc error 73。不是以忽略 PETSc 錯誤或放寬數值門檻處理。

runtime 現在追蹤未結束的插入階段；下次組裝先完成該階段，再清空矩陣／
向量。測試驗證同一 trial 的 rollback／retry 及 abort 後建立新 trial，
都不殘留失敗批次之前的矩陣貢獻，committed state／hash 保持一致。

## 驗證與證據

環境沿用 HPC-00：GCC 11.4、OpenMPI 4.1.2、PETSc 3.15.5 real64／Int32。
此 WSL 工作站顯示 8 個 physical cores、16 個 logical CPUs；同 core 的
配對為 0/1、2/3、…、14/15。建置使用 warning flags，測試保留斷言。

```bash
make -C solvers/cpu compliant_channel_fsi_test compliant_channel_fsi_openmp_test \
  immersed_transient_flow_openmp_test parallel_immersed_volume_test \
  parallel_immersed_volume_serial_test \
  element_assembly_execution_test element_assembly_execution_openmp_test \
  PETSC_DIR=/usr/lib/petscdir/petsc3.15/x86_64-linux-gnu-real
python3 -m unittest discover -s scripts/tests -p 'test_hpc*.py'
```

目前已確認：

- 最終 expanded／compact 回歸：無 OpenMP build 各 1 thread，共兩個配置；
  OpenMP build 各 1／2／4 threads，共六個配置，全部退出碼 0。每個配置
  均檢查實際積分、殘差／Jacobian action、worker 參與、第二批失敗、
  rollback／retry、abort／new trial 與 caller 限制。完整日誌和 RSS 在
  `runtime-final-serial/`、`runtime-final-openmp/`；彙整見 `runtime-final-tests.json`。

- 27-cell 實際 expanded 體積積分的 1／2／4 threads：殘差與 Jacobian action
  完全一致，body-force callback 保留在主執行緒，worker 確實參與計算。
  第二批錯誤、rollback／retry、abort／new trial 及外來 caller 拒絕均通過。
  日誌：`outputs/hpc02/volume/runtime-test3.log`，退出碼 0。
- serial／OpenMP 執行設定測試，各在 MPI FUNNELED 與 SINGLE 初始化執行；
  四次退出碼均為 0，驗證不合法設定與 thread support 拒絕。
  證據：`outputs/hpc02/volume/execution-tests.json`。
- HPC Python 47 個測試通過，包括 CPU 清單、過量配置、實際 team size／
  batch 證據，以及隔離矩陣首次排除、缺漏／失敗拒絕、時間及記憶體門檻；
  日誌 `outputs/hpc02/volume/python-tests-final.log`。
- `build3.log` 建置 serial FSI、OpenMP FSI、固定暫態與實際元素測試通過，
  無編譯 warning。較早 `build.log` 暴露既有固定暫態測試缺少
  `PrescribedSurfaceMotion.hpp` 的明確 include，已修正測試來源。

`runtime-test.log` 是測試 fixture 時間設定錯誤的失敗嘗試；
`runtime-test2.log` 是實際 pending-insertion 回復問題的失敗證據。
`fsi-omp1/` 已主動終止舊候選，collector 記錄 failed；其取消理由在
`cancelled-for-fix.json`，不納入任何成功或效能統計。

## 完整 FSI 參考場比較

1／2／4／8 threads 均已取得 accepted collection，與 HPC-00D 的
`outputs/hpc00/reference-states/fsi-reference/` 比較，各九個物理量的
相對 L2 均為 0。完整原生收斂、force／moment、moving mass、wall leakage
及 continuity gates 均通過；四次 strong iterations 與原參考一致。
執行設定不改變 fixture input identity。下表是同時執行其他正確性測試時
的單次觀測，不用於正式 speedup 結論。

| threads | process wall (s) | peak RSS (bytes) | 比較證據 |
|---:|---:|---:|---|
| 1 | 701.660275 | 75,579,392 | `fsi-omp1-comparison.json` |
| 2 | 442.800829 | 72,986,624 | `fsi-omp2-comparison.json` |
| 4 | 287.577750 | 77,352,960 | `fsi-omp4-comparison.json` |
| 8 | 247.165873 | 79,482,880 | `fsi-omp8-comparison.json` |

所有原始 profile、輸入／來源／binary 雜湊、20 次實際流場 team 與 batch
報告保留於對應的 `fsi-ompN-fixed/collection.json`。4-thread 的 assembly
為 276.876294 s，solver setup 0.768414 s，linear solve 0.008903 s；
這些是各 exclusive phase，並非新增 worker 時間後重複相加。
profile 的 28 次 assembly 包含流場以外的工作；`PretensionedMembrane` 與
`FluidSurfaceTraction` 也使用該 phase，不能把 28 當成流場 team 報告次數。

`source-and-binaries.tar.gz` 封存 143 個檔案，並逐一核對與 accepted
4-thread collection 的 source hashes 一致；另含最終元素／資源測試來源、
執行檔及 Makefile。archive manifest 為 `source-archive.json`。`accepted-evidence.json` 再次核對
四份 collection、36 個物理場、80 份流場 team 報告、8 個實際元素配置與
archive，並綁定 114 個 artifact 雜湊；這份 accepted 僅指本段正確性範圍。

新增的 `scripts/hpc_openmp_matrix.py` 逐次執行完整 fixture 與九場比較，
配置間交替、排除 repetition 0，至少三次才彙整時間。每次更新的
`matrix.json` 先寫暫存再替換；數值通過與能否宣稱 speedup 分別記錄。
它無法排除其他使用者工作，因此正式量測仍需確認資源空閒或取得 allocation。

## 固定幾何 Newton 與接續量測

`fixed-newton/` 的原生 depth-4、兩個 flow controllers 與 pressure gauge
案例，在 1 thread 下到達 1,200 s 時限，wrapper 退出碼 124；沒有完成的
收斂報告或有效 RSS，所以這項仍未通過。完整退出紀錄在 `rank-0/run.json`。

相同案例、16 次 Newton iteration cap、既有 KSP／非線性／線性
殘差門檻的 4-thread 驗證已結束，退出碼 134。測試來源只新增環境控制的 phase profile、
階段標記與即時輸出；未更動模型或驗收條件。建置
`build-fixed-profile.log` 通過且無編譯 warning。`fixed-newton4/inputs.json`
記錄 binary／來源 hashes 與資源設定，對應檔案另行封存在該目錄的
`source-and-binary.tar.gz`；前述 143-file archive 保留較早測試版本。

`fixed-newton4/rank-0/run.json` 記錄 wall time 1,016.982153 s、peak RSS
107,245,568 bytes。七次 Newton 更新後非線性殘差為 `4.1071489592466233e-11`，
低於 `1.0945723152615691e-9` 門檻，但最後一次線性求解的真實相對殘差
為 `2.5883782137720832e-4`，未通過原生 `1e-10` gate；KSP 回報正的
convergence reason 並不足以證明真實殘差合格。此結果屬於失敗證據。

`outputs/hpc02/volume/validation-sequence.json` 已記錄 `fixed_newton_failed`；
後續 `isolated-matrix/` 沒有啟動。舊工具 session 7941 在接續查詢時已不存在，
`session-state.json` 也已改正，不能繼續把它當作執行中的工作。
相同原生 fixture 的 Jacobian／RHS 已擷取。兩種 Gram–Schmidt 候選仍失敗，
FGMRES 的真實相對殘差為 `2.0549559498475232e-15`，通過相同 `1e-10`
門檻。runtime 已採用 FGMRES 並更新 solver configuration hash。完整
Newton 已退出 0，七次線性求解的真實相對殘差皆小於 `2.06e-15`，所有
原生 gates 通過；1／2／4／8-thread FSI 的九場回歸也均已通過，見
[線性求解修正](HPC_02_LINEAR_SOLVE_PROGRESS.md)。1／4-thread 隔離矩陣
已各完成首次與三個正式重複；八份 collection、72 個場比較與來源雜湊
已重新核對。4-thread 速度比 2.300008、RSS ratio 1.021657 均符合
門檻，見 [FSI 效能驗收](HPC_02_FSI_PERFORMANCE.md)。

## 完成核對

1. 修正後 FSI 的 1／2／4／8-thread 場回歸與固定幾何 Newton gate 已通過；
   無 OpenMP 的 `compliant_channel_fsi_test` 也已重建並通過完整九場比較。
   來源、binary、建置命令與原始證據均已封存。
2. 隔離重複量測已完成，端到端中位數改善與 RSS ratio 均通過既定門檻；
   最終完成狀態須一併記錄無 OpenMP 完整場回歸的證據。
3. 貼體 MPI／OpenMP 組裝與固定四核心比較已完成，見
   [混合組裝進度](HPC_02_HYBRID_PROGRESS.md)。保留與 FSI 的獨立證據；
   浸入式 FSI 目前仍不支援 MPI。貼體相同 rank 數、單 thread 基準已完成，
   OpenMP 每 rank 與合計 RSS 增量比值皆不超過 1.009，通過既定預算。
   固定核心數比較的單 rank 記憶體取捨仍完整保留。

最終逐項核對記錄在 `outputs/hpc02/volume/completion-audit.json`，
涵蓋 worker／caller 機制、完整物理場、無 OpenMP、固定核心數比較及 RSS。

完整 collector 命令與執行資源說明見
[HPC_BENCHMARKS.md](../HPC_BENCHMARKS.md#complete-immersed-and-fsi-reference-states)。
本次未改變數值核心公式、`.ntiga` 或 coupling ports／commit 契約。
