# HPC-03C：準靜態 immersed graph 的 MPI 接入

日期：2026-09-09。基準 `d68c19f` 加本批修改。
狀態：完成既有準靜態 graph 的 MPI 接入；HPC-03C 仍部分完成，固定幾何
backward-Euler history／暫態求解尚待接續，未移除該模式的限制。

## 實作

[分散式 static runtime](../../solvers/cpu/include/ImmersedStaticDistributedRuntime.hpp)
新增 collective `SetPortControlValue`。先共同驗證 idle／finite／既有 port ID，
再比較完整 ID／value，最後只更新 doubles。單 rank 無效或有效但不同的值
不會部分發布。允許依序供給多個 flow targets；完整 net-flow compatibility
在 Assemble 前檢查，保留既有序列 graph 的輸入方式。`AbortPrepared` 取消
準備狀態，後續 `FinalizeCommit` 不會發布它。

[ThreeDImmersedDistributedFlowDomain.hpp](../../include/ThreeDImmersedDistributedFlowDomain.hpp)
實作相同 `CoupledDomainRuntime` hydraulic port／step 介面。Domain metadata、
step clock 與 port inputs 共同驗證一致性；會分配記憶體的 input／measurement／
prepared map 先在本地 candidate 完成，再共同決定是否發布。每次 collective
runtime 呼叫前完成本地 preflight，不將單 rank 本地例外帶進另一 rank 的求解。

候選失敗或 coupling rollback 會回復 committed vector 與 committed port controls。
Prepare 成功才允許非拋例外的 Finalize 交換 backend state、port measurements、
controls 與 accepted clock。已取消或已關閉的 backend 不會讓 adapter 發布
prepared metadata。`CommittedOwnedBackendState` 明確只回傳本 rank rows；
port-state getters 是本地讀取已全域歸約的量測。

[ImmersedFlowCase.hpp](../../include/ImmersedFlowCase.hpp) 拆成純本地 `Preflight`
與 collective `InitializeDistributed`。前者保留 schema、contained assets、
geometry catalogs 與 label binding 檢查；後者透過 collective allocation 建立
owned runtime 與 adapter，兩者都成功後才移交所有權。多 rank 不建立序列
runtime。舊 `Load(..., mpi_size)` 保留序列入口及 size 1 要求；MPI caller 使用
預檢／初始化介面。`Runtime()` 只適用於序列 owner，對未初始化或 MPI owner
明確拒絕；MPI owner 使用 `DistributedRuntime()`。

正式 [graph entry](../../solvers/coupling/src/iga_1d_3d_bifurcation.cpp) 在共同本地
preflight 結束後，依一致 domain order 初始化 MPI backend。單 rank 保留序列
參考路徑；多 rank 執行真正分散式 Newton。成功輸出前共同 Close，geometry／
metadata 保持可讀。`IGA_PROFILE=1` 可輸出各 rank 的 owned rows／stencils／
required-state rows，供核對實際並行範圍。既有 entry 未啟用 `PhaseProfile`，
因此此診斷亦直接辨識環境設定；沒有宣稱它已補齊 native entry 的所有階段計時。

此 graph 模型仍是 `dt=0` 的 steady solve，只把 coupling time 當作準靜態負載
取樣時間。沒有增加前一步速度慣性項、移動幾何或 transport；既有這些模式
的 schema／runtime 拒絕仍有效。外部 case、SI ports 與 output CSV 介面保持。

## 驗收

本批共 24 個 MPI 作業、55 份 rank reports。21 個作業退出 0，另外三個作業
在指定提交前故障處退出非零且不發布 output directory，全部符合預期且無 timeout。
所有 affected executables 以 `-Wall -Wextra -Wpedantic` 建置，無 compiler warnings；
Python syntax 與 whitespace checks 通過。

- 正式 operator 的 controls 模式：weighted 1／2／4 ranks，拒絕單 rank unknown
  ID、NaN、不同有效值／ID，拒絕未完整供給的 net-flow targets，供齊後與更新
  相同控制的序列 residual／Jacobian action 及 diagnostics 比較，gate `1e-6`。
- Static Newton flow 四 ranks：active trial 拒絕 port 更新；取消 prepared state
  不發布且可重新 prepare／commit；沿用 HPC-03A 的場解、守恆與故障回復 gates。
- 新 graph-domain C++ 測試：`4×1×1` 有效 cells、112 active nodes、451 rows，
  兩個 flow controllers 與 gauge。1／2／4 ranks 及獨立 1+2 groups 通過。
  包括 metadata／step／input 不一致、候選更新後失敗與同物件重試、prepare
  失敗、rollback、prepared abort、兩步提交，以及 Close 後不發布 prepared
  metadata。逐場 relative L2 最大 `6.59917e-10`，門檻 `1e-6`；port／守恆比較通過。
- 原有 `immersed_case_factory_test` 全部接受／拒絕 assertions 通過，56.515 s，
  peak RSS 39,469,056 bytes。該序列測試內容未改動。

正式 CLI 使用從 `immersed_aneurysm_chain` 複製的 source／sink 1D 設定，將 immersed
surface 換成單位立方體、四個有效 cells 與 depth 2 quadrature，執行兩個 graph
steps。這是小型 regression fixture，不是原有 depth 3 aneurysm closure 重驗。
CLI 配置是一個 flow inlet、一個 pressure outlet，無 gauge，共 449 rows。

| Graph mode | 兩步 coupling iterations（各 ranks 相同） | 最大 port error／gate |
|---|---|---:|
| explicit | 1、1 | 7.84387e-6 |
| fixed | 26、1 | 1.22505e-4 |
| aitken | 6、1 | 1.22539e-4 |

每種模式均執行 1／2／4 ranks、2-rank 第 2 步 precommit failure，以及隨後成功的
新程序重試。新程序重試與 C++ adapter 的同物件 retry 是不同驗收。Accepted ports
逐項比較 time／area／outward flow／mean pressure，gate 為 `1e-12+1e-6*abs(reference)`；
表格是 error 除以此 gate，非 relative L2。Flow interface normalized residual
必須 `<=1e-10`；fixed／aitken 的 normalized pressure residual 必須 `<=1e-6`。
Explicit 保留原有單次 sweep，不宣稱壓力 fixed point 已收斂。

兩 ranks 的 owned rows 為 224／225；四 ranks 為 112／112／112／113。每個 rank
有積分 stencil，required-state buffer 小於全域 vector。Native CLI 最大 rank
peak RSS 為 53,878,784 bytes。這是數值與交易驗收；不以作業 wall time 宣稱
scaling。Static Newton 記錄 assembly／linear solve 時間；native CLI 此批主要
核對 ownership、CSV 與每 rank 資源，沒有完整 solver-phase timing benchmark。
未使用 CUDA，也沒有跨節點證據。

## 重現與接續

```bash
make -C solvers/cpu immersed_distributed_physics_test \
  immersed_distributed_static_flow_test immersed_distributed_domain_test \
  PETSC_DIR=/usr/lib/petscdir/petsc3.15/x86_64-linux-gnu-real
make -C solvers/coupling iga_multidomain_flow iga_1d_3d_bifurcation \
  immersed_case_factory_test PETSC_DIR=/usr/lib/petscdir/petsc3.15/x86_64-linux-gnu-real
python3 scripts/hpc_immersed_assembly_regression.py --kind physics \
  --physics-mode controls --partition weighted --output-dir /path/to/new-controls
python3 scripts/hpc_immersed_domain_regression.py --split \
  --output-dir /path/to/new-domain-results
python3 scripts/hpc_immersed_graph_regression.py --execution fixed \
  --output-dir /path/to/new-graph-results
```

Graph controller 另支援 `explicit`／`aitken`。每次使用新輸出目錄。工作站為
GCC 11.4、Open MPI 4.1.2、PETSc 3.15.5 real64/int32／MUMPS，OMP／BLAS 1。

Ignored evidence：`outputs/hpc03/graph-controls/`。`physical`、`newton`、
`final-domain`、`serial-factory`、`final-graph-*` 是正式結果；`domain` 是增加
Close 後 metadata 防護前的測試，`graph-fixed` 是修正 ownership 診斷條件前的
探索結果。後者計算退出 0，但 controller 因缺少 ownership observation 判為失敗，
不混入正式驗收。`acceptance.json` 與 `commit.json` 保存來源、binary／input／output
hashes 與提交核對。

HPC-03C 接續固定幾何 backward-Euler 的 owned committed velocity history、全域
守恆與 port／trial semantics，再驗證對應入口。HPC-03D 另處理移動 active set、
ownership／halo 更新與 history extension。完整清單仍有 27 項未勾選。
