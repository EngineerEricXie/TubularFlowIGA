# HPC-03A：正式靜態物理 operator 與全域診斷

日期：2026-09-09。基準 HEAD `ed129fd` 加本批工作樹。
狀態：本報告記錄 operator 組裝階段；後續 Newton／守恆與 HPC-03A 完成驗收見
[靜態 runtime 報告](HPC_03A_STATIC_RUNTIME_PROGRESS.md)。

## 實作

[ImmersedStaticFlowSetup.hpp](../../solvers/cpu/include/ImmersedStaticFlowSetup.hpp)
抽出既有序列 runtime 的 options、diagnostics、幾何／catalog preflight、active
node compaction、port ordering 及 scalar row 編號。序列 runtime 私有繼承此共用
setup；數值 kernels、Newton、KSP 與交易狀態流程沿用既有實作。Setup 沒有
PETSc matrix／vector 或 MPI 資源，不會為分散式執行建立完整序列 solver。

[ImmersedStaticDistributedOperator.hpp](../../solvers/cpu/include/ImmersedStaticDistributedOperator.hpp)
建立正式 cell／face／constraint topology，再使用既有 distributed assembly：

- 正值體積的 active cells 依順序分配，active node rows 以完整四場節點分配。
  Ghost face 使用 minus cell 的 owner；controller target 由最後一個 rank 加一次。
  Gauge 與 flow constraint 依 cell 分片，其 scalar owner 不需收集完整場。
- 保留體積 conservative mixed form、全表面 mixed trace、選定 Nitsche wall、
  ghost penalty、pressure／mean-normal-traction／flow-rate ports 的原有 kernels。
  Diagnostic assembly switches 保留 wall delta 與 scalar topology。
- Gauge 權重只在 cell owner 建構時預算一次；數值組裝只讀 required-state halo。
  壓力量測與流量 constraint 使用既有 port measurement 定義，以全域 reduction
  統計 area、outward flow、pressure、normal traction、velocity squared 與 multiplier。
- Port 量測保留序列路徑「走訪所有 surface rules」的範圍。沒有體積貢獻的
  cell 若仍有 port surface 與可用 active nodes，使用只量測、不加物理係數的
  stencil；缺少 active node 則與序列 Gather 一樣拒絕。這種 stencil 仍在 declared
  dense graph 中預留結構，未宣稱 graph 是數學非零項的最小集合。
- 全域 wall 點數／逐 label 點數、cell／face 計數與 physical constant-pressure
  defect 可與序列比較。Gauge row 不納入 physical pressure defect；所有 appended
  scalar 的實際矩陣列均檢查結構零對角，沒有新增 regularization。
- 建構先協調本地 preflight，再比對幾何／quadrature／物理與控制參數。Evaluator
  是本地函數，caller 須在同群提供相同物理定義，不能在 evaluator 中呼叫 MPI。
- 組裝失敗時沿用 backend 的共同 stash drain 與重試；diagnostics 只在所有
  reduction、有限值及結構檢查通過後發布，state 不被組裝修改。Communicator
  與 catalogs 借用至 Close／析構結束，所有成員共同進入建構、Assemble、Close。

Symbolic metadata、active map 與 geometry catalogs 仍複製；persistent numerical
matrix／vectors 與 halo 分散持有。本批沒有 nonlinear distributed KSP、line search、
conservation 或 trial commit／rollback 的完整接入。原有 case／graph MPI size 1
限制與所有外部資料格式維持，不能用本次組裝驗收代表整個 MPI runtime 已完成。

## 驗收範圍

[物理測試](../../solvers/cpu/tests/test_immersed_distributed_physics.cpp) 現在直接
呼叫正式 operator，刪除測試內重複的物理 orchestration。序列 reference 仍使用
既有 runtime；比對完整 residual 與 Jacobian action 的 relative L2 gate `1e-6`，
包括 controller／gauge rows。另比對拓撲、wall counts、port means／flow gates、
pressure measure／defect，並要求 exact declared preallocation、零動態配置。

預設為 27 cells／54 ghost faces／compact depth 2、非零四場與 body force。
Pressure 模式是一個 flow controller 加 pressure outlet；traction 模式使用
mean normal traction 與 pressure、無 flow scalar／gauge。Closed 模式將所有
表面設成 wall；wall-only 關閉 volume／ghost／gauge 的數值貢獻。
Padded 模式有 125 個背景 cells、27 個正值 active cells；expanded 使用展開式
quadrature。Faults 包含單 rank 無效／有效但不同的參數、已有本地積分後的
body-force exception、診斷未發布與 state 未變，以及同物件重試的完整數值比對。

最終物理 executable 的 14 個 MPI 作業、28 份 rank reports 全部退出 0、無 timeout。
各模式最差數值誤差如下；所有全域診斷與 exact preallocation gates 通過。

| 模式 | Ranks | 最大 residual relative L2 | 最大 Jacobian-action relative L2 |
|---|---|---:|---:|
| flow | 1／2／4 | 2.80302e-16 | 7.27974e-16 |
| pressure | 2 | 2.32526e-16 | 7.21201e-16 |
| traction | 2 | 2.343e-16 | 7.31205e-16 |
| closed | 2 | 2.59926e-16 | 6.60885e-16 |
| wall-only | 2 | 1.3011e-16 | 6.97969e-16 |
| faults | 1／2／4 | 2.80189e-16 | 7.27974e-16 |
| padded | 1／2 | 2.48372e-16 | 1.42318e-15 |
| expanded | 1／2 | 2.47667e-16 | 7.84432e-16 |

工作站 GCC 11.4、Open MPI 4.1.2、PETSc 3.15.5 real64/int32。三個 targets
以 `-Wall -Wextra -Wpedantic` 建置，沒有 compiler warnings。MPI 作業固定
OMP／BLAS 為 1；同時最多兩個獨立小案例，以下不是 dedicated scaling benchmark。

| Flow ranks | 最大 assembly s | 各 rank peak RSS bytes |
|---:|---:|---|
| 1 | 7.784927 | 64667648 |
| 2 | 4.167577 | 78020608、69324800 |
| 4 | 2.243229 | 74055680、74723328、70950912、65974272 |

組裝時間包含 halo、局部 kernels、PETSc insertion／assembly 及全域診斷，
不包含 KSP solve。RSS 包含複製的幾何與每 rank 的序列 reference；沒有 CUDA
配置或跨節點測試，不能推論 production 的記憶體或多節點 scaling。

既有 `immersed_static_flow_test` 完整深度 4／深度 2 回歸退出 0，耗時 2121.368 s，peak RSS 197603328 bytes；`immersed_flow_port_test` 退出 0，耗時 114.017 s。
第一次完整 static 回歸達到 600 秒時限而中止，未當作通過；未改動測試內容，
以 2400 秒時限及 OMP／BLAS 1 重跑後完成。超時紀錄保留在 `serial-static/`，
正式結果在 `serial-static-long/`。這兩個既有測試的數值／交易 assertions 均啟用。

## 重現與證據

```bash
make -C solvers/cpu immersed_distributed_physics_test immersed_static_flow_test \
  immersed_flow_port_test PETSC_DIR=/usr/lib/petscdir/petsc3.15/x86_64-linux-gnu-real
python3 scripts/hpc_immersed_assembly_regression.py --kind physics \
  --physics-mode flow --ranks 1 2 4 --output-dir /path/to/new-flow-results
python3 scripts/hpc_immersed_assembly_regression.py --kind physics \
  --physics-mode faults --ranks 1 2 4 --output-dir /path/to/new-fault-results
```

`--physics-mode` 另支援 `pressure`、`traction`、`closed`、`wall-only`、`padded`、
`expanded`。每個輸出目錄必須是新路徑。控制器保存每 rank stdout／stderr、
command、binary hash、退出狀態、timeout 與 peak RSS。

Ignored evidence 在 `outputs/hpc03/static-operator/`。`final/*/summary.json`
對應最終物理 executable；根目錄各模式保留補齊 measurement stencils 前的
探索驗收，舊 binary 留於 `pre-measurement/`，不混為最終 binary 的結果。
`acceptance.json` 彙整正式證據與 source／binary hashes。

## 接續工作

接入 distributed KSP／owned Newton updates、全域 block norms／line search、
守恆與 transactional state，驗證真正 1／2／4-rank 求解與失敗回復，才完成
HPC-03A。之後依原清單加入工作量分區、固定幾何暫態／graph，以及 moving
geometry ownership／history。完整清單仍有 29 項未勾選。
