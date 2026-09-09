# HPC-03A：分散式靜態 Newton 求解與守恆

日期：2026-09-09。基準 `bc5e871` 加本批修改。

## 實作與介面

[ImmersedStaticDistributedRuntime.hpp](../../solvers/cpu/include/ImmersedStaticDistributedRuntime.hpp)
接上正式 distributed operator，使用 owned PETSc vectors 保存 trial、committed、
prepared、Newton update 與線性殘差。數值積分仍由 cell／face 唯一 owner 執行，
所需場值使用 halo；正式 runtime 不建立完整序列求解器或收集完整場。
Geometry catalogs、symbolic topology 與 active map 仍複製。

Newton 使用全域 residual／block norms、flow-controller 門檻與 backtracking。
KSP 預設為 GMRES、右預條件、unpreconditioned norm，以及單／多 rank 一致的
MUMPS LU；有效 backend 會在求解前檢查。既有 `immersed_static_` prefix 可覆寫
KSP／PC。LU shift 只用於預條件器，不修改物理矩陣；另記錄實際 `J delta - rhs`
殘差、KSP reason／iterations 與組裝／線性求解時間。

所有 MPI 操作在 caller 提供的 communicator 上共同進入。Catalogs 與 communicator
借用至 Close／析構結束；evaluator 不得呼叫 MPI，各 rank 須提供相同物理定義。
`SetCommittedOwnedState` 接受本 rank owned rows，先共同驗證大小與有限值。
`SolveTrial` 失敗會共同回復 committed vector。`PrepareCommit` 成功才允許
非拋例外的 `FinalizeCommit` 交換 handles；Rollback／Close 丟棄未發布的準備狀態。
`State`／`CommittedState` 是供唯讀檢查的 borrowed handles。

`ConservationDiagnostics` 每次更新 required-state halo，再唯一積分體積散度與
逐 boundary label 的表面外向流量，歸約 total／wall／port 流量；不依賴前一次
Assemble 的場值。一般 `Diagnostics` 的物理量代表最近成功組裝，不因 Rollback
自動重新組裝；需要當前場守恆時應呼叫 `ConservationDiagnostics`。

此介面支援 C++ 靜態求解。既有 case／graph／暫態 CLI 的單 rank 限制保留，
接入與放行屬 HPC-03C；移動 active set 屬 HPC-03D。外部格式與序列路徑不變。

## 驗收設計

[回歸測試](../../solvers/cpu/tests/test_immersed_distributed_static_flow.cpp)
在每 rank 建立既有 `COMM_SELF` 序列 reference，僅供比較。主要幾何為 27 個
active cells、54 個 ghost faces、depth 2 compact quadrature：

- closed：全壁面、body force `(0.3,-0.2,0.1)` 與 pressure gauge。
- flow：inlet／outlet 目標 `-1e-6`／`1e-6` m³/s、wall 與 gauge。
- pressure：inlet flow `-1e-6` m³/s、outlet pressure `1e-3` Pa，無 gauge。
- empty-work：與背景重合的單位立方體、`2×1×1` cells／321 rows，四 ranks 中
  rank 1／3 有 owned rows，但沒有積分 stencil 或 halo；rank 2 注入 candidate
  failure，rank 3 注入 prepare failure。

序列與 MPI 使用相同 `nonlinear_absolute_tolerance=1e-14`、
`nonlinear_relative_tolerance=1e-11`、`lu_pivot_shift=1e-12`。
預設小流量容差曾在首步 Newton 後停止，速度相對差約 `1.9e-6`；本次收緊
比較案例的求解容差，維持驗收誤差門檻。原生 LU 的探索結果曾無法通過實際
線性殘差門檻，因此正式 MPI 預設採 MUMPS；不是僅更換 KSP norm 就解決。

完整場及各非零場 relative L2 gate 為 `1e-6`。Closed／empty-work 的解析速度
與 gauge multiplier 為零，兩者改用 absolute L2 gate `1e-10`，避免接近零的
reference norm 造成無意義的相對誤差。各 Newton step 的實際線性 relative
residual 不得超過 `1e-6`，KSP reason 必須正值，flow controller 必須通過門檻。
守恆 aggregate 與逐 label 結果對序列 gate 是 `1e-12 + 1e-6*abs(reference)`；
這是離散路徑一致性檢查，不表示 depth 2 積分已精確滿足連續散度定理。

故障驗收涵蓋已做非零 Newton 更新後的單 rank evaluator exception、回復零
committed state、同物件重試、單 rank prepare failure、成功兩階段提交、單 rank
錯誤長度／NaN input 不改變狀態、prepared rollback、Close 後不發布 prepared
state，以及 Close 後明確拒絕 Rollback。另在不 Assemble 的情況輸入零場，要求
守恆診斷立即變零，以檢查 halo 新鮮度。三個 world ranks 分成獨立 1+2 groups
驗證 communicator 隔離。

組裝 backend 額外覆蓋三種配置：均分、全部工作／列在 rank 0、工作全在
rank 0 而列全在最後一個 rank。包含零列 owner 的跨 rank stash 插入與失敗重試。

## 重現

```bash
make -C solvers/cpu immersed_distributed_static_flow_test \
  immersed_distributed_assembly_test immersed_distributed_physics_test \
  PETSC_DIR=/usr/lib/petscdir/petsc3.15/x86_64-linux-gnu-real
python3 scripts/hpc_immersed_static_regression.py --mode closed --split \
  --output-dir /path/to/new-closed-results
python3 scripts/hpc_immersed_static_regression.py --mode flow \
  --output-dir /path/to/new-flow-results
python3 scripts/hpc_immersed_static_regression.py --mode pressure \
  --output-dir /path/to/new-pressure-results
python3 scripts/hpc_immersed_static_regression.py --mode empty-work --ranks 4 \
  --output-dir /path/to/new-empty-results
python3 scripts/hpc_immersed_assembly_regression.py --kind unit --split \
  --output-dir /path/to/new-unit-results
python3 scripts/hpc_immersed_assembly_regression.py --kind physics \
  --output-dir /path/to/new-physics-results
```

控制器限制每 rank 630 s／launcher 660 s，固定 OMP／BLAS 1，記錄 binary hash、
每 rank stdout／stderr、退出狀態、timeout、peak RSS 與時間。輸出目錄必須為新路徑。
工作站使用 GCC 11.4、Open MPI 4.1.2、PETSc 3.15.5 real64/int32／MUMPS。

## 正式結果

以下結果來自 `outputs/hpc03/static-runtime/acceptance/`。四個靜態模式共 11 個 MPI 作業、28 份 rank reports，全部退出 0、無 timeout；
所有逐場、controller、守恆與故障回復 gates 通過。

| 模式 | Ranks | 最大完整場 relative L2 | 最大速度 absolute L2 | 最大壓力 relative L2 |
|---|---|---:|---:|---:|
| closed | 1、2、4、3 (split 1+2) | 4.61891e-14 | 8.1868e-15 | 4.60211e-14 |
| flow | 1、2、4 | 3.33575e-11 | 1.25884e-18 | 3.87249e-14 |
| pressure | 1、2、4 | 9.96051e-14 | 8.08476e-18 | 3.80006e-15 |
| empty-work | 4 | 3.38179e-15 | 1.33814e-15 | 3.19938e-15 |

Flow／pressure 非零速度場的最大 relative L2 分別為 `4.38393e-14`, `2.81552e-13`；非零 scalar 場分別為 `2.38838e-10`, `1.4341e-12`。

| Flow ranks | 最大累計 assembly s | 最大累計 linear s | 最大 rank peak RSS bytes |
|---:|---:|---:|---:|
| 1 | 63.153456 | 0.303692 | 86835200 |
| 2 | 24.310380 | 0.215349 | 98463744 |
| 4 | 13.310130 | 0.216320 | 85381120 |

時間是 runtime 建立後至 accepted trial 的累計成功組裝／線性求解，包含初始
Assemble 與故障 trial 已成功完成的階段；失敗組裝未發布 diagnostics，其耗時
不納入該累計。Wrapper wall time 包括序列 reference、故障、重試及後續交易檢查。
RSS 包括各 rank 的序列 reference；同時最多兩個 Newton 案例，另有短組裝回歸。
這些數據用於驗收紀錄，不能當作 dedicated scaling benchmark，也不代表 CUDA
配置或跨節點效能。

Backend unit 的 1／2／4／split 1+2 共 4 個作業、10 份 rank reports，以及正式
physics operator 的 1／2／4 共 3 個作業、7 份 rank reports 均通過。建置啟用
`-Wall -Wextra -Wpedantic`，無 compiler warnings；Python syntax 與 diff whitespace
檢查通過。既有序列 runtime／setup 未在本批改動，沿用前批完整 static／port
回歸，見 [operator 報告](HPC_03A_STATIC_OPERATOR_PROGRESS.md)。

`acceptance.json` 記錄本批 source／binary hashes 與六組正式 summaries。
`one`、`two`、`review-one`、`true-residual`、`mumps`、`final` 等根目錄保留不同
探索版本與失敗紀錄；不混入正式驗收。對應舊 binaries 亦留在 ignored evidence。

HPC-03A 的靜態 owned rows／halo、正式 Newton、全域守恆及失敗回復已完成。
接續 HPC-03B 工作量分區、HPC-03C 暫態／graph、HPC-03D 移動幾何；HPC-09
跨節點仍須獨立驗收。整份清單尚有 28 項未勾選。
