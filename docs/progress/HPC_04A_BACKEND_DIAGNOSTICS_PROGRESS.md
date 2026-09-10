# HPC-04A 後端與診斷驗證

- 狀態：本機限定範圍通過，HPC-04A 整體仍進行中。
- 日期：2026-09-09（美東）。Runtime 基準 revision `01aa4c8`；本次新增測試、
  harness 與文件，沒有改 production solver、數值預設或資料格式。
- 環境：TsungYehLab 工作站、GCC 11.4、Open MPI 4.1.2、PETSc 3.15.5
  real FP64／32-bit index；`OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1`。
  這是本機功能驗證，沒有跨節點或 scaling 宣稱。

## Factor 能力

`petsc_factor_matrix_test` 在 world 3／4 ranks，再分成 1+2／1+3 子群，
建立每 rank 兩列的 SPD AIJ 矩陣與已知全一解。逐一比較 `MatGetFactorAvailable`
與實際 runtime guard；可用組合須真正求解、回報正收斂原因及 L2 < `1e-10`，
不可用組合須在 preflight 拒絕，不能冒充已求解。

48 個唯一的 rank／backend／PC 組合中，20 個可用、28 個按預期拒絕。
包含子群重複觀測時共 72 筆，32 次成功求解、40 次預期拒絕。

| Backend | 1 rank LU／Cholesky | 2、3、4 ranks LU／Cholesky |
|---|---|---|
| petsc | 可用／可用 | 不可用／不可用 |
| mumps | 可用／可用 | 可用／可用 |
| superlu_dist | 可用／可用 | 可用／可用 |
| superlu | 可用／不可用 | 不可用／不可用 |
| umfpack | 可用／不可用 | 不可用／不可用 |
| unavailable_backend | 不可用／不可用 | 不可用／不可用 |

這張表記錄本機 build 註冊的 factor 介面，不推論 backend 內部實際分解算法，
也不表示這些選項適用於 Navier–Stokes saddle-point 系統。

## 巢狀與 native viewer

`petsc_solver_options_test -test_solver_view` 在 3／4 ranks 及各自子群驗證
block-Jacobi、additive／Schur fieldsplit 的實際子 KSP prefix、子求解器收斂原因及
LU／MUMPS view。共 98 筆解向量紀錄，最大 L2 `5.43896e-16`；另有未開 viewer
的 3-rank 42 筆紀錄通過且沒有 KSP view。原測試預設仍不印 viewer。
尚未把多層 GAMG 的所有 smoother／coarse 配置納入此矩陣。

`scripts/hpc_solver_views.py` 在同一 frozen native binary 上比較 nonlinear_aq、
implicit_1d_pde、species、immersed static／transient 各自的 baseline 與 viewer，
共 10 個正向 3-rank 作業。另有一個 viewer 無法開啟目的檔的負向作業：退出 1，
無 accepted output，沒有 MPI_ABORT。KSP／SNES view 顯示實際 domain prefix；
immersed family alias 也正確映射到 domain。

1D／species accepted files 逐 byte 相同。Immersed 的既有數值 gate 保留：
port 差異 ≤ `1e-12 + 1e-6 * abs(reference)`、flow residual ≤ `1e-10`、
pressure residual ≤ `1e-6`，步數／時間／iteration count 一致，edge 值有限，
檔案集合及非 CSV metadata 一致。最大 port 差異僅占允許誤差的
`1.8118967444608749e-06`；總共 28 個檔案逐 byte 相同。
此新增 viewer 驗證沒有輸出或比較完整體積場，也沒有另測 checkpoint。

初版 harness 要求 immersed CSV 逐 byte 相同而停止；隨後重跑相同無 viewer
baseline 也有五個 CSV 的末位差異（最大 port 差異占 gate `2.2854249358137396e-06`）。
因此 final harness 使用既有數值 gate，沒有放寬它；初版紀錄不列通過。

## 重現與證據

```bash
make -C solvers/cpu petsc_factor_matrix_test petsc_solver_options_test \
  PETSC_DIR=/usr/lib/petscdir/petsc3.15/x86_64-linux-gnu-real
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 \
  mpiexec -np 4 solvers/cpu/petsc_factor_matrix_test
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 \
  mpiexec -np 4 solvers/cpu/petsc_solver_options_test -test_solver_view
python3 scripts/hpc_solver_views.py --binary NATIVE_BINARY --output-dir NEW_OUTPUT
```

本機 `outputs/hpc04/backend-audit/acceptance.json` 彙整 47 份成功 rank report：
`matrix-final-v2`、`matrix-four-v2`、`view-final`、`view-four`、`default-final`、
`native-views-v2`。`audit.py` 驗證退出碼、timeout、RSS、factor 能力、解誤差及 view
內容，並保存測試來源 SHA256。Native binary 為
`outputs/hpc04/immersed/native-v4-binary`（revision `01aa4c8`）；harness 保存 binary、
builder 與輸入 hashes，防止混用。每 rank `run.json` 保存 argv、環境、wall time、RSS
及 log hashes。微型 factor／options fixture 不量測組裝／求解效能，因此分階段效能
比較為 N/A；native 原有 log 的時間也只作功能作業觀測，並行執行不可視為效能基準。

## 剩餘工作

舊版 `iga_transport` options 隔離已完成，見 [legacy transport 驗證](HPC_04A_LEGACY_OPTIONS_PROGRESS.md)；完整 moving 回歸在 rigid wall trace gate 失敗，尚待修復。
多層診斷及 HPC-04B／C 的預條件器、網格與大小案例評估仍待完成。
