# HPC-01C：貼體 3D trial 與快照回復

日期：2026-09-08。狀態：**下列 runtime 邊界通過驗收；HPC-01C 仍部分完成。**
接續 [3D port 報告](HPC_01C_THREE_D_PORT_PROGRESS.md)。

## 實作

flow `BeginStep` 先共同驗證 phase 與控制值，準備 boundary／traction／outlet
候選快照，再比對步數、時間與 Newton／mass 容許值。PETSc state／history
複製成功後才發布本地快照與 trial phase。`SolveTrial` 的 phase、配置存在性與
outlet 複製準備也先協調，避免一個 rank 還沒有配置而其他 rank 已進入求解。

flow `Advance` 的本地 boundary configuration 更新現在獨立協調；失敗時所有
rank 都進入既有快照回復路徑。`RestoreCommittedSnapshot` 先準備本地候選資料，
檢查 state／history 複製回傳碼，再 swap 發布本地快照，保留 trial-ready／
committed 的既有呼叫語義。

flow 與 transport 的 `RollbackTrial`、`PrepareCommitStep` 與 `BeginStep`
前置驗證加入協調。兩者 `AbortStep` 在 committed 提早返回前，先檢查所有
成員的 phase 是否一致；若呼叫者只在部分 rank finalize，現在會共同拒絕
Abort，不會讓其餘成員等待 VecCopy。`FinalizeCommitStep` 仍為本地 noexcept
提交，仍要求呼叫端在共同 prepare 成功後對每個成員提交。

新增 `RequireCollectivePetscSuccess`：PETSc 操作先在 callback 外完成，才協調
已回傳的錯誤碼。接入上述 flow／transport 快照複製及 transport 的矩陣清空、
右端項運算、KSP 設定／求解與 state 發布。transport 的 warm-start 分支先確認
全群一致；convergence reason 的本地檢查也先協調。此 helper 不處理仍卡在
PETSc／MPI 內部的程序，沒有把 collective 本身包進 local-work callback。

outlet trial 的本地初始資料、boundary 更新、表面積分與模型 evaluation 也
分階段協調，並確認 outlet 數量、label 順序與收斂決策一致。模型的 relax／
commit 先在候選 vector 上完成，整群成功才發布。積分使用 RAII read view。
flow `GatherRequiredVelocity` 與 transport `GatherRequiredState` 同樣加入本地
配置／read-view 協調；其 scatter 仍在 callback 外。

flow adapter 的 BeginStep 在原生 BeginStep 成功後才更新 step／清除輸入。
staged adapter 的 BeginStep、prepare commit、三種 rollback 與 transport solve
phase 驗證，也先完成本地錯誤協調再進入已協調的 runtime。其餘 staged 輸入、
transport boundary materialization 與 accounting 路徑仍未全部處理。

沒有更改物理公式、檔案格式、port 單位、耦合方式或既有求解器預設。
新增候選複製與 collective 有成本，尚未量測大型案例的效能影響。

## 故障與數值驗收

新增 [test_three_d_trial_failure.cpp](../../solvers/coupling/tests/test_three_d_trial_failure.cpp)，
使用既有單 cubic element、64-control-point fixture。flow 使用非零入口、
backward Euler；transport 為單位立方體內初值 2、均勻 source 1、無通量邊界，
`dt=0.01` 的解析濃度為 `2.01`。另使用真實 RC outlet 模型驗證求解後失敗。

| 測試 | 結果 |
|---|---|
| Flow BeginStep | 單 rank 非法 Newton 次數、各 rank 合法但不同的時間被拒絕；已提交場與 phase 不變 |
| Flow Advance | 單 rank 不支援的 profile 共同拒絕，整群回到 committed，場逐位元等於原快照 |
| Flow SolveTrial preflight | 只有一個 rank 未提供 trial configuration；共同拒絕後仍為 trial-ready，補上配置可重試 |
| 不合法 prepare／rollback | flow 與 transport 都在進入原生回復或發布前共同拒絕 |
| PETSc 回傳碼 | 僅一個 rank 對本地 `COMM_SELF` KSP 設定非法 rtol，取得真實錯誤碼；helper 將該錯誤協調給整群 |
| 分歧 Abort | flow／transport 各讓一個 rank 提早 finalize，再共同 Abort；兩者均在提早返回／複製前拒絕。補完其餘 rank 的本地 finalize 後，共同 Abort 正常保持 committed 場 |
| Transport 前置輸入 | 單 rank velocity node 清單不符；場與 trial-open phase 不變 |
| Transport 內部 boundary 失敗 | 單 rank 使用 scalar 不支援的 boundary kind；trial 失敗後禁止 prepare，Abort 回到已提交場與步數 |
| Outlet evaluation | 僅一個 rank 中、由測試持有的 model kind 被改成不支援值，真實 flow 求解完成後共同拒絕；候選模型不發布，禁止 prepare，Abort 精確還原流場 |
| 正常重試 | flow／scalar rollback 精確恢復快照；再次求解符合原始結果；RC model 修正後能重新求解、收斂並提交 |
| 數值比較 | 解析濃度、非零 inlet flow、單／多 rank 的 port area／flow／pressure、scalar mass、outlet flow／pressure／capacitor pressure 通過 |

3-rank 群組與 2-rank 子群各 **17 個案例**，單 rank 子群略過三個跨 rank 差異
案例，執行 14 個，合計 **48 個案例**。所有錯誤案例都要求全群收到相同階段
診斷；整個 MPI 作業有 180 秒 timeout，正常退出 0。不同子群執行不同數量的
collective，並各自比較獨立 `COMM_SELF` 參考。

場與各物理量分別使用既有 relative `1e-6`／零參考 absolute `1e-12`。
快照 copy／Abort 比較保持逐位元相等，不將重新求解的浮點誤差混入快照比較。
flow Newton relative `1e-8`、absolute RMS `1e-12`、mass imbalance `1e-6`；
RC outlet 保留既有 fixed-point 容許值與迭代上限。

另重跑：

- 既有 3D port failure test：20 種模式，在 3／1／2-rank groups 共 50 個案例通過。
- flow graph 的獨立 1+2 groups：最大 relative L2 `7.45082e-15`，通過。
- species graph 的獨立 1+2 groups：最大 relative L2 `3.35727e-15`，通過原有守恆 gate。
- `vca_3d_runtime_test`：使用原有 PETSc 預設，原有 lifecycle／staged assertions 通過。
- `vca_3d_smoke_test`：MUMPS，1／2-rank flow／pressure／reservoir 比較、同 rank
  checkpoint／restart 的 flow 與 species state 精確比較、oxygenator 案例均通過。

前次強制 MUMPS 的額外 runtime unit 逐位元斷言問題保留在
[前次報告](HPC_01C_THREE_D_PORT_PROGRESS.md)，此次沒有改動或放寬該斷言，
也沒有重跑該額外配置。

## 重現與證據

環境：WSL 工作站、GCC 11.4、Open MPI 4.1.2、PETSc 3.15.5 real/double、
32-bit PetscInt、MUMPS。HEAD `ee8da2ab528849814b6d12c8186476b5f7f59124`
加未提交工作樹，精確 source／binary 身分由 inventory 記錄。

```bash
make -C solvers/coupling petsc three_d_trial_failure_test three_d_port_failure_test \
  multidomain_subcommunicator_test \
  PETSC_DIR=/usr/lib/petscdir/petsc3.15/x86_64-linux-gnu-real -j2
make -C solvers/cpu iga_navier_stokes vca_3d_runtime_test vca_3d_smoke_test \
  PETSC_DIR=/usr/lib/petscdir/petsc3.15/x86_64-linux-gnu-real -j2
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 \
  PETSC_OPTIONS='-ksp_type preonly -pc_type lu -pc_factor_mat_solver_type mumps' \
  timeout 180 mpiexec --oversubscribe -np 3 \
  solvers/coupling/three_d_trial_failure_test /tmp/three-d-trial-fresh-output
```

最後一個輸出目錄必須尚不存在。沒有新的 compiler warning；sandbox 內 compiler
wrapper 的 socket probe 提示不影響建置，MPI runtime 在允許本地 sockets 的環境執行。

本地 ignored evidence 為 `outputs/hpc01/three-d-trial/`：

- `verified/faults-summary.json`、`verified/faults.log`：48 個案例的最終測試版本、
  命令、timeout、退出碼、source／binary／log SHA256。
- `regression/mpi-summary.json` 與 logs：port、flow graph、species graph。
- `vca/summary.json` 與 logs／保留的 smoke 產物：原生 runtime 與續跑回歸。
- `evidence-summary.json`、建置 logs、`inventory.json`：證據核對與環境身分。

`attempt-1` 是新 scalar fixture 的 source term 未指定 trial field，`attempt-2`
則漏了必要 boundary definitions；兩者在參考 fixture 準備時拒絕，已修正。
`attempt-3` 通過 48 個案例，最後再增加未發布 outlet 候選值的精確檢查；
第一次 `final/` 檢查誤用了要求 committed phase 的 mutable outlet getter；
已改用 const getter 檢查失敗 trial 的唯讀狀態，權威結果以 `verified/` 為準。沒有更動數值容許值。

這是失敗處理與數值回歸，不是效能或跨節點驗收。未量測此版本的每 rank RSS、
GPU allocation 或隔離重複效能樣本；solver logs 的 assembly／linear 時間僅保留
作診斷。沒有新的 GPU、多節點或 OpenMP 加速宣告。

## 仍待補齊

- constructor 後段、初始化、全域 gather 與 checkpoint I/O；flow Newton 內仍有
  未檢查的 PETSc collective 回傳碼。
- flow reference-flow／summary、transport mass／source 等其餘 reduction 的本地邊界。
- runtime 的完整物理配置／外部資產一致性，並非僅本次控制值、outlet count／label
  與分支決策的一致性。rank 內 PETSc 設定也仍需由 embedding preflight 保證。
- staged adapter 的輸入發布、transport boundary materialization、物種 accounting，
  executor bookkeeping 與其餘 CLI。
- MPI 內部失敗或程序消失後的恢復不在本地例外協調的保證範圍，持久化續跑另依 HPC-05。

上述剩餘事項未完成，因此 HPC-01C 與整份清單皆保持未完成。

後續已補強 staged 輸入／邊界／accounting 與 transport mass／source 積分，
72 個故障案例及 flow／species／VCA 回歸見
[3D staged 進度](HPC_01C_THREE_D_STAGED_PROGRESS.md)。其他剩餘邊界仍須完成。
