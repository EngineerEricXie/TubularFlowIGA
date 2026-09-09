# HPC-01C：1D 時間步與 combined trial 的局部錯誤協調

日期：2026-09-08。狀態：**本報告範圍通過；HPC-01C 整體仍部分完成。**

## 實作

原生 `iga_1d` 現在協調 initial diagnostics、每步 inlet／BeginStep、commit
preparation 與 step diagnostics 的本地例外。所有 rank 完成 PrepareCommitStep
後才呼叫 noexcept FinalizeCommitStep，避免部分 rank 已提交而其他 rank
仍因局部前置錯誤離開。原本 output 與 checkpoint 邊界保留。

`OneDFlowRuntime` 新增可選的 `FailureAgreement` 建構參數。combined
`SolveTrial` 的準備、accounting、substep 準備、本地 flow solve、transport
與完成階段先捕捉 `exception_ptr`，再交給 embedding 層協調。
work lambda 不先轉成 `std::function`，避免在進入協調前因轉換配置記憶體。
implicit advance 可能包含 collective，明確保留在 local stage 外。

CLI 與 generic／bifurcation graph 的 combined 1D runtime 已注入借用
communicator 的協調 callback。未提供 callback 時仍為純 C++17 行為，
core／runtime 測試不新增 MPI 或 PETSc 相依性。呼叫者仍須讓同一 group
走相同階段與 substep 數；協調 callback 不會自行修復不同的執行分支。

失敗 trial 會清除 successful-trial flag，包括完成階段的共同拒絕；因此不能
在收到錯誤後仍 PrepareCommitStep。受控故障測試證明 Abort 後可精確回復
flow、species、outlets、動態 radius、step 與 time，再依原有方法重試。
不更改物理公式、時間積分、輸出／checkpoint 格式或既有場數值門檻。

## 驗收與命令

環境沿用 WSL、Intel i9-14900KF、GCC 11.4、Open MPI 4.1.2、PETSc 3.15.5
real/double、32-bit PetscInt、MUMPS。MPI 每 rank OMP／OpenBLAS threads 各 1；
參考求解使用 preonly／LU／MUMPS。HEAD
`ee8da2ab528849814b6d12c8186476b5f7f59124` 加未提交工作樹，精確來源以 inventory 為準。

```bash
make -C solvers/one_d iga_1d one_d_trial_failure_test core-test \
  PETSC_DIR=/usr/lib/petscdir/petsc3.15/x86_64-linux-gnu-real -j2
make -C solvers/coupling petsc multidomain_subcommunicator_test \
  PETSC_DIR=/usr/lib/petscdir/petsc3.15/x86_64-linux-gnu-real -j2
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 \
  PETSC_OPTIONS='-ksp_type preonly -pc_type lu -pc_factor_mat_solver_type mumps' \
  timeout 60 mpiexec --oversubscribe -np 3 \
  solvers/one_d/one_d_trial_failure_test examples/one_d/multispecies_physiology
python3 scripts/hpc_one_d_cli_regression.py --output-dir /tmp/hpc-trial-cli-new
python3 scripts/hpc_one_d_checkpoint_regression.py --output-dir /tmp/hpc-trial-checkpoint-new
```

| 驗收 | 結果 |
|---|---|
| C++ core／coupling／runtime tests | 全部通過，保留無 MPI 相依的預設路徑 |
| trial failure unit | world 與獨立 1+2 groups 均通過六種 local-phase 故障、缺 inlet、缺 implicit callback |
| 提交／回復／重試 | 每種故障皆不可提交；Abort packed 狀態逐位元還原；可重試案例通過 COMM_SELF 參考比較 |
| collective 呼叫數 | 第一個 substep transport 失敗後僅一次 implicit advance；其他 rank 未進入第二次；completion 故障已有三次 advance |
| CLI 21 cases | 八個正常單／三 rank 執行與十三個故障案例全部通過 |
| checkpoint 25 cases | 完整剛體／顯式／六物種續跑、十四種故障與 legacy binary 相容性全部重跑通過 |
| closed-loop CLI | 原生 10-step case 單／三 rank 均通過，flow／species／derived 欄位誤差為 0 |
| graph subcommunicators | 最終接入版本的 flow／species 1+2 groups 均通過，無 world collective 串擾 |

trial unit 的八種模式在 world 執行一次；單 rank group 再執行一次、兩 rank
group 再執行兩次。missing-advance 是永久缺少 callback 的前置錯誤，只驗證
拒絕與 Abort；其餘七種模式會解除故障並重試。retry 的 packed-vector relative
L2 門檻 `1e-6`、零參考 absolute L2 `1e-12` 保留；完整 CLI 回歸另逐欄位檢查。

新增原生 `step-inlet` 案例使用相同 JSON、合法 periodic-table CSV 與相同
初始 oxygen 濃度；只有 rank 1 在第 1 步讀到負濃度。原有 inlet contract
拒絕該值，三 rank 與 launcher 都退出 1。第 0 步 VTP 已存在，第 1 步場及
成功 summary 不存在，證明故障發生於初始化後。harness 同時保留輸入 hashes、
各 rank 真實退出碼、timeout 狀態、RSS 與 log hashes。

四種 CLI 路徑（剛體、顯式、隱式、完整六物種）的單／三 rank 場比較誤差均為 0。
checkpoint 回歸仍維持非時間欄位誤差 0；clock 規則沿用
[既有 checkpoint 報告](HPC_01C_ONE_D_CHECKPOINT_PROGRESS.md)，沒有放寬門檻。
flow／species graph 對各自串行參考的最大 relative L2 分別
`7.45082e-15`／`3.35727e-15`。

CLI batch 最大單 rank wall `0.414158 s`、RSS `40,955,904 bytes`；checkpoint
batch 分別 `0.414208 s`、`41,132,032 bytes`。這些包含啟動與輸出成本，
不是獨立效能重複測試，不能宣稱加速。新增階段協調的成本仍需在 HPC-00C／09D 量測。

## 證據與剩餘範圍

ignored `outputs/hpc01/one-d-trial/` 保存 `cli/summary.json`、
`checkpoint/summary.json`、`trial-unit.json`／log、`closed-loop.json`、
`subcommunicator-final.json`、build logs、`evidence-summary.json` 與 `inventory.json`。
初次重構的編譯失敗 log 保留；權威建置為 final logs，重建與 core tests 通過，
無 compiler warning。graph 注入 callback 前的群組結果另行保留，最終以
`subcommunicator-final.json` 的 binary hash 為準。

後續 staged hydraulic／transport 與 OpenMP worker 例外補強及新版本驗收，見
[staged 進度](HPC_01C_ONE_D_STAGED_PROGRESS.md)。以下是本報告版本當時的剩餘範圍。

**當時仍未涵蓋** implicit PETSc solver 內部的建構／組裝／callback 錯誤、
`SolveHydraulicTrial`／`SolveTransportTrial` 分開執行的 staged 路徑、graph
其他 port／bookkeeping 邊界及其他 CPU／coupling CLI。此 callback 只協調
combined trial 的本地階段；不得宣稱它能救回已卡在 collective 裡的 rank。
沒有 process-loss／OOM-killer 恢復或新跨節點證據。完整持久化與其他 HPC
階段仍依 [待辦清單](../WORKSTATION_HPC_TODO.md) 繼續推進。
