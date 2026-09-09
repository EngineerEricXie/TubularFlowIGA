# HPC-01C：1D staged trial 與 OpenMP 物種例外

日期：2026-09-08。狀態：**本報告範圍通過；HPC-01C 整體仍部分完成。**

## 實作與行為

`OneDFlowRuntime` 的 staged hydraulic preparation、substep preparation、
本地 flow solve、frame capture、completion 與兩種 staged rollback，現在
使用既有 `FailureAgreement` 協調局部例外。implicit advance 可能包含 MPI
collective，保留在本地 work callback 外；此變更尚未涵蓋其 PETSc 內部。
hydraulic completion 被其他 rank 拒絕時也會清除 successful-trial flag，
避免失敗的流場被拿來進行 scalar replay。

staged transport 在暫存狀態中重播所有 hydraulic frames。只有整個程序群
完成協調，才以不拋出例外的 swap／move 發布物種、inlet 與 phase。
某 rank 的物種更新失敗時，所有 rank 都保留已完成的 hydraulic 狀態及原有
scalar image；解除故障後可以重試 transport，不必重算 flow。
這需要額外一份暫存物種狀態，記憶體成本仍需納入 HPC-00C／06A 的大型案例量測。

實際 NaN inlet 故障揭露 `AdvanceOneDTransport` 的既有 OpenMP 物種迴圈
會讓例外逸出平行區域並直接 terminate，即使 `OMP_NUM_THREADS=1` 也如此。
現在每個 species worker 在自己的 exception slot 捕捉例外，離開平行區域後
按物種順序重新拋出第一個錯誤，再交給 runtime 的程序群協調。worker 不呼叫 MPI。
未啟用 OpenMP 的路徑保留，不增加核心 runtime 對 MPI／PETSc 的相依性。

不更改流場／物種公式、coupling scheme、既有輸入驗證或輸出／checkpoint
格式。staged transport 原有不支援 concentration-driven vasodilation 的限制
保留；測試只在記憶體中的配置關閉該功能，沒有修改原生 fixture。

## 驗收

環境：WSL、Intel i9-14900KF、GCC 11.4、Open MPI 4.1.2、PETSc 3.15.5
real/double、32-bit PetscInt、MUMPS。HEAD
`ee8da2ab528849814b6d12c8186476b5f7f59124` 加未提交工作樹；精確來源見 inventory。
MPI unit 使用 3 ranks、每 rank 1 或 4 個 OpenMP threads，`OMP_DYNAMIC=FALSE`，
OpenBLAS threads 為 1。graph 回歸每 rank 1 thread。

```bash
make -C solvers/one_d core-test iga_1d one_d_trial_failure_test \
  PETSC_DIR=/usr/lib/petscdir/petsc3.15/x86_64-linux-gnu-real -j2
make -C solvers/coupling petsc multidomain_subcommunicator_test \
  PETSC_DIR=/usr/lib/petscdir/petsc3.15/x86_64-linux-gnu-real -j2
OMP_NUM_THREADS=4 OMP_DYNAMIC=FALSE OPENBLAS_NUM_THREADS=1 \
  PETSC_OPTIONS='-ksp_type preonly -pc_type lu -pc_factor_mat_solver_type mumps' \
  timeout 60 mpiexec --oversubscribe -np 3 \
  solvers/one_d/one_d_trial_failure_test examples/one_d/multispecies_physiology
```

同一 unit 另以 `OMP_NUM_THREADS=1` 執行。兩次皆由實際 parallel region
計算 thread 數並輸出三個 rank 的成功紀錄；不是只記錄環境變數。

| 驗收 | 結果 |
|---|---|
| 預設 C++ core／coupling／runtime tests | 全部通過 |
| 不含 `-fopenmp` 的 C++17 runtime test | 編譯與執行均通過，無 MPI／PETSc 相依 |
| 既有 combined 八種故障 | world 及獨立 1+2 groups，1／4 threads 均通過 |
| staged 七種故障 | 五個 hydraulic 邊界、transport 完成後共同拒絕、單 rank 真實 NaN oxygen，均通過 |
| hydraulic 失敗 | 驗證提前停止 implicit advance、拒絕 scalar replay；Abort 精確回復 packed snapshot，解除故障可重試 |
| scalar 失敗 | packed flow／scalar、root concentration、phase 保持不變；重試不增加 hydraulic solve 次數 |
| scalar rollback／retry | 精確回到 hydraulic snapshot，再次 replay 產生相同 packed 狀態 |
| staged 數值 | area、flow、pressure、outlet 量與每個 species 分別比較 COMM_SELF 參考，沿用 relative L2 `1e-6`／零參考 absolute L2 `1e-12` |
| 最終 flow graph 1+2 groups | 最大 relative L2 `7.45082e-15`，通過 |
| 最終 species graph 1+2 groups | 最大 relative L2 `3.35727e-15`，通過既有物種守恆門檻 |

每個 fault mode 在 world 執行一次；獨立單 rank group 執行一次、兩 rank
group 執行兩次，避免相同次數意外掩蓋 world collective。staged trial 每次
包含三個 configured substeps。兩個 graph 群組使用既有不同 coupling／routing
與 PETSc options，保留獨立 communicator 的驗收方式。

## 證據與限制

ignored `outputs/hpc01/one-d-staged/` 保存：

- `unit-summary.json`、`unit-threads-1.log`、`unit-threads-4.log`。
- `serial-summary.json`／`serial-runtime.log`，含無 OpenMP 的完整編譯命令。
- `subcommunicator-final.json`、`flow-groups-final.log`、`species-groups-final.log`
  與 `run_final_graphs.py`，含確切命令、環境、退出碼與 binary／log hashes。
- 最終 `openmp`、`coupling-final`、`serial` build logs，無 compiler warning／error；
  sandbox 下 MPI wrapper 的 socket 訊息與實際 C++ 編譯診斷分開判讀。
- `evidence-summary.json` 與 `inventory-final.json`，供核對來源、環境、binary 及產物。
  初版 `inventory.json` 保留；inventory 工具拒絕覆寫既有紀錄。

初次機械修改的編譯失敗、負的有限濃度未觸發既有 validation 的測試失敗、
以及修正前 NaN 導致 OpenMP abort 的 log 均保留。負濃度 fixture 改為 NaN，
沒有為測試改變低階 validation 規則。`subcommunicator.json` 是 OpenMP 修正
前的歷史結果，最終驗收使用 `subcommunicator-final.json` 與相符 binary hash。

本次沒有重跑完整 CLI／checkpoint batch、CUDA 或跨節點 scaling；
[前次 trial 報告](HPC_01C_ONE_D_TRIAL_PROGRESS.md) 的結果屬其記錄版本。
沒有量測加速、組裝／求解時間或峰值記憶體的新主張。

後續 implicit PETSc 的呼叫／組裝／callback 補強與新版本驗收，見
[implicit 進度](HPC_01C_ONE_D_IMPLICIT_PROGRESS.md)。以下保留本報告版本當時的剩餘範圍。

當時仍需補 implicit PETSc 建構／組裝／callback、adapter port／BeginStep／
bookkeeping 與其他 CPU／coupling CLI 錯誤邊界。staged rollback 的錯誤協調
不代表記憶體配置失敗下的完整原子回復。所有 group 成員仍須執行相同階段及
substep 數；本地階段 callback 不能救回已卡在 collective 的 rank，亦不支援
process-loss／OOM-killer 恢復。這次沒有新增或量測 OpenMP 組裝加速，
因此不勾選 HPC-02A。完整清單依 [WORKSTATION_HPC_TODO.md](../WORKSTATION_HPC_TODO.md) 繼續追蹤。
