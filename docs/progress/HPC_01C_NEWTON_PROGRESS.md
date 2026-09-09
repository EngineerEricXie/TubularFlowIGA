# HPC-01C：Newton 求解與 PETSc 操作邊界

日期：2026-09-08。狀態：**本批求解邊界通過驗收，HPC-01C 仍部分完成。**
接續 [流場診斷報告](HPC_01C_FLOW_DIAGNOSTICS_PROGRESS.md)。

## 實作

[TransientFlowRuntime.hpp](../../solvers/cpu/include/TransientFlowRuntime.hpp)
在 Newton 迴圈中協調 MatZeroEntries、VecSet、state／history halo scatter、
MatZeroRows、residual／update VecNorm、KSPSetOperators 與 VecAXPY 的已回傳狀態。
KSPSetUp、KSPSetUpOnBlocks、KSPSolve 現在分開檢查，避免某個 rank 在 setup
失敗後提早返回，其他 rank 已進入下一個 PETSc 操作。setup／linear solve 的
原有 profiling 區段保留，MPI/PETSc 操作皆在 local-work callback 外。

Newton 開始時確認 transient 分支一致；有限殘差、容許值、全域守恆診斷及
converged 決策在返回或進入線性求解前協調。KSP reason、iterations、residual
的本地查詢與診斷文字也先協調，負 reason 或非有限結果不會發布 update。
trial iteration counter 在 VecAXPY 成功後才增加，並檢查累加上限。
原有 KSP reason／殘差／守恆錯誤細節保留在共同診斷中。

收斂與每次更新後的 root logging 也加入本地協調；當 ostream 設定為拋出錯誤時，
其他成員會一起進入原有 trial 失敗／Abort 或 Advance 自動回復路徑。
沒有新增檔案格式、改變物理模型或放寬 Newton／mass 容許值。

[OwnedRowAssembler.hpp](../../solvers/cpu/include/OwnedRowAssembler.hpp) 增加帶
communicator 的 Mat／Vec Assemble overload，每個 Begin／End 分別檢查回傳碼。
flow／transport runtime 改用此 overload，舊 overload 的呼叫者仍保留原行為。
[TransientTransportRuntime.hpp](../../solvers/cpu/include/TransientTransportRuntime.hpp)
另補上兩個 MatZeroRows 與 required-state scatter 的回傳碼協調。

這些檢查處理各成員已回傳的錯誤；不宣稱可以救回卡在 MPI/PETSc 內部、程序
消失或無法完成的通訊。初始化、全域 GatherState、viewer/checkpoint 等路徑
仍未全部處理，不能由這批改動推論所有 runtime 操作均已安全。

## 故障與數值驗收

延伸 [test_three_d_trial_failure.cpp](../../solvers/coupling/tests/test_three_d_trial_failure.cpp)，
沿用 64-control-point、單 cubic element、非零 backward Euler fixture。
3、1、2-rank groups 各 26、23、26 個案例，共 **75 個**，包含原有 48 個案例。
單 rank 仍略過三個跨 rank controls／phase 差異案例。

| 新增驗收 | 實際行為 |
|---|---|
| 單 rank pressure coefficient 為 NaN | velocity history 保持有效，組裝完成後的 conservation result 共同拒絕；Advance 回到 committed。測試再修復原本就有 NaN 的輸入，其他係數不變 |
| MAX_NEWTON=1 | 真實更新後耗盡非線性迭代次數；Advance 自動精確回復 committed state |
| root 輸出失敗 | 僅 rank 0 的 throwing streambuf 在更新後輸出時丟例外；全群收到 iteration logging 診斷，Advance 精確回復，ostream 隨後恢復 |
| KSP 負 reason | Richardson／PC none、max_it=1，取得未收斂 reason；禁止 prepare commit，沒有發布 update 或成功 iteration count，Abort 恢復快照 |
| KSP 真實非零回傳碼 | Richardson 大步長引發發散，搭配 error_if_not_converged；KSPSolve 的錯誤在 flow linear solve 階段共同拒絕，prepare／Abort 與 Advance 自動回復均通過 |

兩個刻意失敗的 KSP runtime 分別驗證 SolveTrial、拒絕 prepare、再次 Advance
失敗後的回復。測試只在建構這些 runtime 時暫時設定 options，隨後恢復原值；
並未把放寬容許值或繞過線性求解當成重試成功。原有正常 flow runtime 的
solve／rollback／retry／commit，以及 RC outlet 求解仍通過。

每個錯誤要求全群收到相同診斷，所有 MPI 作業有 180 秒 timeout。最終 trial
測試退出 0；各獨立 communicator 與各自 COMM_SELF 解比較。沿用 relative
`1e-6`、零參考 absolute `1e-12`；快照回復逐位元比較。flow Newton relative
`1e-8`、absolute RMS `1e-12`、mass imbalance `1e-6`，未改動數值 gate。

## 同版回歸與計時完整性

- staged failure：72 個案例，包含反向物種傳輸／收支。
- port／diagnostics failure：64 個案例，包含非有限量與查詢重試。
- flow graph 獨立 1+2-rank groups：最大 relative L2 `7.45082e-15`。
- species graph 獨立 1+2-rank groups：最大 relative L2 `3.35727e-15`。
- VCA runtime：原有 PETSc 預設下 lifecycle／staged assertions 通過。
- VCA smoke：MUMPS，1／2-rank 場比較、reservoir、同 rank checkpoint/restart
  state 精確比較及 oxygenator 通過。

另以 `IGA_PROFILE=1` 重跑原生 VCA smoke，10 份 CLI rank profile 的 status
皆為 0，含 2-rank 案例。setup／linear solve 均有呼叫紀錄；所有 phase 時間有限、
非負，exclusive 不超過 inclusive，exclusive 總和加 unscoped 與 elapsed 相符。
此檢查驗證 phase stack／時間分解，沒有量測新增協調的獨立效能收益、RSS 或 GPU
allocation，也沒有跨節點或 scaling 宣告。

## 首輪測試修正

首輪測試把 max_it=1 與 error_if_not_converged 合用，預期 KSPSolve 回傳非零。
在實際 PETSc 3.15.5 上取得 status 0、reason -3；production 正確在 linear
convergence 階段共同拒絕，但測試期待不同階段，因此整個 MPI 測試退出 1。
原始 log、binary hash 與測試來源已保留，沒有將該次紀錄算作通過。

最終錯誤碼案例改為 Richardson scale=`1e6`、divtol=`1.01`、max_it=2，
以實際發散觸發非零回傳碼；max_it=1 案例仍獨立驗證負 reason 的拒絕。
選項意義依 [PETSc Richardson 文件](https://petsc.org/release/manualpages/KSP/KSPRICHARDSON/)
與 [KSP tolerances 文件](https://petsc.org/release/manualpages/KSP/KSPSetTolerances/)
核對；特定版本的錯誤回傳行為以本地實測為準。

## 重現與證據

環境：WSL 工作站、GCC 11.4、Open MPI 4.1.2、PETSc 3.15.5 real/double、
32-bit PetscInt、MUMPS。HEAD `ee8da2ab528849814b6d12c8186476b5f7f59124`
加未提交工作樹；source／binary 身分另由 inventory 記錄。

```bash
make -C solvers/coupling petsc three_d_trial_failure_test three_d_staged_failure_test \
  three_d_port_failure_test multidomain_subcommunicator_test \
  PETSC_DIR=/usr/lib/petscdir/petsc3.15/x86_64-linux-gnu-real -j2
make -C solvers/cpu iga_navier_stokes vca_3d_runtime_test vca_3d_smoke_test \
  PETSC_DIR=/usr/lib/petscdir/petsc3.15/x86_64-linux-gnu-real -j2
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 \
  PETSC_OPTIONS='-ksp_type preonly -pc_type lu -pc_factor_mat_solver_type mumps' \
  timeout 180 mpiexec --oversubscribe -np 3 \
  solvers/coupling/three_d_trial_failure_test /tmp/newton-fresh-output
```

案例目錄必須尚不存在。建置通過且沒有新增 compiler warning；MPI runtime 使用
允許本地 sockets 的環境。未改動 production geometry，沒有重跑 mesh-test。

本地 ignored evidence 位於 `outputs/hpc01/newton/`：

- `verified/mpi-summary.json` 與 logs／案例：最終五項 MPI 回歸，包含 75 個 trial 案例。
- `vca/summary.json` 與 logs／案例：兩項原生 VCA 回歸。
- `profile-vca/summary.json`、`profile-vca/phase-validation.json`：開啟 profiling 的 smoke 與 10 份 rank 時間分解檢查。
- `regression/`、`attempt-1-trial-source.cpp`、`attempt-1-notes.json`：首輪未通過的測試與原因。
- `evidence-summary.json`、`profile-evidence.json`、`build/`、`inventory.json`：成功標記、檔案 hashes、建置與環境。

未逐一注入每個 Mat/Vec/Scatter 操作的內部失敗；這些回傳碼邊界具有正常路徑
回歸，直接故障證據以上表為限。HPC-01C 尚需 constructor／初始化、全域 gather、
checkpoint I/O、完整配置／外部資產一致性及其他 adapter／executor／CLI 邊界。
整份 38 個子任務的 goal 保持進行中。

後續 flow 初始化與 transport 全域 gather 的補強及同版驗收，見
[初始化進度](HPC_01C_INITIALIZATION_PROGRESS.md)；constructor 與 checkpoint 等邊界仍待完成。
