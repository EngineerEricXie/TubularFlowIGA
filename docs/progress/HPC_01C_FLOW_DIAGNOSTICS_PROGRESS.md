# HPC-01C：流場參考流量與摘要查詢

日期：2026-09-08。狀態：**這兩項診斷邊界通過驗收，HPC-01C 仍部分完成。**
接續 [3D staged 報告](HPC_01C_THREE_D_STAGED_PROGRESS.md)。

## 實作

[TransientFlowRuntime.hpp](../../solvers/cpu/include/TransientFlowRuntime.hpp)
的 ReferenceBoundaryFlow 在積分前確認全群查詢相同 label，再協調本地 nodal
buffer、connectivity indexing、quadrature 與積分錯誤。boundary velocity 使用
有邊界檢查的 indexing；MPI reduction 留在 callback 外，結果必須有限。
仍保留查詢不存在 label 時回傳零的既有語義。

Summary 先協調 local Vec 長度與 ownership、RAII read view 及有限場值檢查，
成功歸還陣列後才呼叫 VecNorm 和 MPI reductions。VecNorm 的回傳碼由既有
helper 協調，最後共同拒絕非有限的 state／velocity／pressure norms。
有限輸入累加造成平方和 overflow 也會拒絕，不發布無效摘要。

兩個方法仍是整個 runtime communicator 的共同查詢。保留原有 norm 計算、
owned element／node 歸屬及 linear iteration counter，沒有改變解、檔案格式、
求解器預設或 port 方向。可寫 State() 的使用者仍必須遵守 PETSc 的陣列操作契約；
本地例外協調不保證修復使用者造成的 PETSc 內部狀態不一致。

## 測試

延伸 [test_three_d_port_failure.cpp](../../solvers/coupling/tests/test_three_d_port_failure.cpp)：

| 新增模式 | 結果與重試 |
|---|---|
| 最後一個 rank 查詢不同合法 label | 在積分前共同拒絕；單 rank 不執行此模式 |
| owner rank 的 connectivity 改為 -1 | 本地積分階段共同拒絕；修正後 reference inlet 與原值相同 |
| owner rank 的 reference velocity 含 NaN | reduction 結果共同拒絕；恢復後查詢通過 |
| 最後一個 rank 的 state 含 NaN | read view 範圍內檢查失敗，全群在 VecNorm 前拒絕；陣列仍可重新借用並恢復 |
| 最後一個 rank 的有限 state 值為 double max | 平方和／norm 溢位被拒絕；恢復後摘要通過 |

測試在獨立 3、1、2-rank groups 執行，每群各 25、14、25 個案例，共 **64 個**，
包含原有 50 個 port／adapter 錯誤案例。各群有獨立 COMM_SELF 參考；全群錯誤
訊息一致，正常 retry、非零流場求解及 rollback／commit 仍通過。

正常摘要另比對解析離散 norms。64 個控制點的 velocity=(1,2,3)，
pressure=5+x+2y+3z，x/y/z 分別取 0、1/3、2/3、1，因此 velocity 平方和為
`64*14`，pressure 平方和為 `64*(64+14*5/36)`，完整 state 為兩者之和。
在故障前及每次恢復後皆檢查，iteration counter 保持零。
沿用 relative `1e-6`、零參考 absolute `1e-12`，未放寬數值標準。

同一版 production header 重建後，以下回歸也全部通過：

- staged failure：72 個案例，含物種流向反轉與 accounting。
- trial failure：48 個案例，含 outlet 與快照回復。
- flow graph 1+2-rank groups：最大 relative L2 `7.45082e-15`。
- species graph 1+2-rank groups：最大 relative L2 `3.35727e-15`。
- VCA runtime：原有 PETSc 預設，lifecycle／staged assertions。
- VCA smoke：MUMPS，1／2-rank 比較、reservoir、同 rank checkpoint/restart
  state 精確比較及 oxygenator。

上述七項作業皆在各自 180 秒 timeout 內退出 0，成功標記及 binary／log hashes
均已核對。沒有新的 compiler warning。

## 首輪測試修正與限制

第一版故障注入只在目標 rank 呼叫 VecGetArray／VecRestoreArray，其他 rank
保留先前 norm cache 狀態，作業在 180 秒由 timeout 結束、退出 124。
這違反可寫陣列操作的 logically collective 契約，與 norm cache 分歧造成
collective 不一致的機制相符；未取得各 rank stack，故不宣稱已直接定位卡住的指令。

修正為全群借用／歸還陣列，只有目標 rank 改值。舊作業確認 terminal 後才重新
執行，最終 64 個案例通過。未改動 production norm 演算法來適應錯誤測試。
參考 [PETSc VecGetArray](https://petsc.org/release/manualpages/Vec/VecGetArray/)
與 [VecRestoreArray](https://petsc.org/release/manualpages/Vec/VecRestoreArray/)
的共同呼叫要求；實際驗收使用 PETSc 3.15.5。

這些檢查不模擬程序消失、MPI 內部故障或所有 PETSc 後端的 read-view 失敗。
新增 local-size guard 有正常路徑覆蓋，沒有破壞 PETSc vector layout 來注入長度錯誤。
沒有新的效能、RSS、GPU、多核心加速或跨節點證據。

## 重現與證據

環境：WSL 工作站、GCC 11.4、Open MPI 4.1.2、PETSc 3.15.5 real/double、
32-bit PetscInt、MUMPS。HEAD `ee8da2ab528849814b6d12c8186476b5f7f59124`
加未提交工作樹。這是小型故障／數值回歸，沒有更改 production geometry。

```bash
make -C solvers/coupling petsc three_d_port_failure_test three_d_staged_failure_test \
  three_d_trial_failure_test multidomain_subcommunicator_test \
  PETSC_DIR=/usr/lib/petscdir/petsc3.15/x86_64-linux-gnu-real -j2
make -C solvers/cpu iga_navier_stokes vca_3d_runtime_test vca_3d_smoke_test \
  PETSC_DIR=/usr/lib/petscdir/petsc3.15/x86_64-linux-gnu-real -j2
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 \
  PETSC_OPTIONS='-ksp_type preonly -pc_type lu -pc_factor_mat_solver_type mumps' \
  timeout 180 mpiexec --oversubscribe -np 3 \
  solvers/coupling/three_d_port_failure_test /tmp/flow-diagnostics-fresh-output
```

案例目錄必須尚不存在，MPI runtime 需要允許本地 sockets。
本地 ignored 產物位於 `outputs/hpc01/flow-diagnostics/`：

- `regression/mpi-summary.json`、logs 與案例：五項 MPI 作業的命令、環境、退出碼與 hashes。
- `vca/summary.json` 與 logs／案例：兩項 VCA 回歸。
- `attempt-1.log`、`attempt-1-source.cpp`、`attempt-1-identity.json`：首輪錯誤測試及 124 退出紀錄。
- `evidence-summary.json`、`build/`、`inventory.json`：成功標記核對、輸入／輸出身分及 source／環境。

## 尚待完成

HPC-01C 仍需 constructor／初始化、全域 gather／checkpoint I/O、flow Newton
回傳碼、完整配置／外部資產一致性，以及其他 adapter／executor／CLI 邊界。
本批未更新任何子任務為完成；整份 goal 仍保留全部 38 個子任務。

後續 Newton 的 PETSc 操作、收斂與 logging 邊界已補強，75 個 trial 案例與
同版回歸見 [Newton 進度](HPC_01C_NEWTON_PROGRESS.md)；其他剩餘工作仍未完成。
