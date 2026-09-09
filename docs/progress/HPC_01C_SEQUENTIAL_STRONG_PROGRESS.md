# HPC-01C：Sequential strong 的錯誤、收斂與回復協調

接續 [sequential 初始化](HPC_01C_SEQUENTIAL_INITIALIZATION_PROGRESS.md)。
本批補齊 fixed／Aitken 強耦合迴圈的共同 outcome；HPC-01C 仍需完成支援入口
覆蓋稽核，整份清單保持進行中。

## 實作

`RunSequentialFlow` 現在協調 Begin、輸入準備、三個 runtime 的 SolveTrial、
3D port 測量、本地 residual／relaxation、precommit、PrepareCommit 與
accepted-step bookkeeping。每個本地階段的例外先形成共同結果，才進入
下一個 runtime 操作；收斂使用所借用 communicator 的全群 AND 決策。

新增 `SequentialRuntimeStage` 只包住一個原生操作，協調該操作返回或丟出的
結果。操作內部的 MPI 安全性仍由各 runtime 的既有協議負責；不能把多個
collective 或先行配置工作放進這個 callback。`RuntimeConstructionStage`
則只執行本地準備、驗證與診斷。

失敗時依序嘗試 downstream、3D、upstream 的 AbortStep，保留 primary
exception 與三個 abort exception，完成全部清理嘗試後才組合診斷。
各 runtime 的 PrepareCommit 成功並共同確認後，才連續執行 noexcept
FinalizeCommit。Postcommit 錯誤共同退出，保留已提交的物理狀態。

物理公式、收斂容許值、fixed／Aitken 更新、CSV／JSON schema 與 `.ntiga`
介面未變更；增加的同步成本尚未作效能量測。

## 驗收證據

證據目錄：`outputs/hpc01/sequential-strong/`（ignored）。

| 檢查 | 結果 |
|---|---|
| 原生故障與健康重試 | fixed／Aitken 各 27 個故障目標，world 3 ranks 及獨立 1／2-rank 子群，共 162 次故障、162 次重試 |
| 物理狀態回復 | 每個失敗案例比較兩個完整 1D flow state 及各 rank 的本地 3D 向量；precommit 比較進入 step 前的快照，postcommit 比較已接受快照；並檢查 trial phase 已關閉 |
| 清理與診斷 | 除尚未開啟 trial 的 ready 故障外，每次都嘗試全部三個 abort；額外同時注入三個 abort outcome 錯誤，檢查 primary 與各清理診斷均存在 |
| PETSc ownership | 每組累計觀察 648 個非空 object references，在外部參照釋放前確認 runner ownership 已解除 |
| 全群收斂 | 六組測試都在最後一個 rank 否決首次達標的 trial，確認所有 ranks 延後提交、accepted iteration 相同，CSV 記錄正確的全群決策 |
| 健康重試相容性 | 每次重新建立 runtimes、使用同一 communicator 重跑三步；324 份 CSV 與各組健康 baseline 逐位元組相同 |
| 修改前後相容性 | 1／2／3 ranks、fixed／Aitken，共 12 份 CSV 與初始化批次保存的 production binary 逐位元組相同 |
| 原有完整 smoke | 更新檢查器後重新建置並通過；包含 explicit／graph、subcycling、fixed／Aitken 等價與既有拒絕案例 |

`native-final.log` 保存六組通過紀錄，`native-final/` 保存每次故障診斷、
baseline、retry 與 veto 輸出；`comparison.json`／`comparison-runs.json`
保存前後比較與原始命令。所有故障要求 runner 返回 1；健康重試與整個原生
測試返回 0，timeout 不算通過。這些是小型正確性測試，不是 scaling 證據。
`smoke-final.log` 保存完整 smoke 通過紀錄；production、native test 與
smoke build 無 compiler warning。`source-final.json`／archive、
`after-binaries/` 與 `acceptance.json` 保存來源、執行檔與證據身分。

第一次執行 `native.log` 在完成 world/fixed 的 27 個故障／重試後，因檢查器
未考慮注入的全群否決而失敗。現行檢查器只允許明確指定的第一步 iteration
被否決，且要求該 trial 本來達標、有後續 iteration、serialized flag 為 0；
其他數值、工作量與收斂檢查不變。原始失敗輸出保留供追溯。

## 重現命令

本機 WSL、GCC 11.4、Open MPI 4.1.2、PETSc 3.15.5 real64／Int32、MUMPS。
原生測試使用 OMP／BLAS 各 1、core mapping／binding：

```bash
make -C solvers/coupling iga_1d_3d_explicit sequential_strong_failure_test \
  explicit_coupling_smoke_test \
  PETSC_DIR=/usr/lib/petscdir/petsc3.15/x86_64-linux-gnu-real
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 \
  PETSC_OPTIONS='-ksp_type preonly -pc_type lu -pc_factor_mat_solver_type mumps -ksp_rtol 1e-12' \
  timeout --kill-after=5s 600s mpiexec --map-by core --bind-to core -np 3 \
  solvers/coupling/sequential_strong_failure_test /absolute/path/to/new-output
```

原有 smoke 可從 `solvers/coupling/` 執行 `./explicit_coupling_smoke_test`。
設定 `TUBULARFLOWIGA_KEEP_TEST_OUTPUT=1` 與指向新目錄的 `TMPDIR` 可保留案例。

## 剩餘工作

- 完成支援入口的錯誤邊界稽核。已定位 `iga_solve.cpp` 的 Bezier VTKHDF
  初始化仍使用 root `std::exception` catch 加整數廣播，且 catch 中先輸出
  診斷；需改用本地共同階段並補驗收。
- 核對 standalone FSI ParaView exporter 的單程序執行限制與輸出建立順序；
  目前未將它當作支援分散式 FSI 的證據。
- 此批注入的是指定本地階段／原生操作完成後的例外，不是逐 allocator 點的
  窮盡測試，也沒有新增程序失聯恢復或 committed-step 原地回復。
- HPC-01D、持久化協議、浸入式／FSI 分散計算及實際跨節點驗收繼續依原清單。
