# HPC-04A 多層求解器診斷

- 狀態：本機多層 options／view 驗證通過，HPC-04A 整體仍進行中。
- 日期：2026-09-09（美東）；基準 revision `953e53a`。
- 只擴充 `petsc_solver_options_test` 的選用測試，沒有更動 production headers、
  solver 預設、物理模型或資料格式。

`-test_solver_multilevel` 建立每 rank 128 rows、兩個 component 的嚴格對角優勢
SPD AIJ 矩陣；每 component 為 16 × (4 × ranks) 網格，相鄰 node 權重 -1、
對角 6、component 間權重 -0.1。已知全一解提供獨立 L2 判準。
外層為 additive fieldsplit，第一個子區塊 GAMG，第二個 LU／MUMPS。
此模型只驗證階層配置，不能代表 flow saddle-point 系統或 AMG 效能。

在 KSP solve 真正建立階層後，逐層查詢 down smoother 的 Richardson／Jacobi
及完整 level prefix，查詢 coarse GMRES／Jacobi、preconditioned norm 及 prefix；
檢查原始 options database 的使用紀錄。診斷不自行提前建立 smoother。
未加 flag 時保留原小矩陣測試與輸出。

## 本機結果

環境為工作站 GCC 11.4、Open MPI 4.1.2、PETSc 3.15.5 real64/int32；
OMP／OpenBLAS threads=1，CPU test 未編入 OpenMP。World 3／4 ranks 各再分成
1+2／1+3 的 communicator。7 個 ranks 在兩次 solve 階段產生 14 筆多層紀錄：
3-rank world 與 singleton／2-rank 子群為 3 層，4-rank world 為 4 層、
其 singleton／3-rank 子群為 3 層。含原有測試共 112 筆場比較，最大 L2
`1.36214e-12`，通過原門檻 `<1e-10`，正收斂原因及非零 iterations。

另有無 flag 的 3-rank 回歸，42 筆解紀錄、最大 L2 `5.43896e-16`，沒有 KSP view。
合計 10 份成功 rank report，最大單 rank host RSS 43,982,848 bytes。
這是微型配置測試；組裝／求解效能、PDE 守恆及 CUDA allocation 為 N/A，
不從這個模型宣稱擴展性或更改 production solver 策略。

第一次 `world-three` 的 coarse GMRES 沿用了 `NONE` norm，外層回報
`DIVERGED_PC_FAILED`／`SUBPC_ERROR`，測試以 MPI_ABORT 退出 2，未列通過。
明確設置 `fieldsplit_first_mg_coarse_ksp_norm_type preconditioned` 後，
`world-three-v2` 通過；最終再增加 getter 與來源使用紀錄斷言後重建，
`final-three`／`final-four`／`default-three` 全部退出 0。
此失敗說明只指定 KSP type 不保證其餘繼承配置適用。

```bash
make -C solvers/cpu petsc_solver_options_test \
  PETSC_DIR=/usr/lib/petscdir/petsc3.15/x86_64-linux-gnu-real
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 \
mpiexec -np 4 solvers/cpu/petsc_solver_options_test \
  -test_solver_multilevel -test_solver_view
```

本機證據為 `outputs/hpc04/multilevel/acceptance.json`、`audit.py`、各 rank 的
`run.json`／stdout／stderr；保存最終測試來源與 binary SHA256。
原 backend/view suite 的來源紀錄仍屬其原 revision，沒有覆寫為本次結果。

下一步是完整 moving 回歸及 HPC-04B／C：在實際 PDE、網格與 rank 配置上，
維持既有數值與守恆門檻，評估候選預條件器的組裝／求解時間、記憶體與失敗組合。
