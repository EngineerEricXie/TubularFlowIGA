# HPC-04A：PETSc family options 共用層

基準 revision：`e332973f5693ef4ae0a19548a950abbbac278443`。
本紀錄只驗收 `PetscSolverOptions` 共用層；immersed runtime 整合另行驗證，HPC-04A
尚未完成。測試環境為本機 GCC 11.4、OpenMPI 4.1.2、PETSc 3.15.5 real64／int32。

Constructor 最後新增可選 `inherited_prefix` 與 `inherit_unprefixed`，既有 caller
預設仍繼承 root options。優先順序為 runtime defaults、可選 root 基線、family、
完整 domain prefix。這讓只接受 family options 的 runtime 保留原有相容性。

選項 key 依 PETSc 的 ASCII 大小寫規則正規化後比較，避免大寫 explicit domain
option 被小寫 family alias 覆蓋。Value bytes 不變，來源 key 的 usage flag 仍精確
對應。Prefix、fallback policy 與 effective snapshot 均在 runtime communicator
內取得一致；舊的 private database／error-handler scope 與重試行為保留。

```bash
make -C solvers/cpu petsc_solver_options_test \
  PETSC_DIR=/usr/lib/petscdir/petsc3.15/x86_64-linux-gnu-real
mkdir -p NEW_REPORT
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 mpiexec -np 3 \
  python3 scripts/hpc_rank_run.py --output-dir NEW_REPORT --expected-ranks 3 \
  --timeout 180 -- solvers/cpu/petsc_solver_options_test
```

最終共用層驗證見 `outputs/hpc04/immersed/unit-v4/`，3 份 rank report 均 exit 0。
同一 executable 執行 world 3 ranks，再分成獨立 1／2-rank 群組，涵蓋：

- root／domain 繼承、snapshot 不受來源改寫影響、exact value bytes 與 usage flags。
- family fallback、忽略 root KSP、大小寫 family／domain override、nested PC。
- 無效 prefix、rank-local prefix／source／fallback-policy 不一致的 collective 拒絕。
- block-Jacobi 與 additive／Schur fieldsplit 的 child options，實際 PETSc 錯誤、
  C++ exception、scope 重入拒絕，以及 error handler 還原後健康重試。

SPD 測試保留 exact solution 的絕對 L2 `< 1e-10` 與正 convergence reason 門檻。
這是選項語義與錯誤處理測試，不是大型求解器效能證據；RSS 與每次 iterations／error
保存在 rank logs。此提交不改數值 kernel、場或 checkpoint 格式。
