# HPC-01C：Registry 所有權與 graph 建構協調

日期：2026-09-08（本機 EDT）。HPC-01C 及整份清單仍未完成。
接續 [species executor](HPC_01C_SPECIES_EXECUTOR_PROGRESS.md)。

## 問題與修正

`DomainRuntimeRegistry` 原先以值接收 runtime owners，在檢查 graph metadata
前即取得所有權。單一 rank 驗證失敗會提前析構 runtime，其他 ranks 卻可能
繼續建構與求解。以修改前 header、三個 fake runtime 及 3 ranks 重現：
`failed_ranks=1 runtimes_destroyed_before_agreement=3`。這證明提前析構，
沒有把 fake runtime 重現宣稱為原生 PETSc deadlock。

現在 constructor 以 rvalue reference 借用 caller owners，先驗證並建立本地
索引，經可選的 outcome callback 同意後，才以 noexcept move 接管 vector。
任何驗證或索引配置失敗都保留 caller 的原始 runtime 指標，允許同一批物件
修正輸入後重試。純 C++ caller 可省略 callback，沒有引入 MPI 依賴。
新的 `CreateCollectiveDomainRuntimeRegistry` 使用借用的 communicator，
分別協調 callback 配置、registry storage 與 index 建構。

原生 graph 及 sequential explicit caller 使用該 factory，並將本地 adapter
catalog 配置加入建構協調。Graph 的 boundary label／coverage、初始 port
catalog、transport preflight／node agreement 也加入本地錯誤邊界。
含 collective 的 reference-flow measurement 及 runtime construction 保持在
local callback 外。Immersed owner 在預留容量後才移轉，先建立 audit pointer。
沒有修改物理公式、容許誤差、CSV 或資料庫格式。

## 驗收與證據

本批 ignored 證據位於 `outputs/hpc01/registry/`。

| 驗收 | 結果 |
|---|---|
| Registry ownership test | world 3 ranks、獨立 1／2-rank 子群各 9 組，共 27 次故障及同一批 owner 重試 |
| 原生 graph construction | flow、species、0D 三路徑在上述三種 communicator 配置，共 69 次故障及健康重試 |
| 非空 PETSc 清理觀察 | 各配置 flow／0D 各累計觀察 168 個 object references，species 384 個；每次失敗／重試均要求 count > 0，釋放外部參照前確認 runtime 所有權已解除 |
| 新舊數值相容性 | 三路徑 × 三種 rank 數，57 份 CSV 逐位元組相同 |
| Sequential smoke | 原有 explicit／strong-fixed／Aitken 等案例通過，保留輸出 |
| 純 C++ contracts | `make -C solvers/coupling test` 通過 |

Ownership test 涵蓋 kind、domain、ports、count、null metadata，以及 factory
的同步、storage、index 故障。每次失敗都確認 runtime 尚未析構且指標不變；
重試成功後僅析構一次。每個 tracked runtime 持有真正的 MPI Vec，另留 PETSc
參照以檢查 registry 清理。原生測試採 stage hook 注入 caught `bad_alloc`，
確認共同失敗、未建立 output 目錄，再在同一 communicator 重新建構健康案例。
這不是 allocator 每個呼叫位置的窮盡測試，也不代表程序失聯後可以恢復。

環境：本機 WSL、GCC 11.4、Open MPI 4.1.2、PETSc 3.15.5 real64／Int32，
OMP／BLAS 各 1。Native 新舊比較都使用 GMRES、LU／MUMPS、`ksp_rtol=1e-12`。
Sequential smoke 使用其既有每案例設定。沒有放寬原有數值門檻。

```bash
make -C solvers/coupling petsc test registry_failure_test \
  graph_registry_failure_test explicit_coupling_smoke_test \
  PETSC_DIR=/usr/lib/petscdir/petsc3.15/x86_64-linux-gnu-real
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 \
  timeout --kill-after=5s 120s mpiexec --map-by core --bind-to core -np 3 \
  solvers/coupling/registry_failure_test
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 \
  PETSC_OPTIONS='-ksp_type gmres -pc_type lu -pc_factor_mat_solver_type mumps -ksp_rtol 1e-12' \
  timeout --kill-after=5s 240s mpiexec --map-by core --bind-to core -np 3 \
  solvers/coupling/graph_registry_failure_test /absolute/path/to/new-output-directory
```

`native-build-final.log`、`native-final.log`、`unit.log`、`runs.json`、
`comparison.json` 保存驗收紀錄。初次 unit build 缺少明確的 `petscvec.h`
include，已修正；失敗 log 保留，最終 build 無 compiler warning。
`source-final.json`／tar 與 `after-binaries/` 保存本批來源與執行檔。

## 剩餘工作

Sequential strong 迴圈、CSV／manifest close 與支援入口覆蓋稽核仍待完成。
本批沒有重跑全套 immersed smoke；其前一個完整受測版本見 pressure-executor
報告。沒有新增 GPU、跨節點、效能或 committed-step 回復驗收。
清單保持 9／38，完整 checkpoint 發布與恢復協議由 HPC-05 繼續追蹤。
