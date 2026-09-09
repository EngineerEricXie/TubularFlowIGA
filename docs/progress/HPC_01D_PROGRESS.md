# HPC-01D：執行資源與 PETSc 後端檢查

日期：2026-09-08。狀態：本批通過，HPC-01D 部分完成；整份清單範圍不變。

## 實作範圍

[ExecutionResources.hpp](../../solvers/cpu/include/ExecutionResources.hpp) 提供
communicator collective preflight，接入原生 `iga_navier_stokes`、`iga_solve`、
`iga_1d` 及 `RunMultidomainFlow` 入口。Embedding caller 可在數值工作之前，
由 MPI 初始化執行緒呼叫 `RequireExecutionResources(comm, &stream)`。

- 要求 real/double PETSc，並比較 communicator 內的 PetscInt／PetscScalar 位寬。
- 驗證 OMP_NUM_THREADS 正整數清單、OMP_THREAD_LIMIT，及 OpenBLAS／MKL／BLIS
  thread 數；BLAS 的 0 保留為 library default，未設定記為 -1。
- 查詢實際 MPI thread support；編譯了 OpenMP 且 max threads 大於 1 時，要求
  至少 MPI_THREAD_FUNNELED。MPI 已由 caller 初始化時，不嘗試重新初始化或升級。
- 不修改環境、OpenMP／BLAS policy、PETSc options 或 CPU 綁定。不要求各 rank
  使用相同合法執行緒數，摘要保留 min／max。
- CPU flow／transport 的 database 開啟、參數與 partition／rank 檢查改為本地
  準備後協調，會在讀取 mesh 前拒絕不相容分區。Flow 檢查四個 fields 的 row
  容量；transport 的動態 field 容量仍由既有 OwnedRowAssembler 檢查。

`execution_resources` 摘要包含 ranks、PETSc 位寬、MPI thread level 範圍、
OpenMP 是否編譯與 max threads 範圍、BLAS thread requests，以及 root build 的
MUMPS／Hypre flags。MPI level 數值對應使用中 MPI 的 SINGLE／FUNNELED／
SERIALIZED／MULTIPLE constants。Max threads 是 runtime 上限，BLAS requests
不是實測 workers；摘要不證明 CPU／NUMA 綁定或加速。

PETSc 的預設 MPI thread request 與 caller 已初始化 MPI 的情形不同，因此直接查詢
provided level，參見 [PETSc MPI thread 說明](https://petsc.org/release/manualpages/Sys/PETSC_MPI_THREAD_REQUIRED/)。

`RequireFactorBackend` 以實際 Mat 查詢 factor 支援，先共同比較 backend 名稱與
factor kind。`RequireKspFactorBackend` 查詢有效 PC／backend，接到貼體 flow、
in-process transport 的線性求解前，以及 1D nonlinear SNES options 解析後。
只有明確選定的 LU／Cholesky backend 進入 availability 檢查；保留未指定 backend
時由 PETSc 選擇的行為，也不覆蓋 user options。
查詢 API 的本地性見 [MatGetFactorAvailable](https://petsc.org/release/manualpages/Mat/MatGetFactorAvailable/)。

## 驗收證據

| 檢查 | 結果 |
|---|---|
| Resource unit，OpenMP 1／4 threads 與無 OpenMP build | 各 50 個案例，合計 150；world3 與獨立 1+2 groups |
| 實際 MPI_THREAD_SINGLE | OpenMP 4 threads 共同拒絕；OpenMP 1 thread 與無 OpenMP build 正常接受 |
| 原生 resource CLI | 16 項通過，詳見下方 authoritative records |
| 1D 四種 implicit formulation groups | 1+2 ranks 比對各自 COMM_SELF；最大誤差 `4.84544e-16` |
| 1D implicit failure／retry | 既有三群組 solve／callback／vector 故障通過 |
| 1D 原生 CLI | 原有 21 項數值／錯誤案例通過，包含 explicit／implicit／species |
| 3D staged／trial／port | 原有 104／93／64 個案例通過 |
| Flow／species graph groups | 最大 relative L2 `7.45082e-15`／`3.35727e-15` |
| VCA runtime／smoke | 原有預設 runtime；MUMPS smoke 之非零流場、reservoir、1／2 ranks 及同 rank 續跑精確比較通過 |

[Unit](../../solvers/cpu/tests/test_execution_resources.cpp) 在最後一個 rank 設定
非法 thread 環境，涵蓋零、負值、空字串、溢位、錯誤清單與尾端字元；要求共同
diagnostic 並在恢復環境後再次接受。另測合法巢狀清單、不同 BLAS requests、
rank／partition 不符、零 fields、row overflow、root output stream 失敗，
以及不存在、空名稱、跨 rank 不一致和有效 MUMPS 的實際 matrix backend 查詢。

MPI_THREAD_SINGLE 案例在 PetscInitialize 之前改變所需 level，並確認 MPI
確實回報 SINGLE；沒有只修改測試傳入的假 capability 數字。未要求較低 level
後一定有較低 provided level 的通用假設，本機測試明確檢查了實際結果。

[CLI regression](../../scripts/hpc_resources_regression.py) 在四個入口分別注入
單 rank 的非法 OMP／BLAS 設定；CPU flow／transport 另測 partition mismatch
及單 rank 缺失 database，且刻意給不存在的 case 目錄，證明先拒絕 database。
正向案例包含兩 rank 原生 flow／transport，以及 1D 的 1／2-thread 不同配置，
摘要須保留 `omp_max_threads_min=1 omp_max_threads_max=2`。
非零 VCA 案例用 MPI AIJ 搭配不支援的 PETSc 原生 LU，確認共同以 availability
diagnostic 退出 1。正常 MUMPS 配置由 VCA smoke 驗收。

沿用 relative `1e-6`、零參考 absolute `1e-12`，快照／同 rank checkpoint 維持
精確比較。正常 MPI／VCA suite 的 13+2 項作業均在 180 秒內退出 0；CLI
各項由 90 秒 timeout 保護，成功／拒絕退出碼按案例核對。

## 測試修正紀錄

最初 native 測試將 `KSPPREONLY` 套用到獨立 transport。該 CLI 使用非零初始
猜值，因此 PETSc 拒絕此組合；保留失敗 log，正向 transport 測試改為 GMRES＋
MUMPS，沒有放寬數值門檻或改 solver production defaults。

原 adapter fixture 在原生 flow CLI 的 wall tracing 下為零流場，修改 outlet
traction 仍未進入線性求解，因此兩次 backend 拒絕測試未達預期。改用原生 VCA
smoke 的非零案例後 availability 拒絕通過；這兩次零解不作為 backend 驗收證據。
原生 flow 的 quiescent smoke 僅證明入口可用，非零數值仍由 VCA／group tests 證明。

16 項 CLI authoritative records 由 `native-verified/` 的前 14 個通過案例與
`native-final-vca/` 的最後 2 項組成。執行檔 SHA256 相同；修正案例後只重跑剩餘
兩項。`native-authoritative.json` 列出各項來源，失敗 attempt 不改寫成通過。

## 重現

環境：WSL 工作站，GCC 11.4、Open MPI 4.1.2、PETSc 3.15.5 real/double、
32-bit PetscInt、MUMPS。HEAD `ee8da2ab528849814b6d12c8186476b5f7f59124`
加未提交工作樹。受影響 CPU、coupling、1D binaries 重建通過，無新 compiler warning。

```bash
make -C solvers/cpu execution_resources_test execution_resources_serial_test \
  iga_navier_stokes iga_solve vca_3d_runtime_test vca_3d_smoke_test \
  PETSC_DIR=/usr/lib/petscdir/petsc3.15/x86_64-linux-gnu-real -j2
make -C solvers/one_d iga_1d one_d_subcommunicator_test one_d_implicit_failure_test \
  PETSC_DIR=/usr/lib/petscdir/petsc3.15/x86_64-linux-gnu-real -j2
make -C solvers/coupling petsc three_d_staged_failure_test three_d_trial_failure_test \
  three_d_port_failure_test multidomain_subcommunicator_test \
  PETSC_DIR=/usr/lib/petscdir/petsc3.15/x86_64-linux-gnu-real -j2
OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=1 timeout 180 mpiexec --oversubscribe -np 3 \
  solvers/cpu/execution_resources_test
OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=1 timeout 180 mpiexec --oversubscribe -np 3 \
  solvers/cpu/execution_resources_test single
```

Native regression 的 `--case-dir` 使用 staged test 產物中兩 rank 的 `1-1/`，
`--flow-case-dir` 使用 [checkpoint 讀取報告](HPC_01C_CHECKPOINT_READ_PROGRESS.md#重現與證據)
所述保留的 `tubularflowiga-vca-3d-smoke/`。兩者皆需先產生，再執行：

```bash
python3 scripts/hpc_resources_regression.py \
  --case-dir /path/to/staged-output/1-1 \
  --flow-case-dir /path/to/tubularflowiga-vca-3d-smoke \
  --output-dir /tmp/resources-fresh-output
```

Ignored evidence 位於 `outputs/hpc01/resources/`：MPI／VCA／1D CLI summaries、
native attempts／authoritative records、build logs、evidence summary、inventory
保存命令、環境、退出碼、目前 binary 與 log hashes。先前各階段報告保留原 revision
證據，不改成這次重跑的結果。

## 剩餘範圍

後續 CUDA 啟動、版本摘要及索引容量檢查見
[CUDA 資源驗收](HPC_01D_CUDA_PROGRESS.md)；本報告的原始量測維持不變。

HPC-01D 尚需其餘 embedding／工具路徑的完整能力矩陣、scheduler allocation
與綁定配置核對，以及 HDF5 等輸出後端的細部能力檢查。預條件器內部 sub-KSP、
PETSc 自動挑選的 backend 與完整 options 支援矩陣仍依 HPC-04 補齊。
純 C++ 1D、浸入式／FSI 的 embedding callers 尚未自動套用此 MPI 啟動摘要。

沒有執行新的跨節點、GPU、64-bit／complex PETSc、效能或擴展性驗收；沒有以
宏定義、max threads 或 mpiexec 可啟動宣稱加速。Geometry 未變，未重跑 mesh-test。
其餘 HPC-01C 邊界與整份 38 項清單繼續追蹤，沒有縮減 goal 範圍。
