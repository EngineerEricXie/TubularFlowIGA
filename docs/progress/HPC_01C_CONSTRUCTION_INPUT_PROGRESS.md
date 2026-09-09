# HPC-01C：建構呼叫端配置與參數準備

日期：2026-09-08。狀態：已實作，277 個原生建構故障／重試案例、
trial／graph／bifurcation CLI 與 VCA CLI 續跑回歸通過。
HPC-01C 的其他錯誤邊界與外部資產一致性尚未完成。

## 修正範圍

接續 [PETSc 建構清理](HPC_01C_CONSTRUCTION_PROGRESS.md)，處理在建構子
開始協調前，`make_unique` 配置物件或按值傳入容器可能只在一個 rank
拋出例外的問題。其他 rank 此時可能已進入 OwnedRowAssembler 的 collective。

- `AllocateCollectiveRuntime` 先在本地階段配置原始物件空間，再共同確認；
  全員成功後才呼叫會進入 MPI 的建構子。未發布的空間由 RAII 釋放，
  成功物件交由原有 `unique_ptr` 管理。適用於使用一般 host alignment、
  沒有自訂 operator new/delete 的 runtime。
- 貼體 flow 的邊界、標籤、速度、wall basis 與 outlet model 改以 const reference
  接收，內部複製放在本地錯誤協調階段；runtime 仍持有自己的資料。
- transport 的 compiled system、configuration 與 flow 的初始 trial configuration
  透過空 `optional` 在本地階段建構，再以編譯期確認不拋例外的 move 發布。
  quadrature map 的複製也在 transport 本地準備階段進行。
- 正式 graph 的 flow／transport 與 CPU VCA transport 使用新配置工廠。
  stack 建構的 CPU／explicit flow 使用相同的參數準備修正。

呼叫者仍須先協調會配置記憶體的引數運算；工廠不會保護進入函式之前
就已執行的運算。建構子也仍須自己協調內部操作，不能把整個 collective
建構子包進 `CollectiveLocalStage` 的本地 callback。

## 驗收

擴充 `solvers/coupling/tests/test_runtime_construction_failure.cpp`，在最後一個
rank 的配置前、配置後、部分 flow 輸入複製後，以及 configuration／compiled
system 準備後注入 `std::bad_alloc`。各案例確認共同診斷、尚未取得 PETSc
物件，並接續正常建構與數值重試；保留原有 PETSc reference-count 清理檢查。

| communicator ranks | 正常 flow／transport | 空列 flow／transport | fieldsplit／pattern | 總數 |
|---:|---:|---:|---:|---:|
| 3 | 23／21 | 23／21 | 5／1 | 94 |
| 1 | 22／20 | 22／20 | 5／0 | 89 |
| 2 | 23／21 | 23／21 | 5／1 | 94 |

合計 **277 個案例**，包括 world 三 ranks 與獨立一／二-rank groups。
正常 fixture 比較 COMM_SELF 參考的速度、壓力與物種場，沿用 relative L2
`1e-6`／零參考 absolute L2 `1e-12`；常數物種 source 解維持 `2 + dt`
的 `1e-10` 門檻。空列 fixture 只驗證建構及初始場；不宣稱流場求解驗收。

環境沿用 GCC 11.4、OpenMPI 4.1.2、PETSc 3.15.5 real64／Int32 與 MUMPS。
OMP／BLAS threads 為 1；MPI 使用 `--map-by core --bind-to core`。
原生測試退出 0，stderr 為空，外層 timeout 180 秒，未放寬既有數值標準。

```bash
make -C solvers/coupling runtime_construction_failure_test \
  PETSC_DIR=/usr/lib/petscdir/petsc3.15/x86_64-linux-gnu-real
env OMP_NUM_THREADS=1 OMP_THREAD_LIMIT=1 OPENBLAS_NUM_THREADS=1 \
  MKL_NUM_THREADS=1 BLIS_NUM_THREADS=1 IGA_ASSEMBLY_THREADS=1 \
  IGA_ASSEMBLY_BATCH_SIZE=1 PETSC_OPTIONS= \
  timeout --kill-after=5s 180s mpiexec --map-by core --bind-to core -np 3 \
  solvers/coupling/runtime_construction_failure_test /tmp/construction-input-fresh
```

正式 CPU `petsc`、`iga_navier_stokes_openmp`、coupling `petsc` 和受影響的
既有回歸已重建，沒有 compiler warning。同版 trial 回歸通過 93 個故障案例，
獨立 flow／species graph 群組最大 relative L2 分別為 `7.45082e-15`／
`3.35727e-15`，bifurcation CLI 也通過。VCA CLI 的流場、傳輸、守恆、
單／雙 rank checkpoint 與不中斷／續跑比較通過，退出 0；完整輸出保存於
`vca-cli/`，另以 `vca-cli/result.json` 記錄 launch、環境與 binary hashes。

## 保留的失敗與證據

第一版 helper 要求 Value 的預設建構不拋例外，編譯檢查拒絕 CompiledLinearSystem。
修正為空 optional 在協調階段 emplace；沒有刪除安全門檻，仍要求發布用 move
不拋例外。原始編譯失敗日誌保留為 `build-default-guard-failure.log`。

第一次既有 trial 回歸未帶入原驗收所需的 LU/MUMPS 選項，發生參考量比較失敗。
封存的修改前 binary 在相同設定下重現相同失敗；新版本使用原驗收文件的
`-ksp_type preonly -pc_type lu -pc_factor_mat_solver_type mumps` 後通過。
這不證明預設迭代求解器在該 fixture 通過同一門檻；求解器矩陣另依 HPC-04 驗收。

證據根目錄：`outputs/hpc01/construction-inputs/`。
`source.json`／`source.tar.gz` 保存 251 份同版來源，`change.patch` 記錄相對前次
正式建構驗收的六檔修正。`native/evidence.json` 保存案例計數、source／binary／
log hashes；`regressions/` 保存初次失敗，`default-solver-baseline/` 保存修改前
重現，`regressions-configured/` 保存正確設定的四項回歸。
`acceptance.json` 彙整最終重新核對的來源、建置、原生結果與 CLI 證據。

這批測試是正確性證據，與編譯並行執行的時間不作效能統計。
不涵蓋程序遺失、OOM killer 或無法返回的 MPI／PETSc collective。
其餘 runtime／adapter／executor 邊界與外部資產身分仍須依清單接續完成。
