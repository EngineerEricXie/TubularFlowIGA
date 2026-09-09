# HPC-01C：1D checkpoint、VTK 與資源摘要串流

日期：2026-09-08。狀態：本批指定 helper 修正與驗收完成；HPC-01C 尚未完成。
本機證據位於 `outputs/hpc01/text-helpers/`，其中包含修改前來源與執行檔、
執行命令、退出狀態、比較結果、故障變體及修改後來源快照。

## 問題與修正

`ostringstream` 的格式化寫入可能把記憶體配置例外轉成串流錯誤狀態。
若未啟用例外或檢查狀態，接著呼叫 `str()` 可能取得部分文字，讓呼叫端
誤認為 metadata、fingerprint、檔名或摘要已成功產生。

- `OneDNetworkFingerprint` 與 `SerializeOneDCheckpointMetadata` 啟用
  `badbit | failbit` 例外。保留原始 hash seed、字元順序與 17 位精度。
- `ReadOneDCheckpoint` 改用共用 `ReadCheckedText`，完整讀取成功後才解析
  metadata。原先的 `contents << input.rdbuf()` 只檢查來源的 `badbit`，
  無法涵蓋目的串流失敗。本次沒有改動 PETSc state 格式或讀寫協議。
- `VtuStepPath` 與 `hdf_detail::ArraySchema` 啟用相同串流例外。
- 資源摘要的純文字格式化移至 `execution_resources_detail::FormatReport`，
  啟用串流例外，讓它能獨立接受配置失敗測試。實際報告仍在原有共同階段
  由 communicator rank 0 寫出，並保留目的 ostream 的狀態檢查。

這些錯誤會交由既有 MPI 邊界協調：

| 路徑 | 本地工作所在的共同階段 |
|---|---|
| 1D fingerprint 建立／驗證 | `1d checkpoint preparation`／`1d restart validation` |
| 1D metadata 序列化／讀取 | `1d checkpoint write preparation`／`1d checkpoint metadata read` |
| CPU VTU 檔名 | `flow output preparation`／`transport output preparation` |
| CPU VTKHDF array schema | `flow field output`／`transport field output` 的 root callback |
| 資源摘要 | `execution resource report` |
| CUDA 視覺化 | 既有單程序錯誤處理 |

callback 內仍只有本地操作；沒有以外層共同階段包住 MPI collective。
失敗時可能留下先前寫出的部分檔案，原子發布與完整恢復協議仍由 HPC-05 追蹤。

## 驗收證據

| 驗收 | 結果與證據 |
|---|---|
| 五個真實 helper 配置失敗掃描 | 每組依成功呼叫量測的配置次數逐一注入；fingerprint 7、metadata 13、schema 8、檔名 9、資源摘要 2 次。world 3 ranks 與 split 1／2 ranks 合計 117 次故障、117 次健康重試通過；`helper-faults.log` |
| 反向驗證測試有效性 | 五個隔離 header 變體各移除一處例外 mask，全部以退出碼 1 回報 `helper swallowed allocation failure`；`mutant-results.json` |
| 原生資源與後端檢查 | world／split 的 17＋16＋17 項檢查通過，含 root 壞輸出串流、設定與 backend 拒絕；`resources.log` |
| 原生 1D PETSc | 各隱式 formulation、metadata roundtrip、分組 checkpoint 讀寫與狀態檢查通過；`one-d-native.log` |
| 1D checkpoint CLI 矩陣 | 25 個案例、75 份 rank reports 通過，含 rigid／explicit／species 續跑、舊分散式 VecView 檔案相容及 metadata／state／輸出錯誤；`checkpoint-matrix/summary.json` |
| 修改前後 CPU／1D 比較 | 12 次 CPU flow／transport 模擬及 6 次 1D 模擬，190 份輸出逐位元相同；兩組 transport VTKHDF 各六個資料集比較通過。原始 argv 與檔案紀錄在 `comparison-runs.json`、`comparisons.json` |
| VTU／VTKHDF／完整文字單元測試 | 皆退出 0；`vtk-unit.log`、`vtkhdf-unit.log`、`checked-text-unit.log` |
| CUDA 建置與原生比較 | CUDA 12.6、sm89 建置成功。修改前後兩次傳輸模擬的 14 份輸出逐位元相同；四份場的相對 L2 為 0，與本輪 CPU 最終場相對 L2 為 `2.397514096418484e-6`，通過既有 `1e-5` 門檻；`cuda-final-comparison.json`、`cuda-output-bytes.json` |

helper 故障測試先建好輸入，再在 communicator 最後一個 rank 的本地呼叫中
注入一次 `bad_alloc`；所有 rank 都必須觀察到共同失敗。例外診斷與 MPI
agreement 時關閉注入，重試結果必須與健康基準完全一致。這是 helper 層的
逐配置測試；原生呼叫端則另以資源／checkpoint／視覺化測試驗證，未宣稱
逐一在所有 CLI 的每個配置位置注入故障。

新舊輸出比較排除含時間量測的 `summary.json`，保留 CSV、VTU／VTP／PVD、
checkpoint JSON／PETSc state 與場資料等其餘輸出。HDF 比較器檢查六個
資料集的形狀、有限性與數值，沿用相對 `1e-6`／零參考絕對 `1e-12` 門檻。
這些比較不宣稱涵蓋所有 HDF metadata。CUDA 門檻讀取既有
`benchmarks/hpc_baselines.json`，沒有改動標準。

## 重現

本機 GCC 11.4、Open MPI 4.1.2、PETSc 3.15.5 real64／Int32／MUMPS；
CUDA 12.6、RTX 4080 SUPER。OMP／BLAS 固定為 1。CPU transport 使用
GMRES＋LU/MUMPS，flow 使用 preonly＋LU/MUMPS，`ksp_rtol=1e-12`；
1D 使用其既有 preonly＋LU/MUMPS 測試設定。

```bash
make -C solvers/cpu text-helper-failure-test \
  PETSC_DIR=/usr/lib/petscdir/petsc3.15/x86_64-linux-gnu-real
make -C solvers/one_d iga_1d one_d_petsc_test one_d_checkpoint_format_test
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 \
  PETSC_OPTIONS='-ksp_type preonly -pc_type lu -pc_factor_mat_solver_type mumps' \
  timeout --kill-after=5s 90s mpiexec -np 3 solvers/one_d/one_d_petsc_test
python3 scripts/hpc_one_d_checkpoint_regression.py --output-dir /path/to/new-output
conda run -n tubularflow-cuda make cuda CUDA_ARCHS=89
```

執行於本機計算資源，沒有跨節點或排程器驗收，也不作效能提升宣告。
CPU、1D 與 CUDA 的本批建置皆退出 0，沒有 compiler warning。

## 剩餘工作

- 本批結束時尚未檢查的 `OneDOutputWriter::Write` 檔名串流，已在後續
  [原生檔名與 runtime 簽章驗收](HPC_01C_STREAM_BOUNDARY_PROGRESS.md) 補強。
  後續結果不屬於本報告原始 source archive。
- 對支援入口完成完整錯誤邊界清單，逐項核對本地工作、collective 呼叫順序、
  成功摘要與清理；本次指定 helper 通過不代表整份入口稽核完成。
- HPC-01D 的工具／CUDA 能力矩陣及後續 HPC-03～09 保持原範圍。
  整體清單仍為 9／38 項勾選，整份 goal 保持進行中。
