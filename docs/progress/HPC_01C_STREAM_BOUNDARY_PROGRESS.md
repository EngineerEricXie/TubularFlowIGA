# HPC-01C：原生檔名與 runtime 簽章錯誤邊界

日期：2026-09-08。本批實作與列明驗收完成；HPC-01C 尚未完成。
證據根目錄：`outputs/hpc01/stream-boundaries/`。

## 改動與目的

中間 `ostringstream` 若把格式化例外轉成 `badbit`，呼叫端仍可能取得部分
字串。本輪修正九個串流，啟用 `badbit | failbit` 例外：

- `OneDOutputWriter::Write` 的 VTP 檔名。
- `TransientFlowRuntime` 的初始化簽章、BeginStep 控制簽章、port 量測
  數值簽章，以及線性求解失敗診斷所用 `PreciseNumber`。
- `TransientTransportRuntime` 的 GatherState 與 IntegrateFields 簽章。
- `ThreeDBodyFittedFlowTransportDomainAdapter` 的 PlanSignature／InputSignature。

MPI 簽章仍在既有本地準備階段中建立；只有完整建立成功，才進入後續
agreement／gather／reduction。沒有改動物理模型、格式、數值精度、原有
signed flow 語義或 communicator 所有權。

新增 [MPI 錯誤邊界索引](../architecture/MPI_FAILURE_BOUNDARIES.md)，列出支援
入口、runtime／adapter 邊界、驗收索引與 F01～F06 剩餘工作。它是接續稽核
的工作表，不是所有入口均已完成的宣告。

## 真實串流故障

測試專用 `StringStreamFailure.hpp` 安裝 `num_put` locale facet，只在真正
`ostringstream` 的數字格式化中拋錯。涵蓋 `runtime_error`、`bad_alloc` 與
非標準例外；不更動 MPI／PETSc allocation 或既有 file stream。測試固定
OMP／BLAS 為 1，於主執行緒安裝及還原 locale。

### 1D 原生 CLI

`one_d_filename_failure_test` 直接呼叫重新命名的真實 CLI main，保留其 MPI
初始化、輸出階段與退出碼。只對 `profile_1d_` 前綴的指定 step 注入。

- 修改後：1／3 ranks，各測 step 0 與 step 1 的三種例外。12 個故障作業
  全部退出 1，分別在 `1d initial output`／`1d step output` 回報 rank 0。
  沒有 `summary.json`，沒有名為 `profile_1d_` 的截斷檔名。
- 每個故障後用新作業及新目錄重試，共 12 次重試通過。144 份非時間摘要
  輸出與健康基準逐位元相同。加上兩個健康基準，共 26 個案例、52 份
  rank reports；見 `filename-final/summary.json`。
- 修改前反例：隔離載入修正前 `OneDOutput.hpp`、保留同一原生 CLI 與測試
  facet。12 個注入作業皆退出 0，寫出名為 `profile_1d_` 的檔案，並在
  summary 宣告 converged。加兩個健康基準，共 14 個案例、28 份 rank reports；
  見 `filename-old/summary.json`。這些退出 0 是明確記錄的錯誤行為。

重試驗證不宣稱可以重用已部分寫出的同一 writer；故障可能已寫出部分 CSV，
原子發布與崩潰恢復仍由 HPC-05 追蹤。

### MPI runtime／adapter

每個入口對 communicator 最後一個 rank 注入三種例外。world 3 ranks 與
split 1／2 ranks 均執行，全部 rank 必須收到相同的準備階段診斷，不能僅以
後續 agreement 拒絕截斷簽章代替。健康呼叫與既有狀態／串行比較接續驗證。

| 真實呼叫 | 預期共同階段 | 本輪新增故障數 |
|---|---|---:|
| Flow InitializeState | `flow initialization preparation` | 9 |
| Flow BeginStep | `flow begin step preparation` | 9 |
| Flow MeasurePorts | `flow port measurement preparation` | 9 |
| Transport GatherState | `transport global state preparation` | 9 |
| Transport TotalMass | `transport mass preparation` | 9 |
| Staged adapter BeginStep | `3d staged begin preparation` | 9 |
| Staged adapter SetPortInput | `3d staged input` | 9 |
| Staged adapter SetTransportConcentration | `3d staged concentration` | 9 |

合計 72 次新增故障，連同原有案例，trial／port／staged 測試分別記錄
111／73／149 項檢查，共 333 項。最終 `trial-final.log`、`ports-final.log`、
`staged-final.log` 全部退出 0，包含原有解析量、正值／有限性、狀態回復、
物種守恆、流向反轉及 serial/group 數值 gate。

`IntegrateFields` 的 source preparation 使用同一受檢查串流；本輪未另外在
每一種 source term 格式化位置注入。`PreciseNumber` 是原本就會拋錯的 KSP
失敗分支之診斷補強，不算進上表八個原生注入入口。

## 找到並修正的既有測試問題

第一次 port／staged 原生回歸失敗，進一步診斷指出 `flow reference result`
與 `transport mass integration` 的預期故障沒有發生。原因是 runtime 建構
已改為複製輸入，但這兩個舊測試仍保留 caller 資料的指標，誤以為它指向
runtime 內部 storage。它們改動的數值沒有進入被測演算法。

本輪在測試宏下增加 `ReferenceVelocityForTesting`／`VolumeRuleForTesting`，
讓既有測試直接改動 runtime 擁有的資料，再還原原值。production build
不提供這兩個介面。仍使用原先的 NaN reference velocity 與負積分權重，
沒有取消案例、改變預期錯誤階段或放寬數值門檻。

早期 `ports.log`／`staged.log`、加入明確階段後的 `*-diagnostic.log` 均保留
退出 1；最後以 `*-final.log` 為通過證據。這也說明不能直接把早期報告的
測試計數視為目前工作樹的有效覆蓋。

1D controller 初次錯把固定三 rank 的 measurement wrapper 用於單 rank，
在啟動求解器前即失敗。已讓單 rank 使用既有 `hpc_rank_run.py`，三 rank
沿用帶 report rendezvous 的 wrapper；失敗目錄 `filename-matrix/` 保留，
修正後在新目錄 `filename-final/` 驗收。

## 相容性與整合

- 本輪 CPU／1D 修改前後共 18 次模擬，190 份輸出逐位元相同；兩組 transport
  VTKHDF 的六個資料集比較通過。見 `comparison-runs.json`／`comparisons.json`。
  比較排除含量測時間的 summary；HDF 使用既有相對 `1e-6`／零參考絕對
  `1e-12`，不宣稱比較所有 metadata。
- 0D／flow／species graph 在 1／2／3 ranks 共 18 次修改前後模擬，57 份
  CSV 逐位元相同；`graph-runs.json`／`graph-comparisons.json`。
- 完整 `make -C solvers/coupling petsc-test` 通過，包括 sequential
  fixed／Aitken 的物理相容性與原有守恆 gate；`sequential-smoke.log`。
- CPU、1D、coupling 最終建置通過，無 compiler warning。1D `core-test` 通過。
  本輪未更動 CUDA 使用的求解或輸出 helper，未重跑 GPU 或完整 FSI exporter。

修改前 binary 是開始本批時工作區中的執行檔，個別 hash 保存於證據清單。
`before-source.tar.gz` 是上一批文字 helper 的來源快照；未在本批從該 archive
重建全部 baseline binary，不能只用該 archive 推定每個舊 binary 的來源身分。
1D 舊 writer 反例的 header 與實際編譯命令則另存於 `old-filename/`。

## 重現命令

本機 GCC 11.4、Open MPI 4.1.2、PETSc 3.15.5 real64／Int32／MUMPS。
runtime 測試使用 GMRES＋LU/MUMPS、`ksp_rtol=1e-12`，OMP／BLAS 為 1。
這些是小型正確性測試，不是效能量測或跨節點驗收。

```bash
make -C solvers/one_d one_d_filename_failure_test iga_1d core-test
python3 scripts/hpc_one_d_filename_regression.py --output-dir /path/to/new-output
make -C solvers/coupling three_d_trial_failure_test three_d_port_failure_test \
  three_d_staged_failure_test
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 \
  PETSC_OPTIONS='-ksp_type gmres -pc_type lu -pc_factor_mat_solver_type mumps -ksp_rtol 1e-12' \
  timeout --kill-after=5s 180s mpiexec -np 3 \
  solvers/coupling/three_d_staged_failure_test /path/to/new-staged-output
```

其餘兩個 runtime 測試使用相同環境與各自新的 output directory。完整 argv、
退出碼與原始 logs 在證據根目錄的同名 JSON／log。

後續依 [F01～F06](../architecture/MPI_FAILURE_BOUNDARIES.md#可執行的剩餘項目)
接續；主清單仍為 9／38 項完成，整份 goal 保持進行中。
