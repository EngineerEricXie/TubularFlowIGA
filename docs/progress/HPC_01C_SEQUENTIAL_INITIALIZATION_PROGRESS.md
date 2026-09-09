# HPC-01C：Sequential 初始化與輸入協調

日期：2026-09-08（本機 EDT）。接續 [sequential 關檔驗收](HPC_01C_SEQUENTIAL_OUTPUT_PROGRESS.md)。
HPC-01C 與整份 goal 仍未完成。

## 問題與修正

Sequential runner 原先在本地解析、建立 1D runtime、配置 port 與初始狀態時，
可能直接從單一 rank 跳到外層 catch，其他 ranks 卻進入下一個 3D collective。
Graph preflight 也在協調前配置錯誤字串，且沒有檢查實際使用的輸入內容一致。
原有 1D runtime 支援 local trial outcome callback，但 sequential caller 尚未提供。

本批將主程式拆出 `RunSequentialFlow(argc, argv, communicator)`，由 CLI 保持
PETSc 初始化與結束。新 API 借用 communicator，允許同一 PETSc session 中
使用獨立子群；所有 runtime 均在返回前釋放。CLI 與原有數值參數維持不變。

新增協調範圍：

- 解析參數並先確認 graph／positional 模式一致；graph 路徑、manifest 內容、
  definition 與 domain assets 在共同階段完成。
- 比較 stop step、Newton／strong／Aitken 控制及原有 precommit 注入設定；
  比較有效 PETSc options，排除已另外驗證的 application arguments 與輸入位置。
- 對配置、資料庫、mesh、velocity、1D geometry，以及實際使用的 inlet／3D
  boundary periodic tables 建立邏輯名稱 catalog，再比對內容 fingerprint。
  相同內容可位於不同的 rank-local 路徑。輸入必須在執行期間保持不變。
- 本地配置解析、拓樸檢查、1D owners 與 3D host inputs 建構、scalar preflight、
  port identifiers／graph／validation、初始 1D 狀態與 3D boundary input、
  measurement port 配置、lagged pressure 及 history 容量準備。
- 透過既有 collective factory 配置 3D runtime storage；1D runtime 接上
  `FailureAgreement`，使其 trial／substep 本地結果在進入下一個原生操作前協調。

`SequentialLocalValue` 僅包本地工作，結果保留到共同 outcome 後，再以
noexcept move 返回。這些 callback 不建立或析構分散式 PETSc owners；
3D runtime constructor、reference-flow measurement、InitializeState 與
MeasurePorts 等 collective 操作維持在 callback 外。
沒有改動物理公式、容許值、CSV／JSON 格式或 METIS／`.ntiga` 介面。

## 驗收

Ignored 證據位於 `outputs/hpc01/sequential-initialization/`。

| 驗收 | 結果與範圍 |
|---|---|
| 原生階段故障／重試 | positional 19 個 stage、graph 21 個 stage，在 world 3 ranks 及獨立 1／2-rank 子群各執行一次，共 120 次故障與 120 次健康重試 |
| 有效 PETSc 選項差異 | 初始化後僅在最後一個 rank 修改 KSP、prefixed PC 或 unused option；兩種入口、2／3 ranks 共 12 次共同拒絕及 12 次恢復選項後重試 |
| 清理驗收 | 後段故障與每次 stage 健康重試要求觀察到非空 PETSc object 集合；外部參照釋放前逐一確認 runtime ownership 已解除。各 communicator 的 positional／graph 累計觀察 348／372 個 references |
| 真實 CLI 輸入矩陣 | 32 個案例、60 份 rank reports；18 次輸入／控制錯誤共同退出，14 次健康或初始化正規化案例成功 |
| CLI 相容性 | 10 次新版健康執行的 CSV 與對應舊版基準相同，包含本地副本及未使用波形檔案 |
| 原生多 rank 相容性 | positional explicit 與 graph explicit／fixed／Aitken，1／2／3 ranks 的 18 份 CSV 與修改前逐位元組相同 |
| 原有 smoke | 全套 sequential smoke 通過，涵蓋 positional／graph、subcycling、fixed／Aitken 等價與原有拒絕案例 |

Stage 故障包括 caught `bad_alloc`，以及一次非標準例外；失敗後所有 ranks
返回 1 且不建立 output 目錄，再使用同一 communicator 重建健康案例。
前段在 PETSc 建立前失敗，容許觀察到零個物件；後段與健康重試不容許空集合。
這是指定階段的故障協議驗證，不是逐個 allocator 呼叫點的窮盡測試。

CLI 矩陣包含 graph／三個 domain 配置、database、mesh、velocity、兩個
1D geometry 與兩種 selected table 的內容差異；selected table 的缺檔、
directory、FIFO；非法 argument、input mode、stop step、Newton 控制差異。
內容相同的本地副本可正常執行；未使用的 malformed table 不會被額外解析。
所有預期失敗都必須有限時間返回；timeout 不能算通過。

初次 CLI controller 將 rank 間 `PETSC_OPTIONS` 環境差異預期為拒絕，但這套
PETSc 在初始化時將相關條目統一，實際執行成功。保留初次 `input.log` 與
`input-cases/`，最終矩陣將此記錄為初始化行為，並增加初始化後的真正 options
database 差異測試。沒有把未生效的環境注入當成有效選項拒絕的證據。

## 重現方式

環境：本機 WSL、GCC 11.4、Open MPI 4.1.2、PETSc 3.15.5 real64／Int32。
OMP／BLAS 各 1；原生比較使用 `preonly`、LU／MUMPS、`ksp_rtol=1e-12`。
原有 smoke 沿用自身的每案例設定。Production 及 test build 無新增 compiler warning。

```bash
make -C solvers/coupling iga_1d_3d_explicit sequential_initialization_failure_test \
  PETSC_DIR=/usr/lib/petscdir/petsc3.15/x86_64-linux-gnu-real
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 \
  PETSC_OPTIONS='-ksp_type preonly -pc_type lu -pc_factor_mat_solver_type mumps -ksp_rtol 1e-12' \
  timeout --kill-after=5s 240s mpiexec --map-by core --bind-to core -np 3 \
  solvers/coupling/sequential_initialization_failure_test /absolute/path/to/new-output
python3 scripts/hpc_sequential_input_regression.py \
  --case-dir /path/to/retained-explicit-smoke-fixture \
  --reference-binary /path/to/before/iga_1d_3d_explicit \
  --output-dir /path/to/new-input-results
```

CLI fixture 的來源與保留方式見上一批關檔報告。`native-final.log`、
`input-final/summary.json`、`comparison.json` 及 `smoke.log` 保存最終驗收；
`source-final.json`／tar、`after-binaries/` 與 `acceptance.json` 保存對應版本。

## 剩餘工作

Sequential strong-fixed／Aitken 迴圈本身仍需補上 Begin／input／port result、
本地 residual／relaxation 結果、全群收斂、precommit／prepare 與 postcommit
bookkeeping 的協調。Abort 必須嘗試所有 runtime 後才組合診斷，避免單一
1D abort 錯誤讓 peers 進入不同的 3D collective。
本批健康 strong 回歸不代表上述故障邊界已完成。

接著仍需完成支援入口覆蓋稽核、HPC-01D 及後續各階段。清單保持 9／38。
沒有新增程序失聯恢復、committed-step 原地回復、GPU、多節點或擴展性驗收，
新增初始化同步與讀取的成本尚未量測。
