# HPC-04A：1D implicit／SNES 的獨立求解器選項

日期：2026-09-09。基準 commit `672cfc87594528524d8770c8604db00643c5932f` 加本批修改。
**部分完成**：四種 1D implicit 方法、正式 `iga_1d` CLI 及 native graph 已接入；
body-fitted standalone、immersed／moving／FSI、完整 nested diagnostics 與後端矩陣仍待完成。
整份 TODO 維持 38 項範圍、14 項已完成。操作見 [SOLVER_OPTIONS.md](../SOLVER_OPTIONS.md)。

## 實作

`OneDPetscSolverContext` 每個 runtime 持有一份不可變的 `PetscSolverOptions`。native
runner 在本地 1D 建構階段之外先完成 context 的 collective agreement，再由原 callback
借用 context；CLI 也在 runtime 建構前建立 owner。原自訂 advance callback 簽名不變。
省略 context 的直接 Solve／Advance API 在該次 advance 建立空 prefix snapshot，保留
既有未加前綴呼叫方式。context 及 communicator 必須活到呼叫結束，且使用相同群組。

`pressure_network`、`linearized_aq` 的 KSP，以及 `nonlinear_aq`、`implicit_1d_pde` 的
SNES／KSP 都掛入 private options。非線性 initial linearized guess 共用 runtime prefix；
SNES 的線性求解器也接受同一 `..._ksp_*` override。全域選項仍作共同基線，顯式前綴
優先；多 rank 非線性預設 preonly／LU／MUMPS 保留，所有線性路徑也檢查實際 factor
backend 能力。每次 SetFromOptions／Solve 均使用上一批 returning-handler scope，
callback 原始 C++ 錯誤可傳回協調邊界，不跨 PETSc C callback 拋出。

正式入口在 accepted-step 輸出邊界寫 `one_d_solver_configuration`：prefix、step、
有效 KSP／PC／backend、rtol／atol、最後線性 iterations／reason；有 SNES 時另記其
type、iterations／reason。這不是 substeps／trials 的迭代總和。rigid／explicit 1D
沒有 PETSc solve，不產生這筆診斷。

CLI 原本只略過 PETSc option key，會把 `fgmres` 等 value 當成非法位置參數；本批
同步略過其 optional value，PETSc 在初始化時完成真正解析。`PETSC_OPTIONS` 環境
變數仍可使用。配置、場、checkpoint 與 `.ntiga` 格式不變；native source identity
更新為 `1cc091f78e85c779345c5553eb4e056012dfbe108e35bf84bd79e95eaa5ce994`。

## 驗收環境與門檻

本機 `TsungYehLab`，GCC 11.4、Open MPI 4.1.2、PETSc 3.15.5 real64／int32、MUMPS；
OpenMP／OpenBLAS 各 1 thread，允許 CPU 0–15、未固定逐核 binding。沒有 cluster
allocation／跨節點／CUDA 實測。本批屬配置與相容性驗收，不宣稱效能提升。

CLI fixture 由 `examples/one_d/compliant_bifurcation` 複製，改為 implicit_petsc、4 steps、
每步輸出，其餘模型、波形與 RCR 保留。四種 formulation 分別比較：封存舊 CLI、
新 CLI 共同 preonly／LU／MUMPS、前綴 FGMRES（非線性另選 Newton trust region）。
預設輸出 bytes 精確一致；summary 只排除 setup／solve／output 秒數及 RSS。override
逐 CSV column 的 L2 相對差須不超過 `1e-6`（零 reference 用 absolute `1e-12`），
所有數值有限；VTK topology／attributes 精確、numeric DataArray 用 `1e-12+1e-6|reference|`；
其他 metadata bytes 精確。summary 物理量用相同 absolute／relative 門檻、結構及
convergence 狀態一致。這些門檻在比較前固定。

真正 1D–3D graph 由既有 flow fixture 產生，各 1D domain 改 compliant／implicit，
linear wall 的 Young modulus `1e9`、thickness ratio `0.1`。四種 formulation 各跑
reference／default／override；source 使用 FGMRES、branch_a 使用 GMRES，非線性
source 使用 Newton trust region。原 native Newton、coupling 及 mass gates 不變。
另從第 3 步的新 MPI 作業恢復到第 6 步，完整 history 與 field payload 必須精確一致。

## 結果與證據

原始證據位於 `outputs/hpc04/one-d/`；case bytes、argv、binary hashes、逐 rank
stdout／stderr／RSS 及機器可讀結果均保留。主要驗收：

| 範圍 | 結果 |
|---|---|
| `cli-suite-v2`，3 ranks | 22 作業：12 正向、10 預期 exit 1；四種方法預設輸出精確一致；override 最大 CSV 相對 L2 `8.061e-10` |
| `cli-rank1` | 18 作業：12 正向、6 預期 exit 1；單 rank 不把 PETSc 自帶 LU 當成不支援；最大 CSV 相對 L2 `8.061e-10` |
| `cli-ranks4` | 22 作業：12 正向、10 預期 exit 1；最大 CSV 相對 L2 `8.061e-10` |
| `native-implicit-final` reference／default／override | 12 作業通過，四種 formulation 的預設場及 history bytes 一致；override 全場相對 L2 最大 `3.160e-14`，實際 KSP／SNES 與設定相符 |
| `native-implicit-final` save／resume | 8 個新作業通過；四種方法 save-3／resume-6 的輸出與場分片精確一致 |
| `native-checkpoint` | 63 作業通過：41 正向、22 預期失敗；既有完整 graph checkpoint、signal、corruption／identity 拒絕與故障恢復不變 |
| `body-prefix-regression` | 前一批 14 作業通過：12 正向、2 預期失敗；貼體 flow／transport 的配置隔離與原場門檻保留 |
| `unit-mpi` | world 3 及 split 1／2 ranks；四種方法的兩組 KSP／SNES override、原場門檻、未知 KSP／SNES、多 rank LU backend 拒絕及健康重試通過 |
| `core-mpi`、`failure-mpi`、`groups-v2-flow`／`groups-v2-species` | 最終建置的 1D core／checkpoint、單 rank 故障與健康重試、native graph 獨立 1／2-rank 子群組回歸通過 |
| `sanitized-mpi` | 同一套隔離與錯誤重試測試以 ASan／UBSan 通過，無 sanitizer 診斷 |

最終集合共 165 個 MPI 作業（115 正向、50 預期失敗），正向 333 份 rank reports。
`acceptance.json` 保存來源／證據 SHA256 與作業索引；`audit_final.py` 可重新核對。
所有最終正向 stderr 為空、無 timeout；編譯沒有新增 warning。

場／continuity 比較不代表這個短波形 fixture 已通過新的物理精度驗收。舊與新程式的
`maximum_sampled_relative_continuity_residual` 都約為：pressure network `0.00353`、
linearized AQ `0.00355`、nonlinear AQ `0.00355`、implicit PDE `0.12718`。
此診斷由 sampled inlet−outlet−storage 差計算；本批保留這些值並比較差異，沒有提高
允許值或宣稱它們低於另一個尚未設定的 conservation gate。需要更小 dt／網格收斂的
物理驗收仍屬後續工作。

3 ranks CLI override 的成本如下，來自 summary 及最高 rank RSS；既有 1D CLI 的
`solve_seconds` 包含組裝與求解，尚無兩者各自的 phase，因此不能將此列當作純 KSP
時間。通信時間亦未獨立量測。此限制保留給 HPC-04C；本批不作 performance 結論。

| 方法 | setup s | solve（含組裝）s | output s | peak rank RSS MiB |
|---|---:|---:|---:|---:|
| pressure network | 0.000605 | 0.004501 | 0.000639 | 38.41 |
| linearized AQ | 0.000631 | 0.005439 | 0.000420 | 39.10 |
| nonlinear AQ | 0.001172 | 0.016357 | 0.000494 | 40.30 |
| implicit PDE | 0.000613 | 0.012110 | 0.000586 | 40.18 |

## 重現與開發中失敗

```bash
export PETSC_DIR=/usr/lib/petscdir/petsc3.15/x86_64-linux-gnu-real
make -C solvers/one_d -j2 iga_1d one_d_solver_options_test one_d_implicit_failure_test PETSC_DIR="$PETSC_DIR"
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 mpiexec -np 3 \
  solvers/one_d/one_d_solver_options_test examples/one_d/rigid_straight/skeleton_initial.swc
python3 scripts/hpc_one_d_solver_prefixes.py --binary solvers/one_d/iga_1d \
  --reference outputs/hpc04/one-d/iga_1d-reference --output-root NEW_ROOT --ranks 3
```

參考 binary 在修改前封存，SHA256 隨 acceptance 保存；重現跨版本 byte 比較需保留
此 binary。native fixture 生成、配置、比較及續跑由
`scripts/hpc_one_d_graph_prefixes.py --binary outputs/hpc04/one-d/native-v1
--reference outputs/hpc04/options/native-handler --output-root NEW_GRAPH_ROOT --ranks 3`
重現。
ASan／UBSan 編譯完整 argv 在 `sanitized-command.json`；使用
`-O1 -g -fsanitize=address,undefined -fno-omit-frame-pointer`，
`ASAN_OPTIONS=detect_leaks=0:halt_on_error=1`、`UBSAN_OPTIONS=halt_on_error=1:print_stacktrace=1`。
關閉第三方 leak scanning，故不宣稱完整第三方 memory leak audit。

早期 CLI probe 揭露 option value 的 parser 問題，已修正；最初驗收腳本使用錯誤
formulation 拼字 `implicit_pde`，改回設定檔原本的 `implicit_1d_pde`，未改 parser
規格或放寬 gate。既有 failure test 對未知 KSP 的階段名稱從 `KSPSetFromOptions`
更新為精確的 `1d solver options`；狀態未發布、metadata 未變、健康重試等斷言保留。
子群組 flow 測試首次多傳一個不支援的參數而顯示 usage，依原介面重跑。
早期失敗日誌不列入最終通過數。

下一步為其餘 body-fitted standalone、immersed／moving／FSI 的 prefix 整合，以及
完整 nested diagnostics／後端矩陣；HPC-04B/C 還需要候選預條件器與中／大案例成本及
數值驗收。整份目標仍持續，未以本批小案例取代後續必要驗收。
