# 1D 與貼體 graph 的 PETSc 求解器選項

`iga_multidomain_flow`／`iga_1d_3d_bifurcation` 中，每個貼體 3D domain 的 flow 與
transport，以及每個 implicit 1D domain，現在都有獨立 options prefix。
全域 PETSc 選項仍作共同基線，domain 選項優先：

```bash
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
export PETSC_OPTIONS='-ksp_type preonly -pc_type lu -pc_factor_mat_solver_type mumps'
mpiexec -np 3 solvers/coupling/iga_multidomain_flow \
  --graph-case CASE --output-dir NEW_OUTPUT \
  -domain_island_a_flow_ksp_type fgmres \
  -domain_island_b_flow_ksp_type gmres
```

範例假定 case 中兩個貼體 domain 分別叫 `island_a`、`island_b`；保留同一 LU／MUMPS
預條件器，但分別選擇 FGMRES／GMRES。這是配置方式，不是效能推薦。新增配置仍須
通過原有 Newton、守恆與場誤差驗收。

| 目標 | Prefix 範例 |
|---|---|
| `junction` 的流場 | `domain_junction_flow_` |
| `junction` 的物種傳輸 | `domain_junction_transport_` |
| `island_a` 的速度 fieldsplit 子 KSP | `domain_island_a_flow_fieldsplit_0_` |
| `island_a` 的壓力 fieldsplit 子 KSP | `domain_island_a_flow_fieldsplit_1_` |
| `junction` flow 的 block-Jacobi 子 KSP | `domain_junction_flow_sub_` |

例如 `-domain_junction_transport_ksp_rtol 1e-10` 僅改該 transport 的 relative tolerance。
`-domain_island_a_flow_fieldsplit_0_pc_type gamg` 只有實際 PC 選用 fieldsplit 才生效；
指定前綴不會自行改 PC 類型。PETSc 的 prefix 規則見
[KSPSetOptionsPrefix](https://petsc.org/release/manualpages/KSP/KSPSetOptionsPrefix/)。

1–64 個小寫 ASCII 字母、數字、底線組成的 domain ID 使用表中的可讀形式。
其他 ID 用 `domainh_<完整 domain ID 的 SHA256>_<role>_`，避免名稱碰撞及超長 prefix。
每個 accepted step 的 stdout 會列出 `solver_configuration` JSON，記錄實際 domain、
role、prefix、KSP／PC、factor backend、rtol／atol、最大迭代數與最後一次線性求解的
iterations／reason。`last_iterations` 是最後一次 KSP 的值，不是整個 Newton 或所有
trial 的累計；warm start 已收斂時可能為 0。尚未選定的 factor backend 記為 `auto`。

## 優先順序與相容性

1. Runtime 的既有預設：小型貼體 flow 為 FGMRES／block-Jacobi；達既有節點門檻時
   使用原 Schur fieldsplit／GAMG 設定；transport 維持 GMRES／block-Jacobi。
2. 既有未加 domain prefix 的 PETSc 選項，例如 `-ksp_type`、`-pc_type`、
   `-fieldsplit_0_pc_type`，覆蓋相應預設。
3. 該 domain／role 的完整 prefix 選項再覆蓋共同基線。

每個 runtime 建立精確的 options snapshot；空值、flag、含空白或 option-like 文字的
value 都直接複製，不重新解析 `PetscOptionsGetAll` 的顯示字串。不同 runtime 的 snapshot
與 options prefix 先在其 communicator 內取得一致。來源值不會被某個 domain 的 override
改寫；使用旗標會傳回來源，使 `-options_left` 能辨認實際被子求解器使用的選項。

PETSc 3.15 實測中，某些子 KSP 沿用 prefix，卻查詢目前預設的 options database。
因此 SetFromOptions／SetUp／Solve 使用暫時的 private database scope，返回後立即
還原原 stack，再協調錯誤。此 scope 也暫設 returning error handler，避免無效 KSP
類型觸發預設 handler 中止其他 ranks；返回、例外與失敗都還原 caller handler。
fieldsplits 必須在 prefix 設定後建立。此流程與現有 PETSc
呼叫一樣，在 MPI 初始化執行緒同步執行；不提供同程序內多個 host threads 同時呼叫
PETSc 的支援。分屬不同 MPI 程序的獨立 communicator 已驗證。

Embedding caller 可在 `TransientFlowRuntime`／`TransientTransportRuntime` 建構參數
最後提供 `solver_options_prefix`。預設空字串保留舊的未加前綴介面，包括舊 flow
constructor 的全域 fieldsplit defaults 行為；新的 scoped graph 使用 private defaults。
`SolverConfiguration()` 是 runtime 開啟時的本地查詢；外層須自行提供群組錯誤協調。

## 1D implicit 與 SNES

Graph 的 1D domain 使用相同 `domain_<id>_flow_` 規則，例如
`-domain_source_flow_ksp_type fgmres`。獨立 `iga_1d` CLI 則以選定的 equation-system
名稱產生 prefix；系統 `blood_flow_1d` 使用 `domain_blood_flow_1d_flow_`：

```bash
mpiexec -np 3 solvers/one_d/iga_1d IMPLICIT_CASE --output-dir NEW_OUTPUT \
  -domain_blood_flow_1d_flow_ksp_type fgmres \
  -domain_blood_flow_1d_flow_snes_type newtontr
```

這個範例需要 compliant／implicit_petsc case；只有 `nonlinear_aq`、`implicit_1d_pde`
會使用 SNES。`pressure_network`、`linearized_aq` 使用線性 KSP。非線性初始猜測的
linearized solve 與 SNES 的 KSP 共用該 runtime 的 prefix；例如 `..._ksp_type` 同時
影響二者，`..._snes_type` 只影響 SNES。原非線性預設仍是 preonly／LU，multi-rank
預設 MUMPS；explicit override 仍須通過實際 factor backend 能力檢查。

`one_d_solver_configuration` JSON 記錄 prefix、step、有效 KSP／PC／backend、tolerances、
最後一次線性 iterations／reason；非線性另記 SNES type、iterations／reason。數值 gate
通過後才由 CLI／graph 的 accepted-step 輸出邊界寫出；不是所有 substeps 或 trials 的
累計。rigid／explicit 1D 沒有 PETSc solve，因此不產生這筆診斷。

Embedding caller 可建立 `OneDPetscSolverContext`，在 `AdvanceImplicitOneD` 或個別
Solve API 最後傳入指標。context 須使用相同 communicator，並活到整次呼叫完成；
context 建構包含 collective agreement，不得置於僅允許本地工作的準備階段。正式 CLI／
graph 各 runtime 持有一份 immutable snapshot。省略 context 時保留無前綴介面，
在當次 advance 建立暫時 snapshot。使用者提供的自訂 callback 簽名不變。

驗收見 [1D 進度](progress/HPC_04A_ONE_D_OPTIONS_PROGRESS.md)。其他 body-fitted standalone
CLI、immersed／moving／FSI 路徑，以及完整 nested solver 診斷與後端矩陣仍待完成。
貼體基礎驗收見 [原報告](progress/HPC_04A_BODY_FITTED_OPTIONS_PROGRESS.md)。
既有 `.ntiga`、場輸出與 checkpoint payload 格式不變；新的 source identity 會讓舊建置的
checkpoint 明確不相容，續跑須保留相同建置與數值選項。
