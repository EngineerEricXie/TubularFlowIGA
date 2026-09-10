# PETSc 求解器選項

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

## Immersed static、transient 與 moving／FSI

Immersed graph domain 同樣使用 `domain_<id>_flow_`。相容性優先順序為：
既有 runtime defaults → `immersed_static_` 或 `immersed_transient_` family 選項 →
該 domain 的完整 prefix 選項。既有 static／distributed 路徑沒有繼承未加 prefix 的
`-ksp_type`、`-pc_type`；此行為保留。Serial transient 現在也接受 transient family
選項。例：

```bash
mpiexec -np 2 solvers/coupling/iga_multidomain_flow \
  --graph-case CASE --output-dir NEW_OUTPUT \
  -immersed_transient_ksp_type gmres \
  -domain_immersed_flow_ksp_type fgmres
```

此例的 `immersed` domain 使用 FGMRES；其他 transient immersed domain 繼承 GMRES。
Static family 仍使用 `immersed_static_`。未指定選項時，serial static 與 distributed
runtime 保留 GMRES／LU；serial transient／moving 保留 FGMRES／LU。多 rank LU 保留
既有 MUMPS 選擇；不表示任意 PC 都適合這些 saddle-point 系統。

Embedding caller 可設定 static／transient options 最後的 `solver_options_prefix`；
空字串使用原 family prefix。Moving options 的 `flow` 成員提供相同設定，FSI wrapper
在沒有明確 prefix 時依 fluid domain ID 產生。Moving geometry 的所有 committed／trial
epoch 共用建立時的 immutable snapshot；後續改動全域選項不會改掉同一 runtime 的策略。
Serial transient 的可選 shared options owner 必須使用同一 communicator、相符 prefix，
並維持到所有附掛 KSP 銷毀；一般 caller 使用預設 owner 即可。

Transient input hash 現在包含 snapshot 中的 prefixed options：fixed input identity
版本為 `ImmersedTransientInput/v6`，moving 為 `ImmersedTransientMovingInput/v3`。
這項身分更新不改場資料格式。不同 solver 選項的 trial 不再沿用同一輸入身分。
Native graph 會輸出 immersed `solver_configuration`；embedding runtime 也提供本地
`SolverConfiguration()` 查詢。

## Body-fitted standalone 與 VCA

`iga_navier_stokes` 依選定的 Navier–Stokes system 名稱配置 `..._flow_`；舊 boundary
configuration 沒有 system 名稱時用 `domain_flow_flow_`。VCA 的 transport system 另用
`..._transport_`。`iga_solve` 也使用選定 transport system 的名稱，例如：

```bash
PETSC_OPTIONS='-ksp_type gmres -pc_type lu -pc_factor_mat_solver_type mumps' \
mpiexec -np 2 solvers/cpu/iga_solve DATABASE.ntiga CASE_DIR \
  --system transport --output result.txt \
  -domain_transport_transport_ksp_type fgmres
```

兩個 CLI 均接受單個 `-` 開頭的 PETSc key 與 optional value；既有 `--` application
options 與 positional 介面保留。請先列出資料庫／case 位置參數，再列 PETSc options。
PETSc 自行解讀數值的語義；CLI 不另外保證所有負 tolerance 都會被本機 PETSc 拒絕。

`solver_configuration prefix=... ksp=... pc=... factor_backend=... step=... iterations=... reason=...`
記錄每步實際設定與最後一次 KSP 結果。一般 transport CLI 在後續步使用非零 initial
guess；保留其既有 GMRES 策略，PREONLY 不接受這個 warm start。VCA transport 的既有
呼叫方式不變。小案例驗收與限制見 [standalone 進度](progress/HPC_04A_STANDALONE_OPTIONS_PROGRESS.md)。

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

驗收見 [1D 進度](progress/HPC_04A_ONE_D_OPTIONS_PROGRESS.md)。immersed／moving／FSI 進度見 [immersed 報告](progress/HPC_04A_IMMERSED_OPTIONS_PROGRESS.md)。
貼體基礎驗收見 [原報告](progress/HPC_04A_BODY_FITTED_OPTIONS_PROGRESS.md)。
既有 `.ntiga`、場輸出與 checkpoint payload 格式不變；新的 source identity 會讓舊建置的
checkpoint 明確不相容，續跑須保留相同建置與數值選項。


## 求解器診斷

`-domain_junction_flow_ksp_view` 顯示該 domain 的實際 KSP、PC 與已建立的子求解器；
`-domain_source_flow_snes_view` 顯示 implicit nonlinear 1D 的 SNES。
Immersed 也可使用 family alias，例如 `-immersed_transient_ksp_view`，輸出仍顯示
實際的 domain prefix。`-domain_junction_flow_sub_ksp_converged_reason` 記錄
block-Jacobi 子 KSP 的收斂原因；只有實際建立該子求解器時才有輸出。
Fieldsplit 的子名稱依 runtime 使用的名稱配置，例如 `fieldsplit_0_`／`fieldsplit_1_`。

```bash
mpiexec -np 3 solvers/coupling/iga_multidomain_flow \
  --graph-case CASE --output-dir NEW_OUTPUT \
  -domain_junction_flow_ksp_view \
  -domain_junction_flow_sub_ksp_converged_reason \
  -domain_source_flow_snes_view
```

此範例需要含 `junction` 與 nonlinear implicit `source` 的 graph。
Viewer 失敗會回傳錯誤，須檢查 process exit status。驗證與本機 factor capability
矩陣見 [診斷報告](progress/HPC_04A_BACKEND_DIAGNOSTICS_PROGRESS.md)；該矩陣只代表
目前 PETSc build 對測試 AIJ 矩陣註冊的介面，不代表各 PDE 的效能或適用性。


## Legacy transport CLI

`iga_transport` 使用 `domain_neuron_transport_transport_`，未加 prefix 的選項仍為
共同基線。位置參數介面保留；使用環境變數傳入 PETSc 選項，例如：

```bash
PETSC_OPTIONS='-domain_neuron_transport_transport_ksp_type fgmres -domain_neuron_transport_transport_ksp_rtol 1e-12' \
mpiexec -np 2 solvers/cpu/iga_transport DATABASE.ntiga CASE_DIR 2 result.txt
```

既有 GMRES／block-Jacobi 預設不變，後續步仍使用非零初始猜測。新的每步
`solver_configuration` 記錄有效 prefix、KSP／PC、backend、iterations／reason。
範例的較嚴格 tolerance 是候選配置驗收所得，並非普遍效能建議；詳見
[legacy 驗收](progress/HPC_04A_LEGACY_OPTIONS_PROGRESS.md)。
