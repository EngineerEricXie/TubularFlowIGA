# Bifurcation FSI：PSC Bridges-2 時間效率改善與執行交接計畫

建立日期：2026-09-14（PSC 美東時間；UTC 2026-09-15）。讀者：接手實作的 gpt-sol／Codex 與專案維護者。
工作目錄：`/ocean/projects/mch260002p/thsieh1/TubularFlowIGA`。
本文件是待執行計畫，不代表其中的新優化已實作或通過測試。日期以執行環境為準；接手時記錄實際 UTC 時間。

## 0. 任務、優先順序與完成定義

使用者希望最大化時間效率，並避免浪費 PSC CPU。主要目標是縮短「提交工作到取得可驗證 FSI 結果」的時間，包含排隊、編譯、求解、輸出及失敗重做。CPU-hours 是必須一起報告的成本，不以最省核心犧牲大量等待時間，也不因節點有 128 核就直接申請整台。

執行順序：**凍結可信基準 → 便宜的 batch 設定篩選／呼叫計數 → 去除重複完整組裝 → 分離殘差／Jacobian → 優化最熱 kernel 與安全快取 → 合適的核心配置 → 完整心跳驗證**。

所有優化先保持原 PDE、材料、邊界條件、時間步、空間積分與收斂門檻。兩項成果分開交付：

1. `performance_verified`：在同樣問題與驗證準則下，證明時間縮短且未引入數值回歸。
2. `formal_physics_acceptance`：另以正式守恆、網格／時間收斂要求判定。本次展示基準未達此項，不能因變快就勾選。

預定第一里程碑是取得可重現的完整代表性 FSI 加速，開發目標為完整 20 步耗時至少減少 30%，進一步以 2 倍整體加速為研究目標。這些是目標而非承諾；沒有達到時如實記錄測量與限制，不放寬誤差來湊目標。每個有效優化都可獨立保留，不需等所有研究項目成功。

這次使用者要求的是寫計畫。接手收到「執行此計畫」後，按照任務卡自主完成已授權範圍；普通讀取、局部修正與 Slurm 小測試不反覆請求確認。遵守當時的沙箱、作業資源與使用者最新指示。本計畫不要求多代理，gpt-sol 可依序執行。

## 1. 已知事實：不要從零重新猜測

### 1.1 完整案例與硬體

- 原優化心跳作業：`45982619`，節點 `r257`，RM-shared，1 MPI process、8 OpenMP threads，15200 MB 記憶體配額，無 GPU。
- 求解 elapsed：13603.010422705 秒，約 3 小時 46 分 43 秒；Slurm 整體 3:48:41 包含建置／後處理。
- RM 硬體依 PSC 文件是雙 AMD EPYC 7742、每台 128 實體核心；此作業只分配 8 核。
- 背景 grid `16 × 16 × 5 = 1280`；初期 assembly 日誌有 268 usable cells。有效單元數需逐 epoch 量測，不假設永遠相同。
- 固定幾何 benchmark 有 4672 代數 rows。它是該狀態的流體系統大小，不是所有時間步固定 DOF，也不是輸出點數。
- 血管表面 812 vertices／1620 triangles；膜 806 nodes／1582 triangles、32 clamps。
- 20 步、dt=0.05 s、總時長 1 s；壓力 `5 + 15*(1-cos(2*pi*t/1s))` Pa、出口 0 Pa；rho=1060、mu=0.0035。
- small-displacement pretensioned membrane；實際雙向耦合。網格、材料與膜參數以原 driver／run.json／source snapshot 為準，不重抄後猜值。
- cut depth=1，空 cut-cell rescue depth=3；顯示倍率 100。倍率只影響顯示座標，不能拿放大後位移當物理解。
- 展示模式仍檢查 finite、coupling、small-displacement 等；17/20 步未通過原守恆門檻，`formal_acceptance=false`。

### 1.2 成本分解（exclusive wall seconds）

| 階段 | 秒 | 占求解時間 | 呼叫數 |
|---|---:|---:|---:|
| assembly | 8900.323421 | 65.43% | 695 |
| geometry | 2258.880152 | 16.61% | 240 |
| solver_setup | 1469.729249 | 10.80% | 229 |
| coupling 內其他工作 | 862.851873 | 6.34% | 20 |
| linear_solve | 22.703021 | 0.17% | 308 |
| communication | 0.004641 | <0.01% | 1414 |
| unscoped | 88.518065 | 0.65% | — |

共 79 外層耦合迭代，平均每步 3.95 次；平均每步 34.75 次 assembly 呼叫。phase 呼叫數可能包含不同路徑與巢狀 scope，**695 不應在尚未追蹤前解讀成 695 次完全相同的流體矩陣建置，240 也不等於 240 次完整幾何重建**。

`PetscPhaseProfile.hpp` 把 `KSPSetUp` 與 `KSPSetUpOnBlocks` 放進 solver_setup；數值分解／預條件器設定可能在此。因此不能說「所有線性代數僅 0.17%」。不要相加 inclusive assembly、geometry、coupling 時間。

### 1.3 已通過的效能變更與證據

已完成變更：volume、conservative trace、Nitsche wall 在 worker 計算；PETSc scatter／診斷發布仍由 caller 按原 cell 順序進行。

- `45981392`：expanded/compact × 1/2/4 threads，殘差／Jacobian action exact equality，含 wall fault 後 retry／abort 檢查。
- `45981715`：同一 Y 幾何，舊版 8 threads 組裝中位數 28.543797581 s、新版 8 threads 15.982731132 s，新版 1 thread 58.190114988 s；算子雜湊一致、matrix mallocs=0。
- 新版的 1→8 組裝 speedup=3.64081173，parallel efficiency=45.5101%。這不是全程 CPU utilization。
- 同作業的小型 27-cell FSI 單步，baseline 473.56193825 s、candidate 445.392055263 s，單組比較約省 5.95%；不是充分重複實驗，不能宣稱完整心跳也快同樣比例。
- 固定幾何三次 samples 在同一程序內，只有兩種不同 state；不可稱為三次獨立完整 FSI。
- 不要再實作一次已完成的 volume/trace/wall worker 平行化。

### 1.4 目前輸出狀態與來源路徑

基準 evidence root：`benchmarks/heartbeat-optimized-45982619/`。

| 檔案／目錄 | 用途 |
|---|---|
| `solver.stdout`、`solver.stderr` | 完整求解與累積 profile；取最後一筆 profile |
| `results/run.json`、`history.csv`、`coupling.csv` | 問題參數、逐步診斷、外層收斂 |
| `results/step-*/`、`fsi.pvd`、`wall-display.pvd` | 真實已接受時間步的 VTK 輸出 |
| `hashes.txt`、`base-source-snapshot.txt` | 原 binary/header/surface provenance |
| `ImmersedTransientFlowRuntime.hpp` | 已量測的優化 header |
| `export-format-repair.json` | 輸出格式補正前後 SHA256 |
| `results/heartbeat-verification.json` | 膨脹／收縮、數值欄位檢查 |
| `results/heartbeat-vtk-verification.json` | VTK 9.4.1 native reader 檢查 |
| `heartbeat-paraview.tar.gz`、`archive.sha256` | 已打包 ParaView 結果 |

原作業被 Slurm 標記 FAILED 是後處理 metadata 格式缺漏，不是 20 步求解失敗。修復作業 `46006012` 已 COMPLETED、exit 0、1:52；補上 `display_displacement_scale` 的 `NumberOfTuples=1`，不改數值。原錯誤紀錄保留。接手時核對以上檔案後，不再重跑此修復。

其他來源：

- `benchmarks/heartbeat_assembly_comparison.json`
- `benchmarks/heartbeat-performance-45981715.tar.gz`：old/candidate header、benchmark.cpp、integration wrapper、原始 logs。
- `benchmarks/heartbeat-regression-45981392.tar.gz`
- `benchmarks/heartbeat-optimized-launch.sbatch`
- `docs/progress/BIFURCATION_FSI_ASSEMBLY_PERFORMANCE.md`：歷史進度，末尾 pending 描述已被本文件的已完成 evidence 更新。

工作樹已有多項未提交修改與未追蹤檔案。**禁止 git reset --hard、git clean、用 HEAD 直接覆蓋工作樹，或將所有未提交內容當成本次變更。** 原 snapshot 在 `.nsvms_diagnostics`，曾遇共享檔案系統問題；先用上述可讀 archives，不反覆掃描大目錄。

## 2. 修改地圖與不能破壞的契約

路徑以 repository root 為基準；用函式名稱搜尋，不依賴本文件行號。

| 檔案 | 關注函式／責任 | 計畫 |
|---|---|---|
| `solvers/cpu/include/ImmersedTransientFlowRuntime.hpp` | `Assemble`、`SolveTrial`、`InvalidateSolved`、`ResetAttemptWork`、`CreatePetsc`、ports/gauge | P1–P4 |
| `solvers/cpu/include/MovingImmersedTransientFlowFsiRuntime.hpp` | `SolveFluidTrial` 的 BeginTrial→Assemble→SolveTrial | P1 |
| `solvers/cpu/include/ParallelElementBatch.hpp` | `ForEachElementBatch`，caller prepare/consume、worker compute | P2、P4 |
| `solvers/cpu/include/ElementAssemblyExecution.hpp` | thread 與 batch env、MPI caller gate | P0、P2 |
| `solvers/cpu/include/MovingImmersedTransientFlowRuntime.hpp` | `BeginTrial` 建 geometry/Epoch、extension/seed/map | P4、P6 |
| `solvers/cpu/include/MovingCutGeometry.hpp` | geometry/catalog build、幾何與 publication identity | P4、P6 |
| `solvers/cpu/include/ImmersedNitscheWall.hpp` | `BuildImmersedNitscheWallElementImpl`、material-aware overload | P3、P4 |
| `solvers/cpu/include/NavierStokesElement.hpp` | volume 殘差、Jacobian、積分 kernel | P3、P4 |
| `solvers/cpu/include/CutCellGhostPenalty.hpp` | ghost face local operator | P3、P4 |
| `include/PhaseProfile.hpp`、`solvers/cpu/include/PetscPhaseProfile.hpp` | 計時與 scope 語義 | P0 |
| `solvers/cpu/include/StrongFluidStructureCoupling.hpp` | Aitken、rollback、acceptance | P7（條件式） |
| `solvers/cpu/src/bifurcation_fsi.cpp` | 案例 CLI、輸出與收斂 policy | 只增加必要 profiling/export；不改物理預設 |
| `solvers/cpu/Makefile`、`solvers/cpu/tests/` | 建置與回歸 | 全階段 |

契約：MPI 為 FUNNELED；所有 MPI/PETSc 呼叫仍在初始化 MPI 的 caller thread。worker 使用 immutable input、private scratch；不可直接共用 PETSc Vec/Mat、寫共用診斷或做未同步 lazy cache。原本 ordered scatter 保持不變，除非另外證明重排的數值與效益。

所有更動維持 BeginTrial→Solve→PrepareCommit→FinalizeCommit／Rollback／Abort 的交易語義。不得用快取繞過 finite 檢查、stamp、material identity、moving-map、port/gauge 或守恆診斷。失敗後必須能恢復、重試，不能发布半完成結果。

目前 heartbeat 使用 `MatCreateSeqAIJ`、`VecCreateSeq`、`PETSC_COMM_SELF` 的這條 runtime。**不能把 `mpiexec -np 1` 改成 `-np 8` 當成分散式加速，否則可能只是重複計算同個問題。** repository 有 distributed runtime 不代表 heartbeat driver 已連上它。

## 3. P0：凍結基準、建立可重現的小型量測入口

依賴：無。第一個實作任務。先完成此階段，才改演算法。

### P0-A：來源與環境

1. 讀最新使用者指示、AGENTS.md（若存在）、本計畫；查 `git status --short`、`git diff --stat`、當前 commit。
2. 建立本輪 evidence 目錄：`benchmarks/fsi-time-opt/<run-id>/`。建立 `STATE.md` 與 `state.json`，格式見第 10 節。
3. 先比較 archived tested source 與目前工作樹，列出除 VTK writer 以外的差異；若有其他數值變更，先驗證等價或用已驗證 snapshot 建 B0，不能直接拿未核對的 HEAD／dirty tree 當基準。基準 B0 是「已通過的平行組裝版本 + VTK metadata 修正」。修正 XML 的 writer 不影響物理；仍需記錄差異。不要把舊 pre-optimization header 誤當目前 B0。
4. 以 `git ls-files --cached --others --exclude-standard` 取得檔案清單，再明確選入 `include/`、必要 `solvers/`、examples、Makefile 及驗證腳本，保存目前已修改及必要 untracked 來源。不要只用 `git archive HEAD`（會漏現有修改）；不要複製整個含大型 evidence 的 repo。
5. 保存每檔 SHA256、source archive SHA256、binary SHA256、編譯命令、PETSc configuration、`module list`、`mpicxx --version`、`ldd`、`lscpu`、`numactl --hardware`（有則記錄）、Slurm JobID/NodeList/AllocTRES、OMP/BLAS env。不輸出憑證或整份任意環境變數。
6. 原 Makefile 已是 `-O3 -std=c++17`。不要把「開 O3」當新改善；不要先開 fast-math。`-march`/LTO 等獨立最後測，不與數學修改混在一起。
7. baseline 與 candidate 從各自 immutable snapshot 在 node-local 建置；禁止執行中替換 binary 或 header。大型 header 測試第一次最多 `make -j2`，按 MaxRSS 再調整。

### P0-B：細項計數與 profiling

新增低成本診斷，預設關閉、只由一個 profiling 選項開啟；候選名稱 `IGA_PROFILE_DETAIL=1` **尚未實作**。保留現有 coarse profile schema；另寫 `assembly-detail.jsonl`，避免改壞既有報告 parser。

每次 assembly 記錄 step、coupling iteration、Newton iteration、呼叫原因（initial/diagnostic/line-search/next-Newton/other）、geometry epoch、state generation、residual generation、Jacobian generation、active DOF、inside/cut cells、volume/surface quadrature points、ghost faces、線搜尋 damping。

細項至少包括：prepare/gather、volume、trace、wall、ports、ghost、scatter、Mat/Vec assembly、port measurement/gauge、hash/diagnostic；geometry 另拆 classification/spatial index、volume rule、surface rule、ghost catalog、layout/preallocation、extension/seed/map。setup 若可用 PETSc `-log_view` 分解 symbolic/numeric factorization，先核對本地 PETSc 與 wrapper 是否接受。

**計時規則**：現有 PhaseProfile 是 thread_local。worker 的工作秒數不能直接加成 wall elapsed，也不能共享同一 PhaseProfile。使用 per-item 或 per-thread counter，barrier 後 caller 聚合；分别輸出 worker-work-seconds、batch wall seconds、max worker time、prepare/consume wall。不要每個 quadrature point 計時，先以 element/kernel 為單位。診斷版本與無診斷版本做短配對檢查；目標 overhead <3%，超過就採樣或減少欄位。

### P0-C：代表性 workload，避免只測最容易的一步

建立可重跑測試 driver，從已存 archive 的 benchmark/integration driver 整理，不重新寫另一套方程。

- W0：現有 tiny operator／failure regression，秒至分鐘級。
- W1：相同 Y 網格的固定幾何 assembly，zero + 非零 state；三次不是三種 state，報告要清楚。
- W2：相同案例的完整 moving FSI 前 2 步；正式門檻仍能通過的區段可作 regression。不得用全週期 verifier 要求這兩步必須已收縮。
- W3：升壓、峰值、降壓代表狀態，優先使用「完整 accepted solver state」的可重現 replay。

**禁止從 VTU 直接猜 solver restart**：VTU 是取樣點與展示資料，未必有全部 IGA coefficients、history、port multipliers、committed material state、Aitken/lifecycle。若目前沒有合法 checkpoint/replay，先確認現有序列化可用性；沒有就只新增 bounded benchmark snapshot/restore，包含 layout IDs、coefficients、time/index、history/extension、ports/gauge、結構狀態與 immutable kinematics。從相同 committed state 讓 B0/Candidate 都跑相同剩餘步；cold/warm coupling predictor 要一致。

Replay 必須先通過「連續跑 vs 保存再讀回」一致性。未建立可靠 replay 前，用 W2 篩選，最後完整 20 步驗证；不要為每個候選重跑整週期來取得峰值。

**P0 完成門檻**：B0 可重建、W0/W1 可重跑、細項時間能核對 coarse elapsed、列出至少前三個 hot path、已有不可混淆的 source/binary/input identities。未通過時不要同時改 P1–P7。

## 4. 第一批修改：低風險減少工作量

### P1-A：刪除 FSI adapter 的重複前置組裝

位置：`MovingImmersedTransientFlowFsiRuntime::SolveFluidTrial()`。

現状：`moving_.BeginTrial(...)` → `moving_.Assemble()` → `moving_.SolveTrial()`；後者第一輪再次 `Assemble()`。

修改步驟：

1. 用 P0 call tags 證明兩次組裝的 state、geometry、history、ports 等輸入完全一致，且中間無需消費第一份診斷。
2. 移除 adapter 的前置 `moving_.Assemble()`，讓 SolveTrial 負責初始組裝。不要刪除整個 public Assemble API，其他測試／diagnostic caller 仍需要它。
3. 檢查 fault injection 時機、`attempt_assembly_count` 與輸出 stamp。更新只因呼叫數減少而改變的計數預期，不放寬結果驗證。
4. 用 regression 確認 initial zero residual、第一個 assembly throw、rollback/retry、port/gauge、traction publication 全部正常。

理論上每個成功流體 trial 至少少一個前置 assembly；是否等於本案例的 79 次須看 call tags，不能直接保證省 79/695 的時間。此卡是第一個 production code 變更，獨立 diff／證據。

### P1-B：重用線搜尋已接受狀態的完整組裝

位置：`ImmersedTransientFlowRuntime::SolveTrial()`。

現狀：線搜尋在 candidate state 已建好 R/J；accepted 但未達收斂時跳出 while，下輪 for 又對同狀態建 R/J。

先做**SolveTrial 內部區域有效性狀態**，不要一開始就做跨 epoch 全域 cache。第一輪建一次；接受 candidate 後，若所有輸入不变且持有完整 R/J，下輪直接使用；state 恢復、port/control 改變、BeginTrial、rollback、exception 均使其無效。

- 緩存有效與 `diagnostics_.converged` 是不同語義，不能混用。
- 下一輪仍重新執行原本的 residual/block/controller convergence 檢查；保留 initial norm 的定義及 Newton iteration 計數。
- 不用 `memcmp(state)` 每次全量掃描來代替正確生命周期；使用內部 generation/ready flag，對外 state setter 都需覆蓋失效路徑。
- 測正常一步收斂、多步收斂、第一次 line search 被拒／接受較小 damping、最終 backtracking failure、NaN、retry、exact R/J action。

**P1 gate**：同輸入的 operator bitwise 相同，物理解與收斂診斷一致（允許已說明的純工作計數改變），full assembly 次數確實下降，W2 無退化。只有正確性但無節時也先保留測量，不宣稱完成加速。

### P2：利用現有 batch 設定提升 worker 利用率

可與 P1 前後各量一次，先做不改數學的設定篩選。

`ElementAssemblyExecution` 已支援 `IGA_ASSEMBLY_BATCH_SIZE`，預設 `min(threads,8)`。`ForEachElementBatch` 每批 prepare → parallel dynamic compute → join → consume。8 threads/8 slots 在複雜度不均時幾乎沒有餘下工作可平衡；>8 threads 預設仍只有 8 slots。

1. 固定 8 threads，用 W1 比較 batch=8、16、32；若 32 有效、RSS 允許再測 64。一次只改 batch。
2. 記錄每批 worker completion spread、prepare/consume 比例、resident result bytes、MaxRSS；選最小且接近最快的 batch。
3. 考慮按 cost 排程 compute，但結果按原 cell index consume；不得更改全域求和順序。先使用既有 dynamic scheduling，再考慮新排程。
4. 若 prepare/consume 占比明顯，下一步才做 thread-local scratch 重用或 persistent team。`omp single` 的執行 thread 不保證是 MPI 初始化 caller，**不可把含 PETSc/MPI 的 consume 隨意搬入 single 區**。
5. batch 越大同時保留 dense volume/trace/wall 越多；設定 RSS 上限，超過配額 70% 停止擴大。不能為一點 batch 收益盲目占整台節點。

**P2 gate**：expanded/compact、serial/OpenMP、worker throw→retry、ordered scatter exactness 通過；至少 W1 可重現改善、W2 不退化。batch 選擇寫入 profile/manifest，不藏成未記錄環境依賴。

## 5. 第二批修改：主要熱點與資料重用

### P3：殘差／Jacobian 分離，減少線搜尋的矩陣工作

依賴：P0、P1 通過。此項較大，分成 3 個可回退提交／patch。

**P3-A 結構整理，先不改算術**：設計內部 `AssemblyRequest { ResidualOnly, ResidualAndJacobian }`（提案 API，尚不存在）。抽出共用 quadrature、forcing、history、wall/trace/ghost/port/gauge 運算。不可在新 residual path 另抄一套弱式。

**P3-B 真正 ResidualOnly**：不得只是在散佈時不寫 Mat；kernel 必須跳過 dense Jacobian 建構與其他純 J 工作。保留所有會影響 R、finite、wall penalty、controller、mass/flow 判斷的量。trace、Nitsche、ghost、ports、gauge 都要覆蓋，包含非零 wall velocity 與非零 multipliers。

**P3-C 接入 line search**：需要新方向時用 full R/J；每個試探 damping 只算 R；accepted 後若直接收斂，不再建 J；尚需下一輪 Newton 才建 J。P1 的 full-cache flag 不能讓 residual-only 被誤認為持有新 Jacobian。用兩個 generation 表示 R 和 J，只有與當前 state/epoch 相符者可用。接受後是否保存共用 quadrature intermediates 由 profile 決定。

驗證：

- full path 的 R 與 residual-only 在相同 state 應 bitwise 一致（首选保持加總順序）；如果編译优化导致差异，先查原因，不直接放寬。
- 中心差分方向導數檢查 Jv 與 residual 差分，使用多個 epsilon、觀察合理的截斷／roundoff 行為；避開 nonsmooth branch，另測 branch 邊界。不用單一 epsilon 的寬鬆 pass 掩蓋错误。
- 多 damping backtracking、zero initial residual、nonfinite wall、ghost/pressure gauge/port-controller、moving history、abort/retry。
- 分別記 `residual_only_calls`、`full_assembly_calls`、各自時間、matrix mallocs；期待 full 次數下降，不能只改計數。
- 比較 W2 完整結果與峰值／降壓 replay（若可用），決定是否進入 full cycle。

若 coupling/controller 假設 Jacobian 仍屬於當前 state，先修正明確的依賴介面；禁止將其悄悄停用。若此項收益低於實作成本，先保存結構整理與證據，轉 P4。

### P4：由 P0 hottest kernel 決定的安全快取與組裝加速

每次只選占比最高的一項，避免同時實作下表所有 cache。

| 可考慮 cache | 合法重用範圍 | 必須失效的情況 |
|---|---|---|
| 背景 element/connectivity、基底/reference derivatives | 同 grid/order；值還要綁積分點 | grid/order 或 quadrature points 改變 |
| 已算的物理點 basis/gradient/weights | 同一 geometry epoch/rule | target geometry、rule、坐標映射變更 |
| geometry-only wall penalty 幾何部分 | 僅先證明不依賴速度／history／dt 的部分 | gamma、mu/rho/dt、rule、normal、ghost binding 等依賴變更 |
| 固定 ghost linear block／index lists | 同 geometry、layout、viscosity 與 penalty options | cut/ghost 拓樸、DOF、物理係數變更 |
| scatter row/column indices | 完全相同 layout/ports/gauge | 活躍 cell/DOF、connectivity、appended scalar rows 改變 |
| worker scratch vectors/matrices | 同 thread 所有權；每次 reset 正確 | 尺寸變更、exception 後重用 |

在 geometry epoch 建立 immutable cache，再交給 workers；避免 worker 首次寫共享 cache。cache key 列出每個依賴，不以「位移很小」判斷幾何沒變，也不以 row count 相同判斷 pattern 相同。

若 profile 顯示大量 scalar `MatSetValue/VecSetValue`，先批量 gathers/scatters、預存 index array，仍保留原項目順序與 scalar diagonal。現有矩陣已預分配、mallocs=0；「加 preallocation」不是新收益。PETSc COO/blocked insertion 需確認本地版本、零對角與 port/gauge pattern，再以獨立候選測試，不直接換儲存形式。

若 ghost 或 ports 超過總時間 5%，可比照 element 的 prepare/worker-local compute/ordered-caller-consume 平行化。不要讓微小序列區的改造擠掉 volume/wall 大熱點。

**P4 gate**：cache hits/misses/bytes 可見；同 epoch 重用、不同 epoch 失效、係數／port 變更、cut topology 改變、fault injection 後 retry 全通過。保存 disabled-cache 對照入口，測 matched W1/W2。只在證明熱點後改，不预设是哪種特徵值／積分公式最慢。

### P5：求解器設定與預條件器重用

依賴：P1/P3 對最終 matrix 更新頻率已穩定。此階段對準 10.8% setup，不是 0.17% solve。

1. 核對實際 `KSP/PC`、LU backend、pivot shift、symbolic/numeric factorization 時間。源碼預設 FGMRES + LU，但執行選項可能覆蓋；日誌是最終依據。
2. 先試完全相同 pattern 的 symbolic ordering/fill reuse；讀本地 PETSc API/版本，記錄設定。
3. 再試將前一次因子當近似 preconditioner，使用**當前** Jacobian 解當前 Newton 方程。注意「preconditioner reuse」與「直接用舊 LU 當精確新解」不同，後者不允許。
4. 對 KSP 次數增加、true linear residual 不達標、stagnation 自動重建一次；保留重建次數與成本。不能關掉 true residual 檢查。
5. 先限同 geometry epoch；跨 epoch 重用需 P6 的 layout/pattern 契約，不能帶著 stale handles 越過 new Epoch。

SNES lag 選項是參考設計；目前自寫 Newton 並不會自動吃 `-snes_lag_jacobian`。避免提交只有 option 沒有作用的「優化」。不優先換 AMG/MUMPS/GPU；先量這個約幾千 rows 的實際 factorization。

### P6：geometry／runtime 重用（條件式，不列為首批必要大改）

`MovingImmersedTransientFlowRuntime::BeginTrial` 每次建立 target geometry、新 Epoch/PETSc runtime、extension、scalar history、seed/map。240 geometry scope 次數不能當成 240 次完整建置。

先拆時間，復用背景 reference 資料與候選 spatial search；只有在有效 geometry/layout 完全匹配時重用 index/preallocation。小位移仍可能改 cut topology、volume、normal，必須重算依賴這些量的積分與守恆項。不能藉由固定壁面／凍結 cut geometry 冒充等價加速。

不直接跨 runtime 搬 Mat/Vec/KSP ownership。若確有收益，設計獨立 immutable pattern cache，先測 old/new active-set 各類轉換及 fault-stage cleanup；保留原完整建置 fallback。完整物理 geometry reuse 屬較高風險研究，需要精度比較，排在 P1–P5 後。

### P7：耦合加速、warm start 與進階架構（有證據才開）

目前已有 Aitken，每步約 4 次；先把每次 fluid trial 做便宜。

- 若優化後 coupling 次數仍主導，可比較跨步 predictor 或 IQN-ILS。明確保存／清除 history，rollback 不污染下次，內層和外層容差依量綱分別評估。
- Warm start 必須將上一個 coupling trial 的速度／壓力合法映射到新 geometry，且保持時間離散的 committed history 不變；不要把 trial 解改成前一時間步真值。
- 不縮減步數、不提高 dt、不放寬 conservation/coupling，作為效能比較的捷徑。
- Matrix-free/partial assembly、GPU、完整 MPI heartbeat 是後續專案；需要 kernel 支援與預條件器、資料布局，不能只設旗標。若 P0–P6 已達目標，不自動展開這些大改。

## 6. PSC 執行策略：同時優化排隊時間與作業時間

### 6.1 配額與節點規則

- Login node 只做讀檔、編輯、git、Slurm 協調與小型摘要。編譯大量 C++、solver、完整 VTK 檢查、壓縮大量資料都放 Slurm compute allocation。
- 預設 account `mch260002p`、partition `RM-shared`、單節點、1 MPI process。不要永久指定 r257/r253；過去節點能讀檔不代表未來只能用該節點。
- 接手先查 `squeue -u "$USER"`、`scontrol show partition RM-shared`、`sacct`。不取消不屬於此計畫的作業。
- 小型編譯／回歸先 4 cores、7600M、30–60 min；算例 benchmark 常用 8 cores、15200M、30–60 min。這是起始規劃，依測量時間/RSS 調整。不把單步已需 10 分鐘的工作塞進 5 分鐘時限。
- 記憶體使用 `1900M × allocated cores` 作安全起點，避免歷史上 `2G` 換算超出每核心限制的問題。當前查到 DefMemPerCPU=1900；提交前仍核對最新 policy/QOS，不把舊 2000 MB 限制當永遠不變。
- 1-thread 測試若仍占 8-core allocation，Slurm 計費仍是 8 cores。這在同節點公平 scaling 短測可接受，但必須把 allocated-core-hours 和 algorithmic-thread-hours 分開；不能報成只用了 1-core 配額。
- 基本細項計時用程式內 timer，無需整台機器。PSC 官方規定 `-C PERF` 不可用在 RM-shared；若確需硬體 counters，另做短 RM 專用 profiling 作業並記錄整節點費用。不要因 r257 顯示 PERF feature 就假定 shared job 可使用。
- 同時最多一個主要 benchmark／模擬作業；有相依的候選用 `afterok`，先小篩選再選優勝者，避免 array 一次放出多個長跑。讀文件／改下一個獨立小 patch 可在等待時做。
- 目前沒有 GPU kernel、沒有已驗證 MPI heartbeat adapter。此計畫前半段不申請 GPU 或跨節點；16 cores 以上測量只有在前一級顯示收益且 batch 足夠時才開。

### 6.2 Node-local staging 與保存證據

PSC 的 `$LOCAL` 是 allocation 存活期間的節點本機工作區；作業结束後會被清除。使用它建置、跑 solver 與寫詳細 logs，再把必要結果搬回 Ocean。來源：[PSC User Guide，LOCAL 與 PERF 規則](https://www.psc.edu/resources/bridges-2/user-guide/)。

- 用**一個小型 source archive**搬入 `$LOCAL/fsi-time-opt-$SLURM_JOB_ID`，避免 compiler 從共享 FS 反覆讀幾百個 headers。先驗 archive hash，再解壓；本次工作不得依賴舊 allocation 的 `/local`。
- 已有測試曾在共享 FS header read 停住、CPU 接近 0，出現 transport endpoint 問題。這是環境／I/O failure，不要算進 kernel 速度、不改數值程式解決。
- 若 PETSc library/header 所在 Ocean 路徑也有異常，先量最小必要檔案讀取；需要 staging 時只複製相關建置依賴並保留 linkage，不能盲目複製整個 PETSc tree。
- binary、核心 source、input、runner 在提交前凍結；輸出使用唯一 JobID 目錄，絕不覆寫 B0 evidence。
- EXIT trap 盡力保存 status、stdout/stderr、manifest、timing；失敗 log 也要複製。保留原 solver exit code；copy-back 失敗另外標記，不偽裝成 solver 成功。
- 長模擬每個 accepted step 用暫存檔＋rename 發布 step outputs/manifest 到 Ocean；量測 copy 時間。不在 solver 正寫 VTU 時複製半個 step。
- 避免每次 assembly 印一整行與頻繁全目錄同步；每個 accepted step 更新一次摘要，詳細記錄留本機。
- 心跳 driver 的可用 checkpoint/restart 能力要先檢查，不能照搬舊 Graph 的 SIGUSR1/requeue 承諾。沒有 restart 就按量測估計 walltime 並留緩衝；signals/EXIT trap 無法挽救 SIGKILL 前未保存的記憶體狀態。
- 正式 solver、verification、render/package 狀態分開記錄，必要時不同 job；後處理失敗只重跑後處理。沿用 46006012 的教訓，不能因 XML metadata 問題重跑 4 小時 solver。

若 sandbox shell 回報確切的 `bwrap: Creating new namespace failed ... ENOSPC`，按使用者環境規則改用必要的 escalated lightweight command，不重試 sandbox、不歸因 Python/PETSc。若升權被自動審查拒絕，保留原因並遵守當時工具規則；本計畫不繞過權限。

### 6.3 Slurm 基礎範本（待接手建立 runner 後才提交）

以下是**範本而非已建好的可用 benchmark**。接手要先實作 `benchmarks/fsi-time-opt/run_phase.sh`，它接收 case list、寫 machine-readable 狀態、遇驗證失敗回傳非零，再使用此骨架。不要把不存在的 runner 當作已通過的驗證。

```bash
#!/bin/bash
#SBATCH -J fsi-time-opt
#SBATCH -A mch260002p
#SBATCH -p RM-shared
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --cpus-per-task=8
#SBATCH --mem=15200M
#SBATCH -t 01:00:00
#SBATCH -o benchmarks/fsi-time-opt/slurm-%j.out
#SBATCH -e benchmarks/fsi-time-opt/slurm-%j.err
set -euo pipefail
module load anaconda3 openmpi/4.0.5-gcc10.2.0
export PETSC_DIR=/ocean/projects/mch260002p/thsieh1/TubularFlowIGA-hpc-acceptance/petsc
export PETSC_ARCH=arch-bridges2-acceptance CXX=mpicxx
export OMP_NUM_THREADS="$SLURM_CPUS_PER_TASK"
export OMP_THREAD_LIMIT="$SLURM_CPUS_PER_TASK"
export OMP_DYNAMIC=FALSE OMP_PROC_BIND=close OMP_PLACES=cores
export OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 BLIS_NUM_THREADS=1
export IGA_ASSEMBLY_THREADS="$OMP_NUM_THREADS"
export IGA_ASSEMBLY_BATCH_SIZE=16   # 僅候選，先與 8 比較，非已證實最佳
export IGA_PROFILE=1
: "${LOCAL:?Run inside a PSC compute allocation with LOCAL available}"
# runner 負責 frozen archive staging、hash、build、tests、timing、EXIT 保存。
# 提交前建立 Slurm output 的 parent directory；不得在 job 啟動後才建。
# bash benchmarks/fsi-time-opt/run_phase.sh --case-list ...
```

沿用既有 OpenMPI 啟動方式；以下只表示單次啟動形式，實際 binary/參數必須來自 case manifest：

```bash
/usr/bin/time -v -o timing.txt \
  mpiexec -np 1 --map-by slot:PE="$OMP_NUM_THREADS" --bind-to core \
  ./candidate-binary ... > solver.stdout 2> solver.stderr
```

測試 1/2/4/8 threads 时，要同時設定 OMP_NUM_THREADS、IGA_ASSEMBLY_THREADS、MPI PE binding，並確認 log team_size。不在外面再包一層 `srun` 造成雙重 launcher。記錄 Slurm cpuset、affinity、OMP placement；`close`/`spread` 只有 profile 顯示 NUMA/帶寬問題時短測，不主觀認定哪個一定較快。

### 6.4 如何選核心數：優先時間，但不盲目擴張

先用 P1/P2/P3 後較便宜版本短測 1/2/4/8 threads；固定同一 binary、workload、batch 政策明確記錄。報告兩組：固定 batch 的純 thread scaling，以及每個 thread 數最佳合理 batch 的調校結果，不能混稱單一 scaling curve。

對每種配置記錄：

- `T_solver`、`T_job`、queue wait、stage/build/export 時間。
- `allocated_core_hours = AllocCPUS * ElapsedRaw / 3600`。
- `speedup(p)=T(1)/T(p)`、`efficiency(p)=speedup(p)/p`，限定同一 workload/code。
- CPU utilization 用可靠 job/step CPU accounting，不能用上述 efficiency 代替。MPI launcher 的 GNU time 資料未必涵蓋所有子程序，先核對計帳邊界。
- MaxRSS、matrix mallocs、KSP/outer iterations、numerical pass/fail。

預設工程決策（可按新證據調整，需記錄理由）：選在正確性通過且耗時接近最快者中的較少核心配置，例如距最快 5% 內取核心較少者。8→16 的探測只有在 4→8 至少省約 20% 時才進行，並確保 batch≥threads；16→32 同理。若加倍核心只省 <10%，不再擴大。這些是停止無效測試的篩選規則，不是聲稱使用者已設定固定 SU 預算。

觀測到排隊時間更長時，把 queue wait 纳入結果取得時間；`squeue --start` 僅預估，不當確定承諾。不要把硬體利用率接近 100% 當最終目標；做更少工作、較早完成，可能比所有核心一直忙更好。

## 7. 驗證矩陣與效能比較方法

### 7.1 測試分層，避免每次都跑全套長測

| 層級 | 內容 | 何時跑 |
|---|---|---|
| V0 | 編譯、batch/exception/caller 約束、目標 unit regression | 每張修改卡 |
| V1 | 相同 R/J action、residual-only、cache 失效、多種 storage/threads | 改組裝或快取後 |
| V2 | 小型 compliant-channel 完整 FSI + moving lifecycle tests | 通過 V1 的候選 |
| V3 | W2 Y-vessel 2 步與 W3（有合法 replay 才跑） | 最佳少數候選 |
| V4 | 原問題完整 20 步 + CSV／VTK 驗證 + archive | 確定勝出版本一次；新失敗／實質更動才重跑 |

已存在的 Makefile targets（接手先 `rg` 再用，以免變更後失效）：

- `parallel_element_batch_test`、`parallel_element_batch_openmp_test`
- `element_assembly_execution_openmp_test`
- `parallel_immersed_volume_test`、`parallel_immersed_volume_serial_test`
- `immersed_transient_flow_test`、`immersed_transient_flow_openmp_test`
- `moving_immersed_transient_flow_test`
- `moving_immersed_flow_snapshot_test`、`moving_immersed_flow_snapshot_publisher_test`
- `immersed_conservation_roundoff_test`
- `test_compliant_channel_fsi.cpp`、`test_strong_fluid_structure_coupling.cpp` 的 target 用 Makefile 搜尋確認，不憑名字猜。

建置 target 不等於測試已執行。先查 `main` 的 MPI 初始化：部分現有 tests 只有 PetscInitialize，直接設 OMP>1 可能拿到 MPI_THREAD_SINGLE。參考 `45981715` archive 的 wrapper，在 PetscInitialize 之前 `MPI_Init_thread(..., MPI_THREAD_FUNNELED, &provided)`，驗證 provided，並在適當 ownership 下 finalize。或正式補上 OpenMP test 的初始化。**不要放寬 ElementAssemblyExecution 的 thread gate 來讓測試過關。** `-UNDEBUG`／assert 是否真的開啟也要記錄。

測試命令放 case manifest，避免 shell 大段串接難以歸因。不要一次 `make all` 或跑所有 historical distributed acceptance：先完成受影響路徑所需測試，再於要交付時跑既定必要整合檢查。

### 7.2 正確性門檻

**完全等價變更（P1、P2、geometry-only cache）**：同 platform/build settings/輸入下，R 與多個 Jv 的 bitwise equality 作首要門檻。逐步解／診斷應一致；因取消多餘組裝改變的 counters、timing 或證據 hash 欄位要逐項解釋。不要把整份 diagnostic JSON 的 hash 當物理 state hash。

**改求解路徑（P3 或 P5 造成不同迭代路徑）**：保持原 convergence thresholds，另外比較 coefficients、流量、壓力、結構位移、速度、牽引與守恆量。先讓 B0 自身重跑估計 deterministic/roundoff 差異，再凍結比較公差；不看 candidate 結果後才放寬。

可用作初始嚴格比較的設計是 `||xC-xB|| <= atol_quantity + rtol*max(||xB||, characteristic_scale)`，每種量獨立定義 SI scale 與公差。先以 rtol=1e-8 作 algebraic consistency 調查起點；它不是對所有模型保證適用的物理精度門檻。給近零殘差套 relative error 沒意義，不以 1 Pa、1 m 等任意大絕對值掩蓋差異。若需要非 roundoff 級差異，必須附原因、原收斂準則與誤差敏感度，標為「數值等價容差內」，不能稱 bitwise。

- Fluid 用相同 global DOF IDs 的 coefficients，或事先定義共同物理取樣點；不要把兩個不同 geometry 的 VTU 點序號直接相减。
- Membrane 按 global_node_id、reference coordinates 對齊，比較實際位移/速度/traction。100 倍展示坐標不作 accuracy 指標。
- 正式小型案例保留原守恆門檻。展示基準已失敗的守恆值也需比較，不能以「本來就可忽略」容許新增偏差。
- NaN/Inf、ghost coverage、empty rule rescue、time/index、ports/gauge、material/geometry identity 不可略過。
- 結構 BE velocity 與 displacement difference/dt 一致；每個 accepted step 的 coupling residual 仍達原門檻。

### 7.3 公平效能量測

1. 同一 allocation 依序跑 B0/Candidate，固定 thread/binding/input/build flags；不要同時在同節點跑兩個候選互相干擾。
2. W1 每版先暖機一次，再做至少 3 個配對樣本，交替 AB/BA 或固定 seed 的隨機順序。W2 先單對篩選，只對勝出者追加到 3 對。非獨立程序／非獨立 state 要明確註記。
3. 報 median、min/max 或 IQR、每對 speedup、編譯與 stage 是否排除。完整結果取得時間另外報，不能只留最漂亮 kernel 數字。
4. 差異 <5% 或與樣本波動相當：標為不確定；只在此差異會影響選擇時補到最多 5 對。不無限制增加樣本。
5. 如果 candidate 比基準慢 >10% 且正確，先查 profile，停止更大 workload。若硬體/FS 異常造成樣本無效，標記環境失敗並保留原始資料，不把它刪掉假裝沒發生。
6. 每次只改一個可歸因項；最佳版可以累積已通過改善，但同時保留 B0 與上一個通過版本兩種對照。
7. 完整 20 步先與歷史 baseline 13603 s 比較，但它在另一時間/node load，屬歷史參考。要聲稱精確全程 speedup，需配對 full B0/Candidate；只有短測無法代表時才追加此昂貴驗證，不預設重跑多輪 4 小時。

### 7.4 每階段最少交付證據

`manifest.json`（來源/配置）、`case.json`（物理/數值）、`timing.json`（wall/cpu/phase/排程）、`correctness.json`（每項檢查）、stdout/stderr、`summary.md`。每個 record 有 variant、test_id、job_id、node、timestamps；status 至少分 `passed`、`failed_numerical`、`failed_environment`、`inconclusive`，不能將沒有資料默認 passed。

## 8. 限定成本的執行排程與轉向條件

以下是順序／工作量的起始安排，不是硬體速度保證。gpt-sol 依實際測量更新預估，避免時限過長排隊或過短重跑。

| 階段 | 工程工作 | 第一輪計算範圍 | 通過後 |
|---|---|---|---|
| R0 | P0 baseline、call tags、runner | 4-core build + W0/W1；不跑完整週期 | 找 hot path |
| R1 | P1-A | V0/V1 + 一組 W2 | 保存最快可信版 |
| R2 | P1-B、P2 分別比較 | W1 每版 3 對，最佳版 W2 | 決定 residual 分離收益 |
| R3 | P3 | 小型 operator → moving FSI → W2 | 過關才升級 baseline |
| R4 | P4 只做最熱 cache/kernel | 每次一項、同上筛选 | setup 若占比仍高則 P5 |
| R5 | P5、必要 P6 | 先短測，再挑勝出者 | 不必要則跳過且記原因 |
| R6 | thread/batch 最佳配置 | W1 1/2/4/8；W2 前兩名 | 選 production configuration |
| R7 | 全週期驗證／VTK/package | 20 步最佳版一次 | 最終報告與交接 |

若 P1/P2 已取得足夠收益，可以先跑一次完整候選交付，再評估後續；若仍需 P3/P4，避免每張卡跑完整週期。每卡提交前估計 allocated-core-hours；短測沒有方向時停下分析，不用更大作業「碰碰運氣」。排队时可以准备不依赖量测结果的文档/测试，但不能假定排队候选已通过并升级正式binary。

預估收益用 Amdahl 作情境而非承諾：只把 assembly 成本減半、其他不變，13603 s → 約 9153 s（2h33m）；減至四分之一 → 約 6928 s（1h55m）。去掉呼叫數與加速單次 kernel 有相互影響，不能把兩項獨立加速因子直接相乘。硬體 thread scaling 也不能套到 geometry/setup 的序列部分。

## 9. 常見失敗與具體處理

| 現象 | 優先檢查 | 處理 |
|---|---|---|
| >8 cores 仍不快 | batch slots、team_size、MPI binding、串行占比 | 先調 batch／確認 active workers，不增加 nodes |
| 編譯 wall 很長、CPU 幾乎 0 | Ocean header read、transport endpoint | node-local stage、保存環境診斷；不改 solver |
| MPI_THREAD_FUNNELED gate failure | test main 初始化與 provided | 修 wrapper/init，保持安全 gate |
| cache 之後結果只在第二步錯 | epoch/history/layout/ports 失效規則 | 關閉該 cache，比對依賴，新增最小跨步回歸 |
| line search 找不到下降 | residual-only 是否漏 wall/ghost/port/gauge 或 stale state | 比較同 state full R；回退 P3，不放寬 convergence |
| 少組裝後 traction/stamp 錯 | 前置 Assemble 是否有隱藏的資料發布依賴 | 分離明確計算與發布責任，保留交易性 |
| factor reuse 增加 KSP／失敗 | 新 matrix、舊 PC、true residual | 重建 PC 並重試一次，計入時間 |
| numerical pass 但加速不穩 | shared node noise、warmup/order、FS | 有限補測，標 inconclusive，不挑最好一次 |
| solver 成功、VTK/壓縮失敗 | metadata、native VTK、磁碟／路徑 | 只重跑 export/verification；不重算 |
| wall-time 接近 | accepted-step 存檔／真正 restart 是否可用 | 保存 evidence；無 restart 不假稱可 requeue 接續 |

更動觸及過多 lifecycle/contracts、兩次局部修正仍無法得到清楚一致性時，回到上一個已通過 variant，只撤回本次 patch；記錄未解問題並先完成其他有收益工作。不要用 reset/clean 清掉使用者原本修改。

## 10. 給 gpt-sol 的接續規約與狀態檔

每輪只處理一張主要卡。先讀 `STATE.md`，以 evidence 決定下一步；不因對話壓縮就重做已通過試驗。新使用者訊息以最新指示為準。

建立於 `benchmarks/fsi-time-opt/<run-id>/state.json` 的建議 schema：

```json
{
  "schema_version": 1,
  "plan": "docs/plans/BIFURCATION_FSI_PSC_TIME_OPTIMIZATION_PLAN.md",
  "baseline_manifest": null,
  "current_task": "P0-A",
  "accepted_variant": null,
  "tasks": {},
  "jobs": [],
  "next_action": "Freeze current source and inspect baseline artifacts",
  "pending_questions": [],
  "formal_physics_acceptance": false
}
```

每張卡記 `status: pending|running|passed|failed|skipped`、source hash、delta patch、job IDs、gate checks、measured speedup、RSS、reason、next_action。`STATE.md` 限制為容易讀的摘要：已完成、執行中、下一個精確命令或待實作函式、不能做的已知錯誤。

執行中監看至少每次有新 accepted step／stage/failure 時整理一次；對話保持簡短但不只重複「仍在算」。job 不變時不要每幾秒查一次 Slurm。已具備自主完成條件時繼續，只有依賴使用者尚未提供的實質資訊才停下。

完整交付必須包括：

- [ ] immutable B0 + candidate source/binary/input identities 與各自 patch。
- [ ] P1/P2/P3/P4 等每項是否執行、通過或跳過的理由；不需宣稱所有研究分支成功。
- [ ] 至少一個代表性完整 FSI 的公平前後比較，明確區分小型與全週期。
- [ ] thread/batch/core-hours/RSS/排队與求解時間表，推薦 PSC submission 配置。
- [ ] full 20 步結果、coupling/finite/物理差異檢查、原守恆失敗資訊。
- [ ] native VTK 與膨脹／收縮 verifier 通過、可下載的 ParaView archive 及 hash。
- [ ] `docs/progress/BIFURCATION_FSI_TIME_OPTIMIZATION_RESULTS.md`（待建立）記錄結果、限制、未解熱點。
- [ ] 尚有未提交使用者修改時，只整理本次明確 diff；是否 commit/push 依最新授權，不為了交付污染原工作。

## 11. 參考資料與適用邊界

查閱於 2026-09-14（PSC 當地日期）。以下用於設計參考，不是同案例速度實測。

1. [PSC Bridges-2 User Guide](https://www.psc.edu/resources/bridges-2/user-guide/)：RM/RM-shared、LOCAL、Slurm、PERF 和 accounting；實際限制提交前查 scheduler。
2. [PETSc SNESSetLagJacobian](https://petsc.org/main/manualpages/SNES/SNESSetLagJacobian/) 與 [SNESSetLagPreconditioner](https://petsc.org/main/manualpages/SNES/SNESSetLagPreconditioner/)：R/J/PC 更新頻率分離的參考；目前自寫 Newton 需實作相應機制。
3. [MFEM Performance and Partial Assembly](https://mfem.org/performance/)：高階張量結構與 matrix-free／partial assembly 的長期方向。cut-cell 非規則 quadrature 不直接享有同樣結構。
4. [preCICE Acceleration Configuration](https://precice.org/configuration-acceleration)：Aitken 與 IQN-ILS/IMVJ；目前已有 Aitken，不能標成尚缺功能。
5. [COMSOL Fully Coupled vs Segregated](https://www.comsol.com/support/knowledgebase/1258)：不同耦合分組的 tradeoff；不以商業品牌推定必然較快。[COMSOL nonlinear solver 說明](https://doc.comsol.com/6.4/doc/com.comsol.help.heat/heat_ug_modeling.06.37.html) 中 Jacobian 更新討論是一般策略參考，該頁介面為熱傳，不是血管 FSI 預設。
6. [SimVascular svFSI](https://github.com/SimVascular/svFSI)：血流 FSI 外部 benchmark 候選；若另建比較，先匹配膜／固體模型、邊界、空間時間精度、收斂與硬體。先不為評估速度重寫案例或安裝大型替代求解器。

目前沒有經同案例證實的 COMSOL／svFSI 相對加速數字。本計畫先把本專案的可歸因成本降低；外部完整比較是後續獨立工作。
