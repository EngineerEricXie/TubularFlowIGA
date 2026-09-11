# HPC-01C：入口與清理稽核

日期：2026-09-09；2026-09-11 完成最終簽核。基準 `d9eaea9` 加本批 runtime
cleanup 修改。此表記錄逐步核對的呼叫鏈；目前 HEAD 的完整結果見
[F06 最終簽核](HPC_01C_F06_COMPLETION_REPORT.md)。

2026-09-09 接續以 `f9c7223` 核對下列五個 main 的控制流程，並補強 memory
report 的 terminal close；此項驗收見 [memory close 報告](HPC_01C_MEMORY_CLOSE_PROGRESS.md)。

## 已核對的 main 控制流程

下表核對到 main 呼叫與既有 stage/helper 的交界。Runtime、adapter、executor
內部仍依各自報告及下方待簽核項追蹤，不以此表替代其完整內部呼叫圖。

| 入口 | main 控制流程與分支一致性 | 既有驗收依據 |
|---|---|---|
| CPU flow | `flow arguments` 建立含 Newton／步數／輸出／restart／容許值的 controls → exact agreement → PETSc options → database fingerprint／partition → asset catalog／waveforms → boundary input → collective runtime constructor。建構後 VCA／outlet／traction face loops 使用已同意的配置與全域 face counts；沒有本地 face 數決定是否進入下一次 reduction。Initial/restart → BeginStep／SolveTrial／CommitStep → optional VCA transport／ports／budget／circuit → output／checkpoint → final output／Close／summary／profile | [輸入](HPC_01C_FLOW_INPUT_PROGRESS.md)、[步進](HPC_01C_FLOW_STEP_PROGRESS.md)、[輸出](HPC_01C_FLOW_OUTPUT_PROGRESS.md)、[stdout](HPC_01C_SOLVER_STDOUT_PROGRESS.md)、[cleanup](HPC_01C_RUNTIME_CLEANUP_PROGRESS.md) |
| Configured transport | tracking 的初始化前例外保存在 exception_ptr，初始化後共同 rethrow。Controls 包含 memory／velocity／output／checkpoint 分支；database 與實際 fields 分別驗證容量。Assets／waveforms／snapshot manifest 先同意；每步 selection → 選中 snapshot 比對 → interpolation → local assembly → Mat/Vec assembly／KSP → convergence stage → swap → output／checkpoint。記憶體 Record 自帶 measurement／gather／output 協調。Final norm／PETSc destroys → memory Close → summary／profile | [CLI](HPC_01C_TRANSPORT_CLI_PROGRESS.md)、[VTKHDF 初始化](HPC_01C_TRANSPORT_VISUALIZATION_INIT_PROGRESS.md)、[stdout](HPC_01C_SOLVER_STDOUT_PROGRESS.md)、[memory close](HPC_01C_MEMORY_CLOSE_PROGRESS.md) |
| Native 1D | arguments／controls → effective options → configuration text agreement → network／forcing assets → 本地 runtime 建構；callback 只儲存，不在 local constructor stage 執行 implicit solve。`--check` 在共同 output 後早退，尚未建立 PETSc implicit owners。普通分支 initial／restart candidate／validation → writer setup → BeginStep 本地準備 → SolveTrial 內部協調 → PrepareCommit 本地準備 → noexcept FinalizeCommit → diagnostics／output／checkpoint。`OneDFingerprint(config_text)` 只有 uint64 hash 運算，不配置字串 | [CLI](HPC_01C_ONE_D_CLI_PROGRESS.md)、[checkpoint](HPC_01C_ONE_D_CHECKPOINT_PROGRESS.md)、[stdout](HPC_01C_COUPLING_STDOUT_PROGRESS.md) |
| Mesh check | resource preflight → database regular-file／fingerprint／partition／owner checks → fingerprint agreement → owned elements input → geometry reduction → checked result logging → global exit status。空 owned list 可進入 geometry reduction；quality failure 返回 2，輸入或 logging failure 返回 1 | [工具 assets](HPC_01C_TOOL_ASSET_PROGRESS.md)、[工具 stdout](HPC_01C_SERIAL_TOOL_PROGRESS.md) |
| Assembly smoke | resource／returning handler → input／field count／database → field／database／PETSc options agreement → OwnedRowAssembler → CreateMatrix → local insertion → collective assembly → returned MatMissingDiagonal／MatGetInfo → checked destroy → checked summary／global status。沒有把 Mat assembly 包在 local callback | [工具](HPC_01CD_TOOLS_PROGRESS.md)、[工具 assets](HPC_01C_TOOL_ASSET_PROGRESS.md)、[工具 stdout](HPC_01C_SERIAL_TOOL_PROGRESS.md) |

上述歷史報告中的失敗觀察仍保留，尤其 legacy flow 的近零壓力跨 rank relative
gate 未通過，不能由本次入口稽核改寫為通過。配置與資產的 exact agreement
依既有「執行期間內容不變」契約，不提供檔案鎖。MPI 初始化／內部失聯不在
程序存活的錯誤協調保證內。

## 此批已逐項核對的終止呼叫鏈

| 入口／owner | 建構與正常返回 | 例外退棧／清理 | 本批結論 |
|---|---|---|---|
| CPU flow main | `OwnedRowAssembler` → stack `TransientFlowRuntime` → 可選 VCA transport → step／final gather／VTKHDF close → transport Close → flow Close → completion／profile → PetscFinalize | runtime constructor 的 catch 與 destructor 使用同一固定順序；明確 Close 後不重複 destroy | 新增 terminal cleanup 協調；CLI、VCA 與 runtime unit 驗收見 cleanup 報告 |
| Graph runner | 本地準備 native owners → `AllocateCollectiveRuntime` → borrowed adapters／registry → executor accepted histories → 各 3D transport／flow Close → CSV／completion manifest → return → main PetscFinalize | registry owns adapters；貼體 adapters 不在 destructor 呼叫被借用 runtime；native owners 析構在 communicator 釋放前 | 新增 Close，且置於 completion manifest 前 |
| Sequential runner | 1D callbacks／3D runtime → explicit 或 fixed／Aitken loop → accepted histories → 3D Close → CSV／manifest → return → main PetscFinalize | callback outcomes 與內部 runtime stages 分開；strong loop 的 abort 順序沿用既有協議 | 新增 Close，且置於 completion manifest 前 |
| Configured transport main | runtime 物件組直接持有 Mat／Vec／KSP；最後 gather／writer close／VecNorm 後逐一 CheckPetsc destroy，memory report 明確 Close 後才列印 summary | `TransportPetscObjects` 於共同失敗後作 noexcept 後備清理；memory report destructor 只作本地後備 | memory close 的 12 個原生作業及 world3／split1+2 通過 |
| Legacy transport main | `RunLegacyTransport` → output gather objects Close → transport objects Close → summary → return → PetscFinalize | 已有 owner 後備清理，main 協調最終 status | 本批未修改；保留既有 legacy 故障證據 |
| Mesh check | 本地 database／owner／element checks → `InspectGeometry` reductions → result／global status → PetscFinalize | 無 PETSc Mat／Vec owner；database／elements 為本地 C++ 資源 | 空 rank 合法；quality 不佳返回 2，輸入例外返回 1 |
| Assembly smoke | assembler → Mat owner → assembly／info → checked MatDestroy → summary／global status → PetscFinalize | Mat owner 在共同失敗後作後備 destroy | 已有 success 前明確 MatDestroy |
| Native 1D | 本地 configuration／network runtime，只保存 implicit callback；每次 implicit advance 使用局部 PETSc owner 並 Close | `--check` 在任何 implicit advance 前共同印摘要及 finalize；此時 runtime 沒有持有 PETSc solver／vector | `--check` 早退沒有在 finalize 後析構 PETSc 物件；四種 implicit owner Close 沿用既有證據 |

## 共用 owner 與錯誤協議

- `CollectiveLocalStage` 在本地 callback 返回或例外後，用固定大小 stack diagnostic
  做 Allreduce／Bcast；成功者與失敗者完成同一協議後才返回或拋錯。 callback
  不得包含可能與 peer 次序不同的 PETSc／MPI collective。
- `AllocateCollectiveRuntime` 在進入 collective constructor 前協調 host allocation；
  construction input 的配置與複製使用 `PrepareRuntimeConstructionInput`／
  `RuntimeConstructionStage`。建構失敗不得在某一 rank 單獨擁有 distributed owner。
- `OwnedRowAssembler` 持有本地資料容器；成功建立的 Mat／Vec 交給 caller。
  Flow／transport runtime 的失敗建構清理仍由 retained PETSc reference 測試驗證。
- `DomainRuntimeRegistry` 完成本地 index validation 並同步 outcome 後，才以
  noexcept move 接管 adapters。registry 不替代 adapters 所借用 solver 的內部協調。
- `PetscGatherObjects`、checkpoint candidate／local fd、1D PETSc owner 的成功
  路徑已有明確 Close／destroy；destructor 是例外路徑的後備，不作成功標記。

本批新增的 [終止清理契約](../architecture/RUNTIME_CLEANUP.md) 明確列出
各物件順序、重複 Close、程序存活假設與 PETSc 內部故障限制。

## F06 最終簽核

1. 八個 MPI production 入口的參數早退、configuration／external asset、時間步、
   runtime／adapter／executor、writer、checkpoint、Close 與 completion 已逐項對到
   共同階段及故障證據。
2. 93 個含 MPI／PETSc API 的 `CollectiveLocalStage` callback 已掃描；沒有在
   rank-local callback 內發現 MPI／PETSc collective 或 distributed destroy。
   七種主要 owner 的明確 Close、固定釋放順序與 destructor 後備亦已核對。
3. 目前 HEAD 的 world3、split1+2 protocol／executor／construction／cleanup tests
   及 112 個原生 CLI 作業通過；完整數量、來源 hash 與限制見最終簽核報告。

CUDA／FSI exporter 與前處理仍保留各自的單程序限制；本批不宣稱這些入口
支援跨 rank runtime。HPC-01D 能力矩陣及 HPC-03～09 的必要驗收保持原範圍。
