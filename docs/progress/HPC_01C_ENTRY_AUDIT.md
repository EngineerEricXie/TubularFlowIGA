# HPC-01C：入口與清理稽核

日期：2026-09-09。基準 `d9eaea9` 加本批 runtime cleanup 修改。
狀態：F06 進行中；此表記錄已核對的呼叫鏈，未把歷史測試自動視為本版實測。

## 此批已逐項核對的終止呼叫鏈

| 入口／owner | 建構與正常返回 | 例外退棧／清理 | 本批結論 |
|---|---|---|---|
| CPU flow main | `OwnedRowAssembler` → stack `TransientFlowRuntime` → 可選 VCA transport → step／final gather／VTKHDF close → transport Close → flow Close → completion／profile → PetscFinalize | runtime constructor 的 catch 與 destructor 使用同一固定順序；明確 Close 後不重複 destroy | 新增 terminal cleanup 協調；CLI、VCA 與 runtime unit 驗收見 cleanup 報告 |
| Graph runner | 本地準備 native owners → `AllocateCollectiveRuntime` → borrowed adapters／registry → executor accepted histories → 各 3D transport／flow Close → CSV／completion manifest → return → main PetscFinalize | registry owns adapters；貼體 adapters 不在 destructor 呼叫被借用 runtime；native owners 析構在 communicator 釋放前 | 新增 Close，且置於 completion manifest 前 |
| Sequential runner | 1D callbacks／3D runtime → explicit 或 fixed／Aitken loop → accepted histories → 3D Close → CSV／manifest → return → main PetscFinalize | callback outcomes 與內部 runtime stages 分開；strong loop 的 abort 順序沿用既有協議 | 新增 Close，且置於 completion manifest 前 |
| Configured transport main | runtime 物件組直接持有 Mat／Vec／KSP；最後 gather／writer close／VecNorm 後逐一 CheckPetsc destroy，再列印 summary | `TransportPetscObjects` 於共同失敗後作 noexcept 後備清理 | 已有明確成功路徑清理；本批未修改該 CLI |
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

## 尚未簽核的 F06 範圍

1. 把各入口的參數早退、configuration／external asset、時間步控制、
   runtime／adapter／executor 與 writer 呼叫，逐項對到現有 stage 及故障證據；
   本表只完成終止路徑，不以「已看過 main」替代完整呼叫鏈。
2. 核對共用 helper 的 error-handler stack／cleanup fallback 與呼叫端的假設，
   區分可控制的本地／返回錯誤和 MPI／PETSc 內部程序故障；目前不同 helper
   的 Close 在第一個共同失敗後由 destructor 清理剩餘物件，與新 runtime 的
   全部釋放後協調不同，不能混稱同一實作。
3. 依上述核對結果補最後的必要原生／subcommunicator 整合驗收；既有 F01～F05
   已驗證項目可引用精確 revision／binary 證據，不必機械式重跑所有歷史測試。

CUDA／FSI exporter 與前處理仍保留各自的單程序限制；本批不宣稱這些入口
支援跨 rank runtime。HPC-01D 能力矩陣及 HPC-03～09 的必要驗收保持原範圍。
