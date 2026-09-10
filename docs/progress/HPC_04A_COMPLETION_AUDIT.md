# HPC-04A 可配置求解器完成稽核

狀態：HPC-04A 已完成；HPC-04B／C 與整份工作站／HPC 清單仍進行中。
日期：2026-09-10。稽核 repository revision：`d3bb714`。

本項要求是各 PETSc runtime 的獨立選項、既有預設相容性、有效 solver 診斷與
多 rank LU 後端檢查。以下逐項連結實作及實際驗收，不以跨節點 scaling 或
分散式 FSI 的尚未完成代替本項判定；那些要求仍保留在其原子任務。

| 要求 | 實作／證據 | 判定 |
|---|---|---|
| 各 runtime 獨立 prefix | Graph flow／transport／1D 依 domain 與 role；immersed static／transient family 與 domain override；moving epoch 共用 immutable snapshot、FSI wrapper 按 fluid domain ID；standalone／legacy CLI 保留 embedding 相容介面。見 [操作契約](../SOLVER_OPTIONS.md)、[貼體](HPC_04A_BODY_FITTED_OPTIONS_PROGRESS.md)、[1D](HPC_04A_ONE_D_OPTIONS_PROGRESS.md)、[immersed](HPC_04A_IMMERSED_OPTIONS_PROGRESS.md)。 | 通過 |
| 保留已驗證預設 | Default／reference bytes、override field／history、checkpoint、健康重試與 family precedence 均有各自驗證；貼體 132 作業、standalone 27 作業、immersed 最終 options 集合 81 作業。原 unfitted rigid 失敗保留，另有 fitting 修復後完整回歸通過。 | 通過 |
| 記錄有效 KSP／PC 與迭代 | Accepted-step solver_configuration／one_d_solver_configuration，區分最後一次 KSP 與累計迭代；SNES type／reason 與 child viewer、GAMG smoother／coarse prefix 均驗證。見 [診斷](HPC_04A_BACKEND_DIAGNOSTICS_PROGRESS.md)、[多層](HPC_04A_MULTILEVEL_OPTIONS_PROGRESS.md)。 | 通過 |
| 多 rank LU 能力，包括 1D nonlinear MUMPS | 1／2／3／4-rank、48 個 backend／PC 組合；可用後端核對正 reason／數值，不可用者拒絕。1D 四 formulation、SNES 與 native CLI 驗證，單 rank 不誤拒 PETSc LU；invalid backend 不印成功。 | 通過 |
| 共同 options 層覆蓋 | 搜尋 production CPU／1D 的 KSPCreate／SNESCreate 與 SetFromOptions，共 10 處配置呼叫均經 options Call；private database／returning handler／群組一致性及失敗 cleanup 有測試。 | 通過 |

## 來源與證據完整性

`outputs/hpc04/completion-audit-v1/` 保存兩次核對：

- `evidence.json`：貼體、1D、standalone 原 manifests 列出的 10867／8601／539
  個檔案，合計 20007 個 SHA256 重算一致，無遺失或修改。
- `acceptance.json`／`audit.py`：較晚 20 組 suite、211 份 rank reports 的
  stdout／stderr hashes、timeout／退出狀態、預期負向作業、後端能力與多層
  數值門檻、完整 FSI／moving 結果及 manifest hashes。與前述 suite 可能重疊，
  不相加宣稱獨立作業數。

七個核心 options／runtime 檔案與 `01aa4c8` 逐位元一致：PetscSolverOptions、
TransientFlowRuntime、TransientTransportRuntime、ImmersedStaticFlowRuntime、
ImmersedDistributedNewtonRuntime、ImmersedTransientFlowRuntime、OneDImplicit。
多層測試來源與 binary 亦與其通過 receipt 相同。後續 CLI 記憶體／PVTU 修改由
其實際 PDE、checkpoint、writer 故障與 ParaView 驗證覆蓋；沒有更動這些核心選項。

完整 moving frozen `c4fee00` 的來源另核對 162 檔，rigid gap 5.52939e-17，原
門檻 5.68434e-14；該結果不改寫成舊 `01aa4c8` 的通過紀錄，也不宣稱所有
executables 都在目前 HEAD 重跑。後續 fitting 計數／identity 修正保留 focused
驗證與獨立 provenance，未改本項 solver prefix／snapshot 契約。

本項性能／RSS 觀測分別保存在各連結報告；本次重算 hash 與來源稽核本身沒有
新增求解效能數字。環境仍為工作站 PETSc 3.15.5 real64/int32、Open MPI 4.1.2、
GCC 11.4，CUDA allocation 對本 PETSc options 稽核不適用。

## 保留的限制與下一步

HPC-04B／C 仍須其他網格／rank／pressure／immersed 條件與無干擾成本比較。
同 source revision 的跨節點 allocation 屬 HPC-09，分散式 FSI 屬 HPC-07；
本項沒有用單機 MPI 取代其驗收。既有預設求解策略不因候選測試通過而變更。
