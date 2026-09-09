# 耦合 checkpoint 狀態契約

狀態：HPC-05A 已完成的介面與狀態稽核；**本文件定義待實作的契約，並不表示
graph 已可 checkpoint／restart**。程式基準為
`70b1715992912871eebb4a6458634bedb1b9a542`，2026-09-09。
實作與驗收分別由 [HPC-05B／C／D](../WORKSTATION_HPC_TODO.md) 承接。
更新：05B 已交付 [bundle v1 發布／載入元件](COUPLED_CHECKPOINT_BUNDLE.md)；
以下 domain payload／restore 與完整 graph 續跑仍由 05C／D 實作。
05C 已補上 [0D 與 pressure／donor state 的首批介面](../progress/HPC_05C_ACCEPTED_STATE_PROGRESS.md)；
下文的程式缺口盤點以開頭的原始基準 revision 為準，最新交付見該進度報告。

## 1. 唯一可保存的邊界

第一版只保存所有 domain **完成 accepted macro-step** 的狀態。設已接受步數為
`N >= 1`，必須依序完成：所有 trial 成功、守恆檢查與 before-commit callback、
所有 domain prepare、所有 domain finalize、executor donor 發布、runner 的
accepted-result／下一步壓力猜值更新，以及 accepted clock 更新，才可擷取 checkpoint。
這是 [native graph runner](../../solvers/coupling/src/iga_1d_3d_bifurcation.cpp)
每次迴圈尾端、下一個 `BeginStep` 前的新同步點。

不在 Newton iteration、hydraulic／transport 兩階段之間、FSI iteration、prepare
與 finalize 之間或 before-commit callback 中寫 checkpoint。後者可能在
`CollectiveLocalStage` 的 local lambda 內執行，不能加入 MPI collective。
FSI 必須等流體、結構與兩側 surface publication 全部完成配對 finalize。
第一版不保存尚未接受的 trial，也不增加 multirate。

所有參與 rank 必須同意 domain catalog、epoch、accepted clock 及 idle 狀態，
才可依相同順序執行擷取和發布。signal／排程時間預警只設定待保存旗標；正常控制流程
在下一個上述邊界處理旗標。signal handler 不執行 PETSc、MPI 或檔案 I/O。
若下一步失敗或作業被終止，新的作業只能恢復最後已完整發布的 epoch。

## 2. 時鐘、身分與資料表示

| 欄位 | 契約與驗證 |
|---|---|
| `accepted_macro_steps` | 整數 `N`；所有整數轉換檢查範圍與溢位；不由時間反推 |
| `accepted_time_s` | 保存 runner 實際逐步 `EndTime()` 累加的 binary64 值；不以 `N*dt` 重算 |
| `macro_dt_s`、設定終止步數 | 與配置一致；第一版固定原有時間表，不允許載入時靜默改 dt 或總步數 |
| 下一個 graph context | `step_index=N`、`start_time_s=accepted_time_s`；第一個 graph step 的 index 為 0 |
| 0D clock | `CommittedStepIndex=N-1`、`CommittedStepCount=N`、`CommittedTime=accepted_time_s` |
| 浸入式暫態 clock | backend `index=N`、`time_s=accepted_time_s`；trial inactive；不能用 graph 的前一步 index 取代 |
| 1D clock | 各 domain 保存自己的 `completed_step`、`physical_time`、`internal_substeps`；依原有整數 subcycling 驗證與 macro-step 對應，不能一律設成 `N` |
| 3D transport clock | 保存 `Steps()` 真值；每個 macro-step 一次 transport 的現有 graph 為 `N`，未來多 system 必須各自驗證 |
| epoch 身分 | 所有 metadata、shard、歷史片段均屬於同一唯一 generation，並記錄前一個完成 generation；不能只靠相同時間或檔名判斷 |
| 物理／配置身分 | graph domain／edge／port IDs、方向、物種／field 順序、單位、模型參數、coupling method／tolerance／routing controls、所有外部輸入內容身分 |
| 分布身分 | domain communicator 成員的穩定 rank 映射、rank 數、owned row／node／cell／surface ID 覆蓋與 layout／partition 身分；不保存 MPI handle |
| 執行相容性 | payload schema、runtime kind／版本、real／complex、scalar 與 index 寬度、byte order、PETSc 與數值選項及程式 revision；第一版不承諾不同 backend／build 的自動轉換 |

外部輸入包括 `.ntiga`、分區、1D network、表面、FSI reference／patch map、
boundary waveform、replay CSV、速度場來源與其時間表。使用內容 SHA-256 和位元組數，
路徑只作定位資訊；檔名相同不代表相容。原有 FNV fingerprint 可以留作 legacy 載入
檢查，但不能取代新 bundle 的完整 payload checksum。有效 PETSc options 也會影響
續跑，必須納入相容性記錄；輸出目錄／checkpoint 頻率等不影響數值的設定另行記錄。

浮點值須有限且可往返 binary64；文字使用足夠有效位數，二進位明定編碼。
optional 欄位保存 presence，空物種 map 與缺失 map 不混用；不把 absent pressure
改成零。所有 map 以穩定 ID 為 key，所有序列明定順序及長度；載入先檢查重複、缺漏、
非法 enum、非有限值、長度與上限，再配置大型記憶體。SHA 校驗後仍須做物理及身分驗證。

## 3. 跨步耦合狀態：明確保存／重建決策

| 狀態 | 決策與程式依據 |
|---|---|
| 下一步壓力猜值 | **保存** runner 的 `pressure[edge_id]`，它在 accepted result 後更新為最後 iteration 的 `measured_pressure_pa`。不能改用原始 `initial_pressure_pa` 或最後 applied pressure |
| Pressure-flow Aitken | **重建**。[PressureFlowComponentExecutor::Advance](../../include/PressureFlowComponentExecutor.hpp) 每次呼叫建立 local `previous_residual`，並以 controls 的 relaxation 初始化 omega；沒有跨 macro-step 的殘差向量要載入 |
| Species hydraulic Aitken | **重建**。[SpeciesPressureFlowComponentExecutor::SolveHydraulics](../../include/SpeciesPressureFlowComponentExecutor.hpp) 同樣每一步建立 local residual／omega；保存 method 與所有 controls，下一步按原演算法初始化 |
| Species donor hysteresis | **保存** executor 的 `CommittedDonorOwnership()`，key 為 `(edge_id, species_id)`，value 為 `First` 或 `Second`。它只在所有 domain finalize 後 swap 發布；近零流量由 [ResolveSpeciesDonor](../../include/SpeciesCoupling.hpp) 讀取上一個 donor，缺值會失敗；不可由零流量重新猜方向 |
| Species routes／transport order | **重建**：下一個 accepted hydraulic trial 與已恢復 donor 決定新 routes／order；保留 graph 邊的 first／second 順序。載入 donor 必須覆蓋所有有 species 的 edge pairs，不接受未知 key |
| Strong FSI Aitken | **重建**。[StrongFluidStructureCoupling::Execute](../../solvers/cpu/include/StrongFluidStructureCoupling.hpp) 每次開頭 `aitken_.Reset()`。下一步 predictor 為已提交 scalar displacement 加 `dt*velocity`，套用既有 clamp 與 immutable reference normals；故這些 state／model 都必須恢復 |
| Coupling iteration 歷史 | **保存 accepted output prefix**，包括每步迭代次數、殘差、applied／measured pressure、relaxation、ports、species amounts／donors、0D 帳目與 FSI 接受身分。它不是下一步的 Aitken warm state，但屬於完整歷史驗收 |

最後一列採用版本化的 accepted records／不可變歷史片段，manifest 指定其範圍、
筆數和 checksum。新作業輸出由 prefix 加新 accepted records 組成，每步恰好一次；
不得依賴原作業仍存活的 vector、未驗證的舊 CSV，或只輸出 restart 之後的尾段。
時間、ID、有效性旗標與已接受紀錄必須精確保留；後續新求解數值按既定容許值比較。
wall time、MPI 時間與 RSS 是每次執行的量測，分開記錄，不能偽造為連續作業的時間。

## 4. 各 runtime 的 accepted payload

下表的「保存」可由 shard 直接承載，或引用 bundle 內有 checksum 的不可變資產。
任何必要種類缺少 exporter／loader 時，該 graph 的 checkpoint 必須明確拒絕；
不可以略過該 domain 後發布一份看似完整的 manifest。

### 4.1 0D 與 outlet

[ZeroDFlowDomain](../../include/ZeroDFlowDomain.hpp) 保存 `stored_pressure_pa`、
上述三個 committed clock／counter、optional committed `PortState` 與
`ZeroDFlowStepAccounting`。帳目包含初／末儲存量、source、distal sink、graph-port
積分量與 residual；不是只保存 pressure。模型、role、port 身分及 model digest
須一致。`base_state_`、input、trial、staged records 在 idle 恢復為 inactive。
目前 constructor 的 initial pressure／time 不能恢復 index、count 與 publication，
05C 需加入完整、驗證後的恢復入口。

[OutletCheckpoint](../../include/OutletCheckpoint.hpp) 已保存 outlet label／kind、
flow、pressure、capacitor pressure。新 payload 延續這些狀態，並驗證全部 resistance、
capacitance、reference／distal pressure 與 boundary mapping 的配置身分。
不得重設 RCR capacitor pressure，也不能把 graph port pressure 猜值當成其替代品。

### 4.2 1D hydraulic 與 transport

[OneDRuntime](../../solvers/one_d/include/OneDRuntime.hpp) 的 Ready phase 是 accepted
邊界；`CloseTrialState()` 已清除 hydraulic frames 與 rollback snapshot。

- 保存完整 [OneDFlowState](../../solvers/one_d/include/OneDFlow.hpp)：cell area、flow、
  pressure、node pressure、segment flow、各 outlet 動態值、inlet flow、completed
  step、physical time、internal substeps。node／segment／cell 順序與 network 綁定。
- 保存各 [OneDSpeciesState](../../solvers/one_d/include/OneDTransport.hpp) 的 concentration、
  root／outlet native flux、`boundary_flux_valid`、step initial mass、root／outlet／source
  amounts 與 `step_accounting_valid`。另保存動態 `inlet_value`、`inlet_waveform`：
  `ApplyOneDCoupledInlet` 會更新入口值並清除該物種的 waveform，後續入口省略該物種時
  仍沿用這些變動。definition、wall kind／value／coefficient／exterior value 由驗證過的
  配置恢復。保存完整 `LastInlet()`，包含 optional
  flow／pressure、species、temperature、hematocrit 與 metadata；Ready-phase port 查詢
  會讀取這些 accepted 邊界資料。
- 動態 vasodilation 改變 network segment 的 `radius0`、`area0`、`resistance`；保存
  accepted radius，按既有公式由 radius、length、viscosity 重建後兩者並驗證。
  `baseline_radius0` 屬於不可變模型。`ApplyOneDVasodilation` 的 configuration 雖為 const，
  `ApplyOneDCoupledInlet` 會改寫 physiology 與 coupling perfusate oxygen 的
  hematocrit／hemoglobin（共四個 scalar）；這些必須保存。後續入口可以不帶 Hct，
  僅保存 `LastInlet()` 不足以重建先前的血液狀態。不必序列化整份 rollback snapshot。
- 既有 integer configured substeps 與 adaptive internal substeps 保留原語義。
  hydraulic frames、trial inlet、trial outlet overrides、transport routing／rollback
  cache 在下一次 Begin／Solve 重建，不保存未接受的 substep。

[OneDCheckpoint v2](../../solvers/one_d/include/OneDCheckpoint.hpp) 已涵蓋 hydraulic
陣列、outlet 三個動態值、radius 與濃度，但沒有完整 `LastInlet()`／物種通量帳目。
`RestoreCommittedState(flow, transports, network)` 也未恢復 `last_inlet_`。
05C 必須補齊這些欄位及驗證，不能把 standalone checkpoint 直接包成 graph checkpoint。
staged species adapter 目前拒絕 dynamic vasodilation；restart 不擴大它的物理支援範圍。

05C 已補上 `OneDFlowCheckpointState`、Ready accepted capture／fresh-candidate restore，
以及 [metadata／field codec](../../solvers/one_d/include/OneDAcceptedCheckpoint.hpp)。
只支援本契約的固定 macro dt：保存 accepted macro count、最後 start／dt，並核對
configured step count 的整數倍關係。provider 必須傳入涵蓋完整配置、selected system、
network 與外部輸入的已驗證 SHA-256 identity；省略 identity 的 legacy owner 不能使用
新 checkpoint API。metadata 與 fields 必須同屬一個完整 bundle，通過所有分片校驗、
consumer Finish、runtime validation 與群組 agreement 後才發布候選 owners。
詳見 [1D 狀態與串流驗收](../progress/HPC_05C_ONE_D_STATE_PROGRESS.md)。

### 4.3 貼體 3D hydraulic 與 transport

[TransientFlowRuntime](../../solvers/cpu/include/TransientFlowRuntime.hpp) 在 Committed
phase 保存 `state_` 的全部 `(u,v,w,p)` owned rows、outlet 動態值、accepted resolved
boundaries／pressure tractions，以及累積 linear iteration count。保存 graph clock，
不能從 `trial_step_` 的 stale storage 判斷 accepted index。配置的 waveform 及 port
mapping 仍須驗證，下一步 trial controls 由 adapter 按新時間重新 materialize。

Backward Euler 的下一次 `BeginStep` 會把 `state_` 複製到 `committed_state_` 與
`previous_`。因此 accepted 邊界只需要一份當前完整場；前一 trial 殘留的 `previous_`
不是額外的跨步積分歷史。載入後以 accepted field 初始化必要 snapshot／ghost，
不沿用新建 runtime 的零場。steady graph 也保存 accepted field 作 Newton warm start。

[TransientTransportRuntime](../../solvers/cpu/include/TransientTransportRuntime.hpp)
保存所有 ordered fields 的 `current_` 及真實 `steps_`。`previous_` 在這個 class 是
**Mat**，每次 SolveTrial 重新組裝，不能誤認為需要另一份 concentration history。
`committed_` 由當前場重建，next／rhs／forcing／KSP／matrix／scatter 重建。
現有 `ReadState()` 只把 `steps_` 設為 1；它足以打開原有 warm-start 分支，仍不符合
graph clock／完整診斷契約，05C 需提供帶正確 counter 的 validated restore。

[FlowCheckpoint v1](../../include/FlowCheckpoint.hpp) 與
[TransportCheckpoint v1](../../include/TransportCheckpoint.hpp) 的 PETSc vector
可以是新 shard 的承載方式，但舊 metadata＋state 兩檔不具備跨 domain 原子發布。
相同 rank 數也須核對 `.ntiga` 分區、node IDs、field 順序與 ownership，不能只比 vector 長度。

05C 已加入 flow／transport 的 typed accepted capture／fresh candidate restore：
flow 的 graph macro dt 與 steady kernel dt=0 分開保存／驗證；transport 使用真實
accepted step count，舊 ReadState 的 steps=1 行為只留在 legacy API。場按 owned rows
保存，replicated boundaries／outlets 經 SHA agreement；clock、iteration count 與
snapshot／BE history 一起恢復。詳見 [貼體 3D 狀態進度](../progress/HPC_05C_BODY_FITTED_STATE_PROGRESS.md)。
3D 磁碟 codec 與 local bundle producer／loader 已加入：metadata、replicated boundary
arrays 每 domain 保存一次，owned fields 按 bundle world rank 分片，stream buffer 上限
64 KiB。decoder 的形狀取自 verified fresh runtime；全部 payload／SHA 與群組 agreement
通過後才 typed restore。三種模式已由新 MPI 作業恢復，詳見
[3D 分片進度](../progress/HPC_05C_BODY_FITTED_BUNDLE_PROGRESS.md)。全作業 provider
publication 尚未接入 native graph，不能把 typed restore 當作 live graph 的跨 domain 原子替換。

### 4.4 VCA／外部 circuit

[VcaExternalCircuit](../../include/VascularCoupling.hpp) 保存 reservoir volume、所有
species concentrations、temperature、hematocrit，以及 `last_arterial_species_`。
後者由 `InletState()` 設定，並被 `Advance()` 的 device-source report 使用；完整恢復
選擇直接保存，不能用更新後 reservoir 重算上一個 arterial observation。
同一 accepted 邊界的 flow、transport、outlet、reservoir 與 circuit history 必須一起發布。
原有 pump／oxygenator／dialyzer／infusion 模型參數及 oxygen state 定義屬於配置身分；
目前沒有額外裝置積分器狀態。

[VcaCheckpoint v2](../../include/VcaCheckpoint.hpp) 已有 reservoir 與 transport state
連結及 identity，沒有上述完整 circuit cache／history，也沒有多檔完成協議。
[CPU VCA runner](../../solvers/cpu/src/iga_navier_stokes.cpp) 會從恢復濃度重算
`vca_previous_mass` 與 required-node species；這些可重建，歷史 report 則必須保存。
[1D CLI](../../solvers/one_d/src/iga_1d.cpp) 目前明確拒絕 closed-loop checkpoint／restart。
此契約不把 VCA 當成已有 graph adapter；納入前須補齊相應 transaction／restore 實作。
replay／open-loop 的未來輸入以保存時間及經 checksum 驗證的時間序列重讀，不能依賴舊 file cursor。

### 4.5 固定幾何浸入式

[ImmersedDistributedNewtonRuntime](../../solvers/cpu/include/ImmersedDistributedNewtonRuntime.hpp)
保存 owned accepted vector、layout／partition、controller／gauge 身分、port control
values，以及邏輯 commit count。其他 begin／rollback／prepare 等執行診斷作為 accepted
audit prefix 保存；新作業的計時與執行計數另列。Newton candidate、update、prepared
vector、ghost、matrix／KSP 不保存，載入 idle runtime 後重建。

[ImmersedTransientDistributedRuntime](../../solvers/cpu/include/ImmersedTransientDistributedRuntime.hpp)
另存 committed time／index。accepted finalize 會釋放 frozen velocity／body-force inputs；
下一次 BeginTrial 從已恢復場建立 history，按新時間凍結 force。force／motion callback
必須由可驗證且可重建的模型與輸入提供；任意捕捉 process-local mutable state 的 callback
不能聲稱可 restart，除非 producer 額外提供持久化狀態。

現有 `SetCommittedOwnedState` 可設場與 clock，未恢復 graph adapter clocks／accepted
ports／commit count；[transient lifecycle](../../include/ThreeDImmersedTransientDistributedFlowDomain.hpp)
constructor 還要求 fresh index 0。05D 必須一起處理，不能繞過 clock guard。
穩態／暫態共有 adapter 的 committed boundary controls、port publications 和 graph
time／index／count 也屬於 payload；其 trial inputs 與 rollback copies 不屬於 payload。
legacy serial [ThreeDImmersedFlowDomain](../../include/ThreeDImmersedFlowDomain.hpp)
同樣保存 committed controls／measurements／clock／counter 與 accepted field，不能只
替它包裝分散式 adapter 的 schema 而遺漏序列 runtime 本身的狀態。

### 4.6 移動幾何與暫態場

[MovingImmersedTransientFlowRuntime](../../solvers/cpu/include/MovingImmersedTransientFlowRuntime.hpp)
在 idle 保存當前 epoch 的完整
[ImmersedGlobalFlowState](../../solvers/cpu/include/ImmersedTransientState.hpp)：time、index、
geometry identity、stable node IDs、四個場係數、port IDs／multipliers、gauge presence／value。
另存 port controls 與全部 `CommittedConservation()` transition primitives／身分／scales；
這些 accepted 舊→新量測不能只從當前場反推。

幾何 payload 保存 [MaterialSurfaceKinematics](../../solvers/cpu/include/MaterialSurfaceKinematics.hpp)
的 reference material vertices、當前 source vertices／velocities、source triangle
IDs／ordered vertices／labels、canonical-corner provenance、evaluated time、step start／end、
material／topology／content／geometry-epoch identities，並區分 prescribed producer 的
歷史 epoch hash 與 producer-neutral content hash。只保存 canonical 表面或 VTP 不足夠。
prescribed motion 的原始 reference 與 motion model 也必須驗證，以便計算下一個 epoch。

由完整 payload 與 grid／cut／quadrature／ghost／extension options 重建當前 geometry、
layout 與積分 catalogs，重算並核對 geometry identity。下一次 trial 的 velocity／pressure
extension、seed、force snapshot 與 trial map 從當前 accepted epoch 重建；不必保存舊
trial PETSc objects 或所有歷史 epoch 的完整 mesh。

**Publication identity 有額外來源依賴**：
[MovingCutGeometry::HashPublicationState](../../solvers/cpu/include/MovingCutGeometry.hpp)
包含當前 geometry digest、上一個 geometry digest、unchanged-cell count、transition
counts 與逐 cell old／new classifications。須保存這份有限的 predecessor provenance
及 publication digest。05D 的 verified restore builder 必須檢查 cell ID／分類、完整計數、
當前 geometry 一致性並重算 digest；不能用一般 constructor 的 genesis publication
冒充原 publication，也不能提供未驗證的任意 hash override。
目前 `SetCommittedGlobalState()` 會清空 conservation，沒有完整 publication restore，
故尚不足以實作上述續跑。

### 4.7 結構與 FSI publication

[PretensionedMembrane](../../solvers/cpu/include/PretensionedMembrane.hpp) 保存按
owned global surface-node IDs 排列的 scalar displacement 與 velocity。reference mesh、
reference normals／areas、clamped mask、material constants、layout／partition 與 model
identity 須可驗證重建。BE 結構模型下一步使用這兩個 committed fields；沒有另一個
acceleration history。dense factorization／matrix 及 trial owner token／generation capability
在新 runtime 重建，不能反序列化成可使用的舊 trial handle。

[PretensionedMembraneFsiRuntime](../../solvers/cpu/include/PretensionedMembraneFsiRuntime.hpp)
還須保存 committed kinematics 與完整 stamp。
[MovingImmersedTransientFlowFsiRuntime](../../solvers/cpu/include/MovingImmersedTransientFlowFsiRuntime.hpp)
保存上述移動流場、`committed_full_` material kinematics、committed surface traction、
traction diagnostics、composition identity、patch map 身分與完整 surface stamp。
完整 stamp 含時間、步數、accepted coupling iteration、reference／layout／partition
及 producer-state identities；不可為了重啟而重設 iteration 或換成新 trial 的 stamp。

[FsiTrialLifecycle](../../include/FsiDomainRuntime.hpp) 載入時為 Idle，恢復
`has_committed_output_` 與 `committed_output_`；不恢復 expected input、actual／prepared
trial output 或 active context。整個配對 candidate 驗證後一起 publish，不能只恢復
membrane displacement 便讓 fluid 的 surface stamp 停在初始狀態。Strong FSI 的下一步
predictor／Aitken 遵循第 3 節，accepted raw／relaxed／traction identities 留在歷史 prefix。
目前 membrane／strong coordinator 仍限單 partition；分散式恢復在 HPC-07 與 05D 整合。

## 5. 實作邊界與 restore 順序

[CoupledDomainRuntime](../../include/CoupledDomainRuntime.hpp) 現在只有 trial／commit
介面，沒有 checkpoint capability。05B／C 新增的 provider 應明確區分：local metadata
擷取、需要 communicator 的分散式場 I/O、candidate 驗證，以及不丟例外的 publish。
名稱與 wire format 由實作確定；本文件不是已存在的 C++ API 宣告。

1. 讀取完成 manifest，驗證版本、所有檔案 checksum／大小／epoch 與配置相容性。
   不先修改 live runtime 或輸出；05B 定義暫存 shard、sync、完成 manifest 最後發布。
2. 以原 rank 數與 domain membership 建立模型、ownership、geometry、PETSc objects；
   檢查所有穩定 ID 覆蓋，不允許缺片、重複 owned rows 或以零補洞。
3. 在獨立 candidate 中載入全部 domain state、clocks、publications、pressure map、
   donor map 與歷史 prefix。重新計算模型／場／publication identities，驗證有限性、
   正值、clock／counter 關係、donor coverage、FSI 配對與 output prefix 末步。
4. 全 rank 一致確認恢復成功後，才將完整 candidate 設為 accepted；失敗時銷毀
   candidate，保留原 bundle 與輸出。重新建立 trial scratch／ghost，而非執行一個假 step。
5. 第一次續跑從 `N+1` 的人類可讀步號開始，輸出 prefix 與新紀錄不重複、不漏步。

05C 先交付既有 0D／1D／貼體 3D graph 的相同 rank 恢復；05D 分批加入浸入式、移動與
FSI。未支援 runtime／版本／rank 變動要在開始寫 shard 前拒絕。不同 rank 恢復須另做
global-ID ownership 搬移，尤其 `.ntiga` 分區和 surface partition stamp 不可直接沿用。
現有 standalone checkpoint 檔案介面保留；不把舊格式重新標記成新原子 bundle。
CUDA 的 standalone raw-state checkpoint 保留其 backend／format 身分；目前沒有
CUDA graph provider，不能將其位元組當成 PETSc vector 載入。新增 provider 時仍須遵守
相同 accepted field／outlet／species／clock 契約並通過 backend 對照驗收。
[Phase 7 visualization snapshots](../POST_PHASE_7_IO_HARDENING.md) 仍是視覺化格式。

## 6. 後續必須通過的驗收

| 驗收 | 判定與對應任務 |
|---|---|
| 拒絕非 accepted 邊界 | 所有 trial／prepare／部分 domain finalize／species staged／FSI active 狀態拒絕擷取，不發布 manifest（05B–D） |
| 完整歷史續跑 | 獨立 MPI job 載入非零 `N` 後完成；對照不中斷作業的所有 accepted records、0D/outlet state、1D internal counters、3D fields／species；IDs／時鐘／presence 精確一致，後續場使用既定 relative L2 及守恆 gate（05C） |
| Coupling 記憶辨識 | explicit／fixed／Aitken 的兩個以上 accepted steps；pressure initial guess 不能碰巧與 last measured 相同；near-zero species step 接在有確定 donor 的非零 step 後，另含 donor 反轉；比較迭代歷史（05C） |
| History 辨識 | 非零 transient history、RCR capacitor、物種通量帳目及 reservoir composition；FSI 非零 displacement／velocity，moving 有真實 active-cell transition；不能只測零場／stationary 重啟（05C／D） |
| 故障發布與全作業重啟 | commit 前、分片寫入中、manifest 發布前終止；新 job 只恢復最後完整 epoch。截斷、byte corruption、缺片、混 epoch、錯配置／rank／field order 全部拒絕（05B／C） |
| Publication 恢復 | 舊→新 conservation、geometry predecessor provenance、FSI committed stamp 與後續 predictor 重現；禁止新建 genesis identity 冒充（05D） |
| 規模與資源 | rank 1／2／4 與 split communicator；assembly／solve、I/O／sync／通訊時間及每 rank peak RSS 分開記錄；large／cross-node 經 scheduler（05C／D） |

這些是 05B／C／D 的待執行 gates，並非本次 05A 稽核的數值通過紀錄。
05A 的證據、目前缺口與文件檢查見 [驗收報告](../progress/HPC_05A_CHECKPOINT_CONTRACT_REPORT.md)。
