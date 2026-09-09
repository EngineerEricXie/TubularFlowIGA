# 單機與超級計算機開發待辦清單

建立日期：2026-09-07。用途：**分階段開發待辦與進度追蹤，供後續 goal 指定範圍**。

最近更新：2026-09-09。HPC-03C 已補 owned 暫態體積組裝與凍結外力；完整暫態 runtime 仍待接續；
其他既有完成狀態沿用對應報告，未宣稱本批重新驗收整份清單。

快速導覽：[進度總覽](#接續開發的狀態總覽) ·
[優先順序](#建議開發順序) · [goal 範本](#可複製的-goal-範本) ·
[任務依賴](#階段與依賴) · [驗收紀錄範本](#每個子任務的完成紀錄)。

本文件接續已完成的 [Foundational Phases 0–9](progress/FINAL_ROADMAP_REPORT.md)，
目標是讓既有物理模型能在單機有效利用多核心，並逐步支援跨節點 MPI、
長時間執行與可靠續跑。此處的 `HPC-xx` 是新的工作編號，不沿用舊 Phase 編號。

## 目標執行模式

「單機」描述硬體位置，「單程序」描述程序數量，兩者不同；單機也可以使用
多個 MPI rank，單一程序則可以透過 OpenMP 使用多核心。以下是開發目標，
各模組目前的支援程度仍以後文基線及驗收紀錄為準。

| 執行模式 | 用途 | 對應待辦 |
|---|---|---|
| 單機、單 rank、單 thread | 小案例、除錯與數值參考 | HPC-00、HPC-09 |
| 單機、單 rank、多 threads | 加速浸入式／FSI 等局部組裝 | HPC-02 |
| 單機、多 ranks，可搭配 threads | 分散較大問題，比較純 MPI 與混合配置 | HPC-01～04、HPC-06～08 |
| 多節點、多 ranks，可搭配 threads | 超級計算機上的大型、長時間模擬 | HPC-03～09 |
| 單機、單 GPU | 維護既有 CUDA 路徑及 CPU/GPU 數值比較 | HPC-00、HPC-09；多 GPU 另立目標 |

局部組裝使用多核心不代表整個 FSI 求解已平行，也不代表已支援跨節點。
每項驗收需分別記錄數值正確性、實際並行範圍與端到端效能。

## 如何使用本文件開發

可選擇一個尚未完成的子任務，例如 `HPC-01C`，作為單次 goal 的範圍。
若使用者明確指定完成整份清單，則保留全部 38 個子任務的目標，依相依條件
分批實作、測試與記錄，不能以完成其中一項宣告整個 goal 完成。
本文件本身不啟動 goal。

文件保留工作區既有的完成與部分完成紀錄。後續啟動 goal 時，應先核對報告、
目前 revision 與剩餘工作，並沿用使用者已指定的開發範圍；建立或閱讀本清單不代表啟動開發。
基準工具與使用方式見 [HPC_BENCHMARKS.md](HPC_BENCHMARKS.md)。

### 接續開發的狀態總覽

以下依本文件現有核取方塊彙整；本批新增的驗收範圍見相應進度報告。

| 狀態 | 任務 |
|---|---|
| 已勾選完成（11 項） | HPC-00A、HPC-00B、HPC-00C、HPC-00D、HPC-01A、HPC-01B、HPC-02A、HPC-02B、HPC-02C、HPC-03A、HPC-03B |
| 已有部分進度、尚未完成（3 項） | HPC-01C、HPC-01D、HPC-03C |
| 其餘待辦（24 項） | HPC-03D 至 HPC-09；既有程式能力不等於已通過各項驗收 |

合計尚有 27 項未勾選，其中 3 項已有部分進度。此數量依現有紀錄彙整，
後續 goal 應依實際驗收結果更新。

啟動後續 goal 時：

1. 指定子任務編號或明確指定整份清單，並閱讀該項的前置條件及進度報告。
2. 檢查 `git status` 與目前實作，保留既有未提交修改；確認哪些證據適用於目前版本。
3. 列出剩餘工作、驗收命令及所需硬體，再開始實作。只完成部分驗收時，保留未勾選狀態。

## 建議開發順序

以下共 10 個階段、38 個子任務；每個子任務都有固定編號，可直接用於指定 goal。
完成狀態以各階段的核取方塊與驗收報告為準。

優先級：P0 是後續開發的基礎，P1 是單機效率與長時間計算的主要工作，
P2 是分散式耦合與部署驗收。優先級用來選擇下一項工作，實際開工仍須滿足
「階段與依賴」列出的前置條件；HPC-09 的建置與測試配置可隨各階段逐步補齊。

| 任務範圍 | 優先級 | 開發重點 | 子任務數 |
|---|---|---|---:|
| HPC-00A–D | P0 | 建立數值與效能基準 | 4 |
| HPC-01A–D | P0 | Communicator、所有權與錯誤處理 | 4 |
| HPC-02A–C | P1 | 單機 OpenMP 多核心組裝 | 3 |
| HPC-03A–D | P1 | 浸入式／移動流場 MPI | 4 |
| HPC-04A–C | P1 | 線性求解器與預條件器 | 3 |
| HPC-05A–D | P1 | 耦合 checkpoint／restart | 4 |
| HPC-06A–D | P1 | 記憶體、輸出與前處理 | 4 |
| HPC-07A–D | P2 | 分散式 FSI | 4 |
| HPC-08A–C | P2 | 多域程序群與資源配置 | 3 |
| HPC-09A–E | P2 | 建置、排程與跨節點驗收 | 5 |

1. **基準與可靠性：HPC-00／HPC-01。** 先確認數值基準、量測方式、
   communicator 與錯誤處理，讓後續效能改善可以被驗證。
2. **單機多核心：HPC-02。** 針對已量測的組裝瓶頸加入 OpenMP，
   與相同總核心數的純 MPI 比較。
3. **長時間與大型計算：HPC-03／04／05／06。** 按下方相依條件推進
   浸入式 MPI、求解器、續跑與 I/O；既有 graph 的 checkpoint 可提早進行。
4. **分散式耦合與部署：HPC-07／08／09。** 完成 FSI、domain 程序群，
   最後取得實際跨節點驗收證據。

已有部分完成紀錄，後續可先指定 `HPC-01C` 的剩餘工作；
每次 goal 應以當時程式與測試證據核對清單狀態。

### 下一次開發的具體入口

| 子任務 | 優先核對的剩餘工作 | 接續閱讀 |
|---|---|---|
| HPC-01C | F01／F03／F04 與 F02 的 CPU／1D／耦合／四個序列工具已驗證；CUDA stdout 亦已驗證；F05 診斷亦已驗證；接續 F06 完整入口覆蓋 | [MPI 錯誤邊界索引](architecture/MPI_FAILURE_BOUNDARIES.md)、[本輪串行工具與 FSI 驗收](progress/HPC_01C_SERIAL_TOOL_PROGRESS.md) |
| HPC-01D | 剩餘 embedding／工具入口、PETSc 型別、配置與後端能力矩陣；CUDA 啟動檢查已補入 | [資源檢查進度](progress/HPC_01D_PROGRESS.md)、[工具進度](progress/HPC_01CD_TOOLS_PROGRESS.md)、[CUDA 進度](progress/HPC_01D_CUDA_PROGRESS.md) |

上述入口用於定位下一批工作；各子任務仍須滿足下方完整驗收條件。
Checkpoint 的完整發布與恢復協議繼續由 HPC-05 追蹤。
1D／graph／sequential 的 216 個原生作業與四個序列工具的 48 個作業亦已通過，
詳見 [耦合與工具 stdout 驗收](progress/HPC_01C_COUPLING_STDOUT_PROGRESS.md)。
CUDA 的 38 個作業與 27 個場比較亦通過，見
[CUDA stdout 驗收](progress/HPC_01C_CUDA_STDOUT_PROGRESS.md)。
F05 的原始診斷例外保存與恢復亦通過，見
[Bezier 診斷驗收](progress/HPC_01C_BEZIER_DIAGNOSTIC_PROGRESS.md)。

## 可複製的 goal 範本

可複製以下任務描述，將工作編號換成當次目標：

```text
請以 docs/WORKSTATION_HPC_TODO.md 為依據，建立並執行一個 goal：完成 HPC-01C 的剩餘工作。
先閱讀 AGENTS.md、該任務的前置條件及目前相關實作，再進行開發。
保留既有數值與檔案介面，完成該子任務的實作、適當測試與文件。
在實作前確定驗收案例和數值容許值；不得為通過測試而放寬既有標準。
更新本清單狀態，並在 docs/progress/ 下記錄命令、環境、結果與限制。
沒有執行的測試必須明確標示，不可把編譯通過當作數值或平行驗收通過。
若缺少硬體或排程資源，先完成可做的部分並記錄尚待驗證項目，不宣稱全部完成。
本次範圍止於指定子任務，不自動進入下一階段。
```

若要開發整份清單，可改用以下描述：

```text
請以 docs/WORKSTATION_HPC_TODO.md 為依據，建立並執行一個 goal：完成整份清單的剩餘工作。
先核對工作區既有修改、完成報告與驗收證據，再依前置條件分批開發。
遵守 AGENTS.md，保留既有數值與檔案介面，不為通過測試而放寬標準。
每批完成後更新清單與 docs/progress/ 報告，記錄測試命令、環境、結果及限制。
部分交付不代表整個 goal 完成；需要叢集或其他硬體驗收的項目，在取得證據前保持未勾選。
```

勾選規則：只有實作與該子任務的必要驗收均完成，才將 `[ ]` 改為 `[x]`。
部分完成、環境受阻、尚未量測都保持未勾選，並補上報告連結及剩餘工作。
建議報告檔名為 `docs/progress/HPC_00A_REPORT.md`，依任務編號調整。

## 目前基線與設計原則

開始任務時重新查核此表；本表是建立文件時的閱讀結果，不是新執行的效能測試。

| 路徑 | 目前能力 | 本計畫要補的部分 |
|---|---|---|
| CPU 貼體 3D 流場／傳輸 | MPI/PETSc、owned-row 組裝 | 擴展性證據、可配置 communicator、求解與 I/O 改善 |
| 原生 1D | 部分隱式求解使用 MPI；部分顯式／物種更新使用 OpenMP | 釐清各方法的執行模式、減少大型案例的重複資料 |
| 0D／1D／貼體 3D graph | 有 MPI 耦合案例及 trial/commit/rollback | 全域 checkpoint、程序群配置與可靠失敗處理 |
| 浸入式／移動流場 | 靜態 runtime／準靜態 graph 已接 MPI；暫態／移動維持序列 | 固定幾何 velocity history、移動 ownership／halo 與狀態移轉 |
| FSI | 單分區膜結構與流體強耦合 | 分散式介面、全域收斂與回復 |
| Spline | 部分 OpenMP、分塊與串流輸出 | 大型案例的時間／記憶體證據 |
| Mesh／packer | 單程序工具 | 先量測，只有成為瓶頸時才擴充 |
| CUDA | 單程序、單 GPU | 維持數值相容與回歸；多 GPU 另立計畫 |

必須保留的原則：

- 沿用既有求解器、元素核心及 coupling contracts，不整套重寫。
- 保持 `controlmesh.vtk`、`bzmeshinfo.txt`、`.igacache`、legacy extraction、
  METIS 分區、速度場與 `.ntiga` 的既有相容性。新格式必須明確版本化。
- 保留既有 schema 的接受／拒絕行為；執行資源設定不得暗中改變物理模型。
- 保持 SI coupling ports、向外流量為正、既有守恆檢查及 trial/commit 語義。
- MPI 單 rank 是正式支援模式。OpenMP 是可選加速，純 MPI 必須繼續可用。
- 小型 0D 與廉價 1D 工作可以集中或複製執行；不要求每個模組都做 MPI 分割。
- 不將 `mpiexec` 可啟動、多 rank 結果相近、或編譯含 `-fopenmp` 當作加速證據。
- 大型測試在排程器分配的計算資源上執行。保留現有工作目錄修改，
  不提交二進位、案例產物、VTK 結果或 scheduler logs。

## 階段與依賴

| 編號 | 工作 | 前置條件 |
|---|---|---|
| HPC-00 | 效能與數值基準、環境紀錄 | 無 |
| HPC-01 | Communicator、所有權介面與一致失敗處理 | HPC-00A/B |
| HPC-02 | 單機 OpenMP 組裝 | HPC-00；與 HPC-01 的資料介面對齊 |
| HPC-03 | 浸入式／移動流場 MPI | HPC-01；串行基準已建立 |
| HPC-04 | 可擴展線性求解 | HPC-00；MPI 驗收依賴 HPC-03 |
| HPC-05 | 耦合 checkpoint／restart | HPC-00、HPC-01；可先於 HPC-03 完成既有路徑 |
| HPC-06 | 分散式記憶體與 I/O | HPC-00、HPC-01；新路徑依賴 HPC-03 |
| HPC-07 | 分散式 FSI | HPC-03、HPC-04、HPC-05A/B 的持久化契約 |
| HPC-08 | 多域程序群與資源配置 | HPC-01、HPC-03、HPC-05 |
| HPC-09 | 單機／跨節點整合驗收 | 欲宣告支援的前述階段均通過 |

建議先完成 HPC-00，接著 HPC-01／HPC-02；在 HPC-03 推進期間，
及早補 HPC-05 的既有 graph 路徑。不要等所有平行化完成才處理續跑。
表中的獨立性代表技術依賴，不代表必須同時開發。

## HPC-00：建立可比較的基準

入口：[BENCHMARKS.md](BENCHMARKS.md)、[BRIDGES2.md](BRIDGES2.md)、
[check_dependencies.sh](../scripts/check_dependencies.sh)、
[既有階段報告](progress/FINAL_ROADMAP_REPORT.md)。

- [x] **HPC-00A：執行環境與案例清單。** 選定小型貼體血流、傳輸、
  浸入式及 FSI 基準，記錄 revision、輸入雜湊、編譯器、PETSc/MPI/HDF5、
  CPU/GPU、rank/thread 數、綁定方式與求解器選項。
  區分 source fixture、生成資料與權威參考結果。
  已交付工具、四個 fixture 清單與實際環境／輸入紀錄，見
  [HPC-00A 報告](progress/HPC_00A_REPORT.md)。數值與擴展性驗收仍依以下子任務追蹤。
- [x] **HPC-00B：分階段量測。** 記錄幾何更新、積分／組裝、預條件器建立、
  線性求解、耦合迭代、通訊及 I/O；定義計時區間，避免重疊時間被重複相加。
  保留每 rank 時間／RSS、最大值與不平衡程度，CUDA 另記 peak allocation。
  四個選定基準及支援的 CUDA 路徑已產生分階段紀錄，完成數值回歸、
  每 rank RSS／不平衡彙整及 GPU 配置峰值追蹤，見
  [HPC-00B 報告](progress/HPC_00B_REPORT.md)。計時範圍與未歸類工作有明確標示；
  這些並行執行的驗證不是獨立效能重複測試。
- [x] **HPC-00C：單機測試矩陣。** 依實際硬體測試 1、2、4、8 等有效核心數，
  區分純 MPI、純 OpenMP、混合模式及單 GPU。尚未支援的模式標示 N/A。
  分開冷啟動與重複執行，報告重複次數、時間中位數與變動。
  已新增首次／三次重複的 CPU MPI 矩陣工具，保留完整 rank 量測與場比較，
  失敗不納入成功統計。貼體流場與傳輸 1／2／4／8 rank 共 32 次執行通過，
  流場另通過 16 次獨立質量檢查，並量測到 required 元素重複積分與記憶體增長；見
  [HPC-00C 進度](progress/HPC_00C_PROGRESS.md)。首次求解不等於受控冷快取；
  後續亦完成 GPU 流場／傳輸與浸入式共 12 次執行及 4 次獨立 GPU 質量檢查，
  見 [單 GPU／單程序進度](progress/HPC_00C_SERIAL_PROGRESS.md)。FSI 完整四次亦
  已通過；共 48 次選定模式執行的證據齊備，未支援模式標示 N/A。
- [x] **HPC-00D：數值基準與效能門檻。** 為各案例固定相對 L2、Jacobian、
  殘差、質量／物種守恆、FSI 合力／力矩的驗收方式。
  在優化前登記容許值與有意義的效能目標，不能事後挑選有利案例。
  CPU rank 比較、浸入式及 FSI 原生 gate 已執行；修正 CUDA 壁面 trace
  後穩態／暫態比較通過。早期紀錄見
  [進度報告](progress/HPC_00BD_PROGRESS.md)。
  獨立物種收支已通過 69 項解析／錯誤檢查、CPU 1／2／4／8 rank 與 GPU
  十個時間步及五種 CLI 錯誤驗收，見
  [物種收支驗收](progress/HPC_00D_TRANSPORT_BUDGET.md)。浸入式／FSI 各兩份
  完整參考及共 12 個物理量的相對 L2 比較亦已通過，並驗證拒絕錯誤及
  封存來源／binary。四個選定基準的來源、門檻與完成範圍見
  [HPC-00D 報告](progress/HPC_00D_REPORT.md)；不代表後續平行模式已完成。

驗收：機器可讀的量測結果與可重現命令齊全；可分辨求解時間與輸出時間，
也可辨識「沒加速」或「記憶體增加」。本階段不承諾任何加速比例。

## HPC-01：統一平行執行介面

入口：[OwnedRowAssembler.hpp](../solvers/cpu/include/OwnedRowAssembler.hpp)、
[OneDImplicit.hpp](../solvers/one_d/include/OneDImplicit.hpp)、
[DomainRuntimeRegistry.hpp](../include/DomainRuntimeRegistry.hpp)。

- [x] **HPC-01A：Communicator 注入。** 盤點 runtime 內的 `COMM_WORLD`／
  `COMM_SELF`，定義 communicator 的擁有、複製與釋放責任；保留舊 CLI 預設。
  增加兩個互不相交 communicator 的小測試，證明沒有意外的 world collective。
  1D 四種隱式路徑、貼體 3D 流場／傳輸、flow／species graph 已通過獨立
  程序群與串行解比較，完成 borrowed lifetime 稽核；完整 CLI 相容性回歸
  在保留逾時紀錄後以較長時限重跑通過，見
  [HPC-01A 進度](progress/HPC_01A_PROGRESS.md)。
- [x] **HPC-01B：所有權與 halo 介面。** 明確區分全域 ID、本地 ID、owned、
  ghost，以及元素積分／矩陣列／介面量的歸屬；避免把相接元素重複計入全域積分。
  不強制所有後端採用相同資料容器，但要有一致的語義與可檢查條件。
  已實作共用 ID／row 映射、分散式精確覆蓋與元素貢獻檢查，並接入貼體
  runtime 的 halo row 建立。空 rank、遺漏／重複貢獻、表面發布權、真實案例
  halo 更新與數值回歸均通過；浸入式／FSI 分散計算仍由 HPC-03／07 驗收，見
  [所有權契約](architecture/PARALLEL_OWNERSHIP.md)與[進度報告](progress/HPC_01B_PROGRESS.md)。
- [ ] **HPC-01C：一致失敗處理。** 在輸入、配置、組裝與輸出階段建立跨 rank
  錯誤協議，避免一個 rank 丟例外而其他 rank 永久等待。
  對可控制的局部錯誤進行故障注入，並用 timeout 證明可有限時間退出。
  明確區分協同失敗退出與作業系統殺死 rank；本任務不承諾 MPI 程序故障後原地續跑。
  已接入 graph 參數／輸入／domain preflight／輸出、owned-row 預配置，以及
  貼體流場／傳輸元素組裝；三 rank 局部故障、非零作業退出及 rollback 已驗證。
  graph manifest、0D／1D／3D 配置快照、停止步數與 Newton 控制已增加一致性
  檢查；graph 啟動時也會比較可見 PETSc 選項，22 項案例、原生 CLI 的
  rank 設定差異拒絕及獨立 communicator 數值回歸通過。先前的 1005-node
  貼體回歸另有版本化紀錄。其餘 CLI、外部資產一致性、建構後段與耦合步驟
  錯誤邊界仍需補齊，故保持未勾選，見
  [HPC-01C 進度](progress/HPC_01C_PROGRESS.md)與
  [配置／建構補強](progress/HPC_01C_CONFIGURATION_PROGRESS.md)及
  [PETSc 選項驗收](progress/HPC_01C_PETSC_OPTIONS_PROGRESS.md)。原生 1D CLI
  也已補上啟動配置及 root 輸出錯誤協調，18 項原生 CLI 案例通過，見
  [1D CLI 驗收](progress/HPC_01C_ONE_D_CLI_PROGRESS.md)。1D checkpoint 的本地
  PETSc I/O、metadata／state 一致性及 restore 邊界亦已補強，25 項續跑／故障／
  格式案例及獨立群組 unit 通過，見
  [1D checkpoint 驗收](progress/HPC_01C_ONE_D_CHECKPOINT_PROGRESS.md)。CLI 時間步
  與 combined trial 的局部階段已加入協調及 rollback／retry 驗證，21 項 CLI、
  25 項 checkpoint 與最終 graph 群組回歸通過，見
  [1D trial 驗收](progress/HPC_01C_ONE_D_TRIAL_PROGRESS.md)。staged hydraulic／transport
  局部階段及 OpenMP 物種例外亦已補強，七種 staged 故障、1／4 threads、
  無 OpenMP runtime 與最終 graph 群組回歸通過，見
  [staged 驗收](progress/HPC_01C_ONE_D_STAGED_PROGRESS.md)。四種 implicit 方法的本地
  組裝、回傳碼、SNES callbacks 及狀態發布亦已補強，故障／retry、四種方法群組、
  21 項 CLI、25 項 checkpoint 與 graph 回歸通過，見
  [implicit 驗收與限制](progress/HPC_01C_ONE_D_IMPLICIT_PROGRESS.md)。兩種 1D adapter
  的本地操作／求解準備與 staged 輸入發布亦已補強，20 種故障、1／4 threads、
  純 C++、獨立 graph 群組及原生 precommit 失敗退出通過，見
  [1D adapter 驗收](progress/HPC_01C_ONE_D_ADAPTER_PROGRESS.md)。貼體 3D port 的查詢
  準備／積分／結果、跨 rank 查詢描述與 adapter flow 輸入／求解準備亦已補強；
  20 種錯誤模式、解析 port 場、非零流場 retry 及 flow／species graph 群組回歸通過，見
  [3D port 驗收](progress/HPC_01C_THREE_D_PORT_PROGRESS.md)。flow／transport 的主要
  trial／rollback／Abort／prepare、候選快照、required-state 查詢及 outlet evaluation
  也已補強；48 個錯誤案例與 port／graph／VCA 續跑回歸通過，見
  [3D trial 驗收](progress/HPC_01C_THREE_D_TRIAL_PROGRESS.md)。staged 3D 的輸入候選發布、
  plan 一致性、傳輸邊界與物種 accounting，以及原生 transport mass／source
  積分亦已補強；72 個故障案例、流向反轉與 port／trial／graph／VCA 回歸通過，見
  [3D staged 驗收](progress/HPC_01C_THREE_D_STAGED_PROGRESS.md)。flow 參考流量與
  Summary 也已補上本地錯誤及非有限結果協調；64 個 port／診斷案例及同版
  staged／trial／graph／VCA 回歸通過，見
  [流場診斷驗收](progress/HPC_01C_FLOW_DIAGNOSTICS_PROGRESS.md)。Newton 的 assembly／
  scatter／KSP setup／solve／update 回傳碼、收斂決策與 logging，以及 transport
  assembly／boundary rows 亦已補強；75 個 trial 案例、同版 graph／VCA 回歸與
  10 份 rank profile 完整性檢查通過，見
  [Newton 驗收](progress/HPC_01C_NEWTON_PROGRESS.md)。flow 初始化候選發布、resolved
  邊界一致性，以及 transport 全域 gather 亦已補強；93 個 trial、83 個 staged、
  12 個 gather 故障案例及同版 port／graph／VCA 回歸通過，見
  [初始化與 gather 驗收](progress/HPC_01C_INITIALIZATION_PROGRESS.md)。3D checkpoint
  讀取亦已補上本地檔案驗證、候選發布與 metadata 協調；67 個讀取故障案例、
  98 個 staged 案例、14 項原生 CLI 損毀續跑與同版 graph／VCA 回歸通過，見
  [checkpoint 讀取驗收](progress/HPC_01C_CHECKPOINT_READ_PROGRESS.md)。3D checkpoint
  寫入也已補上 owned-row 暫存檔、共同錯誤與單 Vec 發布；65 個 writer 故障、
  104 個 staged 案例、17 項原生輸出故障及同版讀取／graph／VCA 回歸通過，見
  [checkpoint 寫入驗收](progress/HPC_01C_CHECKPOINT_WRITE_PROGRESS.md)。3D adapters 的其餘邊界、
  executor bookkeeping、其他 runtime 邊界與外部資產一致性仍待完成；
  完整持久化協議另依 HPC-05 追蹤。建構期間的 PETSc 清理、boundary catalog 與
  halo 錯誤協調已整合，正式路徑通過 229 個故障／重試案例及耦合 CLI 回歸；
  見 [建構失敗修正進度](progress/HPC_01C_CONSTRUCTION_PROGRESS.md)。貼體 runtime
  呼叫端配置、參數複製與設定初始化亦已補強，277 個建構案例及同版
  trial／graph／bifurcation CLI 及 VCA 續跑回歸通過；其餘邊界保持待辦，見
  [建構輸入與配置進度](progress/HPC_01C_CONSTRUCTION_INPUT_PROGRESS.md)。graph 外部資產
  已加入內容比對，40 個讀取故障／重試、31 項 graph 回歸及原生副本 CLI 通過；
  同時修正 runtime 誤拒絕本地路徑差異的問題。其餘 CLI 的資產檢查仍待補齊，見
  [資產一致性驗收](progress/HPC_01C_ASSET_PROGRESS.md)。原生 1D CLI 也已補上
  network、實際使用的 table／replay 內容比對及設定檔特殊型別拒絕；最終版
  35 項資產觀察、22 項 CLI 與 25 項 checkpoint 回歸通過，其餘 CPU CLI
  與 runtime 邊界仍待完成，見 [1D 資產驗收](progress/HPC_01C_ONE_D_ASSET_PROGRESS.md)。
  獨立 CPU flow 的建構前輸入亦已補強：51 筆輸入／相容性觀察、16 項資源檢查、
  VCA 續跑及 OpenMP 版本回歸通過其指定 gates；legacy 近零壓力的跨 rank
  比較仍保留四筆未通過觀察。其餘 CPU CLI 與建構後段保持待辦，見
  [CPU flow 輸入驗收及限制](progress/HPC_01C_FLOW_INPUT_PROGRESS.md)。
  `iga_solve` 的輸入、速度序列、時間步、memory report、PETSc 狀態與輸出／
  checkpoint 協調亦已補強；81 筆 CLI、106 個 reader 故障／重試、舊 binary
  讀新 checkpoint、實際 HDF5 資料及資源／VCA 回歸通過。其餘 CLI、flow 後段
  與共用 I/O helper 保持待辦，見 [傳輸 CLI 驗收](progress/HPC_01C_TRANSPORT_CLI_PROGRESS.md)。
  共用 VTU／PVD／速度序列／physiology writer 已補上關檔錯誤檢查，12 個故障與
  8 次重試、CPU CLI 及 GPU 小案例場比較通過；其餘 I/O 與呼叫端仍待完成，見
  [文字輸出關檔驗收](progress/HPC_01C_TEXT_OUTPUT_PROGRESS.md)。
  CPU flow 的 gather／場檔／索引與成功摘要順序亦已補強，33 個 writer 故障及
  重試、46 筆原生輸出觀察、HDF5 資料及既有 flow／VCA 回歸通過指定條件；
  建構後其他邊界仍待完成，見 [flow 輸出驗收](progress/HPC_01C_FLOW_OUTPUT_PROGRESS.md)。
  原生 flow／VCA 的每步本地輸入、結果／budget 與 circuit／history 協調也已補上；
  66 個函式故障及健康重試、13 筆 CLI 觀察、最後接受步 checkpoint 比較與
  既有 flow／VCA 回歸通過指定條件。整體 rollback 與其餘邊界仍待完成，見
  [flow 每步驗收](progress/HPC_01C_FLOW_STEP_PROGRESS.md)。
  共用 geometry report／coupling history 與 CPU／CUDA HDF5 明確關檔亦已補強；
  17 個文字故障、5 個 HDF5 故障與重試、原生 MPI／CUDA 關檔注入及資料讀回通過。
  預設 flow 收斂設定下的 CPU/GPU 場差異仍保留為未通過觀察；採既有基準的
  較嚴格設定後比較通過，未放寬數值門檻。其餘 runtime、executor 與清理邊界
  保持待辦，見 [共用 I/O 關檔驗收](progress/HPC_01C_IO_FINALIZATION_PROGRESS.md)。
  Legacy `iga_transport` 的資產／有效 PETSc 選項比對、局部組裝、矩陣／vector
  回傳碼與明確清理亦已補齊；46 個故障及健康重試、43 筆 CLI／82 份 rank
  reports、39 項既有工具回歸通過。舊版接受不一致參數的錯誤結果保留，
  新版在組裝前拒絕；其餘工具與耦合邊界仍待完成，見
  [Legacy transport 驗收](progress/HPC_01C_LEGACY_TRANSPORT_PROGRESS.md)。
  Mesh checker／assembly smoke 亦已加入資料庫內容比對；mesh owner 範圍與
  index／record 不一致會被拒絕，smoke 清理失敗不再先印結果摘要。
  38 筆原生案例／67 份 rank reports 及 39 項既有工具回歸通過；保留舊版
  三種空 geometry coverage 的成功誤報及清理順序重現，見
  [MPI 工具資產驗收](progress/HPC_01C_TOOL_ASSET_PROGRESS.md)。其餘 runtime／
  adapter／executor 與支援入口覆蓋稽核仍待完成。
  Pressure-flow executor 已接入跨 rank outcome／全域收斂協調，60 次故障與
  60 次同 step 健康重試通過；兩個原生 3-rank 案例的 14 份 CSV 與修改前
  完全相同。其完整 multidomain smoke 已退出 0，見
  [executor 進度](progress/HPC_01C_PRESSURE_EXECUTOR_PROGRESS.md)。
  Species executor 亦已加入 outcome／全群 hydraulic 收斂與 donor／transport
  排程協調，107 次故障與 107 次健康重試通過；原生 schema-v6 的 1／2-rank
  案例共 10 份 CSV 與修改前相同，見
  [species executor 進度](progress/HPC_01C_SPECIES_EXECUTOR_PROGRESS.md)。
  Registry 已改為完成共同驗證後才接管 owners；27 次 ownership 故障／重試、
  69 次原生 flow／species／0D graph 故障／重試通過，57 份 CSV 與修改前相同。
  Adapter catalog 與初始 port／transport 建構邊界亦已補強，見
  [registry 進度](progress/HPC_01C_REGISTRY_PROGRESS.md)。
  Sequential CSV／manifest 已增加明確關檔與共同錯誤診斷；51 個原生案例及
  82 份 rank reports 通過預期，22 份 CSV／JSON 與修改前相同，見
  [sequential 關檔驗收](progress/HPC_01C_SEQUENTIAL_OUTPUT_PROGRESS.md)。
  Sequential 初始化與輸入亦已協調，並接上 1D trial outcome；120 次 stage
  故障／重試、12 次有效 PETSc 選項拒絕／重試及 32 個 CLI 案例通過。
  原生 1／2／3-rank 的 18 份 CSV 與修改前相同，見
  [sequential 初始化驗收](progress/HPC_01C_SEQUENTIAL_INITIALIZATION_PROGRESS.md)。
  Sequential strong-fixed／Aitken 亦已協調錯誤、全群收斂與全部 abort outcome；
  162 次故障／重試、6 次全群否決及原有完整 smoke 通過，12 份 CSV 與修改前
  相同，見 [sequential strong 驗收](progress/HPC_01C_SEQUENTIAL_STRONG_PROGRESS.md)。
  Configured transport 的 VTKHDF 初始化亦已改用共同錯誤階段；34 個原生案例、
  51 份 rank reports 通過預期，128 份文字檔相同及 16 次 HDF 資料集比較通過，
  見 [transport 初始化驗收](progress/HPC_01C_TRANSPORT_VISUALIZATION_INIT_PROGRESS.md)。
  共用完整文字讀取與指定序列化亦已補強，15 次原生配置失敗／重試、132 個
  flow／transport CLI 案例與既有 unit／smoke 通過；57 份 graph CSV 相同，
  CUDA 前後場一致並符合原定 CPU/GPU 門檻，見
  [文字串流驗收](progress/HPC_01C_CHECKED_TEXT_PROGRESS.md)。1D checkpoint metadata／fingerprint、
  VTU 檔名、VTKHDF schema 與資源摘要亦已補強，117 次配置失敗／重試、25 個
  checkpoint 案例通過；190 份 CPU／1D 輸出與修改前相同，兩組 HDF 資料集
  比較通過，見 [文字 helper 驗收](progress/HPC_01C_TEXT_HELPER_PROGRESS.md)。
  1D VTP 檔名與 flow／transport／staged adapter 簽章亦已補強，原生格式化例外
  與既有故障測試詳見 [串流邊界驗收](progress/HPC_01C_STREAM_BOUNDARY_PROGRESS.md)。
  支援入口、共同階段及剩餘 F01～F06 已整理於
  [MPI 錯誤邊界索引](architecture/MPI_FAILURE_BOUNDARIES.md)；完整覆蓋尚未完成。
  F01 abort 診斷與 F03 snapshot 的既有修正已重建驗證；packer、Womersley、
  transport budget 及 FSI helper 已補齊讀寫檢查，原生故障／重試與 cache/text
  相容性驗證詳見 [串行工具進度](progress/HPC_01C_SERIAL_TOOL_PROGRESS.md)。
  完整 FSI 兩版皆 4 次迭代收斂，5 份輸出逐位元相同、29 個數值 arrays 相對
  L2=0；CPU／1D／coupling 核心回歸與 84 個 MPI stdout 案例亦通過。
  F02 的 CPU flow／runtime、configured／legacy transport 亦已通過 412 個
  原生作業、618 份 rank reports 與 310 個場比較（relative L2=0），見
  [CPU stdout 驗收](progress/HPC_01C_SOLVER_STDOUT_PROGRESS.md)。
  F02 其他入口與 F05 已完成相應驗收，見上方報告。F06 的終止呼叫鏈已核對，
  runtime Close、144 次故障／重試與最終 112 個原生 CLI 作業通過；
  graph／sequential 清理失敗不發布 completion manifest，見
  [runtime cleanup 驗收](progress/HPC_01C_RUNTIME_CLEANUP_PROGRESS.md)與
  [入口稽核](progress/HPC_01C_ENTRY_AUDIT.md)。F06 其餘呼叫鏈與
  HPC-01D 剩餘能力矩陣繼續追蹤。
  後續已核對五個 main 的控制流程，並修正 memory report 關檔失敗仍回報成功：
  12 個原生作業、10 個場比較及 world3／split1+2 關檔與重試通過，見
  [記憶體報告驗收](progress/HPC_01C_MEMORY_CLOSE_PROGRESS.md)。

- [ ] **HPC-01D：執行資源與能力檢查。** 檢查分區數／rank 數、執行緒配置、
  PETSc index/scalar 型別及必要後端能力；提供清楚錯誤資訊與執行摘要。
  已接入 CPU flow／transport、1D 與 graph 啟動檢查及摘要，並查詢明確指定的
  factor backend。150 個 resource unit 案例、實際 MPI thread level、16 項
  原生 CLI 資源案例及既有數值回歸通過；其餘入口／配置與後端矩陣仍待
  完成，故保持未勾選，見 [HPC-01D 進度](progress/HPC_01D_PROGRESS.md)。後續亦已
  接入 mesh check、assembly smoke 與 legacy transport，39 項工具回歸、
  幾何錯誤同步及既有 runtime／VCA 回歸通過，見
  [其餘 MPI 工具驗收](progress/HPC_01CD_TOOLS_PROGRESS.md)；其餘能力矩陣保持待辦。
  Standalone FSI exporter 已在建立輸出前拒絕多 rank；2／3-rank 拒絕與既有
  輸出保留檢查通過，單 rank 完整案例亦已以 4 次耦合迭代收斂，見
  [exporter 入口](progress/HPC_01C_EXPORTER_ENTRY_PROGRESS.md)。
  CUDA 五個入口已加入單程序、thread requests、索引容量與版本摘要檢查；
  84 項 host checks、23 個原生啟動案例及 CPU 共用 parser 的 150 項檢查通過。
  數值與限制詳見 [CUDA 資源驗收](progress/HPC_01D_CUDA_PROGRESS.md)。

驗收：既有單程序與小型 MPI 數值門檻通過；subcommunicator 測試無串擾；
錯誤注入不掛住、不寫出成功標記。不以替換 communicator 常數宣稱分散式完成。

## HPC-02：單機多核心組裝

入口：[ImmersedTransientFlowRuntime.hpp](../solvers/cpu/include/ImmersedTransientFlowRuntime.hpp)、
[NavierStokesElement.hpp](../solvers/cpu/include/NavierStokesElement.hpp)、
[CutCellVolumeQuadrature.hpp](../solvers/cpu/include/CutCellVolumeQuadrature.hpp)。

進度：有界批次執行器已接入浸入式暫態／moving／FSI 的體積積分，
實際 27-cell 的 1／2／4-thread 殘差、Jacobian action 與錯誤回復已通過。
修正線性求解器前的完整 FSI 1／2／4／8-thread 九場比較已通過
（共 36 個比較，相對 L2 皆為 0）。修正後的 FSI 隔離矩陣亦已完成，
八次執行與 72 個場比較通過；4-thread 速度比 2.30、RSS ratio 1.022 均符合
預定門檻，見 [FSI 效能驗收](progress/HPC_02_FSI_PERFORMANCE.md)。
無 OpenMP binary 的修正後九場回歸也已通過。貼體同 rank 單 thread
基準的 12 次執行亦已通過；OpenMP 的每 rank 與合計峰值 RSS 比值
皆不超過 1.009，符合 1.25 預算。HPC-02A／B／C 的必要驗收已完成，
最終證據核對在 `outputs/hpc02/volume/completion-audit.json`。
後續固定幾何測試發現真實線性殘差失敗；FGMRES 修正已通過完整 Newton
與 2／4／8-thread FSI 九場回歸；修正後的 1-thread 完整場回歸也已通過，
1／4-thread 隔離矩陣各完成首次與三個正式重複，見
[線性求解修正](progress/HPC_02_LINEAR_SOLVE_PROGRESS.md)。
貼體 MPI／OpenMP 元素組裝已接入，2／4 ranks 的 1／2-thread 故障回復、
場比較及既有獨立 communicator 回歸通過；固定四核心矩陣的 12 次執行及
與封存單 rank 參考的 24 個場比較也已通過，見
[混合組裝進度](progress/HPC_02_HYBRID_PROGRESS.md)。
見 [體積積分整合進度](progress/HPC_02_VOLUME_PROGRESS.md)；早期 helper 機制證據
保留於 [批次執行器進度](progress/HPC_02_BATCH_PROGRESS.md)。

- [x] **HPC-02A：平行化一個已量測瓶頸。** 優先選元素積分／局部矩陣，
  使用 thread-local scratch、明確工作切分與有上限的批次暫存。
  捕捉工作執行緒錯誤，離開平行區域後統一處理；不得直接跨 OpenMP 區域拋出例外。
- [x] **HPC-02B：安全合併。** 計算與 PETSc 插入分離；不讓多執行緒未經驗證地
  同時修改同一 Mat/Vec。不在 OpenMP worker 呼叫 MPI，除非明確要求、檢查並驗證
  對應的 MPI thread support。
- [x] **HPC-02C：執行配置與比較。** 支援關閉 OpenMP、設定 thread 數與 CPU 綁定，
  管理 BLAS 執行緒避免過量配置；比較固定總核心數的純 MPI 與混合模式。
  後續才依量測擴充幾何分類、surface／ghost 組裝等區段。

驗收：1 thread 回歸通過，2／4／可用上限的結果通過預設數值門檻；
無共享 scratch 污染，峰值記憶體符合預先設定的預算；以端到端時間判斷成效。

## HPC-03：浸入式與移動流場的 MPI

入口：[ImmersedFlowCase.hpp](../include/ImmersedFlowCase.hpp)、
[ImmersedStaticFlowRuntime.hpp](../solvers/cpu/include/ImmersedStaticFlowRuntime.hpp)、
[MovingImmersedTransientFlowRuntime.hpp](../solvers/cpu/include/MovingImmersedTransientFlowRuntime.hpp)。

- [x] **HPC-03A：靜態 owned rows 與 halo。** 分配 active cells／自由度，
  建立精確稀疏配置與必要鄰接資料交換，保留單 rank 比較路徑。
  處理 pressure gauge、port multiplier、ghost penalty 跨分區貢獻與空 rank。
  正式 operator 與 Newton runtime 已完成 owned Mat／Vec、required-state halo、
  唯一積分、全域 convergence／port／wall／pressure／守恆診斷與交易回復。
  27 cells／54 ghost faces 的 closed、flow、pressure 案例在 1／2／4 ranks
  對照既有序列場解通過；零列／零工作配置、split 1+2、候選與 prepare 失敗
  後重試亦通過。C++ 靜態介面已支援，case／graph／暫態入口依 HPC-03C
  另行接入與放行。見 [組裝進度](progress/HPC_03A_ASSEMBLY_PROGRESS.md)、
  [operator 進度](progress/HPC_03A_STATIC_OPERATOR_PROGRESS.md) 與
  [靜態 runtime 驗收](progress/HPC_03A_STATIC_RUNTIME_PROGRESS.md)。
- [x] **HPC-03B：按計算量分區。** 以積分點、切割與穩定化工作量建立初始權重，
  比較加權／未加權分區的最大 rank 耗時與通訊成本；避免只按 cell 數量分配。
  已加入穩定 cell 順序的 weighted contiguous 選項與 checked work model，
  驗證分區最小最大權重、全域物理一致性、Newton／空工作 rank 及故障回復。
  1／2／4 ranks 固定核心比較已完成；此案例的 cell-count 原已達模型最適，
  未觀察到明確加速，故保留預設。見
  [工作量分區與單機比較](progress/HPC_03B_WORK_PARTITION_PROGRESS.md)。
- [ ] **HPC-03C：固定幾何暫態與 graph。** 驗證上一時間步速度、port 測量、
  全域守恆及 trial rollback。只有這些門檻通過後，才移除對應已支援路徑的
  MPI size 1 限制；未完成的模式維持明確拒絕。
  準靜態 graph 的 collective port 更新、owned runtime／adapter 與 case 建立
  已接入；正式 explicit／fixed／Aitken 入口通過 1／2／4 ranks、介面門檻與
  precommit failure，序列 schema 回歸通過。固定幾何已新增 owned velocity
  snapshot／必要 halo、精確 clock 與失敗重試；owned backward-Euler 體積
  組裝已凍結外力，驗證 compact／expanded 與分散式殘差／Jacobian action。
  完整暫態 operator／Newton／graph 尚未接入，故保持未勾選。見
  [準靜態 graph MPI 進度](progress/HPC_03C_STATIC_GRAPH_PROGRESS.md) 與
  [分散式 history 進度](progress/HPC_03C_DISTRIBUTED_HISTORY_PROGRESS.md)、
  [暫態體積組裝進度](progress/HPC_03C_TRANSIENT_VOLUME_PROGRESS.md)。
- [ ] **HPC-03D：移動幾何。** 支援 active set 改變後的 ownership／halo 更新，
  保存穩定 ID、history extension 與提交／回復語義。先做固定背景分區，
  再依負載變化的量測決定是否加入動態重新分區及狀態搬移。

驗收：靜態、暫態、移動案例依序通過 1／2／4 rank 比較、介面守恆及失敗重試；
案例須有足夠元素供分配。多 rank 只能複製完整求解的實作不算完成。
跨節點驗收另記 HPC-09，不把單機 MPI 結果當作跨節點證據。

## HPC-04：線性求解器擴展性

- [ ] **HPC-04A：可配置求解器。** 為各 runtime 建立獨立 PETSc options prefix，
  保留已驗證預設；記錄有效 KSP/PC 與迭代數。明確檢查多 rank LU 所需後端，
  包含現有 1D 非線性路徑的 MUMPS 需求。
- [ ] **HPC-04B：候選預條件器。** 以現有 LU／block-Jacobi／ILU 作基準，
  評估速度壓力分塊、Schur 與適合子區塊的 AMG；顧及 pressure nullspace、
  port constraints、切割小元素與 ghost stabilization，不對整個 saddle-point
  系統盲目套用標量 AMG。
- [ ] **HPC-04C：網格與硬體擴展。** 分別測試建立成本、求解成本、記憶體及
  網格加密／rank 增加後的迭代數。僅在有穩定證據時變更預設求解策略。

驗收：相同數值門檻下完成小／中／大案例比較；報告失敗與負收益組合。
不以 Krylov 迭代較少作為唯一成功條件，也不宣稱 AMG 必然較快。

## HPC-05：耦合 checkpoint 與可靠續跑

入口：[CoupledDomainRuntime.hpp](../include/CoupledDomainRuntime.hpp)、
[FlowCheckpoint.hpp](../include/FlowCheckpoint.hpp)、
[TransportCheckpoint.hpp](../include/TransportCheckpoint.hpp)、
[I/O 待補項目](POST_PHASE_7_IO_HARDENING.md)。

- [ ] **HPC-05A：Checkpoint 狀態契約。** 盤點所有續跑必需狀態：時間、步數、
  流場與歷史、物種、0D／outlet／VCA、結構、移動幾何、donor hysteresis 等。
  對 Aitken／耦合歷史明確選擇保存或按既有演算法重建，不能猜測為空即可。
  第一版只在所有 domain 完成 accepted macro-step 後寫入。
- [ ] **HPC-05B：版本化與發布協議。** 分片先寫暫存，完成校驗及同步後，
  最後發布含版本、案例／配置身分、時間與各分片校驗資訊的完成 manifest。
  載入時拒絕截斷、損毀、缺片、不同 epoch 混用與不相容配置。
- [ ] **HPC-05C：既有 0D／1D／貼體 3D graph 續跑。** 先支援相同 rank 數，
  比較不中斷與中斷後續跑的完整歷史；測試 commit 前、寫分片中與 manifest
  發布前的中斷。整個作業終止後可由新的作業載入最後完整 checkpoint。
- [ ] **HPC-05D：新增 runtime 與重分區續跑。** 在相應 runtime 完成後加入
  浸入式／移動流場與 FSI；最後評估不同 rank 數的恢復。
  此模式需明確搬移 ownership，並處理 `.ntiga` 分區相依性，不能僅改啟動參數。

HPC-05D 可依 runtime 分批交付；其 FSI 部分在 HPC-07A/C 完成後與 HPC-07D
整合驗收。HPC-07 不以整個 HPC-05D 完成為前置，避免循環依賴。

驗收：完整 restart 與不中斷運行在既定數值容許值內一致，失敗發布不會覆蓋
最後可用結果；排程時間即將結束時能在安全時間步邊界保存。
本階段先不引入 multirate，以免同時改變時鐘與持久化語義。

## HPC-06：記憶體與輸出規模化

入口：[OneDImplicit.hpp](../solvers/one_d/include/OneDImplicit.hpp)、
[iga_navier_stokes.cpp](../solvers/cpu/src/iga_navier_stokes.cpp)、
[TemporalVtkHdf.hpp](../include/TemporalVtkHdf.hpp)、[VISUALIZATION.md](VISUALIZATION.md)。

- [ ] **HPC-06A：記憶體稽核。** 盤點 gather-to-all、gather-to-root、完整 mesh
  副本與移動幾何雙份狀態的大小、必要性及生命週期；優先移除大型案例熱路徑的
  不必要複製，保留便宜小模型的簡單實作。
  貼體小案例已有 required 元素分布與 RSS 初步證據：8 rank 的 required
  元素總次數為單 rank 的 5.23 倍，詳見 [HPC-00C 進度](progress/HPC_00C_PROGRESS.md)。
  後續需釐清節點分配與重複積分成本，保留外部 ID／輸出與資料庫相容性。
- [ ] **HPC-06B：分散式視覺化。** 先選擇與 ParaView 相容的分片輸出或
  Parallel HDF5 路徑，分別驗證格式、shared-point／cell 身分、時間索引與
  collective 規則。保留既有序列輸出作相容路徑。
- [ ] **HPC-06C：I/O 控制。** 分別設定場輸出、純量診斷與 checkpoint 頻率；
  測試 rank 增加時的檔案數、metadata 成本與寫入吞吐。
  依結果決定是否引入有限數量的 I/O aggregator，避免每 rank 每步大量小檔案。
- [ ] **HPC-06D：前處理大型案例。** 量測 mesh、spline、packer 的工作量與
  記憶體上限；優先保留串流／分塊設計。只有已證明為瓶頸時才增加其平行化。

驗收：大型場輸出不再要求 rank 0 收集完整場；輸出可被目標 ParaView 版本讀取，
數值與既有輸出一致；記錄所有 rank 的峰值而非只記 rank 0。

## HPC-07：分散式 FSI

入口：[FSI_ARCHITECTURE.md](architecture/FSI_ARCHITECTURE.md)、
[DistributedSurfaceInterface.hpp](../include/DistributedSurfaceInterface.hpp)、
[StrongFluidStructureCoupling.hpp](../solvers/cpu/include/StrongFluidStructureCoupling.hpp)。

- [ ] **HPC-07A：分散式表面介面。** 定義 owned／ghost surface nodes、
  參考面積權重與跨分區三角形貢獻；驗證全域覆蓋、唯一歸屬與物理身分。
  保留目前匹配介面，不同時導入非匹配投影。
- [ ] **HPC-07B：牽引力與結構更新。** 分散式計算／傳遞 nodal force、traction、
  位移與速度；用全域合力、力矩及離散功／能量一致性檢查防止遺漏或重算。
  可先由單一 owner 求解小型膜再分送結果，但須明確記錄其串行成本；
  若膜成為瓶頸，再加入分散式結構矩陣與求解。
- [ ] **HPC-07C：全域強耦合。** Aitken 內積、加權 RMS／最大殘差與收斂決策
  使用全域量；所有 rank 同步接受、回復或拒絕同一次 trial。
- [ ] **HPC-07D：驗證與續跑。** 比較單／多 rank 的位移、速度、牽引力、
  流量、守恆及收斂歷史，測試局部失敗、rollback／retry 與 FSI checkpoint。

驗收：同一 FSI 問題確實分散計算，單機多 rank 與跨節點結果通過數值門檻；
清楚區分「流體已平行、膜集中求解」與「流體／結構皆分散」的完成範圍。
不要求浮點 reduction 在不同配置逐位元一致。

## HPC-08：多域資源配置

- [ ] **HPC-08A：程序群配置。** 為大型 3D domain 配置 communicator；
  小型 0D／1D 可指定 owner。設計跨群 port 資料交換及全域步驟協調，
  不將原有的 `PETSC_COMM_WORLD` 呼叫遺留在子域內。
- [ ] **HPC-08B：合法並行排程。** 依 graph 依賴判斷哪些 domain 可同時求解，
  保留需要前序邊界結果的順序；壓力／物種方向反轉及 staged transport 均需驗證。
  不以平行化為由暗中改成不同 coupling scheme。
- [ ] **HPC-08C：實測資源分配。** 比較單一 communicator 與 domain 分群的
  通訊、閒置、記憶體及端到端時間；保留較簡單配置作為小案例預設。

驗收：至少一個具有兩個足量 3D 區域的案例證明程序分組與 port 交換正確，
不同分組符合數值門檻；checkpoint 包含必要資源映射資訊或明確重建策略。

## HPC-09：部署、CI 與最終擴展性驗收

- [ ] **HPC-09A：可重現建置。** 延伸既有 Makefile／依賴檢查，提供工作站與
  叢集配置；記錄 MPI、PETSc、HDF5 的 ABI／能力及 CUDA 目標架構。
  不把新增另一套 build system 當作必要前提。
- [ ] **HPC-09B：分層測試。** 快速單元測試、小型 1／2／4-rank 數值測試、
  有硬體的 GPU 回歸、排程執行的跨節點與大型測試分開。
  測試有 timeout、非零退出與機器可讀結果；跳過必須有原因。
- [ ] **HPC-09C：排程與執行配置。** 提供匹配 `ntasks`、`cpus-per-task`、
  OpenMP／BLAS threads、CPU／NUMA 綁定的範例；整合輸入 staging、
  安全時間步 checkpoint 及重提交流程。依實際站點政策驗證。
- [ ] **HPC-09D：Strong／weak scaling。** 固定問題增加核心測 strong scaling；
  隨核心數增加工作量測 weak scaling。報告速度比、平行效率、每 rank 工作量、
  solver iterations、記憶體及通訊。不要使用只有兩個元素的 Phase 9 fixture
  宣稱四 rank 或大規模擴展性。
- [ ] **HPC-09E：相容性與文件收尾。** 更新 root README 的能力矩陣、限制與
  可重現命令；保留既有案例、cache/text packing 與 CPU/CUDA 數值介面。
  如尚有階段或硬體測試未完成，明確列出，不將本清單整體標為完成。

驗收：以相同 source revision 在工作站與實際多節點 allocation 完成宣告範圍的
案例，保留環境／job ID／數值／效能證據；清楚標示沒有收益的配置。

## 每個子任務的完成紀錄

報告至少包含：

1. 任務編號、實作範圍、變更的資料／數值／檔案介面。
2. 基準 revision、輸入身分、硬體與軟體、MPI/thread 配置及命令。
3. 預先設定的驗收標準、實際結果、退出碼與收斂原因。
4. 組裝與求解分開的時間、通訊／I/O、host RSS；涉及 CUDA 時另列 peak allocation。
5. 單／多 rank、thread 或 CPU/CUDA 的場相對 L2 與守恆比較。
6. 未執行、未通過或受環境阻擋項目，以及後續工作的前置條件。

僅文件、配置或介面稽核任務可將不適用的數值／效能欄位記為 N/A，
但須說明原因。不可把舊階段的通過紀錄寫成新 revision 的實測結果。

後續 goal 可使用以下報告範本。每次接續時更新「剩餘工作與下一步」，
讓下一次開發能直接定位尚未完成的實作或驗收；任務狀態仍以本清單為準。

```markdown
# HPC-XXY 開發與驗收紀錄

- 狀態：未開始／進行中／待硬體驗收／已完成
- 日期與 revision：
- 任務範圍及前置條件：
- 相關程式與介面變更：
- 驗收案例、輸入身分及預先設定的容許值：
- 執行環境：CPU／GPU、MPI／PETSc、rank／thread、綁定方式
- 建置與測試命令：
- 實際結果：退出碼、數值比較、守恆、效能與記憶體
- 證據位置：日誌、機器可讀結果及產生方式
- 未執行或未通過項目與原因：
- 剩餘工作與下一步：
```

## 另立目標的後續工作

以下不屬於本清單的必要完成條件，避免同時擴大硬體與物理模型範圍：

- 多 GPU／跨節點 GPU 與 GPU FSI。
- Multirate 耦合時鐘與自適應 macro-step。
- 非匹配 FSI、monolithic FSI、ALE／remeshing、瓣膜接觸。
- 完整閉迴路心臟、0D species、臨床參數校準。
- 只有在量測需要時才做的 mesh／packer MPI 重寫。

## 參考資料

- [專案既有架構與完成範圍](progress/FINAL_ROADMAP_REPORT.md)
- [耦合架構](architecture/COUPLING_ARCHITECTURE.md)
- [移動流場架構](architecture/MOVING_DOMAIN_ARCHITECTURE.md)
- [PETSc 效能分析](https://petsc.org/release/manual/profiling/)
- [PETSc 線性求解器與預條件器](https://petsc.org/main/manual/ksp/)
- [PETSc 執行緒說明（3.15 文件；開發時須核對實際版本）](https://www.mcs.anl.gov/petsc/petsc-3.15/docs/miscellaneous/threads.html)
- [Parallel HDF5](https://portal.hdfgroup.org/documentation/hdf5/latest/_intro_par_h_d_f5.html)
