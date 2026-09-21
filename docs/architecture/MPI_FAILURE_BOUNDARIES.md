# MPI 錯誤邊界

本文件列出 production 入口的 MPI 錯誤協調邊界。「已接入」表示可定位
協調程式與相應 regression，不表示 MPI／PETSc 內部程序故障可原地恢復。

## 協調契約

- `CollectiveLocalStage` 只包本地工作；它捕捉本地例外後協調失敗 rank 與訊息。
- PETSc collective 位於本地 callback 外，返回後再由
  `RequireCollectivePetscSuccess` 或相應 helper 協調回傳錯誤。
- runtime／adapter／executor 共用借用的 communicator。同群須以相同順序進入
  相同階段；控制參數與有效輸入的簽章必須先完整建立，再做群組一致性比較。
- `Stage`／returned-outcome wrapper 不能替代被呼叫 runtime 的內部協調。
  只在 runtime 呼叫已返回後捕捉錯誤，無法救回已卡在不同 collective 的 peer。
- 這個協議以程序仍存活且 MPI 可通訊為前提。MPI 初始化、MPI 內部故障、
  SIGKILL／OOM killer 與節點失聯不以此協議宣稱可原地恢復。

## 支援入口

| 入口 | 已定位的邊界與協議 | Regression 範圍 | 支援界線 |
|---|---|---|---|
| [CPU flow CLI](../../solvers/cpu/src/iga_navier_stokes.cpp) | 參數、asset／configuration、boundary、runtime 建構、initialization、step、VCA、checkpoint、field／VTKHDF／completion | 輸入、步進、輸出及故障注入 | main、runtime/helper、明確 Close 與 completion 順序 |
| [Configured transport CLI](../../solvers/cpu/src/iga_solve.cpp) | controls、assets、velocity series、runtime、step、checkpoint、field、VTKHDF init／close | CLI、初始化、輸出及故障注入 | main、共用 helper、memory report 與 PETSc owners |
| [Legacy transport](../../solvers/cpu/src/iga_transport.cpp) | input snapshots、effective PETSc options、本地組裝、返回錯誤、資源清理與 completion | Legacy CLI regression | 與 configured 路徑分開；destructor 只作失敗後備 |
| [Mesh check](../../solvers/cpu/src/iga_mesh_check.cpp) | database／owner 檢查、owned element 讀取、geometry reduction、result logging | 工具 asset regression | 空 rank 保持合法 |
| [Assembly smoke](../../solvers/cpu/src/iga_assembly_smoke.cpp) | input、PETSc options、owned-row 建構／插入、assembly、matrix info／destroy、logging | Assembly tool regression | 已返回的 PETSc 錯誤與內部故障分開 |
| [原生 1D CLI](../../solvers/one_d/src/iga_1d.cpp) | input／assets、runtime、implicit solve、initial／step／final output、checkpoint | CLI、checkpoint 與 stream regression | 新舊 checkpoint 發布仍不是原子 bundle |
| [Multidomain／bifurcation CLI](../../solvers/coupling/src/iga_1d_3d_bifurcation.cpp) | 兩個 executable 進入相同 borrowed-communicator runner；catalog、registry、flow／species executor、history／manifest | Registry、pressure 與 species executor regression | runner/helper 與清理已核對；程序群排程由 deployment layer 處理 |
| [Sequential CLI](../../solvers/coupling/src/iga_1d_3d_explicit.cpp) | borrowed communicator、input／runtime、fixed／Aitken loop、global convergence、全部 abort outcome、輸出／completion | 初始化、strong coupling 與輸出 regression | runner/helper、明確 Close 與 publication 順序 |
| [Immersed graph loader](../../include/ImmersedFlowCase.hpp) | 本地 preflight 後 collective 初始化 steady／fixed-transient owned runtime；模式／clock agreement 與 graph rollback | Case／graph regression | serial `Load` 僅供 size 1；MPI 使用 `InitializeDistributed` |
| [FSI ParaView exporter](../../solvers/cpu/src/phase8_compliant_channel_fsi_paraview.cpp) | 在輸出路徑變更前拒絕多 rank；單 rank 執行流固耦合與 exporter | Exporter entry regression | 分散式 FSI 使用獨立 runtime |
| [CUDA CLI](../../solvers/cuda/src/iga_cuda.cu) | 單程序、單 GPU，既有主例外處理與共用 writer | CUDA CLI 與 stream regression | 支援的多程序 launcher 會明確拒絕 |
| [Mesh](../../preprocessing/mesh/src/main.cpp)／[spline](../../preprocessing/spline/main.cpp) | 串行入口／可選 OpenMP；無跨 rank 協調契約 | 本地失敗與 OpenMP worker regression | 大型前處理應在配置的 compute resource 執行 |
| Packer／inspect／case-check／config-check／flow-validate／transport-validate／Womersley 工具 | 串行 CLI，資料格式相容是 pipeline 的必要條件 | CPU README 與各工具測試 | 不宣稱為 MPI runtime |

## 共用 runtime 與 adapter

| 元件 | 本輪核對結果 | 限制 |
|---|---|---|
| `TransientFlowRuntime` | initialization、BeginStep、MeasurePorts 的簽章在本地準備階段啟用 stream exceptions；下一步才做 agreement。非收斂診斷的數值字串亦啟用檢查 | 完整 Newton／PETSc 呼叫證據見既有報告；本輪新增三個原生格式化故障入口 |
| `TransientTransportRuntime` | GatherState 與 IntegrateFields 在完整簽章建立後才 agreement／gather／reduction | 本輪驗證 gather／mass，source preparation 使用同一受檢查的 stream；不可把它當成每一 source term 的新故障驗收 |
| `ThreeDBodyFittedFlowTransportDomainAdapter` | PlanSignature／InputSignature 受檢查；BeginStep、mixed input、staged concentration 的 candidate 在共同成功後才發布 | 新增三種入口故障；仍沿用既有 hydraulic／transport trial 及 rollback 驗證 |
| `OneDRuntime`／1D adapter | stage outcome 與 implicit runtime 內部協調分開；CLI output 使用共同本地階段 | 本輪 filename 故障不代表可以重用部分寫出的同一 writer；健康重試使用新作業及新目錄 |
| Registry／pressure／species executor | storage 建立、runtime 所有權、stage outcome、global convergence、abort 順序已有測試 | 兩個 abort 診斷 formatter 已驗證 F01；不能僅依外層 wrapper 宣告安全 |
