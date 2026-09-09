# MPI 錯誤邊界與剩餘工作

本文件是 HPC-01C 的入口稽核索引，最近更新於 2026-09-09。
「已接入」表示可定位協調程式與相應測試，不表示該入口的所有錯誤均已驗收。
完成仍須核對本文件的剩餘項目與 [主清單](../WORKSTATION_HPC_TODO.md)。

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

| 入口 | 已定位的邊界與協議 | 驗收索引 | 尚需核對 |
|---|---|---|---|
| [CPU flow CLI](../../solvers/cpu/src/iga_navier_stokes.cpp) | 參數、asset／configuration、boundary、runtime 建構、initialization、step、VCA、checkpoint、field／VTKHDF／completion | [輸入](../progress/HPC_01C_FLOW_INPUT_PROGRESS.md)、[步進](../progress/HPC_01C_FLOW_STEP_PROGRESS.md)、[輸出](../progress/HPC_01C_FLOW_OUTPUT_PROGRESS.md) | F02；將所有呼叫／清理與現有 stage 逐項對照 |
| [Configured transport CLI](../../solvers/cpu/src/iga_solve.cpp) | controls、assets、velocity series、runtime、step、checkpoint、field、VTKHDF init／close | [CLI](../progress/HPC_01C_TRANSPORT_CLI_PROGRESS.md)、[初始化](../progress/HPC_01C_TRANSPORT_VISUALIZATION_INIT_PROGRESS.md) | F02；共用 runtime 最終覆蓋 |
| [Legacy transport](../../solvers/cpu/src/iga_transport.cpp) | input snapshots、effective PETSc options、本地組裝、返回錯誤、資源清理與 completion | [legacy 驗收](../progress/HPC_01C_LEGACY_TRANSPORT_PROGRESS.md) | F02；與 configured 路徑分開核對 |
| [Mesh check](../../solvers/cpu/src/iga_mesh_check.cpp) | database／owner 檢查、owned element 讀取、geometry reduction、result logging | [工具 asset](../progress/HPC_01C_TOOL_ASSET_PROGRESS.md) | F02 的目的串流狀態；空 rank 保持合法 |
| [Assembly smoke](../../solvers/cpu/src/iga_assembly_smoke.cpp) | input、PETSc options、owned-row 建構／插入、assembly、matrix info／destroy、logging | [工具驗收](../progress/HPC_01CD_TOOLS_PROGRESS.md) | F02；已返回的 PETSc 錯誤與內部故障分開 |
| [原生 1D CLI](../../solvers/one_d/src/iga_1d.cpp) | input／assets、runtime、implicit solve、initial／step／final output、checkpoint | [CLI](../progress/HPC_01C_ONE_D_CLI_PROGRESS.md)、[checkpoint](../progress/HPC_01C_ONE_D_CHECKPOINT_PROGRESS.md)、[本輪串流](../progress/HPC_01C_STREAM_BOUNDARY_PROGRESS.md) | F02；新舊 checkpoint 發布仍不是原子 bundle |
| [Multidomain／bifurcation CLI](../../solvers/coupling/src/iga_1d_3d_bifurcation.cpp) | 兩個 executable 進入相同 borrowed-communicator runner；catalog、registry、flow／species executor、history／manifest | [registry](../progress/HPC_01C_REGISTRY_PROGRESS.md)、[pressure executor](../progress/HPC_01C_PRESSURE_EXECUTOR_PROGRESS.md)、[species executor](../progress/HPC_01C_SPECIES_EXECUTOR_PROGRESS.md) | F01／F02；程序群排程另由 HPC-08 完成 |
| [Sequential CLI](../../solvers/coupling/src/iga_1d_3d_explicit.cpp) | borrowed communicator、input／runtime、fixed／Aitken loop、global convergence、全部 abort outcome、輸出／completion | [初始化](../progress/HPC_01C_SEQUENTIAL_INITIALIZATION_PROGRESS.md)、[strong](../progress/HPC_01C_SEQUENTIAL_STRONG_PROGRESS.md)、[輸出](../progress/HPC_01C_SEQUENTIAL_OUTPUT_PROGRESS.md) | F02；完整失敗清理的最終整合核對 |
| [Immersed graph loader](../../include/ImmersedFlowCase.hpp) | 明確要求 MPI size 1；現有 owner／executor 協調本地執行結果 | [建構](../progress/HPC_01C_CONSTRUCTION_PROGRESS.md) | 不視為已支援多 rank；分散式資料／求解由 HPC-03 完成 |
| [FSI ParaView exporter](../../solvers/cpu/src/phase8_compliant_channel_fsi_paraview.cpp) | 在輸出路徑變更前拒絕多 rank；單 rank 執行流固耦合與 exporter | [入口驗收](../progress/HPC_01C_EXPORTER_ENTRY_PROGRESS.md) | F04 的 reader／array serializer；分散式 FSI 由 HPC-07 完成 |
| [CUDA CLI](../../solvers/cuda/src/iga_cuda.cu) | 單程序、單 GPU，既有主例外處理與共用 writer | [文字 helper](../progress/HPC_01C_TEXT_HELPER_PROGRESS.md) | F02／F04 的相應串行路徑；不能把 mpiexec 重複啟動當作 MPI 支援 |
| [Mesh](../../preprocessing/mesh/src/main.cpp)／[spline](../../preprocessing/spline/main.cpp) | 串行入口／可選 OpenMP；無跨 rank 協調契約 | 主清單 HPC-00／02 的相應案例 | 本地失敗／OpenMP worker 傳遞另核對；大型前處理由 HPC-06D 驗收 |
| Packer／inspect／case-check／config-check／flow-validate／transport-validate／Womersley 工具 | 串行 CLI，資料格式相容是 pipeline 的必要條件 | CPU README 與各工具測試 | F04；逐項確認完整讀取、flush／close 與成功回報 |

## 共用 runtime 與 adapter

| 元件 | 本輪核對結果 | 限制 |
|---|---|---|
| `TransientFlowRuntime` | initialization、BeginStep、MeasurePorts 的簽章在本地準備階段啟用 stream exceptions；下一步才做 agreement。非收斂診斷的數值字串亦啟用檢查 | 完整 Newton／PETSc 呼叫證據見既有報告；本輪新增三個原生格式化故障入口 |
| `TransientTransportRuntime` | GatherState 與 IntegrateFields 在完整簽章建立後才 agreement／gather／reduction | 本輪驗證 gather／mass，source preparation 使用同一受檢查的 stream；不可把它當成每一 source term 的新故障驗收 |
| `ThreeDBodyFittedFlowTransportDomainAdapter` | PlanSignature／InputSignature 受檢查；BeginStep、mixed input、staged concentration 的 candidate 在共同成功後才發布 | 新增三種入口故障；仍沿用既有 hydraulic／transport trial 及 rollback 驗證 |
| `OneDRuntime`／1D adapter | stage outcome 與 implicit runtime 內部協調分開；CLI output 使用共同本地階段 | 本輪 filename 故障不代表可以重用部分寫出的同一 writer；健康重試使用新作業及新目錄 |
| Registry／pressure／species executor | storage 建立、runtime 所有權、stage outcome、global convergence、abort 順序已有測試 | 兩個 abort 診斷 formatter 已驗證 F01；不能僅依外層 wrapper 宣告安全 |

## 可執行的剩餘項目

以下是已定位的缺口，不是完整性證明。修正後應更新此表，避免每次只從一個
grep 結果重新開始；新增入口也必須補進上表。

| ID | 已觀察的缺口 | 下一步與完成證據 |
|---|---|---|
| F01（已驗證） | 兩個 AbortAll formatter 已啟用 fail/bad exceptions，所有 abort outcome 返回後才格式化診斷 | 3-rank world 與 split 1／2 ranks 的 24 次 allocation 故障、30 次健康／重試通過，primary error 與全 domain abort 保留；見 [本輪報告](../progress/HPC_01C_SERIAL_TOOL_PROGRESS.md) |
| F02（已驗證） | mesh check／assembly smoke 的 result logging 與 resource summary 已 flush 並通過 84 個原生作業；CPU flow／runtime、configured／legacy transport 亦已補強，412 個作業及 310 個場比較通過；1D／graph／sequential 另通過 216 個作業與 660 個輸出比較，四個序列工具通過 48 個作業；CUDA 另通過 38 個原生作業與 27 個場比較 | [CUDA stdout 報告](../progress/HPC_01C_CUDA_STDOUT_PROGRESS.md)、[耦合與工具 stdout 報告](../progress/HPC_01C_COUPLING_STDOUT_PROGRESS.md) 與 [CPU stdout 報告](../progress/HPC_01C_SOLVER_STDOUT_PROGRESS.md) 含 prefix／stage 對應與限制。區分可控制 stream failure 與 SIGPIPE／程序死亡 |
| F03（已驗證） | Publisher 的檔名、JSON、VTU arrays、Number、Nonce 使用受檢查串流，ReadAll 使用完整讀取 | Publish／retry／collection rebuild 共 219 次故障、222 次重試通過，accepted epoch 不變且重試輸出逐位元一致；見 [本輪報告](../progress/HPC_01C_SERIAL_TOOL_PROGRESS.md) |
| F04（已驗證） | Packer／FSI reader、FSI Values、Womersley filename／close 與 budget JSON 已補強 | 21 個原生 CLI 案例、30 個輸出比較、cache/text packing 與 FSI helper 故障通過；完整 FSI 兩版均 4 次迭代收斂、5 份輸出逐位元相同及 29 個數值 arrays 相對 L2=0，見 [本輪報告](../progress/HPC_01C_SERIAL_TOOL_PROGRESS.md) |
| F05（已驗證） | `BezierVisualization` 的錯誤訊息 formatter 啟用 exceptions；packer formatter 已加 exceptions | 舊實作重現例外被吞掉；修正後三種 formatter 例外及三次健康診斷恢復通過，原有幾何／array tests 通過；見 [Bezier 診斷驗收](../progress/HPC_01C_BEZIER_DIAGNOSTIC_PROGRESS.md)。這原本就是拒絕路徑，不是靜默數值成功 |
| F06 | 入口／constructor／destructor／早退分支的完整呼叫圖與最後一輪整合矩陣仍未逐項簽核 | 已完成終止呼叫鏈核對與 runtime Close 補強：144 個 runtime 故障／重試及最終 112 個原生作業通過，見 [入口稽核](../progress/HPC_01C_ENTRY_AUDIT.md) 與 [cleanup 驗收](../progress/HPC_01C_RUNTIME_CLEANUP_PROGRESS.md)。其餘呼叫鏈與最後整合覆蓋仍需簽核 |

HPC-01D 的型別／配置／後端能力矩陣、HPC-05 的持久化契約與 HPC-09 的
跨節點驗收仍維持原要求。此索引不移除或縮小任何主清單任務。
