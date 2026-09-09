# HPC-01C：貼體 3D staged 輸入與物種收支

日期：2026-09-08。狀態：**下列邊界已驗收，HPC-01C 與整份清單仍未完成。**
接續 [3D trial 報告](HPC_01C_THREE_D_TRIAL_PROGRESS.md)。

## 實作與行為

[staged adapter](../../include/ThreeDBodyFittedFlowTransportDomainAdapter.hpp)
在 BeginStep 開啟原生交易前，先共同驗證 transport phase、步驟與執行計畫。
計畫包含依序的 port、boundary label、方向、required/provided quantities、
物種與 logical/native bindings。初始 mass 查詢若失敗，flow 與 transport
都 Abort 回到原始快照，成功後才發布初始 mass。

SetPortInput 與 SetTransportConcentration 先建立候選 map，完成本地驗證，
再精確比較跨 rank 輸入描述，最後 swap 發布。描述包含 port、時間、各 optional
hydraulic quantity 及物種 map。某個 rank 的錯誤或合法但不同的輸入，現在會
共同拒絕；修正後能重新送入，不會留下其他 rank 已接受的半份輸入。
這兩個方法現在需要 communicator 的所有成員依相同順序呼叫。

SolveTransportTrial 的 waveform materialization 與每個 port 的 boundary
設定分別協調。hydraulic port 查詢留在 local-work callback 外；查詢完成後
才協調濃度與邊界驗證。流向反轉後只向流入端施加遠端濃度，仍保留原有的
流出端拒絕行為。原生求解也留在 callback 外。

GetSpeciesStepAccounting 先共同驗證 trial 成功，再查詢 mass、source 與各
species port。每個 port 只查詢一次，快取發布及 logical/native map 轉換均
協調後才繼續。仍以 backward Euler 的步末 flux/source 乘 dt 計算 amount；
SUPG accounting residual 仍是診斷量，沒有新增精確守恆宣告。

[TransientTransportRuntime](../../solvers/cpu/include/TransientTransportRuntime.hpp)
的兩種 TotalMass 與 SourceIntegrals 共用分階段積分流程：準備欄位與緩衝區、
比較有序的物理描述、本地積分、MPI reduction、結果驗證。欄位必須非空且唯一，
數量不得超過 MPI int 上限；source equation 與有限係數在每個 rank 檢查，
包含沒有 owned element 的 rank。結果中的 NaN 或 overflow 明確拒絕。
全域 state buffer 的長度及 indexing 也檢查；分散式場值本身不做逐 rank 相等比較。

沒有更改物理公式、求解器預設、port 單位、檔案格式或快照語義。
精確描述比較針對應複製一致的 metadata/input，不適用於各 rank 不同的局部場。

## 故障與數值驗收

[test_three_d_staged_failure.cpp](../../solvers/coupling/tests/test_three_d_staged_failure.cpp)
使用既有單 cubic element、64-control-point fixture，非零 backward Euler 流場、
兩個 native species 對應 red/blue，初值分別 1.25/3.5。測試在 compiled system
加入 source 0.05/0.10；單位立方體、dt=0.01 的解析 source amounts 為
0.0005/0.001。這些 source 是測試內配置，不是 fixture JSON 的內容。

| 範圍 | 驗收內容 |
|---|---|
| 開始步驟 | 初始 mass 積分失敗後兩個 runtime 都回到 committed，場與步數精確還原；跨 rank port 順序及 species binding 差異共同拒絕 |
| 輸入發布 | 非法時間、未知物種／port、不完整濃度及合法但不同的輸入拒絕；修正後重送成功 |
| staged 順序 | 缺少流入濃度、提前 accounting、重複 solve、成功後非法 concentration 更新被拒絕；有效的先前狀態保留 |
| 本地積分 | 壞 quadrature、重複欄位、非法 source equation、不同欄位／source 描述、錯誤 state 長度均共同拒絕 |
| 非有限結果 | state NaN、有限 source 係數累加造成 overflow，在 reduction 後拒絕結果發布 |
| 邊界設定 | 單 rank 的重複 scalar BC 在 hydraulic 查詢後失敗，Abort 還原先前 flow／scalar 快照 |
| 流向反轉 | 向新流出端施加遠端濃度被拒絕並可 Abort；只向新流入端送入濃度後成功求解及提交第二步 |
| 數值與快照 | 正向 port flow／pressure、物種場、mass、source、port amounts 與獨立串行解比較；反向比較物種場及 amounts，另檢查 flow 正負；rollback/Abort 快照逐位元相等 |

3-rank 群組與 2-rank 子群各 26 個案例，單 rank 子群跳過六個跨 rank 差異案例，
執行 20 個，共 **72 個故障案例通過**。各子群執行不同數量的 collective，並與
各自的 COMM_SELF 參考比較。錯誤診斷要求全群一致；整個作業在 180 秒 timeout
內退出 0。沿用 relative L2 `1e-6`，零參考 absolute L2 `1e-12`；flow Newton
relative `1e-8`、absolute RMS `1e-12`、mass imbalance `1e-6`。

故障注入僅在測試中短暫改動所持有的 quadrature storage／compiled metadata，
在後續求解前恢復；未增加 production 故障開關。

其餘受影響回歸的退出碼、成功標記、執行檔與 log SHA256 已逐項核對：

- trial failure：48 個案例通過；port failure：50 個案例通過。
- flow graph 的獨立 1+2-rank groups：最大 relative L2 `7.45082e-15`。
- species graph 的獨立 1+2-rank groups：最大 relative L2 `3.35727e-15`，原有收支 gate 通過。
- VCA runtime：原有 PETSc 預設下 lifecycle／staged assertions 通過。
- VCA smoke：MUMPS，1／2-rank 場比較、reservoir、同 rank 的 checkpoint/restart
  state 精確比較及 oxygenator 案例通過。

先前額外強制 MUMPS 的 runtime unit 逐位元斷言問題仍見
[3D port 報告](HPC_01C_THREE_D_PORT_PROGRESS.md)，本批沒有重跑或放寬該配置。

## 重現與證據

環境：WSL 工作站、GCC 11.4、Open MPI 4.1.2、PETSc 3.15.5 real/double、
32-bit PetscInt、MUMPS。HEAD `ee8da2ab528849814b6d12c8186476b5f7f59124`
加未提交工作樹；source／binary 身分由 inventory 與最終 summary 記錄。

```bash
make -C solvers/coupling petsc three_d_staged_failure_test three_d_trial_failure_test \
  three_d_port_failure_test multidomain_subcommunicator_test \
  PETSC_DIR=/usr/lib/petscdir/petsc3.15/x86_64-linux-gnu-real -j2
make -C solvers/cpu iga_navier_stokes vca_3d_runtime_test vca_3d_smoke_test \
  PETSC_DIR=/usr/lib/petscdir/petsc3.15/x86_64-linux-gnu-real -j2
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 \
  PETSC_OPTIONS='-ksp_type preonly -pc_type lu -pc_factor_mat_solver_type mumps' \
  timeout 180 mpiexec --oversubscribe -np 3 \
  solvers/coupling/three_d_staged_failure_test /tmp/three-d-staged-fresh-output
```

最後的案例目錄必須尚不存在。受影響建置通過且沒有 compiler warning；MPI
runtime 使用允許本地 sockets 的執行環境。沒有新的幾何實作，未重跑 mesh-test。

本地 ignored evidence 位於 `outputs/hpc01/three-d-staged/`：

- `final/faults-summary.json` 與 log：最終 72 個案例、命令、環境、timeout、退出碼、
  source／binary／log hashes；本批重新執行確認通過。
- `regression/mpi-summary.json`、`vca/summary.json` 及 logs／案例產物：上述六項回歸。
- `evidence-summary.json`：成功標記與身分核對、輸入／輸出檔案 hashes。
- `build/` 與 `inventory.json`：建置 logs、環境與工作樹身分。

先前 `attempt-1` 是 63-case 版本，`attempt-2` 已加入反轉／非有限结果並通過
72 個案例，均保留。最終版本以 `final/` 的 source／binary hashes 為準。
本批證據整理首輪誤把 port test 的 `modes=... retry=passed` 標記當成 trial
test 的 `cases=... passed`，驗證器拒絕；核對原始測試後修正標記並全部核對通過，
未改動測試或數值 gate。

這是故障與數值驗收；沒有量測新增協調／候選複製的隔離效能或 RSS，沒有新的
OpenMP、GPU 或跨節點宣告。原有 logs 中的 assembly／linear 時間僅作診斷。

## 剩餘工作

- constructor／初始化、全域 gather／checkpoint I/O，以及 flow Newton 未檢查的 PETSc 回傳碼。
- flow ReferenceBoundaryFlow 與 Summary 的本地失敗／reduction 邊界。
- 完整物理配置、外部資產身分與 embedding PETSc options preflight。
- 其餘 adapter／executor bookkeeping／CLI 邊界；本次 plan 比較不代表全部配置一致性。
- 全域持久化由 HPC-05 追蹤；MPI 內部卡住或程序消失後的恢復不由本地例外協調保證。

以上尚未完成，因此 HPC-01C 保持未勾選，完整 goal 保留全部清單範圍。

後續 ReferenceBoundaryFlow／Summary 邊界與 64 個 port／診斷案例已完成，見
[流場診斷進度](HPC_01C_FLOW_DIAGNOSTICS_PROGRESS.md)，其餘工作仍待補齊。
