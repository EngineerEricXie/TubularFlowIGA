# HPC-01C：貼體 3D port 查詢與 adapter 求解準備

日期：2026-09-08。狀態：**以下查詢路徑通過驗收；HPC-01C 整體仍部分完成。**

## 變更與呼叫契約

`TransientFlowRuntime::MeasurePorts` 現在分別協調本地準備、表面積分與結果
建立的例外。MPI scatter／reduction 保留在本地 callback 外；積分使用
`PetscReadArray`，使本地積分失敗時仍能歸還借用的 PETSc array。

在 scatter 與 reduction 前，比對有長度分隔的完整查詢描述：port ID、解析後
的 boundary label、順序、方向、物種名稱與順序、時間，以及有效的 advection／
diffusion 係數。不同 port／species 數量或不同積分意義會共同拒絕，避免相同
MPI reduction 合併了不同量。另檢查 MPI count 容量。分散式 species state
本來就可不同，沒有要求各 rank 的本地 state 相同。

VCA overload 的 metadata 準備與結果轉換、runtime `GetPortState` 的 phase
驗證與結果複製，也接上相同協議。積分公式、SI 單位、方向與既有格式不變。

`ThreeDBodyFittedFlowDomainAdapter` 的 port lookup、輸入驗證與求解準備增加
協調。輸入先更新候選 map，整群成功後才 swap；某 rank 拒絕輸入時，其他
rank 不會留下該次已發布的輸入。準備 callback 內只有本地 boundary configuration
與 override 更新，原生 `SolveTrial` 在 callback 外執行。

staged flow/transport adapter 的 hydraulic 準備先建立 pending inputs，再於
callback 外呼叫已協調的 flow adapter。hydraulic／transport getter 的本地
前置驗證、logical port 準備與物種結果轉換也增加協調。

這些 3D 操作要求 communicator 的每個成員以相同順序參與；3D port 測量原本
就包含 collective，不能只由 root 呼叫。新加入協調的 adapter `SetPortInput`
也遵守此規則，且不可再包進同一 group 的 local-work callback。
`Communicator()` 回傳 runtime 已借用的 handle，沒有複製或轉移生命週期責任。
這與 1D adapter 的本地 getter 契約不同。

## 驗收與限制

新增 [test_three_d_port_failure.cpp](../../solvers/coupling/tests/test_three_d_port_failure.cpp)，
重用既有 bifurcation 的 64-control-point、單 cubic element fixture。測試直接
呼叫實際 runtime，沒有在 production code 新增故障開關。

| 案例 | 驗收 |
|---|---|
| 5 個準備錯誤 | 單一 rank 的非有限時間、空 port、非法 locator、錯誤 species buffer 長度、compiled fields 不符，全部成員得到相同階段診斷 |
| 10 個查詢差異 | port 數量／順序／ID／方向／label、時間、species 名稱／數量、advection 係數與交叉 diffusion 係數不符，在 collective 積分前拒絕 |
| 積分途中錯誤 | 只修改 element owner 的 connectivity，使真實積分在取得 PETSc read view 後拋出 lookup 例外；恢復 connectivity 後可重查 |
| 零面積 | 所有 rank 查詢不存在的 boundary label，reduction 後共同拒絕 |
| adapter 錯誤 | 單 rank 輸入時間錯誤、缺必要輸入、單 rank 未知 port lookup，皆共同拒絕 |
| 輸入發布 | 被拒絕的候選輸入不發布；後續 missing-input 診斷的第一失敗者仍是 rank 0，包含原本輸入合法的 rank |
| 每次測量重試 | 修正後的 area、pressure、flow、各 species concentration／flux 均符合原始參考 |
| 解析場 | 單位立方體、常速度 `(1,2,3)`、壓力 `5+x+2y+3z`、常濃度 `2`／`3`；逐個 port 比較解析 area、flow、pressure 與物種通量 |
| VCA 路徑 | VCA flow／pressure 與相同 label 的 generic measurement 一致 |
| 非零暫態求解 | 使用 packed face label 的 wall trace；明確要求 inlet flow 非零，Newton 收斂、mass imbalance 門檻通過，rollback／重試及提交成功 |
| 程序群與串行參考 | 3-rank group 與獨立 1+2 groups；每群非零流場 port 結果比較獨立 `COMM_SELF` 求解 |

以上共 **20 種模式**。3-rank 與 2-rank 群組各執行 20 個；單 rank 略過 10 個
跨 rank 差異，執行 10 個，合計 50 個案例，不重複計算各 rank。
單 rank 與兩 rank group 執行不同數量的 collective，可暴露意外的 world 操作。
所有量各自比較，容許值沿用 relative `1e-6`，零參考 absolute `1e-12`，
不混合不同物種或物理量的 norm。暫態求解使用 Newton relative `1e-8`、
absolute RMS `1e-12` 與 relative mass imbalance `1e-6`。

既有 `multidomain_subcommunicator_test` 也重跑 flow 與 species graph。兩者都在
獨立 1+2 groups 通過；最大 relative L2 分別為 `7.45082e-15` 與
`3.35727e-15`，各自保留原有守恆與接受迭代檢查。這也驗證 staged adapter
修改後的正常 hydraulic／transport 路徑與 logical species 映射。

建置三個 coupling CLI、新測試、graph subgroup test、CPU `iga_navier_stokes`
及 VCA runtime／smoke test 均成功，沒有編譯器 warning。Open MPI compiler wrapper 在 sandbox 內印出的 socket probe 訊息
不影響編譯；MPI runtime 驗收在可使用本地 sockets 的環境執行。

原有 `vca_3d_runtime_test` 使用預設 PETSc 選項通過，涵蓋既有 adapter／staged
trial、query、rollback 與 transport 狀態斷言。原生 `vca_3d_smoke_test` 使用
MUMPS 通過：1／2 ranks 的 flow／pressure 比較、VCA reservoir 守恆、同 rank
checkpoint／restart 的 flow 與 species state 逐位元相等，以及 oxygenator 案例。
此驗收沒有改動原有測試斷言或門檻。

額外將 runtime unit 強制改用 PREONLY／MUMPS 時，原有第 879 行的 outflow
重新求解逐位元相等斷言失敗。獨立診斷副本保留該斷言，只在失敗前輸出差異：
relative L2 `2.20687e-12`、最大 absolute `1.2695e-11`。這筆額外配置保持失敗
紀錄，沒有換用容許值使斷言通過；其差異成因與各 solver 下的精確相等期待
仍需另外釐清。本報告只宣告原有預設配置的 runtime unit 通過。

初版測試誤用只根據 VTK corner labels 的 wall-trace overload，導致 adapter
測試成為零解；新增非零流量 gate 後確認失敗，再改用與 production graph 相同的
`WallTraceBasis(database, mesh, 0)`。初版結果不作為非零流動驗收。中間 build
缺少 `BoundarySupport.hpp` 的錯誤亦已修正。保留中間紀錄以區分測試版本。

## 可重現命令與證據

環境：WSL 工作站、GCC 11.4、Open MPI 4.1.2、PETSc 3.15.5 real/double、
32-bit PetscInt、MUMPS；精確工具與來源身分見 inventory。基準 HEAD 為
`ee8da2ab528849814b6d12c8186476b5f7f59124` 加未提交工作樹。

```bash
make -C solvers/coupling petsc three_d_port_failure_test multidomain_subcommunicator_test \
  PETSC_DIR=/usr/lib/petscdir/petsc3.15/x86_64-linux-gnu-real -j2
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 \
  PETSC_OPTIONS='-ksp_type preonly -pc_type lu -pc_factor_mat_solver_type mumps' \
  timeout 180 mpiexec --oversubscribe -np 3 \
  solvers/coupling/three_d_port_failure_test /tmp/three-d-port-fresh-output
```

輸出目錄必須尚不存在。graph 測試另以相同環境執行
`multidomain_subcommunicator_test FRESH_OUTPUT`，species 版本再加 `species`。

本地 ignored evidence 位於 `outputs/hpc01/three-d-port/`：

- `final/faults-summary.json`、`final/faults.log`：最終 20-mode 測試、命令、
  timeout、退出碼、source／binary／log SHA256；launcher 正常退出 0。
- `graph-summary.json`、`graph-flow.log`、`graph-species.log`：graph 回歸與
  實際 binary hash，兩個作業皆正常退出 0。
- `hpc-three-d-port-final-build.log`、`hpc-three-d-port-graph-build.log`：建置紀錄。
- `vca-final/summary.json` 與各 log：原有 VCA runtime 的預設求解器回歸，及
  原生 CPU 1／2-rank VCA smoke／restart，兩個作業退出 0。
- `vca/summary.json`：強制 MUMPS 的額外 runtime unit 失敗，保留其逐位元
  斷言拒絕紀錄，不算作通過；原生 smoke 的 MUMPS 結果另見 `vca-final/`。
- `vca-mumps-diagnostic/summary.json`：保留原斷言的診斷副本所量測的差異，
  副本 source／binary 位於 `/tmp/hpc-three-d-vca-mumps-diagnostic*`，身分已記錄。
- `inventory.json`：來源、環境、測試執行檔與摘要身分。
- `faults-1.log` 為零流場版本；`faults-2.log` 是非零 gate 確認失敗的版本；
  `faults-3.log` 已修正但僅含 18-mode 測試。權威結果以 `final/` 為準。

這是數值／失敗處理驗收，沒有新效能、OpenMP 加速、GPU 或跨節點證據。
測試輸出的 assembly／linear 時間僅為診斷，未收集此版本各 rank RSS 或隔離
重複效能樣本；相關效能門檻仍由 HPC-00C/D 追蹤。

## 後續邊界

本報告沒有宣告整個 3D runtime／adapter 已完成一致失敗處理。以下仍待實作或
獨立故障驗收：

- flow `BeginStep`／`SolveTrial` 的原生前置與 snapshot／restore 邊界，constructor
  後段、PETSc collective 的回傳碼，以及 reference flow／summary 等其他 reduction。
- `GatherRequiredVelocity`／transport `GatherRequiredState` 的本地借用與配置失敗；
  transport mass／source reduction 與原生 trial／rollback／I/O 邊界。
- staged adapter 的其餘輸入、transport solve、commit／abort、accounting 與 rollback。
  本次 staged 正常 graph 回歸不等於這些故障路徑已驗證。
- executor bookkeeping、其他 CLI 與外部 mesh／SWC／profile 資產的一致性。

本地例外協調以所有 MPI 成員仍存活且依序參與為前提，不保證 MPI 內部故障或
程序遭作業系統終止後可原地恢復。持久化續跑仍依 HPC-05 驗收。

後續 trial／快照回復、required-state 查詢與 outlet evaluation 的補強及新版
回歸結果，見 [3D trial 進度](HPC_01C_THREE_D_TRIAL_PROGRESS.md)。本文件保留
port 階段當時的驗收與限制。
