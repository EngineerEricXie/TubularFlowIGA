# HPC-01C：初始化候選發布與全域 gather

日期：2026-09-08。狀態：**本批邊界通過驗收，HPC-01C 仍部分完成。**
接續 [Newton 報告](HPC_01C_NEWTON_PROGRESS.md)。

## 初始化

[TransientFlowRuntime.hpp](../../solvers/cpu/include/TransientFlowRuntime.hpp)
的兩種 InitializeState 現在共同進入初始化流程，要求 committed phase。
先在本地候選物件上 resolve 邊界與 pressure tractions，檢查 topology、indexing
及有限值，再精確比較 resolved 初始邊界與呼叫模式的描述。
描述涵蓋 configured／transient／是否提供配置、節點數、有效 velocity／pressure
constraint 及其值、pressure tractions；不把分散式自由場值拿來逐 rank 比較。

暫態初始化使用既有 update／rhs scratch vectors 準備 state 與 history 候選，
協調插入／assembly／copy 完成後，才使用檢查回傳碼的 VecSwap 發布。公開的
State() Vec handle 保持不變，之後才 swap 發布本地邊界候選。非法配置、跨 rank
差異及受控的本地準備失敗，不會先改動其他 rank 的有效邊界或初始場。

BoundaryValue 現在檢查負 row 與 vector indexing。原有 trial 配置更新改用共用
的本地 ResolveConfiguredBoundaries，避免初始化與時間步解析行為分叉。
穩態仍不改動場係數，但 configured 模式可更新 resolved 邊界；未配置模式仍忽略
傳入的配置。初始化沒有新增大型 PETSc vectors，但候選邊界與描述字串有暫時記憶體成本。

這是本地驗證後的候選發布，不是 PETSc 任意內部失敗下的原子交易承諾；若 backend
在發布或物件建立中途失效，仍受既有錯誤恢復能力限制。其他 constructor 操作及
完整物理配置／外部資產身分仍待補齊。

## Gather

新增 [PetscGather.hpp](../../solvers/cpu/include/PetscGather.hpp) 的 GatherAllPetscReal：
在建立 scatter 前共同核對來源 Vec 全域長度與預期 row 數、配置結果 buffer；
逐項檢查 CreateToAll／Begin／End，並協調結果長度、read view 與複製。
正常清理明確檢查 Destroy 回傳碼；RAII owner 負責共同失敗後的已建立暫存物件。
部分 PETSc creation 失敗或程序消失不在此 owner 的恢復保證內。

[TransientTransportRuntime.hpp](../../solvers/cpu/include/TransientTransportRuntime.hpp)
的 GatherState 先檢查欄位數、唯一且非空的名稱、乘法範圍、assembled row 數，
比較全群有序欄位描述，再呼叫 helper。維持全域 node／field 順序與原始數值，
包含 NaN、Infinity 與 signed zero，供診斷讀取。此方法仍在每個 rank 複製完整場；
沒有把它當作 HPC-06 的分散式記憶體改善。

## 驗收

| 測試 | 結果 |
|---|---|
| trial／初始化 | 3、1、2-rank groups 各 33、27、33 個，共 **93 個故障案例** |
| staged／global state | 各 30、23、30 個，共 **83 個故障案例** |
| Petsc gather helper | 六組 shape／communicator 配置各兩個尺寸錯誤，共 **12 個故障案例**，另含正常值、空向量、空 rank 與原始非有限值比較 |
| port／diagnostics | 原有 **64 個案例**通過 |
| flow graph 1+2-rank groups | 最大 relative L2 `7.45082e-15` |
| species graph 1+2-rank groups | 最大 relative L2 `3.35727e-15` |
| VCA runtime | 原有 PETSc 預設下 lifecycle／staged assertions 通過 |
| VCA smoke | MUMPS，1／2-rank 比較、reservoir、同 rank checkpoint/restart state 精確比較與 oxygenator 通過 |

初始化新增故障涵蓋：單 rank 不支援 profile、NaN scale、constraint topology
改變，以及合法但不同的 scale／traction／overload。每次拒絕後檢查 state 與
traction 不變，再用無參數 InitializeState 證明沒有候選邊界洩漏。trial 中呼叫
初始化也共同拒絕，phase 與快照保持原狀。合法 scale 更新後重新初始化回原配置，
確認場與原值逐位元相同、Vec handle 不變。

另驗證 configured steady、unconfigured steady、unconfigured transient 三種
模式，分別保留既有場係數／邊界更新／忽略配置的語義。原有非零 flow solve、
history、rollback／retry 及 outlet 比較仍通過。

staged 測試在成功 trial 後破壞單 rank 的 field count、重複／空名稱，或提供
合法但不同名稱；GatherState 在通訊前共同拒絕，恢復欄位後結果與原全域場精確相同。
全域結果也逐 node／field 比對 required-node 場，確認排序沒有改變。

[test_petsc_gather.cpp](../../solvers/cpu/tests/test_petsc_gather.cpp) 使用 world3、
獨立 1+2 groups，配置 `(ranks, rows)=(3,2)、(3,0)、(1,3)、(1,0)、(2,1)、(2,0)`。
包含 ranks 多於 rows 與全空向量。單 rank 預期尺寸錯誤或 UINT64_MAX 都共同
拒絕，修正後可再次 gather；NaN、Infinity、負零由同群每個成員取回且保留。
所有 array mutation 都由全群共同借用／歸還，只有 owning rank 改值。

全部八項 MPI／VCA 作業在各自 180 秒 timeout 內退出 0。錯誤訊息要求全群一致，
獨立 subgroup 的既有數值回歸使用 COMM_SELF 參考。沿用 relative `1e-6`、
零參考 absolute `1e-12`；初始化／回復與原始 gather 值採精確比較。
沒有改動求解器預設、物理公式、檔案格式或原有數值 gate。

## 重現與證據

環境：WSL 工作站、GCC 11.4、Open MPI 4.1.2、PETSc 3.15.5 real/double、
32-bit PetscInt、MUMPS。HEAD `ee8da2ab528849814b6d12c8186476b5f7f59124`
加未提交工作樹；inventory 記錄 source／binary 身分。

```bash
make -C solvers/coupling petsc three_d_trial_failure_test three_d_staged_failure_test \
  three_d_port_failure_test multidomain_subcommunicator_test \
  PETSC_DIR=/usr/lib/petscdir/petsc3.15/x86_64-linux-gnu-real -j2
make -C solvers/cpu petsc_gather_test iga_navier_stokes vca_3d_runtime_test vca_3d_smoke_test \
  PETSC_DIR=/usr/lib/petscdir/petsc3.15/x86_64-linux-gnu-real -j2
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 \
  timeout 180 mpiexec --oversubscribe -np 3 solvers/cpu/petsc_gather_test
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 \
  PETSC_OPTIONS='-ksp_type preonly -pc_type lu -pc_factor_mat_solver_type mumps' \
  timeout 180 mpiexec --oversubscribe -np 3 \
  solvers/coupling/three_d_trial_failure_test /tmp/initialization-fresh-output
```

案例目錄必須尚不存在。MPI runtime 需要本地 sockets。
最終 CPU／coupling 建置通過且沒有新 compiler warning；首輪 signedness warning
及測試將 vector boundary value 寫成 scalar 的編譯錯誤已修正，原始 logs 保留。

本地 ignored evidence 位於 `outputs/hpc01/initialization/`：

- `verified/mpi-summary.json`、logs／案例：六項 MPI 作業、命令、退出碼及 binary／log hashes。
- `vca/summary.json`、logs／案例：兩項原生回歸。
- `evidence-summary.json`、`build/`：成功標記核對、248 個證據檔 hashes 與最終 build logs。
- `build-attempts/`、`build-attempts.json`：已修正的編譯問題紀錄，不作為通過證據。
- `inventory.json`：目前工作樹、環境、執行檔與 evidence 身分。

沒有新的 performance／RSS、GPU、profiling 或跨節點驗收。未更改 production
geometry，未重跑 mesh-test。Gather 的 post-create read-view/backend 內部錯誤
未逐一注入，故障覆蓋以上述案例為限。

## 剩餘工作

HPC-01C 仍需兩種 runtime 的 constructor 後段、transport 初值建立、viewer／
checkpoint I/O、完整配置與外部資產一致性，以及其餘 adapter／executor／CLI
錯誤邊界。初始化比較只涵蓋 resolved 初始邊界；不能替代全部物理配置 preflight。
完整 checkpoint 發布／恢復契約另由 HPC-05 追蹤。整份 goal 保留全部 38 個子任務。
