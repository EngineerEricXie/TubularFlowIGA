# HPC-03A：分散式浸入式組裝層

日期：2026-09-09。基準 HEAD `0b9df67` 加本批工作樹。
狀態：owned rows／halo／物理 operator 組裝已驗證，HPC-03A 部分完成。

本頁保留 `ed129fd` 的首批組裝層證據。後續正式物理 orchestration、共用 setup
與全域診斷的實作／驗收見 [static operator 進度](HPC_03A_STATIC_OPERATOR_PROGRESS.md)。

## 實作與資料契約

[ImmersedDistributedAssembly.hpp](../../solvers/cpu/include/ImmersedDistributedAssembly.hpp)
持有 distributed PETSc AIJ matrix、state／residual vectors，以及 required-state
VecScatter。Caller 提供同群一致的 row offsets 與 stencil catalog；建構先驗證
範圍、唯一 ID／rows、owner、pattern，再 exact agreement，之後才建立 PETSc
objects。Communicator 借用至明確 Close 或析構結束，不使用 WORLD 替代 caller。

- 每個 stencil 只有一個 owner。數值 callback 只在 owner 執行，以 ADD_VALUES
  將跨 owned-row 貢獻交給 PETSc assembly。沒有在每 rank 複製整個數值組裝。
- Cell volume／mixed trace／Nitsche 採 dense field stencil；ghost penalty 只含
  same-field coupling；gauge／port 採最後一列 scalar 與其餘 rows 的雙向連接。
  Caller 將 gauge／port 依 cell 分片，scalar owner 不需為它收集完整場。
- 每個 rank 只為 owned rows 建立 sorted／unique columns，分別計數 diagonal／
  off-diagonal preallocation。所有 row 的零對角線先保留為結構，包含 saddle-point
  scalar rows；數值組裝不增加 regularization。固定 graph 與 allocation 檢查啟用。
- Required rows 是本地 stencils 使用的 union，保留 global row IDs；scatter 將
  遠端 state 帶到局部連續 buffer。callback 用 StateAt 查找，不能把 halo offset
  當作 global row。空 owned rows／空 stencils／零長度 halo 都可參與 collective。
- Callback 只做本地物理計算。例外先同群協調，再共同 drain Mat／Vec stashes；
  部分完成的 operator 不能由 Matrix／Residual 取得。下一次 Assemble 重置 matrix／
  residual 後重做，state 維持原值。Callback 非有限值、錯誤尺寸或不符 pattern
  的非零項明確拒絕。Close 依固定順序釋放後協調，重複 Close 不重做釋放。

Symbolic catalog 目前仍由每 rank 掃描，建構時會有完整 signature；persistent
numerical storage 為 owned matrix／vectors 與本地 halo。這不是大型 symbolic
拓撲或 geometry catalog 記憶體已規模化的宣稱，後續由 HPC-06 量測。
MPI 內部失聯、部分成功的 PETSc collective object creation、程序死亡不在
可恢復錯誤契約內。此層不建立 time step／checkpoint 的 commit 契約。

## 驗收

[基礎測試](../../solvers/cpu/tests/test_immersed_distributed_assembly.cpp) 使用
63 個 dense／same-field／scalar stencils、70 rows：1／2／4 ranks 及獨立 1+2
群組均對照獨立 global sum；殘差與 MatMult 一致，absolute gate `1e-13`。
每個 communicator 都另測全部 work／rows 在 rank 0 的配置，其餘 rank 完全空。
總 symbolic entries 恆為 1,008，`nz_allocated == nz_used`、`mallocs == 0`。

故障包括：已有 stencil 插入後的單 rank exception、不合法 field coupling、
單 rank 重複 stencil ID、row 越界與有效但不同的 topology。所有同群診斷一致；
integration 故障後的 operator accessor 拒絕，重試再次符合原始數值。另驗證
重複 Close、closed-state 拒絕，以及 1+2 群組互不混用 communicator。

[物理測試](../../solvers/cpu/tests/test_immersed_distributed_physics.cpp) 使用
3×3×3 背景、邊界位於 0.1／0.9 的切割 cube、compact depth-2 quadrature，
含 27 cells、54 ghost faces、非零 body force、非零四場狀態、Nitsche side walls、
兩個總流量相反的 flow-rate ports，以及 pressure gauge。共 867 rows。

分散式 callback 直接呼叫既有 element／wall／ghost／port kernels，且依 cell
積分 gauge；不從序列 runtime 取得局部 block。序列 reference 使用原本
ImmersedStaticFlowRuntime::Assemble。128 個物理 stencils 的全域使用次數
逐個為 1。對照完整 negative residual 及 Jacobian action，包含 scalar rows，
沿用 CPU relative L2 `1e-6`：

| Ranks | Residual relative L2 | Jacobian action relative L2 | 全域結構項目 |
|---:|---:|---:|---:|
| 1 | 1.61922e-16 | 2.19980e-16 | 477363 |
| 2 | 2.31969e-16 | 7.27956e-16 | 477363 |
| 4 | 2.80432e-16 | 7.06624e-16 | 477363 |

每 rank 的精確 preallocation 都通過，沒有新配置。物理與合成測試沒有更動
既有 runtime 或 case／database／checkpoint 格式；既有入口的 MPI size 1
限制維持，沒有以 wrapper 重複執行序列 solver 宣稱已完成分散式 runtime。

## 單次量測與重現

本機 GCC 11.4、Open MPI 4.1.2、PETSc 3.15.5 real64/int32；OMP／BLAS 1。
兩個新 targets 建置無 compiler warnings。正式驗收共 7 個 MPI 作業、17 份
rank reports，全部退出 0、無 timeout，保留 stdout／stderr／argv／binary hashes。

| Ranks | 各 rank owned rows | 各 rank halo rows | 組裝最大 s | 各 rank peak RSS bytes |
|---:|---|---|---:|---|
| 1 | 867 | 867 | 8.204642 | 66387968 |
| 2 | 432、435 | 838、691 | 4.122970 | 78024704、69378048 |
| 4 | 216、216、216、219 | 714、774、690、483 | 2.327930 | 73973760、73474048、70721536、65798144 |

組裝時間只包 required-state scatter、owned physics integration 及 Mat／Vec
assembly；本批沒有線性 solve 或 nonlinear solve 時間。RSS 包含每個 validation
程序的序列 reference 與 replicated geometry，不能當作 production runtime 的
記憶體 scaling。單次正確性作業不是 dedicated strong-scaling benchmark。

```bash
make -C solvers/cpu immersed_distributed_assembly_test immersed_distributed_physics_test \
  PETSC_DIR=/usr/lib/petscdir/petsc3.15/x86_64-linux-gnu-real
python3 scripts/hpc_immersed_assembly_regression.py \
  --output-dir /path/to/new-results --ranks 1 2 4 --split
```

可用 `--kind unit` 或 `--kind physics` 分開跑，控制器建立父目錄並使用
`hpc_rank_run.py --expected-ranks N` 保存逐 rank 證據。全部是工作站小型 MPI
驗收，沒有跨節點／scheduler 結果。Ignored evidence 在
`outputs/hpc03/assembly/`：正式來源為 `matrix/summary.json`（1／2 ranks 與
split1+2）、`unit-four/summary.json`、`physics-4-final/rank-*/run.json`。
`acceptance.json` 統一列出本批來源及 hashes。最早兩次四 rank wrapper 參數／
父目錄錯誤未啟動物理 executable，其 logs 保留，沒有算成物理失敗或通過。

## 接續工作

1. 將目前在物理驗收中建立的 cell／face／constraint topology 與局部 block
   orchestration 接入正式靜態 runtime，保持外部 active-node IDs 及 port ordering。
2. 改為 distributed KSP／owned updates，將 line search、field block norms、
   pressure defect、port flow／traction、守恆與收斂決策改為適当的全域量；驗證
   真正 1／2／4-rank nonlinear solve、rollback／retry 後，才完成 HPC-03A。
3. HPC-03B 增加積分工作量權重及比較；HPC-03C/D 依序加入 previous state、
   graph trial、moving geometry ownership／history。這些都仍是原清單必要工作。
