# HPC-03C：暫態 graph adapter

日期：2026-09-09。基準 `0ee54f0` 加本批修改。
狀態：分散式 C++ graph adapter 已實作；case／正式 graph 暫態入口與整體耦合
驗收尚未完成，HPC-03C 保持未勾選。

## 實作

[ThreeDImmersedTransientDistributedFlowDomain.hpp](../../include/ThreeDImmersedTransientDistributedFlowDomain.hpp)
將固定幾何 transient runtime 的生命週期接到 `CoupledDomainRuntime`。共用
[ThreeDImmersedDistributedDomain.hpp](../../include/ThreeDImmersedDistributedDomain.hpp)
保存既有 port metadata、collective inputs、port measurements、prepare／finalize
與 accepted controls 的交易流程。steady adapter 保留原 C++ 名稱及 dt=0 行為；
不同 lifecycle 名稱納入 collective constructor signature。

adapter 要求初始 backend time 與 constructor time 一致且 index=0、無 active
trial；這是 fresh-run adapter，尚未宣稱 checkpoint 恢復。graph 的 step_index
從 0 起，對應 backend 第一次接受後的 index=1。BeginStep 核對雙方 accepted
clock；所有端口輸入備齊後，SolveTrial 才開始 backend trial 並凍結 previous velocity。

Newton／量測失敗、coupling rollback 與整步 abort 會呼叫 backend `AbortTrial`，
釋放 frozen history，然後恢復 accepted port controls。下一次 coupling trial
重新從相同 accepted field 開始，能重新供給完整邊界資料。backend 的 Newton
`Rollback` 只保留同一 trial 的輸入，因此不適合直接當 coupling rollback。

prepare 成功前先配置量測候選；noexcept finalize 同時發布 backend 場／時鐘與
adapter 的 port measurements／controls。prepare 失敗、取消、重複 finalize
或 backend Close 不得造成 graph metadata 單獨前進。數值 kernel、格式與
native case 的支援範圍沒有在本批改動。

## 驗證

```bash
make -C solvers/cpu immersed_distributed_domain_test \
  immersed_transient_distributed_domain_test CXX=mpicxx \
  PETSC_DIR=/usr/lib/petscdir/petsc3.15/x86_64-linux-gnu-real
python3 scripts/hpc_immersed_domain_regression.py \
  --output-dir outputs/hpc03/transient-domain/static --ranks 2 --split
python3 scripts/hpc_immersed_transient_domain_regression.py \
  --output-dir outputs/hpc03/transient-domain/accepted-source --split
make -C solvers/coupling iga_multidomain_flow CXX=mpicxx \
  PETSC_DIR=/usr/lib/petscdir/petsc3.15/x86_64-linux-gnu-real
```

TsungYehLab、GCC 11.4、OpenMPI 4.1.2、PETSc 3.15.5 real64/int32，OMP／BLAS
threads=1。小型本機 MPI correctness tests；無跨節點或 CUDA。受限編譯中的
`opal_ifinit socket errno=1` 為環境訊息，MPI runtime 以允許 socket 的權限執行。

unit cube、4×1×1 cells、112 nodes、compact quadrature，入口／出口為 balanced
flow controls。兩步 dt=0.125，流量由 ±1e-4 變成 ±0.5e-4。每一步先完成一次
trial、rollback 再重新求解後接受。對照 legacy moving owner 在 stationary
geometry 上的完整 previous velocity continuation；參考只存在測試內。

legacy setter 每次 mutation 都檢查 all-flow compatibility，首輪測試因先改入口
而拒絕中間狀態。修正 reference fixture：在建立 source epoch 時一次帶入完整
邊界設定，按 stable node IDs 重綁 accepted coefficients、controller／gauge
狀態，再解下一步。production setter 與 epoch guard 未改動。失敗證據保留在
`final/`，`reference-fixed/` 通過後再加入非零 accepted field／controls 的精確回復檢查，
最終有效矩陣為 `accepted-source/`。

預設 gate 為 velocity／pressure／scalar block relative L2 ≤ 1e-6；reference
norm < 1e-10 時檢查 absolute L2 < 1e-10。port area／flow／pressure／traction
使用 1e-12 + 1e-6 |reference|，global surface flow／volume divergence 差 ≤ 1e-12。
clock／index／commit count、controls 回復與 unpublished field 必須精確符合。

負面測試涵蓋 metadata／rank agreement、錯誤初始 clock、step mismatch、
非有限或 rank 間不同輸入、缺漏／重複輸入、force freeze failure、線性更新後
candidate failure、prepare failure、coupling rollback、prepared abort、重試、
雙重 finalize 及 Close 後禁止發布。steady adapter 的 2 ranks 與 split 1+2
另外回歸，並重建 native graph executable。

## 證據與後續

最終 transient 4 jobs／10 rank reports 全部 exit 0、無 timeout，最大 field
relative L2 3.709381e-14，低於 1e-6。最大 rank wall time 6.927811 s、peak RSS
60,989,440 bytes；runtime 累計 assembly 最大 4.535736 s，linear solve 最大
0.191809 s。第二步非零 accepted field／controls 的 rollback 精確保留檢查通過。

steady 2 jobs／5 reports 也通過，最大 field relative L2 6.599167e-10。
本批最終共 6 MPI jobs／15 rank reports；native graph 重建成功，無新增 compiler
warning，尚未在本批重跑 native 多域耦合矩陣。

每 rank 時間、RSS、stdout／stderr、controller summaries 與 build logs 保存於
`outputs/hpc03/transient-domain/`。transient completion line 分別記錄 runtime
assembly 與 linear solve 累計時間；wall／RSS 另含 fixture、serial reference、
MPI 啟動。部分測試與編譯並行，這些結果不作獨立 scaling／效能比較。

下一步擴充 `ImmersedFlowCase` 與正式 graph 初始化，只有 fixed stationary
backward-Euler 模式通過多步 graph、介面守恆與 precommit failure 後才放行。
移動幾何、FSI、checkpoint 與跨節點不屬於本次 adapter 的完成宣告。
