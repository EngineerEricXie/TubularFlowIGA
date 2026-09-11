# HPC-01C F06：入口與清理最終簽核

日期：2026-09-10。基準 revision `b74b1ccb08ef06690ae3ddb54c0f0fc99f84e16c`。
本批完成 production MPI 入口、rank-local callback、distributed PETSc owner
與終止發布的最後呼叫鏈稽核；F01–F05 的既有精確證據沿用
[錯誤邊界索引](../architecture/MPI_FAILURE_BOUNDARIES.md)。

## 呼叫邊界結果

逐一核對八個 MPI production 入口：CPU flow、configured／legacy transport、
mesh check、assembly smoke、native 1D、multidomain／bifurcation 與 sequential。
參數早退、配置／外部資產、runtime construction、trial／prepare／finalize、
writer／checkpoint、明確 Close 及 completion 發布順序均對到現有共同階段。
單程序 CUDA、FSI exporter 及前處理入口維持其明確能力限制，不宣稱為 MPI runtime。

機器稽核解析 `include/`、CPU、coupling 與 1D production source 中的
`CollectiveLocalStage` callback。93 個含 MPI／PETSc API 的 callback 只使用：

- communicator／thread／版本／options／object metadata query；
- local Vec array access、ownership query、local `MatSetValues`／`VecSetValues`；
- solver convergence reason／iteration query與本地格式化。

沒有在 rank-local callback 中發現 `MPI_Allreduce/Bcast/Gather/Barrier`、
`Mat/VecAssembly`、`KSPSolve/SNESSolve`、`VecScatterBegin/End`、norm 或 distributed
object destroy。這些 collective 位於 callback 外，回傳後由共同 PETSc outcome
協調。`CollectiveFailure.hpp` 本身是協議實作，未被當成呼叫端 callback 掃描。

七種主要 distributed owner 均有顯式 `Close()`；flow／transport、assembly、
Newton、history 與 extension 的 destructor 只作例外退棧後備。正常 production
路徑在 completion summary／manifest 前完成 Close；每個 Close 嘗試完固定順序的
所有 destroy 後才協調第一個錯誤。Moving runtime 對 trial／retired／committed
epoch 使用同一順序，candidate restore 失敗也會先 collective cleanup。

## 目前 HEAD 驗收

GCC 11.4、Open MPI 4.1.2、PETSc 3.15.5 real64／Int32／MUMPS；本機
OMP／BLAS 各 1 thread。

| 驗收 | 結果 |
|---|---|
| 共同錯誤協議 | 3-rank world 與 split communicator 通過；transport retry relative L2 `2.22045e-16` |
| Pressure executor | world3、split1+2 的 abort formatter allocation faults、20 次 stage faults／retries 及全域 convergence 通過 |
| Species executor | world3、split1+2 的 allocation faults、35／36 次 stage faults／retries、反向 flow 與全域 convergence 通過 |
| Runtime construction／cleanup | world3、split1+2 的 flow／scalar 普通與空 rank、scalable option rollback、每種 owner 12 個 cleanup faults、重複 Close 及 retained references 全部通過 |
| Native CLI 最終矩陣 | 112 個作業全部符合預期：40 個 cleanup failure、72 個前版／健康／新作業 retry；296 份輸出比較、12 個 binary 與 423 個輸入 hash 通過 |

Native 矩陣涵蓋 1／2-rank CPU flow、VCA、flow／species graph、bifurcation，
以及 sequential explicit／fixed／Aitken。故障注入位於真實 PETSc destroy
回傳後；每個故障都不得印出成功摘要，coupling 路徑不得發布 completion 目錄，
後續新作業 retry 的場與其他輸出須符合封存前版。完整結果在
`outputs/hpc01/f06-final/native-all-head/summary.json`，機器可讀來源／呼叫邊界／
限制索引在 `outputs/hpc01/f06-final/completion-audit.json`。

## 完成界線

HPC-01C 對「程序仍存活且 communicator 可通訊」時的輸入、配置、組裝、輸出、
runtime lifecycle 與可控制返回錯誤提供共同失敗語義。MPI 初始化失敗、MPI／PETSc
內部程序故障、SIGKILL、OOM killer 或節點失聯不能原地協調，須依 HPC-05／09
從完整 checkpoint 啟動新作業。這項限制不是未完成的 F06 呼叫鏈。
