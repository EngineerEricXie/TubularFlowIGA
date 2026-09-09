# HPC-03C：固定幾何暫態 case 與正式 graph 驗收

日期：2026-09-09。基準 `e5e23f1` 加本批修改。
狀態：HPC-03C 已完成。完成範圍是 fixed stationary geometry 的 owned MPI
backward-Euler flow runtime、case 與 graph；移動幾何、跨節點、FSI、checkpoint
依其他任務追蹤。整份 HPC TODO 尚未完成。

## case 與入口

[ImmersedFlowCase.hpp](../../include/ImmersedFlowCase.hpp) 現在按 Navier–Stokes
`time_integration` 建立 steady 或 `backward_euler` backend。既有 steady `Load`
仍用序列 runtime；單程序 transient `Load` 使用 COMM_SELF owned runtime。
多 rank 先執行無分散式資源的 local `Preflight`，再在 caller communicator 上
collective `InitializeDistributed`，不在外層 local-failure callback 建立跨 rank
runtime。初始化前核對 time-integration mode 與 transient dt／steps 一致性。

case 的 fixed geometry identity 綁定 canonical surface、grid 與 volume quadrature
depth／rescue／storage；其他 catalog／operator 配置仍由既有 collective signature
核對。geometry catalogs 保持 immutable，wall material velocity 定義為零。
owned Mat／Vec、required-state halo、committed／prepared field 與 previous velocity
沿用已驗證的 transient runtime，production 多 rank 沒有複製完整數值求解。

設定沿用 `simulation_config.json` 與 `immersed_geometry.json`：

- 流場 equation system 使用 `"time_integration": "backward_euler"`。
- domain 與 graph 的 dt／steps 必須完全相同；每步 context 的 dt 也必須符合 case。
- transient `runtime` 可選 `"wall_inertial_gamma0": 0.6`，省略時沿用零值。
  steady 配置不接受此選項，以免靜默忽略。
- 此入口仍拒絕 moving geometry 配置、species／transport、body-fitted mesh／
  outlet-model／reference-profile 假設及不支援的 port controls。

`Load(..., mpi_size)` 是單程序相容介面，仍拒絕 size≠1；多 rank 的正式路徑
使用 `Preflight`／`InitializeDistributed`。這不是以多個 serial Load 複製求解。

[正式 runner](../../solvers/coupling/src/iga_1d_3d_bifurcation.cpp) 保留相同
explicit／fixed／Aitken coupling scheme。每次成功提交後，以 collective
conservation measurement 檢查 backend accepted clock／index／commit count；
`IGA_PROFILE=1` 另輸出 ownership、time integration、每步全域通量與分開的
assembly／linear-solve 累計時間。JSON／CSV 結果檔格式保持相容。

## HPC-03C 完成證據鏈

| 必要條件 | 權威實作與驗收 |
|---|---|
| previous velocity 唯一 owned snapshot、必要 halo、時間與失敗 guard | [history](HPC_03C_DISTRIBUTED_HISTORY_PROGRESS.md)、[frozen volume](HPC_03C_TRANSIENT_VOLUME_PROGRESS.md) |
| backward-Euler volume、wall inertia、ghost、port、gauge 與全域量 | [operator](HPC_03C_TRANSIENT_OPERATOR_PROGRESS.md)，殘差與 Jacobian action 對序列比較 |
| 多步場解、conservation、Newton／prepare／abort／Close、history 確實影響第二步 | [runtime](HPC_03C_TRANSIENT_RUNTIME_PROGRESS.md)，1／2／4 ranks、split、weighted；field relative L2 最大 1.47572e-9 |
| coupling rollback 使用 accepted history，恢復非零 controls／field，場與 clock 一起發布 | [adapter](HPC_03C_TRANSIENT_DOMAIN_PROGRESS.md)，對 legacy stationary moving owner；field relative L2 最大 3.709381e-14 |
| case、正式多域 graph、介面守恆、precommit failure、未支援模式明確拒絕 | 本報告的 native graph／factory／negative matrices |

前述兩份 field 數值引用各自已封存 revision 的測試；本批 native graph 比較
accepted port histories 與全域守恆，並非重新量測完整場 relative L2。各層證據
共同覆蓋要求，不把同一個綠色測試當成所有層的驗收。

## 可重現命令

```bash
make -C solvers/coupling iga_multidomain_flow iga_1d_3d_bifurcation immersed_case_factory_test \
  CXX=mpicxx PETSC_DIR=/usr/lib/petscdir/petsc3.15/x86_64-linux-gnu-real
python3 scripts/hpc_immersed_graph_regression.py \
  --output-dir outputs/hpc03/transient-case/explicit --execution explicit --transient
python3 scripts/hpc_immersed_graph_regression.py \
  --output-dir outputs/hpc03/transient-case/fixed --execution fixed --transient
python3 scripts/hpc_immersed_graph_regression.py \
  --output-dir outputs/hpc03/transient-case/aitken --execution aitken --transient
python3 scripts/hpc_immersed_graph_regression.py \
  --output-dir outputs/hpc03/transient-case/inertial --execution explicit \
  --transient --wall-inertial-gamma0 .6
python3 scripts/hpc_immersed_graph_regression.py \
  --output-dir outputs/hpc03/transient-case/steady-fixed --execution fixed
python3 scripts/hpc_immersed_case_rejection_regression.py \
  --output-dir outputs/hpc03/transient-case/rejections-final
env OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 \
  mpiexec -np 1 solvers/coupling/immersed_case_factory_test
```

控制器建立新的 output directory，從既有 aneurysm-chain source fixture 複製
並生成 unit-cube、4×1×1 cells／112 nodes 的小型 1D→immersed 3D→1D 案例。
不修改原始 reference case。每次執行記錄 binary／input hashes、完整命令、
MPI rank affinity、環境、日誌與資源，執行後核對 binary／inputs 未改變。

TsungYehLab、GCC 11.4、OpenMPI 4.1.2、PETSc 3.15.5 real64/int32、MUMPS；
OMP／BLAS threads=1。測試用允許 MPI socket 的本機權限。初次 `-j2` build 在
factory translation unit 的 GCC RTL/ira pass 發生 ICE；更新 factory test 後的
單獨重建成功，native runner 亦重建成功。沒有據此認定 ICE 原因。最終 build
無 compiler warnings；受限環境 `opal_ifinit socket errno=1` 為環境訊息。

## 數值與故障驗收

各 graph matrix 包含 1、2、4 ranks 成功求解、2-rank 第二步 precommit failure，
以及新的 2-rank 程序重試。這是新程序健康重試，沒有宣稱 checkpoint restart；
同一物件的 rollback／retry 由前述 adapter 測試證明。

固定兩步 dt=0.01；逐項比較所有 domain／port 的 time、area、outward flow、
mean pressure，門檻為 1e-12 + 1e-6 |one-rank reference|。flow interface
normalized residual ≤1e-10，fixed／Aitken 的 pressure residual ≤1e-6。
global surface／wall flow、volume divergence 亦以 1e-12 + 1e-6 |reference|
比較；backend time／index／commit count 必須精確一致。所有成功案例完整輸出
2 step rows、4 edge rows、12 port rows。4 cells 供 4 ranks 分配，owned row
總和等於 global rows，多 rank required rows 嚴格小於 global rows。

| 模式 | 成功步的迭代數 | 最大 port 誤差／容許值 | 最大 rank wall s | peak RSS bytes |
|---|---|---:|---:|---:|
| transient explicit | 1、1 | 7.08909e-7 | 1.319541 | 55,136,256 |
| transient fixed | 26、10 | 1.57989e-6 | 15.206690 | 55,115,776 |
| transient Aitken | 6、3 | 1.91383e-6 | 4.428638 | 54,558,720 |
| transient explicit、wall inertia=0.6 | 1、1 | 1.66441e-6 | 1.266392 | 54,075,392 |
| steady fixed 回歸 | 26、1 | 1.21066e-4 | 16.875175 | 54,034,432 |

各 rank 數的 fixed／Aitken 迭代數一致。transient accepted-step 診斷的最大
runtime 累計 assembly 14.297154 s、linear solve 0.624669 s；最大 residual norm
1.0391721e-10，實際 Newton convergence 仍由既有配置的 absolute／relative
gate 判定。這些 timings 含 coupling retries，非單一 Newton solve 的時間。
failure 第二步在 precommit 結束，其未接受工作的耗時仍包含在 rank wall time。

額外觀察：固定耦合 transient 與 steady 的 immersed port pressure 最大差
1.192927e-3；explicit wall inertia=0.6 與零值的最大差 8.537356e-5。這些是
兩種不同模型的差異觀察，證明入口沒有把模式／參數靜默忽略，不是 parity gate。

6 個 2-rank negative cases 分別測 dt mismatch、steps mismatch、負 wall inertia、
steady 使用 transient-only inertia、moving key、species field。全部由 rank 0
回報預期原因，所有 rank exit 1、無 timeout、無結果發布。第一版檢查器誤要求
每個 rank 都印錯誤，已按 CLI 的 root-only 診斷契約修正；原結果保留在
`rejections/`，最終為 `rejections-final/`。

factory regression 保留既有幾何／配置負面案例，新增 backward-Euler Load、
owned row coverage、dt guard、abort clock、重複初始化與 Close；exit 0，
wall 67.050881 s、peak RSS 66,154,496 bytes。總共 32 個有效 MPI jobs／
68 rank reports 通過其預期 gate，其中 11 jobs 為預期拒絕／precommit failure。

## 證據、限制與後續

`outputs/hpc03/transient-case/` 保留 summaries、每 rank logs／RSS、generated
fixtures、build logs、acceptance manifest 與 archive binaries。部分 correctness
jobs／編譯同時執行；不得以此宣稱效能加速、strong scaling 或跨節點能力。
CUDA 與既有長時間完整 closure suite 未在本批重跑。

已更新 root 能力矩陣、CPU README、coupling architecture 與失敗邊界索引。
HPC-03C 可勾選；HPC-03D 移動 ownership／halo、HPC-04 求解擴展、HPC-05
checkpoint 以及 HPC-09 真實多節點驗收均保持未完成，整體 TODO 為 12／38。
