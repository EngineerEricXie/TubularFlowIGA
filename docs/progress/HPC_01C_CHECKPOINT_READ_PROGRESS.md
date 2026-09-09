# HPC-01C：3D checkpoint 讀取與拒絕失敗

日期：2026-09-08。狀態：本批讀取邊界完成驗收；HPC-01C 保持部分完成。
接續 [初始化與 gather 報告](HPC_01C_INITIALIZATION_PROGRESS.md)。

## 變更與契約

[PetscCheckpointRead.hpp](../../solvers/cpu/include/PetscCheckpointRead.hpp) 提供
collective `ReadPetscCheckpointVector`。先協調路徑、向量大小與 ABI 描述，
各 rank 再以本地檔案操作讀取 header 及自己擁有的連續 rows。使用 PETSc 的
binary read／seek 處理位元組順序；不在本地錯誤 callback 中呼叫 collective viewer。
本地開檔採 nonblocking，隨後要求 regular file，因此 FIFO 不會等待另一端開啟。

讀取要求檔案長度恰為兩個 PetscInt header 加上所有 PetscScalar，驗證 Vec class、
row count 與 owned coefficients 的有限值。缺檔、截斷、尾端多餘資料或任何 rank
的本地錯誤先共同拒絕；全部通過後才組裝候選 Vec 並複製到 target，保留其 handle。
沒有 gather 完整場至 root；有 owned buffer、row indices 與候選 Vec 的暫時成本。

契約是所有 rank 讀取同一份不可變的共享 regular file，格式為目前 PETSc build
的普通、未壓縮、單一 binary Vec。沿用原生 `.state` 格式，忽略相鄰 `.info`
中的 options，避免載入時改變 solver 配置。本批只實測 real/double、32-bit
PetscInt；沒有宣稱跨 ABI、壓縮檔或任意 PETSc viewer 設定的相容性。

[TransientTransportRuntime::ReadState](../../solvers/cpu/include/TransientTransportRuntime.hpp)
先共同要求 committed phase，將成功候選發布至 current／committed snapshots。
沿用既有 `Steps() == 1` warm-start 語義；實際續跑時間仍由外層 metadata 決定。
檔案或前置驗證失敗不改動有效場、phase 與步數。

[iga_navier_stokes.cpp](../../solvers/cpu/src/iga_navier_stokes.cpp) 使用上述讀取器。
Flow metadata 的讀取、驗證、outlet 候選準備及精確描述比較先完成，才載入 flow。
VCA metadata、transport 路徑、reservoir restore 及 restart logging 的本地錯誤
亦加入協調。保留既有物理公式、solver defaults、數值 gate 與 outlet／VCA 格式。

## 驗收結果

| 檢查 | 結果 |
|---|---|
| binary Vec 讀取 unit | 67 個故障案例通過；world 3 ranks 與獨立 1+2 groups |
| transport staged | 98 個案例通過：3／1／2 ranks 分別 35／28／35 |
| flow trial | 93 個案例通過：33／27／33 |
| port／診斷 | 64 個案例通過：25／14／25 |
| generic gather | 12 個案例通過，包含空 rank／全空向量 |
| flow graph groups | 最大 relative L2 `7.45082e-15` |
| species graph groups | 最大 relative L2 `3.35727e-15` |
| VCA runtime | 原有 PETSc 預設下通過 lifecycle／staged assertions |
| VCA smoke | MUMPS；1／2-rank 數值、reservoir、oxygenator 與同 rank 續跑精確比較通過 |
| 原生 CLI 損毀續跑 | 14 項通過；每項退出 1，沒有 restart 成功標記或新 checkpoint／最終場 |

[binary unit](../../solvers/cpu/tests/test_petsc_checkpoint_read.cpp) 用真實 COMM_SELF
VecView 產生原格式參考，再驗證解析係數。測試缺檔、空檔、錯 class／row count、
截斷、尾端多餘資料、非有限值、directory、FIFO，以及最後一個 rank 的不同路徑、
NUL 或空路徑。每次錯誤要求共同 diagnostic、target／handle 不變，之後成功重讀。
六種 `(ranks, rows)` 為 `(3,2)、(3,0)、(1,3)、(1,0)、(2,1)、(2,0)`；
零 rows 不測非有限係數，單 rank 不測跨 rank 路徑差異。
惡意測試 sidecar 的 block size 與 KSP options 不會進入 options database。

[staged test](../../solvers/coupling/tests/test_three_d_staged_failure.cpp) 每個群組新增
五個拒絕案例：初始缺檔、trial 中載入，以及完成兩步後的缺檔、截斷和單 rank NUL。
拒絕後核對場、phase、clock。另保存第一步、完成流向反轉的第二步，確認兩份場不同，
再讀回第一步並 begin／abort，驗證 current 與 committed snapshot 均已更新。
原有非零流場、物種 accounting、rollback／retry 和 serial/group 比較維持通過。

[原生 CLI regression](../../scripts/hpc_checkpoint_read_regression.py) 複製 VCA smoke
保留的 checkpoint，在 1／2 ranks 分別測試缺失／非法 flow metadata、截斷 flow、
缺失／非法 VCA metadata、截斷／多餘尾端 transport。檢查 stage、退出碼與成功產物，
並核對來源 fixture hashes 不變。這些是可控制的檔案故障，不是殺死 MPI rank。

七項 MPI 和兩項原生 VCA 作業均在各自 180 秒時限內退出 0；十四項預期失敗的
CLI 作業均在 90 秒時限內退出 1。數值容許值沿用 relative `1e-6`、零參考
absolute `1e-12`；快照和同 rank checkpoint 比較保持精確相等。

## 重現與證據

環境：本地 WSL 工作站，GCC 11.4、Open MPI 4.1.2、PETSc 3.15.5
real/double、32-bit PetscInt、MUMPS，OMP／OpenBLAS 各 1 thread。
HEAD `ee8da2ab528849814b6d12c8186476b5f7f59124` 加未提交工作樹。
CPU、coupling tests 與 coupling CLI 重建通過，沒有新 compiler warning。

```bash
make -C solvers/cpu iga_navier_stokes petsc_checkpoint_read_test petsc_gather_test \
  vca_3d_runtime_test vca_3d_smoke_test \
  PETSC_DIR=/usr/lib/petscdir/petsc3.15/x86_64-linux-gnu-real -j2
make -C solvers/coupling petsc three_d_staged_failure_test three_d_trial_failure_test \
  three_d_port_failure_test multidomain_subcommunicator_test \
  PETSC_DIR=/usr/lib/petscdir/petsc3.15/x86_64-linux-gnu-real -j2
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 timeout 180 mpiexec --oversubscribe -np 3 \
  solvers/cpu/petsc_checkpoint_read_test /tmp/checkpoint-read-fresh
```

原生故障腳本需要 smoke 保留的 fixture。以下以新的暫存目錄隔離測試資料：

```bash
checkpoint_scratch=$(mktemp -d /tmp/iga-checkpoint.XXXXXX)
(
  cd solvers/cpu
  TMPDIR="$checkpoint_scratch" TUBULARFLOWIGA_KEEP_TEST_OUTPUT=1 \
    OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 \
    PETSC_OPTIONS='-ksp_type preonly -pc_type lu -pc_factor_mat_solver_type mumps' \
    timeout 180 ./vca_3d_smoke_test
)
python3 scripts/hpc_checkpoint_read_regression.py \
  --case-dir "$checkpoint_scratch/tubularflowiga-vca-3d-smoke" \
  --output-dir "$checkpoint_scratch/native-failure"
```

本地 ignored 證據位於 `outputs/hpc01/checkpoint-read/`：`verified/mpi-summary.json`、
`vca/summary.json`、`native-failure/summary.json` 保存命令、環境、退出碼與 hashes。
`evidence-summary.json` 核對成功標記及目前 binary 身分；`build/` 保存建置 logs，
`inventory.json` 記錄 source／binary／環境與證據身分。最初的 unit 預跑也保留在
`unit.log`，驗收數量以 `verified/` 的同版執行為準，不重複計數。

## 剩餘範圍

本批沒有改寫 WriteState／WriteCheckpoint，沒有建立整個 flow＋transport＋reservoir
的原子發布交易。HPC-05 仍需完整 state contract、校驗碼、epoch／manifest、
寫入故障恢復與重分區續跑。同路徑文字相同不代表不同節點的檔案內容相同；
共享不可變檔案是本讀取器的前置條件，沒有新增全檔內容一致性或並行覆寫保護。

沒有承諾任意 PETSc 內部配置／allocation／publish 失敗都能回復，亦沒有注入
rank process crash、網路或檔案系統失去回應。Constructor、其餘 adapter／executor
邊界及完整 CLI／外部資產一致性仍由 HPC-01C 追蹤。

未執行新的效能／RSS、GPU、跨節點或 64-bit／complex PETSc 驗收；本批不宣告
擴展性改善。未變更 geometry，未重跑 mesh-test。完整 38 項 goal 保持未完成。
