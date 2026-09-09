# HPC-05C：0D 與耦合控制狀態保存／恢復進度

日期：2026-09-09。基準 revision：`e4bffb34d94e87a6ac980d6062385f9172c8cfd7`
加本批修改。狀態：**部分完成，HPC-05C 保持未勾選**。本批交付 accepted-state
介面與 payload codec，並驗證新的執行程序可恢復真實 0D 模型；完整 MPI graph
checkpoint／restart CLI 及 1D／3D state providers 在本批尚未完成。後續 1D runtime
與兩分片 codec 的進展另見 [1D 狀態進度](HPC_05C_ONE_D_STATE_PROGRESS.md)。

## 已實作介面

[ZeroDFlowDomainRuntime](../../include/ZeroDFlowDomain.hpp) 新增：

- `CaptureCheckpointState()`：只在已有 accepted step 的 Committed phase 擷取
  domain／model 身分、完整 accepted step context、stored pressure、PortState 與
  storage accounting。初始狀態、TrialReady、TrialSolved、Prepared 都拒絕。
- `RestoreCheckpointState()`：只初始化 fresh、idle candidate runtime。驗證 model／
  domain、clock、port 數量與 presence、RCR／source 物理關係和 storage balance 後，
  以不丟例外的 swaps／scalar assignments 恢復所有 committed 欄位。
  失敗不改變任何 accepted datum，已運行／已恢復的 runtime 不接受重複覆寫。
- 額外保留 `committed_step_`，在 finalize 時才更新。因而下一個 trial 使用不同 dt
  再 abort 後，checkpoint 仍攜帶原 accepted dt／start time。
  `BeginStep` 同時避免 committed index 位於 `INT_MAX` 時的有號溢位。

0D payload 由 [ZeroDFlowCheckpoint.hpp](../../include/ZeroDFlowCheckpoint.hpp)
編碼，版本為 `IGA_ZERO_D_ACCEPTED/1`。帶 bundle epoch 的 parse overload 同時驗證
accepted count、time 與 dt。原本的 flow／transport／VCA／1D 格式及 0D model digest
保持不變；0D kernels 的更新公式沒有改動。

[CheckpointMetadataCodec.hpp](../../include/CheckpointMetadataCodec.hpp) 提供共用
metadata 編碼：little-endian uint64／IEEE binary64、length-prefixed strings、optional
presence、ordered maps、完整 PortState 與 DomainStepContext。單一 metadata payload
上限 16 MiB、string 4096 bytes、map 16384 entries；先檢查長度／範圍再取值，拒絕
trailing bytes、重複／未排序／空 key、非法 presence、非有限值與 step index 溢位。
保留 signed zero 與 absent 值。大型場仍由 05B streaming shard 承載，不塞入此 metadata buffer。

[PressureFlowCheckpointControls.hpp](../../include/PressureFlowCheckpointControls.hpp)
保存 accepted context、下一步壓力 map 與 donor map。建構 helper 使用最後 accepted
iteration 的 **measured pressure**，與 native runner 一致。解析時對照 graph 的完整
edge／species catalog 及 bundle clock；`NextStep()` 直接延續實際 EndTime，不使用
`N*dt`。此 helper 必須由 caller 在 executor 已返回 accepted result 後使用。
全步歷史與 graph／配置 compatibility digest 仍由後續 runner provider 承擔。

[SpeciesPressureFlowComponentExecutor](../../include/SpeciesPressureFlowComponentExecutor.hpp)
新增 idle accepted capture 與 fresh-executor donor restore。restore 要求所有
`(edge,species)` keys 完整、無未知項，enum 僅接受 First／Second。active-step guard
使 before-commit callback 期間無法擷取或恢復 checkpoint；RAII 使成功與失敗都離開
active 狀態。既有 Aitken per-step reset、donor finalize 順序與 trial rollback 不變。

## 預先設定的驗收條件與結果

| 驗收 | 結果 |
|---|---|
| 真實 source reservoir、terminal RCR，含正／反向輸入 | 第 4 步保存，新的 `exec` 程序由 bundle 載入；與同程序重新計算的不中斷參考相比，accepted snapshot 及後續 6 步完整 payload 逐位元一致 |
| 狀態覆蓋 | stored pressure、time／step count／index、完整 port 與 accounting 一起恢復；模型／clock／presence／帳目不一致全部拒絕，失敗後 fresh state 完全不變 |
| lifecycle | 初始與三個 active phases 拒絕 capture；重複 restore 拒絕；以不同 dt 開 trial 再 abort 不污染 accepted context |
| 格式與範圍 | 每個 0D payload 截斷位置、尾端資料、oversized string／map、重複 key、非法 optional、錯 bundle count／time／dt 均拒絕；optimized 與 sanitizer 各通過 614 個拒絕判定 |
| graph controls／donor | 以原始 pressure guess=0、accepted measured pressure=1 驗證保存值；新 executor 恢復後可完成近零 flow step，隨後合法切換為反向 donor；prepare failure 保持上一個 committed donor |
| MPI executor regression | 一個 3-rank job，包含 COMM_WORLD 與 split 的 1／2-rank groups，既有 fault／retry／convergence／reverse／allocation-failure tests 全部通過 |
| production regression | 3-rank 0D–3D–0D case 執行兩步，與封存的 pre-change executable 對照；7 份 CSV、28 筆資料列、204 個 numeric values，最大絕對差 0 |

production 比較的預設門檻為 `abs(error) <= 1e-12 + 1e-6*abs(reference)`；非數值
欄位與 rows／columns 必須一致。實測所有數值完全相同。此處比較 accepted port／edge／
0D／balance histories，不宣稱本批量測完整 3D field relative L2 或已完成 graph restart。
一份空的 `one_d_initialization.csv` 只核對檔案格式，本案例沒有 1D runtime。

## 執行與證據

工作站 Linux x86_64、GCC 11.4、OpenMPI 4.1.2、PETSc 3.15.5 real64/int32；
`OMP_NUM_THREADS=1`、`OPENBLAS_NUM_THREADS=1`，未固定 CPU affinity。
編譯使用 C++17、`-Wall -Wextra -Wpedantic`；沒有新警告。

```bash
make -C solvers/coupling zero_d_checkpoint_test zero_d_flow_test \
  species_pressure_flow_executor_test CXX=g++ -j2
mkdir -p outputs/hpc05/state
solvers/coupling/zero_d_checkpoint_test "$PWD/outputs/hpc05/state/final"
solvers/coupling/zero_d_flow_test
solvers/coupling/species_pressure_flow_executor_test
make -C solvers/coupling iga_multidomain_flow species_collective_failure_test \
  CXX=mpicxx PETSC_DIR=/usr/lib/petscdir/petsc3.15/x86_64-linux-gnu-real -j2
mkdir -p outputs/hpc05/state/species-mpi
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 mpiexec -n 3 \
  python3 scripts/hpc_rank_run.py --output-dir outputs/hpc05/state/species-mpi \
  --expected-ranks 3 --timeout 120 -- solvers/coupling/species_collective_failure_test
```

每次重跑使用新的 output child directory。0D sanitizer 使用 `-O1 -g
-fsanitize=address,undefined -fno-omit-frame-pointer` 重新編譯同一 source，並以
`ASAN_OPTIONS=detect_leaks=1 UBSAN_OPTIONS=halt_on_error=1` 執行；ASan／UBSan／
LeakSanitizer 通過、退出碼 0、stderr 為空。MPI sockets 與 leak detector 使用所需
執行權限；restricted compile log 的 `opal_ifinit errno=1` 未影響編譯或實際 MPI 執行。

production fixture 為
`outputs/hpc01/assets/graph-regression/zero-d-config/fixture`：schema-v5 graph，
explicit coupling、dt=0.1，`zero_junction.ntiga` 為 3 partitions、1 element、64 nodes。
本批使用 `--stop-after-step 2`，無 shared-login-node 大型運算。兩個命令均經
`hpc_rank_run.py` 包裝，完整 argv、環境、退出碼及 source inputs 身分保存在證據中：

```text
solvers/coupling/iga_multidomain_flow --graph-case <fixture> --stop-after-step 2
    --output-dir outputs/hpc05/state/native-current-results
outputs/hpc03/transient-case/binaries/iga_multidomain_flow
    --graph-case <fixture> --stop-after-step 2
    --output-dir outputs/hpc05/state/native-baseline-results
```

`outputs/hpc05/state/` 包含 `final.log`、`sanitized.log/.err`、`species-final.log`、
`build-native.log`、三個 MPI job 的 9 份 rank reports、`native-comparison.json`、
`resources-summary.json`、封存執行檔與 `acceptance.json` 的 SHA-256。
所有最終測試退出碼 0，沒有 timeout；中間的 `run1`／`run2` 不是最終來源的替代證據。

production current job 的每-rank wall time 為 0.414411／0.414342／0.414169 s；
peak RSS 為 41,582,592／41,353,216／41,373,696 bytes。
root log 各 iteration 的 assembly 時間相加為 0.097449 s，linear solve 為 0.00275562 s。
species MPI job 三個 ranks 的 peak RSS 為 32,276,480／32,395,264／32,370,688 bytes。
不將 RSS 相加當成同時的總峰值。0D payload I/O 與通信未另作效能量測；本批不宣稱
性能改進，CUDA allocation／CPU-GPU 比較不適用。

## 接續工作

1. 補齊 1D state／LastInlet／species flux accounting 和動態 radius 的 validated restore。
2. 補齊 body-fitted flow owned fields／outlet／boundary state 與 transport 真實 step counter。
3. 接上完整 graph catalog／compatibility digests、accepted history prefix、MPI receipts／
   failure agreement、CLI checkpoint／restart 與安全 step-boundary 的保存旗標。
4. 以新的 MPI job 驗證 0D／1D／3D graph 不間斷對續跑的全歷史，及 commit 前、寫片中、
   manifest 發布前全作業中斷。immutable bundle 核心的先前測試不能取代此整合驗收。

本批沒有變更 05D 的 immersed／moving／FSI／repartition 範圍，也沒有縮小整份 goal。
