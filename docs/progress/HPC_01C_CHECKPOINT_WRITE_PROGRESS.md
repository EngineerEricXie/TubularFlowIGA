# HPC-01C：3D checkpoint 寫入錯誤協調

日期：2026-09-08。狀態：本批通過驗收；HPC-01C 仍部分完成。
接續 [checkpoint 讀取報告](HPC_01C_CHECKPOINT_READ_PROGRESS.md)。

## 實作

[PetscCheckpointWrite.hpp](../../solvers/cpu/include/PetscCheckpointWrite.hpp)
取代原生 flow 與 in-process transport 的未檢查 PETSc viewer 寫入。
所有 rank 先準備 owned coefficients 的副本，驗證有限值、索引／byte offset
容量及共同 path／layout。Root 在目標檔同一目錄建立唯一暫存檔，設定長度並寫入
原有 Vec header，再傳遞暫存檔路徑。每個 rank 以本地 binary write 寫入自己的
連續區段，檢查 seek、write、fsync、close；全部完成後，root 才 rename 至正式檔。

本地 callback 中沒有 collective viewer。資料副本也避免直接修改 solver Vec：
PETSc binary write 會在原 buffer 上轉換位元組順序，詳見
[PETSc 官方說明](https://petsc.org/release/manualpages/Sys/PetscBinaryWrite/)。
實際 3.15.5 build 的格式相容性另由以下測試確認。

發布前的可控制錯誤共同拋出，root RAII 嘗試移除暫存檔，原 Vec 檔保留。
不 gather 全域場；額外記憶體為 owned scalar 副本及路徑描述。
要求同一共享檔案系統、每個目的路徑只有一個 writer group；獨立 groups 使用
不同路徑。目的地只接受 regular file 或不存在的路徑，拒絕 symlink／FIFO／directory。
新檔由 mkstemp 建立，權限為 owner-only；替換既有 regular file 時保留其 0777
範圍內的 permission bits。沒有產生 `.info` sidecar；普通單 Vec binary 格式不變。

[TransientTransportRuntime::WriteState](../../solvers/cpu/include/TransientTransportRuntime.hpp)
共同要求 committed phase，避免輸出尚未接受的 trial。原生 flow CLI 的 metadata
準備、序列化內容一致性、root 寫入及 logging 亦加入錯誤協調；VCA metadata
採同樣流程。Metadata 寫入明確 close 後檢查 stream，捕捉延遲 flush 的錯誤。

這是單一 Vec 的發布，不是 flow／outlet／transport／reservoir 的完整原子交易。
Metadata 寫入失敗仍可能發生在部分 Vec 已發布之後；HPC-05 仍需完成 manifest、
版本／epoch、checksum、整體提交及最後完整 checkpoint 的恢復契約。

## 驗收

| 測試 | 結果 |
|---|---|
| Writer unit | 65 個故障案例通過；7 種 rank／row 配置 |
| Reader unit | 原有 67 個案例通過 |
| Staged transport | 104 個案例通過：3／1／2 ranks 分別 37／30／37 |
| Flow trial、port、gather | 原有 93、64、12 個案例通過 |
| Flow／species graph groups | 最大 relative L2 分別 `7.45082e-15`／`3.35727e-15` |
| VCA runtime | 原有 PETSc 預設通過 |
| VCA smoke | MUMPS；1／2-rank 場、reservoir、oxygenator 及同 rank 續跑精確比較通過 |
| 原生 CLI output faults | 17 項皆退出 1，正確 stage、無 checkpoint 成功標記／最終場／遺留暫存檔 |
| 新 writer 檔案的損毀續跑 | 原有 14 項皆按預期退出 1 |

[Writer unit](../../solvers/cpu/tests/test_petsc_checkpoint_write.cpp) 測試空／NUL
路徑、單 rank 不同路徑、缺少父目錄、父路徑為 regular file、目的地為 directory／
FIFO／symlink、單 owning rank NaN，以及真實 file-size-limit I/O error。
`RLIMIT_FSIZE` 暫時設為 0，保留 hard limit，暫時忽略 SIGXFSZ 並在呼叫後恢復。
有非空末端 rank 時錯誤發生於非 root 的 data write；空 rank 配置及單 rank
則驗證 root 檔案準備失敗。沒有改 production code 來假造 I/O 回傳碼。

每次使用不同的新係數，確認失敗不改變 source bits，也不覆蓋舊檔；確認沒有
遺留暫存檔、錯誤訊息全群相同，再解除錯誤成功寫入並重讀。另比較原 PETSc
VecView 的完整 binary bytes，並用原 VecLoad 獨立讀取新 writer 結果。
配置為 `(ranks, rows)=(3,4)、(3,2)、(3,0)、(1,3)、(2,3)、(1,0)、(2,0)`，
涵蓋空 rank、全空 Vec 與獨立 1+2 subcommunicators。

[Staged test](../../solvers/coupling/tests/test_three_d_staged_failure.cpp) 新增 trial
期間拒絕 WriteState，以及完成兩步後寫入缺失父目錄的錯誤；驗證 phase、clock、
state 不變。先前的正常寫入、讀回、rollback 與流向反轉驗收維持通過。

[CLI writer regression](../../scripts/hpc_checkpoint_write_regression.py) 在 1／2 ranks
測試缺失父目錄、state directory／FIFO、flow／VCA metadata directory、transport
state directory，並將 metadata 導向 `/dev/full` 以驗證實際 flush／close 失敗。
另以 MPMD launcher 給兩個 ranks 不同 checkpoint 路徑，確認共同拒絕。
這些是作業仍存活時的協同失敗測試，不是 MPI process crash recovery。

八項 MPI 與兩項 VCA 正常作業均在 180 秒時限內退出 0；31 項原生 CLI 預期失敗
作業均在 90 秒時限內退出 1。數值容許值未改：relative `1e-6`、零參考 absolute
`1e-12`；source／rollback 快照及同 rank checkpoint 比較使用精確相等。

## 重現與證據

環境：WSL 工作站，GCC 11.4、Open MPI 4.1.2、PETSc 3.15.5 real/double、
32-bit PetscInt、MUMPS，OMP／OpenBLAS 各 1 thread。
HEAD `ee8da2ab528849814b6d12c8186476b5f7f59124` 加未提交工作樹。
受影響 CPU 與 coupling binaries 重建通過，沒有新 compiler warning。

```bash
make -C solvers/cpu petsc_checkpoint_write_test petsc_checkpoint_read_test \
  petsc_gather_test iga_navier_stokes vca_3d_runtime_test vca_3d_smoke_test \
  PETSC_DIR=/usr/lib/petscdir/petsc3.15/x86_64-linux-gnu-real -j2
make -C solvers/coupling petsc three_d_staged_failure_test three_d_trial_failure_test \
  three_d_port_failure_test multidomain_subcommunicator_test \
  PETSC_DIR=/usr/lib/petscdir/petsc3.15/x86_64-linux-gnu-real -j2
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 timeout --kill-after=5s 180 \
  mpiexec --oversubscribe -np 3 solvers/cpu/petsc_checkpoint_write_test /tmp/checkpoint-write-fresh
```

原生 CLI regression 使用 [讀取報告](HPC_01C_CHECKPOINT_READ_PROGRESS.md#重現與證據)
所列方式產生並保留 smoke fixture，然後執行：

```bash
python3 scripts/hpc_checkpoint_write_regression.py \
  --case-dir "$checkpoint_scratch/tubularflowiga-vca-3d-smoke" \
  --output-dir "$checkpoint_scratch/native-write"
python3 scripts/hpc_checkpoint_read_regression.py \
  --case-dir "$checkpoint_scratch/tubularflowiga-vca-3d-smoke" \
  --output-dir "$checkpoint_scratch/native-read"
```

本地 ignored 證據：`outputs/hpc01/checkpoint-write/`，包含 `verified/mpi-summary.json`、
`vca/summary.json`、`native-write/summary.json`、`native-read/summary.json`、各項 logs／
fixture、`build/`、`evidence-summary.json` 與 `inventory.json`。記錄命令、環境、退出碼、
執行檔與資料 hashes；unit 預跑另保留於 `unit.log`，不重複計入驗收數量。

## 限制與剩餘工作

沒有驗證多節點檔案系統、斷電後 directory entry durability、rename／unlink
失敗、檔案系統卡住、MPI rank 遺失或同一路徑的並行 writer。沒有 directory fsync，
也沒有承諾任意環境故障後都能清除暫存檔。沒有新增 64-bit／complex PETSc、GPU、
效能／RSS／擴展性證據；geometry 未變，未重跑 mesh-test。

HPC-01C 的 constructor、其餘 adapter／executor／CLI 及外部資產一致性仍待完成。
完整持久化協議依 HPC-05 推進；本批不將 HPC-01C 或整份 38 項清單標為完成。
