# HPC-01C：清理診斷、snapshot 與串行工具驗收

日期：2026-09-08。HPC-01C 仍進行中；整份 38 項 goal 不變。
基準 HEAD：`ee8da2ab528849814b6d12c8186476b5f7f59124`，含既有未提交 HPC 修改。
本次入口狀態保存於 `outputs/hpc01/serial-tools-audit/initial-status.json`。

## 實作與核對

本輪開始時 F01 的兩個 AbortAll formatter，以及 F03 snapshot publisher 的
受檢查讀取／序列化，已存在於工作區。重新建置並驗證這些修改，沒有把文件中
過時的缺口當作尚未實作。F02 的 mesh／assembly summary 也已有明確 flush，
其他 CLI 的完整覆蓋仍待補。

本輪新增：

- Packer 的 mesh configuration 使用 `ReadCheckedText`，拒絕開啟失敗及讀取 EIO。
  不改 geometry transform、cache/text 選擇或 `.ntiga` 格式。錯誤訊息 formatter
  啟用 fail/bad exceptions；成功摘要明確 flush。
- Womersley 檔名 formatter 啟用 exceptions；每一 velocity field 在加入 manifest
  前明確 close，manifest close 後才回報成功，stdout 亦明確 flush。
- Transport budget JSON 的兩個中間串流啟用 exceptions，完整 JSON 明確 flush。
  Field reader 在最後一列後區分 EOF 與 I/O failure；守恆與殘差門檻不變。
- FSI exporter 的 ReadAll 使用完整讀取，Values 啟用 exceptions。
  沒有修改流體／結構方程、耦合收斂門檻或 VTU arrays。

## 故障與相容性證據

證據目錄：`outputs/hpc01/serial-tools-audit/`。
GCC 11.4、Open MPI 4.1.2、PETSc 3.15.5 real64／Int32，工作站本地執行。
MPI 測試在 sandbox 外允許建立本地通信 sockets；sandbox 原始啟動因 PMIx
listener 無法建立而失敗，未當作求解器錯誤或數值驗收。

- `pressure-abort.log`／`species-abort.log`：world 3 ranks 與 split 1／2 ranks
  全部退出 0。每個 executor／communicator 逐一注入 4 個 report allocation
  failure，合計 24 次，並有 30 次健康／重試。全部 domain 先完成 abort outcome，
  保留 primary failure，report failure 經共同 stage 回報，健康 retry 與 serial
  完整結果一致。原有 60 次 pressure 與 107 次 species fault/retry 亦通過。
- `snapshot.log`：219 次故障、222 次重試通過。涵蓋 publish、idempotent retry、
  collection rebuild 的 formatter，metrics／VTU read EIO，以及發布中斷。
  前一 accepted epoch 不變，重試後所有檔案與健康參考逐位元相同。
- `native/summary.json`：21 個原生 CLI 案例、30 個逐位元比較通過；含三種
  真實 formatter 例外、Womersley field／manifest close failure、stdout `/dev/full`、
  packer config read EIO。原生 Womersley 八場與 manifest、budget JSON 及 database
  均與本輪保存的修改前 binary 相同。Transport 使用既有 1-rank accepted history，
  原門檻 `rtol=1e-6`、`zero-atol=1e-12`，沒有重跑 transport 求解或 GPU。
- `packing.json`：重新從既有 control mesh 產生 cache 及 legacy text，修改前後
  四次 pack 全部退出 0，四份 `.ntiga` 逐位元相同。
- `fsi-helper.json`：使用真實 exporter 的 ReadAll／Values，三種 formatter failure
  與健康重試通過，完整 VTU 讀取後將 EOF 改成 EIO 必須退出 1；健康重試退出 0。
  該 helper 使用既有 VTU 作輸入，不代表已重跑完整 FSI。
- `old-counterexamples.json`：保留的舊 packer 在 config EOF 改成 EIO 後仍退出 0；
  舊 Womersley 在 field／manifest close 注入失敗後亦退出 0。三個反例皆重現
  靜默成功，對應修正後案例退出 1。這些退出 0 是被測出的缺陷。
- `stdout/summary.json`：mesh check／assembly smoke 的 1／2-rank 資源摘要與
  結果摘要，各注入 short-write、flush、標準例外、allocation、非標準例外。
  84 個原生作業、126 份 rank report、40 個健康重試摘要比較通過；故障 rank 0
  的診斷與所有 rank 的退出 1 均核對，健康 geometry 保持正 Jacobian。
- `core-regression.log`：`make cpu-test one-d-test coupling-test` 全部退出 0，
  包含原有 1D 的 1／2-rank PETSc 測試。`python-tests.log` 的 50 個 benchmark
  工具單元測試通過。新建置未出現 compiler warning。

完整 FSI 前後案例均退出 0、4 次耦合迭代收斂。`fsi-comparison.json` 記錄
5 份輸出逐位元相同，29 個數值 arrays 的相對 L2 皆為 0；含流體／膜的場、
幾何及連接資料。每次含 94,312 個流體點與 9 個膜節點。兩版皆維持
`moving_mass=0.015561258177167596`、`wall_leakage=0.015561258176288955`、
`continuity=1.8222488628470318e-14`，通過原有物理 gate。

新版 wall time 671.437 s、peak RSS 178,040 KiB；保留 baseline 為 670.592 s、
177,916 KiB。這是 I/O 修改的數值／格式相容性驗收，非效能驗收；兩個作業
部分重疊，不能用此時間宣稱 speedup。組裝／求解分階段時間及 CUDA allocation
本輪 N/A：沒有更改或量測平行效能，未執行 GPU。

## 重現命令

```bash
make -C solvers/coupling pressure_flow_collective_failure_test species_collective_failure_test text_close_preload.so
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 timeout --kill-after=5s 100s \
  mpiexec --oversubscribe -np 3 solvers/coupling/pressure_flow_collective_failure_test
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 timeout --kill-after=5s 100s \
  mpiexec --oversubscribe -np 3 solvers/coupling/species_collective_failure_test
make -C solvers/cpu PETSC_DIR=/usr/lib/petscdir/petsc3.15/x86_64-linux-gnu-real \
  iga_pack iga_womersley_reference iga_transport_validate \
  womersley_stream_failure_test transport_budget_stream_failure_test \
  phase8_compliant_channel_fsi_paraview fsi_export_stream_failure_test \
  snapshot_publication_failure_test text_read_preload.so
python3 scripts/hpc_serial_tool_regression.py \
  --baseline-dir outputs/hpc01/serial-tools-audit/before \
  --flow-case outputs/hpc00/matrix/flow-inputs \
  --transport-case outputs/hpc00/matrix/transport-inputs \
  --transport-history outputs/hpc00/transport-budget/validation1/cpu1 \
  --read-preload solvers/cpu/text_read_preload.so \
  --output-dir /path/to/new-output
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 timeout --kill-after=5s 120s \
  solvers/cpu/snapshot_publication_failure_test
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 timeout --kill-after=5s 900s \
  solvers/cpu/phase8_compliant_channel_fsi_paraview /path/to/new-fsi-output
solvers/cpu/fsi_export_stream_failure_test /path/to/new-fsi-output/t1.05/fluid.vtu
# Read-failure negative: must exit 1, with complete-text-stream diagnostic.
LD_PRELOAD="$PWD/solvers/cpu/text_read_preload.so" \
  IGA_TEST_READ_PATH=/absolute/path/to/new-fsi-output/t1.05/fluid.vtu \
  solvers/cpu/fsi_export_stream_failure_test /absolute/path/to/new-fsi-output/t1.05/fluid.vtu
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 make cpu-test one-d-test coupling-test
python3 -m unittest discover -s scripts/tests
```

Baseline binaries 是本輪開始時保存的實際檔案，hash 在 machine-readable report；
不宣稱僅 checkout 上述 HEAD 就能重建它們。新 source 與驗收報告會一起 commit；
執行檔、generated fixtures、VTK、database、logs 不提交。

剩餘工作：接續錯誤邊界索引 F02／F05／F06，並補 HPC-01D
能力矩陣。HPC-03～09 與跨節點驗收仍保留原範圍。
