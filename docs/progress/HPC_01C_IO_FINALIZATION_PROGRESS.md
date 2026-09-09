# HPC-01C：共用 I/O 的明確關檔與失敗協調

日期：2026-09-08。狀態：本批關檔修正與指定驗收已完成；HPC-01C 及整份
清單仍未完成。工作區包含未提交修改，以封存來源、binary 與輸入雜湊辨識版本。

## 行為變更

`WriteBezierGeometryReport` 原先在關檔前檢查 stream；`CouplingHistoryWriter`
原先沒有檢查最後寫入狀態。兩者現在明確 close，再檢查 stream。原先的小檔案
緩衝區寫入失敗會回傳成功；history 的大型串流寫入失敗也會漏報。
健康 JSON 的格式與內容保持一致，沒有改變 checkpoint 或數值介面。

`TemporalVtkHdfWriter::Close()` 檢查 flush、root group close、本 file ID 的
殘留 objects 及 file close。CPU flow／configured transport 在共同失敗階段中
呼叫，CUDA flow／configured transport 亦在最後成功摘要前呼叫。
HDF5 允許仍有開啟物件時延後真正關檔，因此不能只檢查 `H5Fclose` 的返回值；
本地物件檢查使用 `H5F_OBJ_LOCAL`，不要求其他 reader 一併關閉。
此處依據 [HDF5 官方文件](https://support.hdfgroup.org/documentation/hdf5/latest/group___h5_f.html)。

成功 Close 可重複呼叫。開始 Close 後禁止 Append；關檔失敗時保留仍持有的
handle，允許呼叫端再次 Close。destructor 仍只負責 best-effort 清理，避免
unwinding 中拋出第二個例外。原生 CLI 收到關檔錯誤即失敗退出，不自動重跑模擬。

## 驗收證據

證據目錄為 `outputs/hpc01/io-finalization/`，產物保持 ignored。`acceptance.json`
彙整指定 gates 與另外保留的未通過觀察；來源快照見 `source-final.json` 及
`source-final.tar.gz`。本批未修改原有數值容許值。

| 驗收 | 實際結果與範圍 |
|---|---|
| 文字 writer | 6 個 writer、17 個故障、11 次健康重寫通過。新增 geometry／history；geometry 固定小 payload，不冒充大型 stream 測試 |
| 修改前文字版本 | 同一測試搭配封存標頭，重現 3 個漏報，退出 1；保留 `text-before.log` |
| 健康 JSON | geometry／history 修改前後逐位元組相同，另用 JSON reader 解析 |
| HDF5 wrapper 單元測試 | 5 個故障：flush、group close、file close、dataset close 遺留 handle、object count 查詢；全部拒絕並在解除故障後 Close 成功 |
| HDF5 狀態與資料 | 檢查重複 Close、關閉中／關閉後拒絕 Append、file objects 回到原數量；6 次實際讀回 scalar／time，同檔獨立 reader 仍可用 |
| 既有格式測試 | Bezier visualization 與 temporal VTKHDF schema／三個時間步／resume 測試通過；後者同時保留 destructor 相容路徑與 explicit Close 路徑 |
| 原生 CPU 關檔注入 | 10 次執行、14 份 rank reports；新版 flow／transport 1、2 ranks 都共同退出 1，沒有成功摘要；兩個舊版單 rank 程式重現退出 0 且印成功 |
| CPU output／transport 回歸 | flow 46 筆／83 rank reports、68 份場比較；transport 81 筆／154 rank reports、67 份場比較，全部符合各案例預期 |
| CUDA 關檔回歸 | RTX 4080 SUPER 上 flow／transport 共 8 次執行，包含舊／新健康與故障版本；48 份修改前後場比較完全一致，兩個新版故障案例退出 1 且沒有成功摘要 |
| HDF5 原生輸出讀回 | 14 組檔案比較、84 個 dataset 比較通過，涵蓋 CPU rank、故障後清理及 CUDA 修改前後；包含時間、點、connectivity、offsets 與兩個場 |
| 耦合 | 重建 1D CLI／coupling test，1D history 與既有 VCA flow／transport／checkpoint／restart smoke 通過 |

文字故障使用子程序 `RLIMIT_FSIZE=64`，忽略 SIGXFSZ，確認實際 regular file
截斷為 64 bytes；父程序隨後重寫，與健康 bytes 比較。HDF5 單元測試使用
GNU linker `--wrap`，原生 CLI 使用測試專用 `LD_PRELOAD` shim；只對 controller
指定的完整 VTKHDF 路徑回傳一次 H5Fclose 錯誤。production code 沒有故障開關。
測試保留真實 handle 與檔案，不能推論所有真實儲存故障都能恢復。
dataset 遺留案例由測試顯式釋放故意遺留的 handle 後才重試，writer 不會自動
修復 orphan。既有 CLI 場比較包含參考自比，不能把所有比較都視為獨立證據。

CPU transport 場比較最大 relative L2 為 `3.0631271369706667e-12`，flow output
為 `1.8156697926380728e-13`；native close 為 `1.4601410874306183e-12`。
CPU gate 為 relative `1e-6`／零參考 absolute `1e-12`。CUDA 修改前後為 0；
另外 4 個 CPU/GPU transport 場比較最大 `2.397514096418484e-6`，通過 `1e-5`。

## 保留的失敗觀察與較嚴格數值比較

原生關檔 harness 首次誤將 flow 的 `preonly` 套用於 transport，健康傳輸案例
回報 PETSc error 83。這次結果保留在 `native-close/`；修正腳本為原有 GMRES
選項後，完整十個案例在新目錄 `native-close-verified/` 通過。沒有修改求解器
實作或放寬驗收條件來通過該測試。

額外比較預設 flow 收斂條件下的 CPU/GPU 場時，8 份比較中 6 份未達 `1e-5`，
最大 relative L2 為 `0.0024196303841121387`。完整結果保留在
`cpu-gpu-fields.json`，未將它改寫為通過；同 backend 修改前後完全一致。

改用既有 HPC 基準登記的 `--nonlinear-rtol 1e-8 --nonlinear-atol 1e-12
--mass-rtol 1e-6`，CPU 與 GPU 各重新求解一次，8 份場比較全部通過，最大
relative L2 為 `6.900038140369433e-7`。另讀回 6 個 HDF5 datasets 比較通過，
pressure 最大 relative L2 為 `6.835389232930688e-7`。命令／數值／日誌位於
`flow-accuracy.json`、`flow-accuracy/`、`hdf-cpu-gpu-accuracy.json`。
這是較嚴格求解設定下的驗收，不表示預設收斂條件保證同樣的場誤差。

## 重現方式

環境：GCC 11.4、OpenMPI 4.1.2、PETSc 3.15.5 real64／Int32、MUMPS，
CUDA 12.6／SM 89。CPU rank 綁定、OMP／BLAS 與求解器選項記錄在命令和 rank
reports。CPU 編譯器在沙箱內曾輸出 `opal_ifinit` socket 探測訊息，但三個
target 都建置退出 0；其後原生 MPI／GPU 驗證在授權環境執行並通過。
未出現新的 C++／CUDA 編譯警告。

```bash
make -C solvers/cpu text-output-failure-test vtkhdf-close-failure-test
make -C solvers/cpu temporal_vtkhdf_test bezier_visualization_test
./solvers/cpu/temporal_vtkhdf_test
./solvers/cpu/bezier_visualization_test
make -C solvers/cpu iga_solve iga_navier_stokes iga_navier_stokes_openmp \
  PETSC_DIR=/usr/lib/petscdir/petsc3.15/x86_64-linux-gnu-real -j2
make -C solvers/one_d iga_1d one_d_coupling_test \
  PETSC_DIR=/usr/lib/petscdir/petsc3.15/x86_64-linux-gnu-real -j2
conda run -n tubularflow-cuda make cuda CUDA_ARCHS=89
```

先用 `hpc_flow_output_regression.py` 與 `hpc_transport_cli_regression.py` 建立
各自的新輸出目錄。flow output 的舊成功順序重現須使用
`outputs/hpc01/flow-output/iga_navier_stokes-before`；transport 使用本批封存的
`outputs/hpc01/io-finalization/iga_solve-before`。這兩個 harness 的完整命令
與本批 logs／summary 均封存在證據目錄。

原生 HDF 關檔注入在 Linux 執行；將 `EVIDENCE` 設為既有本批證據根目錄，
而 `NEW_RESULT` 指向尚不存在的新結果目錄：

```bash
g++ -std=c++17 -O2 -Wall -Wextra -Wpedantic -shared -fPIC \
  -I/usr/include/hdf5/openmpi -I/usr/lib/x86_64-linux-gnu/openmpi/include \
  solvers/cpu/tests/hdf_close_preload.cpp -ldl -o "$EVIDENCE/hdf-close-preload.so"
python3 scripts/hpc_vtkhdf_close_regression.py \
  --flow-fixtures "$EVIDENCE/flow-output" --transport-fixtures "$EVIDENCE/transport" \
  --reference-flow "$EVIDENCE/iga_navier_stokes-before" \
  --reference-transport "$EVIDENCE/iga_solve-before" \
  --preload "$EVIDENCE/hdf-close-preload.so" --output-dir "$NEW_RESULT"
```

HDF include 路徑須依實際安裝調整。CUDA／較嚴格 flow 命令亦記錄在各自 JSON
與 `run_cuda.py`／`run_flow_accuracy.py`；重跑必須使用新目錄，避免混入舊證據。

## 剩餘工作

HPC-01C 的其他 runtime／配置、adapter／executor 與清理邊界仍須繼續稽核。
本次檢查不能涵蓋所有 transient dataspace／datatype／property-list 的析構
錯誤，也未提供程序失聯後的 MPI 恢復、全步 rollback、原子輸出發布或 fsync。
截斷檔案仍可能存在；完整 checkpoint 發布契約繼續由 HPC-05 追蹤。
沒有新增 ParaView 實際載入、跨節點或效能驗收；小案例重疊執行只用於正確性。
