# HPC-01C：Transport VTKHDF 初始化錯誤協調

日期：2026-09-08（EDT）。本批指定範圍驗收通過，HPC-01C 與整份 goal
仍需完成其餘入口稽核與後續階段。

## 問題與修正

`iga_solve` 的 root-only Bezier mesh／VTKHDF 初始化原本只捕捉
`std::exception`，且在通知 peers 前寫出診斷。非標準例外會跳過這個廣播，
造成 root 異常終止而其他 ranks 等待。

現在整個本地初始化使用 `CollectiveLocalStage`，標準與非標準例外均先形成
共同錯誤結果，再讓各 rank 進入既有清理與退出路徑。Root-only 幾何建構、
路徑檢查、geometry report、writer 建構與初始化診斷都在這個階段內。
沒有將 HDF5 collective 放入本地 callback；這條輸出路徑仍由 root 寫入。
物理參數、輸出格式與數值容許值維持不變。

## 原生 CLI 驗收

`outputs/hpc01/transport-visualization-init/cli/summary.json` 保存 **34 個案例、
51 份 rank reports**，controller 已返回 0。

| 範圍 | 證據 |
|---|---|
| 修改前健康基準 | 1／2 ranks 各一次 |
| 修改後健康執行 | 1／2 ranks 各一次 |
| 修改後故障 | 1／2 ranks 各 7 種：H5Fcreate 返回錯誤、bad_alloc、非標準例外；geometry report 與 HDF 路徑各為 directory／FIFO |
| 健康重試 | 每次故障後用新程序、新輸出目錄重跑，共 14 次 |
| 舊版問題重現 | 單 rank root 返回 134；雙 rank root 返回 134、peer 超時返回 124 |
| 文字場比較 | 16 次新版健康執行的 128 份場／欄位名稱檔，與相同 rank 數的修改前結果逐位元組相同 |
| HDF5 比較 | 16 次健康執行，檢查 time values、points、connectivity、offsets、three_red、three_blue 的 shape、有限性與數值 |

所有新版故障都要求每個 rank 返回 1、沒有 timeout、包含共同階段診斷，且
尚未發布初始場、最終場或 physiology manifest。HDF 建構失敗前可能已寫出
geometry report；本批不宣稱整個輸出目錄具有原子發布性。

HDF 比較沿用 relative L2 `1e-6`、零參考 absolute `1e-12`；文字檔另採更嚴格
的逐位元組比較。此為指定資料集驗收，不宣稱比較所有 HDF metadata 或完成
ParaView GUI 驗收。組裝／求解計時與 RSS 保留在 rank logs，未作效能比較。

測試 interposer 僅匹配指定的完整輸出路徑，於 C++ caller 的 H5Fcreate
邊界、進入 HDF5 C library 前注入錯誤；不讓 C++ 例外跨越 library C stack。
未修改 production binary 的建置旗標或加入 production 故障開關。

## 重現方式

環境：本機 WSL、GCC 11.4、Open MPI 4.1.2、PETSc 3.15.5 real64／Int32，
GMRES＋MUMPS LU、`ksp_rtol=1e-12`、OMP／BLAS 各 1。
Production、preload 與 reader builds 均無 compiler warning。

先保留修改前 binary，並依
[transport CLI 報告](HPC_01C_TRANSPORT_CLI_PROGRESS.md)產生 retained fixtures。
以下使用本機既有 HDF5 配置；在其他機器須使用該站的 header／library：

```bash
make -C solvers/cpu iga_solve \
  PETSC_DIR=/usr/lib/petscdir/petsc3.15/x86_64-linux-gnu-real
g++ -std=c++17 -O2 -Wall -Wextra -Wpedantic -shared -fPIC \
  $(pkg-config --cflags hdf5) solvers/cpu/tests/hdf_create_preload.cpp \
  -ldl -o /path/to/hdf_create_preload.so
g++ -std=c++17 -O2 -Wall -Wextra -Wpedantic $(pkg-config --cflags hdf5) \
  solvers/cpu/tests/compare_transport_vtkhdf.cpp \
  /usr/lib/x86_64-linux-gnu/libhdf5_serial.so.103 -o /path/to/compare_transport_vtkhdf
ulimit -c 0
python3 scripts/hpc_transport_visualization_init_regression.py \
  --fixtures outputs/hpc01/transport-cli/cli-final \
  --reference-binary /path/to/before/iga_solve \
  --preload /path/to/hdf_create_preload.so \
  --hdf-reader /path/to/compare_transport_vtkhdf \
  --output-dir /path/to/new-results
```

每個 child timeout 20 秒，rank report rendezvous 30 秒，launcher timeout
60 秒、kill grace 5 秒。原始退出碼與 timeout 保存在 report，舊版異常案例
沒有被歸類為正常錯誤退出。新舊二進位、輸入及 rank log hashes 均有核對。

## 接續稽核

- Standalone FSI ParaView exporter 原本沒有 MPI size gate；本輪已加上在建立
  輸出前拒絕多 rank 的檢查，驗收狀態另見
  [exporter 入口報告](HPC_01C_EXPORTER_ENTRY_PROGRESS.md)。分散式 FSI 仍由 HPC-07 負責。
- CPU flow／transport 的控制值序列化，以及 1D／graph 的 ReadText、JSON escape
  與 sequential CSV 中間字串，仍需核對 stream failure 是否會被正確轉為共同
  失敗；不能僅因外層有 collective callback 就忽略 stream 自身的錯誤狀態。
- 其餘入口與共用 helper 的覆蓋稽核尚未完成，HPC-01C 保持未勾選。
