# HPC-01C：CPU configured transport 的錯誤協調

日期：2026-09-08。狀態：本批指定驗收通過，HPC-01C 與整份清單仍未完成。

## 實作

`iga_solve` 現在先比較 system、步數、velocity override、輸出／checkpoint／
memory report 分支等執行控制，再比較 PETSc 初始化後可見的 options。
database 先作 regular-file／SHA-256 檢查，接續原有 partition／row capacity
驗證；編譯 transport system 後再以實際 field 數確認 row capacity。

設定、mesh、prescribed velocity 或 velocity manifest，以及邊界實際引用的
periodic tables，均在 parser／初始化之前比對檔案型別和內容。本地設定物件、
labels、boundaries、初值及 coupling pattern 的配置與讀取納入 local stages。
允許不同 rank 本地路徑，但內容須相同，且在整個執行期間保持不變。

snapshot series 在每步選取 lower／upper 後，先協調選取錯誤、比對實際選中
檔案，再讀取及插值。保留未使用的 snapshot／temporal definition 不開檔的行為。
若後續時間超出 manifest 範圍，或後續才用到的檔案有問題，於該步開始求解前
共同退出，先前已接受步的輸出仍保留。waveform 乘上邊界值後若溢位為非有限值，
在 `transport step input` 拒絕，不將非有限邊界送進 KSP。

其他協調範圍包括：

- memory report 的本地建構、量測／root buffer 配置、gather 後的 root 寫入。
  在 MPI 初始化前啟用 allocation tracking 的例外也保存到初始化後共同報告。
  保留原有量測欄位與階段，`Record` 使用 literal stage 名稱以避免呼叫前的
  暫時字串配置。memory report 的 FIFO／目錄在開啟前拒絕。
- 元素矩陣配置、局部積分／插入、assembly begin/end、boundary rows、
  Mat/Vec 操作與每步邊界值配置。既有矩陣、向量及 KSP 由本地 owner 在共同
  unwind 時清理；借用 array 先 restore 再釋放 vector。
- KSP 選項設定、必要 factor backend、setup、block setup、solve 與收斂結果。
  returning error handler 讓已返回的 PETSc 錯誤進入共同處理；setup 的兩個
  操作分別協調，避免在本地錯誤後只有其他 ranks 進入 block setup。
- 輸出路徑配置、gather、root 場檔／VTU 寫入、時間索引 bookkeeping、PVD／
  physiology manifest，以及成功 summary。text fields 與 field-name sidecar
  明確 close 後檢查狀態；正常 PETSc cleanup 成功後才列印最後 summary。

沒有更動方程、quadrature、KSP/PC 預設、initial-guess 規則或數值容許值。
例如 `KSPPREONLY` 與第二步的 nonzero initial guess 原本不相容，現在明確
共同退出 1，而不是修改 initial-guess 行為以使它通過。

## Checkpoint 相容性

metadata 的本地讀取／驗證及序列化內容先協調，之後才讀取 state。
`ReadReplicatedPetscCheckpointVector` 新入口先自行驗證 rank 本地副本內容相同，
再沿用既有 owned-row 讀取、候選 vector、非有限值檢查及最後發布。
原有 `ReadPetscCheckpointVector` 仍要求相同絕對路徑，其接受／拒絕行為保留。
兩個入口共用實作，沒有以未驗證的 caller flag 跳過內容一致性。

寫入改用既有 `WritePetscCheckpointVector`，先完成獨立 state 檔的暫存寫入與
發布，再寫 metadata；metadata 複本須一致，root 寫檔錯誤會共同報告。
既有 PETSc binary Vec 與 JSON 格式不變；修改前 binary 實際讀取新 checkpoint
並續跑通過，包括單／雙 rank 及從單 rank checkpoint 在雙 rank 續跑。

**這仍不是完整的 crash-safe checkpoint bundle。** state 與 metadata 分開發布，
尚未綁定 waveform／完整案例身分；若中間中斷，不能由本批結果推論最後完整
checkpoint 一定可用。完整多檔發布與耦合狀態契約繼續由 HPC-05 驗收。

## 驗收

| 驗收 | 範圍與結果 |
|---|---|
| 新增原生 CLI regression | 81 筆觀察、154 份 rank reports；18 筆健康、63 筆預期失敗，均符合指定 gates |
| 健康場比較 | prescribed、snapshot series、VTKHDF、override、memory、unused definitions、同／不同 rank restart；67 份比較中有 12 份參考自比，55 份為實際候選比較 |
| 原有 shared-path reader | 67 個故障／不變狀態／重試案例通過 |
| 新 replica reader | 39 個故障／不變狀態／重試案例通過；含 world、獨立 1+2 groups、空 vector 與無 owned rows 的 rank |
| 舊執行檔讀新 checkpoint | 3 次續跑與最終場比較通過 |
| HDF5 實際 datasets | 兩次前後比較，各讀時間、Points、Connectivity、Offsets、兩個物種的 PointData |
| 既有資源 regression | 16 項通過 |
| 重建 flow 的 VCA smoke | flow／transport／reservoir、守恆與單／雙 rank 續跑 gates 通過 |

三筆 `*-before` 是修改前參考執行，不算獨立改版驗收。CPU 場最大 relative L2
`3.0631271369706667e-12`，維持 `1e-6` 與零參考 absolute `1e-12` 門檻。
HDF5 的三個時間值與幾何／拓撲 datasets 相同；雙 rank PointData 的 red／blue
relative L2 分別為 `9.7316818706642584e-13`／`8.4068019897486364e-13`。
HDF5 比較使用保存的 C++ reader 實際讀檔，不以文字場通過代替 HDF5 資料比較。

最終 snapshot fixture 的上端速度為下端的 1.5 倍，實際驗證非恆定速度插值；
兩端速度 L2 差為 `0.7901234567901236`。舊 binary 續跑與 HDF5 補充比較
仍使用 `cli-replica/` 的 prescribed／VTKHDF 產物；這兩組 fixture 與正式 binary
在最終 CLI 回歸中未變更。較早的封存另保留於
`source-before-distinct-snapshots.*` 與 `acceptance-before-distinct-snapshots.json`。

63 筆預期失敗涵蓋參數／控制差異、各種資產的有效內容差異、missing／FIFO／
directory、malformed input、第二步才遇到的 snapshot 差異／時間超界、
waveform overflow、KSP 不相容／後端不可用、memory report、初始／每步／
最終場／index／physiology 輸出與 checkpoint metadata 寫入失敗，以及損毀續跑。
內容差異部分使用有效格式的換行差異，驗證精確位元組契約，不宣稱每次差異
都改變物理模型。每 rank 均核對預期階段、原始退出 1、無 timeout、無最後成功
summary；late-input／overflow 另外檢查最後已接受步與下一步輸出的界線。

## 命令、環境與證據

環境：GCC 11.4、OpenMPI 4.1.2、PETSc 3.15.5 real64／Int32、MUMPS。
一般 CLI 回歸以 GMRES＋MUMPS LU，`ksp_rtol=1e-12`；OMP／BLAS 固定 1，
core mapping／binding。child timeout 60 秒、report rendezvous 15 秒、
job timeout 90 秒、kill grace 5 秒；VCA smoke 外層 180 秒。
部分作業重疊，這些小案例是正確性證據，不能作效能或跨節點擴展性宣告。

```bash
make -C solvers/cpu iga_solve petsc_checkpoint_read_test \
  iga_navier_stokes iga_navier_stokes_openmp \
  PETSC_DIR=/usr/lib/petscdir/petsc3.15/x86_64-linux-gnu-real -j2
python3 scripts/hpc_transport_cli_regression.py \
  --case-dir outputs/hpc01/checkpoint-write/verified/staged-failure/1-1 \
  --reference-binary outputs/hpc01/transport-cli/iga_solve-before \
  --output-dir outputs/hpc01/transport-cli/cli-final
env OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 \
  timeout --kill-after=5s 90s mpiexec --map-by core --bind-to core -np 3 \
  solvers/cpu/petsc_checkpoint_read_test outputs/hpc01/transport-cli/reader-fixture
```

以上 output directories 是已完成證據位置；重跑須換新目錄。其他完整 argv／
日誌見 `resources/`、`old-reader-*/`、`vca-smoke/`、`hdf-validation.json`。
HDF reader 的來源在同一證據根目錄 `compare_hdf.cpp`，本機重建命令：

```bash
g++ -std=c++17 -O2 -Wall -Wextra -Wpedantic \
  -I/usr/include/hdf5/openmpi -I/usr/lib/x86_64-linux-gnu/openmpi/include \
  outputs/hpc01/transport-cli/compare_hdf.cpp \
  /usr/lib/x86_64-linux-gnu/libhdf5_serial.so.103 \
  -o outputs/hpc01/transport-cli/compare_hdf
```

`outputs/hpc01/transport-cli/acceptance.json` 彙整最終 scope 與 hashes。
已重新核對所有 CLI rank reports、日誌、輸入及 binary 身分；144 份相關來源
保存在 `source-final.json`／`source-final.tar.gz`，正式 binary 另行保存。
`iga_solve` SHA-256：
`7912a09febc4fca956bab84f82576e2aff15566fc32b3cd7ae3c0c30af010254`。

首輪 build 的縮排 warning 已修正，最終 transport／flow／unit builds 無新增
compiler warning。首輪 harness 把 `output-control` 誤當作檔案破壞模式，在
執行該案例前發生 KeyError；`cli/` 與舊腳本保留。`cli-replica/` 通過較早的
fixture，最終 CLI 驗收以 `cli-final/` 為準。HDF reader 起初使用本機不存在的 serial header 目錄，已按實際
OpenMPI header／既有 serial runtime 配置修正，未改 HDF5 輸出或比較門檻。

## 接續工作與限制

HPC-01C 仍待其他 CPU／explicit CLI 的資產與錯誤邊界、flow 後段、adapter／
executor 的剩餘操作，以及共用 I/O helper 的寫入完成檢查。此批沒有故障注入
覆蓋每一個配置點，也不承諾從程序被殺死、檔案系統失聯或永不返回的 PETSc
collective 中原地恢復。輸入雜湊亦沒有檔案鎖／快照；完整持久化與大型 I/O
量測仍屬 HPC-05／06。既有 legacy flow 壓力比較限制也仍保留於
[flow 輸入報告](HPC_01C_FLOW_INPUT_PROGRESS.md)，未由本次 transport 通過取代。
