# HPC-01C：獨立 CPU 流場入口的輸入協調

日期：2026-09-08。狀態：runtime 建構前的輸入檢查與限定範圍的相容性驗收
通過；legacy 跨 rank 壓力比較仍有未通過觀察，HPC-01C 整體保持未勾選。

## 實作範圍

`iga_navier_stokes` 及其 OpenMP build 現在先協調以下本地操作，再建立包含
PETSc collective 的 runtime。沒有將 runtime 建構子包進 local-stage callback。

- CLI 解析及數值／控制選項的完整字串比對，包含 Newton、停止步數、輸出與
  checkpoint 分支、視覺化模式與三個容許值；允許本地路徑不同。
- 比對 PETSc 初始化後實際可見的 options；消耗 application options 時檢查
  `PetscOptionsClearValue` 的回傳碼。CLI 語法與 solver 預設維持原樣。
- `.ntiga` 先用非阻塞 descriptor 檢查 regular-file 並計算 SHA-256，再開啟
  Database／驗證 partition 和 row capacity，之後協調 database 內容一致性。
- 蒐集 mesh、initial velocity、modern configuration 或 legacy parameters／
  optional case configuration。先比對邏輯目錄，防止各 rank 選到不同設定路徑，
  再驗證檔案型別與內容。
- 在 collective local stage 中建立 `FlowCaseInput` 預設值、讀取拓撲／速度／
  設定、建立 wall trace，然後比對暫態邊界實際引用的 periodic tables。
  未引用的 tables 及 steady 路徑不讀取的 tables 維持原有行為。
- 初始 waveform、outlet／VCA 本地狀態、resolved boundaries、步數與
  visualization 條件檢查及啟動 logging 均先完成錯誤協調。

沿用 [資產 helper 契約](HPC_01C_ASSET_PROGRESS.md)：輸入在整個執行期間
必須不變，各 rank 可有不同的本地副本，但內容須完全相同。沒有檔案鎖或快照，
每 rank 額外完整掃描資產；沒有大型共享檔案系統 I/O 效能證據。
物理公式、檔案格式、solver 預設及數值容許值均未修改。

## 驗收與精確範圍

| 驗收 | 結果 |
|---|---|
| 新增 CPU input CLI harness | 51 筆觀察、94 份 rank reports；輸入契約與指定的相容性 gates 通過 |
| 既有 resource harness | 16 項通過，含缺 database、分區錯誤、後端與執行緒檢查 |
| 原生 VCA smoke | 通過非零流場、transport／reservoir、守恆及單／雙 rank 續跑 gates |
| OpenMP 正式 CLI | 1／2 ranks、每 rank 設定 2 threads，4 份 velocity／pressure 場比較通過 |

新 harness 的 51 筆包含 6 筆修改前參考執行、10 筆修改後健康執行與 35 筆
故障執行。正向路徑包括帶 table 的 flow、VCA、legacy defaults 與 optional
legacy configuration。32 份場比較中 12 份是參考自比，只有其餘 20 份是
修改後的相容性比較。最大 relative L2 為 `1.8156697926380728e-13`，沿用
CPU `1e-6`、零參考 absolute `1e-12`。flow／VCA 與修改前單 rank 比較；
legacy 的相容性比較使用相同 rank 數的修改前 binary，原因見下一節。

35 筆故障涵蓋單 rank 參數錯誤、控制值／是否輸出不同、有效但不同的 database、
configuration、mesh、velocity、table、legacy parameters／configuration，
optional configuration 選取不同，以及各資產的 FIFO、目錄、必要檔案缺失。
所有 ranks 都讀到 malformed configuration、table 或 parameters 時，也共同
拒絕。每份 rank 日誌均核對預期階段，保留退出 1，無 timeout、無結果目錄。
設定與 mesh 的有效內容差異包含純格式差異；此處驗證精確位元組契約，並不
聲稱每一項差異都代表物理模型不同。

OpenMP 案例只有一個元素，僅證明該 build 與配置仍可用，不證明多 worker
加速。上述作業可能重疊執行，全部屬小案例正確性檢查，不是效能或跨節點驗收。

## 保留的失敗與限制

最初的 uniform legacy 案例在單 rank 下 pressure L2 為
`3.9847761006167724e-15`；雙 rank 壓力差的 absolute L2 為
`4.625827780195998e-15`，relative L2 為 `1.1608752068855268`，**未通過**
既有相對門檻。以修改前 binary 重新跑雙 rank 得到完全相同的觀察；新舊版
在同 rank 數的場比較一致。兩種 legacy 模式各保留修改前／後的跨 rank
失敗，共四筆，列在 `acceptance.json` 及各 case 的 `cross_rank_observations`。
沒有放寬容許值、把近零的非零參考改成零，或宣稱 legacy 跨 rank 壓力驗收通過。

另外嘗試的非均勻 legacy 入口在修改前 binary 已未通過原生 Newton／mass gate；
補上 channel wall／interior labels 的版本仍有 mass gate 失敗。這些不是本次
可接受的數值基準，資料、命令與日誌保留於 `cli-profile/`、`cli-channel/`。
後續應在 legacy 數值／求解驗收中釐清此案例限制；本次只交付輸入可靠性與
有證據的相容性範圍，不以其他案例通過取代這些失敗。

首輪測試 fixture 漏填 periodic table 的必要 `units`，在修改前 parser 即遭
拒絕，保留於 `cli/`。另一次試圖只改 rank 1 的 `PETSC_OPTIONS` 環境以觸發
options mismatch，初始化後沒有產生預期的可見差異，反而進入正常求解；
該失敗保留於 `cli-rank-matched/`。最終將這個環境情境作為健康相容性觀察，
不將它列為有效 PETSc options 差異拒絕的證據；共用 options helper 的實際
差異測試仍見 [既有 options 報告](HPC_01C_PETSC_OPTIONS_PROGRESS.md)。

第一次 build 的暫存 aggregate 錯用了 velocity 容器型別，編譯失敗紀錄保留。
修正為既有 `vector<array<double, 3>>` 後，兩個正式 target 均重建成功，
`build-corrected.log` 無 compiler warning；Python 腳本語法與 diff whitespace
檢查通過。

## 重現與證據

環境：GCC 11.4、OpenMPI 4.1.2、PETSc 3.15.5 real64／Int32、MUMPS。
普通 MPI 測試固定 OpenMP／BLAS 為 1，使用 core mapping／binding。
CLI worker timeout 60 秒、report rendezvous 15 秒、job timeout 90 秒、
kill grace 5 秒；VCA smoke 外層 180 秒。所有預期非零退出保留原始值。

```bash
make -C solvers/cpu iga_navier_stokes iga_navier_stokes_openmp \
  PETSC_DIR=/usr/lib/petscdir/petsc3.15/x86_64-linux-gnu-real -j2
python3 scripts/hpc_flow_input_regression.py \
  --case-dir outputs/hpc01/assets/vca-cli/tubularflowiga-vca-3d-smoke \
  --reference-binary outputs/hpc01/flow-inputs/iga_navier_stokes-before \
  --output-dir outputs/hpc01/flow-inputs/cli-final
python3 scripts/hpc_resources_regression.py \
  --case-dir outputs/hpc01/checkpoint-write/verified/staged-failure/1-1 \
  --flow-case-dir outputs/hpc01/assets/vca-cli/tubularflowiga-vca-3d-smoke \
  --output-dir outputs/hpc01/flow-inputs/resources \
  --launcher 'mpiexec --map-by core --bind-to core'
```

output directories 以上述已完成證據為例，重跑須換新目錄。VCA／OpenMP 完整
argv、環境與日誌分別保留於 `vca-smoke/`、`openmp/`。`baseline.json` 與
`before/` 保留修改前來源和 binary 身分；`source-final.json`／tar 保存
142 份相關來源，兩個正式 binary 另行保存。
普通 flow binary SHA-256：
`d47342cb3ea69eccc4119a18da3b5bde982f46d90089600df4797b9e2cf1d06f`。
`acceptance.json` 明列成功範圍與四筆失敗的 legacy 壓力觀察。

後續 `iga_solve` 補強與驗收見 [傳輸 CLI 報告](HPC_01C_TRANSPORT_CLI_PROGRESS.md)。
HPC-01C 仍待其他 CPU CLI 資產檢查，以及 flow 建構後的
port catalogue、VCA identity／history 準備、每步本地輸入、輸出及其他 adapter／
executor 邊界。本次輸入雜湊也沒有建立完整 restart 身分，後者仍由 HPC-05 驗收。
