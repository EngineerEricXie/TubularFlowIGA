# HPC-01C：原生 1D CLI 外部資產一致性

日期：2026-09-08。狀態：本項實作與最終版回歸通過；HPC-01C 整體保持未勾選。

## 變更與行為

`iga_1d` 原先只比對設定文字，各 rank 的相同相對檔名仍可能指向不同 network、
waveform 或 replay。現在沿用 `CollectiveAssetInput.hpp`，先比對 SWC／OBJ
network，再讀取拓撲；建立 runtime 後依實際選定的 flow 與 transport species
蒐集 periodic table，以及啟用的 replay CSV／JSON，於初始狀態及輸出前比對。
邏輯 key 不含本地路徑，允許內容相同的 rank 本地副本。未使用的 temporal
definition 維持不開檔的既有行為。replay parser 的失敗由獨立 local stage 協調。

設定檔本身也先透過非阻塞 descriptor 檢查 regular-file，避免 parser 開啟 FIFO
時等待 writer。所有上述驗證同樣適用於 `--check`。這不修改物理公式、KSP/PC
預設、數值容許值、檔案格式或 checkpoint metadata。

案例輸入必須在執行期間保持不變；雜湊檢查與後續 parser 之間沒有鎖或快照，
不承諾抵抗刻意替換檔案的競態。每 rank 額外完整掃描資產，設定文字也多一次
regular-file／雜湊讀取。共享 helper 使用 Linux／WSL POSIX 介面及固定 64 KiB
buffer；本次沒有大型共享檔案系統的啟動 I/O 效能驗收。

## 最終版驗收

環境沿用工作站 GCC 11.4、OpenMPI 4.1.2、PETSc 3.15.5 real64／Int32 與 MUMPS。
測試固定 OpenMP／BLAS 為一執行緒，使用 core mapping／binding。child timeout
60 秒、job timeout 90 秒、kill grace 5 秒；既有多 rank report rendezvous 為
15 秒。以下測試可能在時間上重疊，屬正確性回歸，不是獨立效能量測。

| 測試 | 觀察數 | rank reports | 結果 |
|---|---:|---:|---|
| `hpc_one_d_asset_regression.py` | 35 | 85 | 通過 |
| `hpc_one_d_cli_regression.py` | 22 | 58 | 通過 |
| `hpc_one_d_checkpoint_regression.py` | 25 | 75 | 通過 |

資產測試包含五種案例的修改前單 rank、修改後單 rank 與三 ranks：flow table、
species table、CSV replay、JSON replay、OBJ network。五筆修改前執行是參考
觀察，不算獨立前後比較。其餘健康解與修改前的完整 hydraulic／species CSV
欄位比較最大誤差為 0；仍沿用 CPU relative L2 `1e-6`、零參考 absolute `1e-12`。

故障涵蓋六種有效但內容不同的 network／table／replay、缺檔、network／設定
目錄、network／設定／table／replay FIFO、所有 rank 相同的 malformed replay。
另測未使用的缺檔 table，以及 `--check` 成功／內容不同／FIFO。故障保留原始
非零退出碼、沒有 timeout，且不建立求解輸出目錄。

既有 CLI 測試中原本只在 rank 1 改變未來 oxygen table 的案例，現在正確地在
啟動期拒絕。因此保留它作 `step-inlet-assets`，另以所有 ranks 相同、初值合法
而未來濃度非法的 table 保留真正的 `step-inlet` 故障：初始輸出存在，step 1
輸出不存在。這是新增覆蓋，沒有把原本時間步故障改成僅測啟動。
checkpoint 回歸保留不同 rank 數續跑、舊格式、損毀與寫入失敗等既有 gates。

```bash
make -C solvers/one_d petsc \
  PETSC_DIR=/usr/lib/petscdir/petsc3.15/x86_64-linux-gnu-real
python3 scripts/hpc_one_d_asset_regression.py \
  --output-dir outputs/hpc01/one-d-assets/assets-final \
  --reference-binary outputs/hpc01/one-d-assets/iga_1d-before
python3 scripts/hpc_one_d_cli_regression.py \
  --output-dir outputs/hpc01/one-d-assets/cli-final \
  --launcher 'mpiexec --map-by core --bind-to core'
python3 scripts/hpc_one_d_checkpoint_regression.py \
  --output-dir outputs/hpc01/one-d-assets/checkpoint-final \
  --launcher 'mpiexec --map-by core --bind-to core'
```

以上 output directories 是已完成證據位置；重跑須使用新的目錄。編譯完成且無新增
warning，三份 Python regression 腳本通過語法檢查。

## 證據與剩餘範圍

`outputs/hpc01/one-d-assets/acceptance.json` 彙整最終 82 筆觀察、218 份 rank
reports。各 suite 保留 command argv、原始日誌、退出碼及 SHA-256，已重讀所有
rank report 並核對日誌雜湊。`source-final.json`／`source-final.tar.gz` 保存
147 份相關來源，`iga_1d-final` 保存最終 binary；`baseline.json`、`before/`
與 `iga_1d-before` 保存修改前 provenance。最終 binary SHA-256：
`0cb8d2f4bc0d004f7edcc86e55e3da38a38ba2e10b6ac9c4f8bb0bb4221d3cb1`。

較早的 `assets/`（33 項）、`cli/`、`checkpoint/` 與
`source-before-config-probe.tar.gz`／`iga_1d-before-config-probe` 保留在加入
設定檔 regular-file probe 前的證據，不將它們冒充最終版驗收。

本次比對僅證明當次 communicator 內資產一致，未將 waveform／replay 身分
綁定既有 checkpoint；完整持久化契約仍由 HPC-05 驗收。其餘獨立 CPU CLI
輸入與 adapter／executor 邊界仍是 HPC-01C 的剩餘工作。
