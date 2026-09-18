# 分岔 FSI：R7 重寫前精選證據

整理日期：2026-09-18。數值工作完成日期：2026-09-17。Slurm Job：`46243889`。

本目錄保存**小型輸入、驗證摘要及來源指紋**，供重寫前基線追溯。原始檔案沒有搬移、修改或刪除；副本與原始路徑的關係記錄在 [provenance.json](provenance.json)。

## 狀態：不是正式物理驗收

- 求解器完成 20 步，solver elapsed 為 6968.953 s。
- 原 R7 的 19 項等價性／輸出／封裝等檢查通過。
- 原有 **17／20 步守恆門檻失敗**仍保留，`formal_physics_acceptance=false`。
- `final/status.json` 與 `final/correctness.json` 的 `passed` 是該工作的執行／既定檢查狀態；不能取代 `final/results/run.json` 的 `visualization_only`。
- 本次只整理資料與核對小型檔案，**未重跑模擬、未執行新物理驗證、未驗證整合後 main**。
- 歷史 13603 s 與 R7 時間是跨 job／跨節點比較，不是配對完整週期 speedup。

詳細结果：[原加速報告](../../progress/BIFURCATION_FSI_TIME_OPTIMIZATION_RESULTS.md)。未來設計：[TODO](../../../todo.md)、[技術矩陣](../../../TECHNOLOGY_MATRIX.md)。

## 保存內容

| 路徑 | 用途 |
|---|---|
| [inputs/surface.txt](inputs/surface.txt) | 實際測試的封閉、帶標籤三角表面，不以重新生成幾何代替原輸入 |
| [final/case.json](final/case.json) | 原 R7 設定、歷史比較檔案雜湊、容許值與已知限制 |
| [final/manifest.json](final/manifest.json) | Job、來源封存、binary、surface、runner、case 指紋及環境定位 |
| [final/correctness.json](final/correctness.json) | 19 項原檢查及物理驗收旗標 |
| [final/status.json](final/status.json) | 原工作終止狀態與時間 |
| [final/timing.json](final/timing.json) | 時間、互斥／包含式 profile、組裝與預條件計數 |
| `final/results/` | history、coupling、run 與 heartbeat／VTK 驗證摘要；**沒有完整場或動畫檔** |
| `historical/` | Job 46243889 當時保存的五份歷史參考小檔，可對照 case 內的凍結 hashes |
| `provenance/` | 原來源清單、封存／surface hashes、來源 base HEAD、B0 manifest |
| [provenance.json](provenance.json) | 23 份逐位元副本的原始路徑、大小與 SHA-256，以及本次整理範圍 |
| [external-artifacts.json](external-artifacts.json) | 大型封存／binary／原命令及外部 evidence 目錄索引；**不是備份完成證明** |
| [SHA256SUMS](SHA256SUMS) | 本目錄所有資料／說明檔的可攜相對路徑 checksum，不包含本檔自身 |

原始 provenance 內的 `/local/...` 可能已不存在；原 PSC 絕對路徑亦不保證在別台機器有效。這些是歷史身分資料，故未改寫；現有位置與副本位置請看本目錄兩份 JSON 索引。

`solver-options.txt` 原件是沒有結尾換行的命令字串，本輪僅將其列入外部索引，不把補換行的版本冒充逐位元副本。`final/case.json` 與原 runner 仍保留參數定義。

## 輕量完整性檢查

從 repository 根目錄執行：

```bash
(cd docs/evidence/bifurcation-fsi-r7 && sha256sum --check SHA256SUMS)
```

這只驗證精選檔案沒有變動，不代表原大型封存仍完整，也不證明數值物理正確。

若要檢查工作樹中的已測來源是否仍與 R7 一致，可從 repository 根目錄執行：

```bash
sha256sum --check --quiet docs/evidence/bifurcation-fsi-r7/provenance/R7Final-source-files.sha256
```

該清單共 442 個項目，包含未必已提交的原始 source／tests／benchmark 工具。新 clone 若缺檔，表示尚未取得所有基線來源；不要因此刪除檢查或將缺檔當作通過。未來正常重寫後不匹配是預期現象，清單不得隨新程式覆寫。

## Git 與外部保存邊界

適合進 Git：本目錄的精選 surface、JSON、CSV、hash 清單與說明。原 job 目錄繼續忽略，不整批 `git add -f`。

留在 Git 外：

- 45 份 source archive，合計 221214720 bytes，約 211 MiB。
- 最終 `heartbeat-paraview.tar.gz`，93174492 bytes，約 89 MiB。
- 原執行檔、原始大量日誌、完整場、所有中間／失敗工作證據。
- 比較 workload archive 與歷史完整結果。

[external-artifacts.json](external-artifacts.json) 保存可定位的現有路徑、檔案大小與**既有** checksum。建立索引時未全面重讀大型資料；這些 hashes 不應被描述為本輪重新驗證。

原 evidence tree 與索引內個別檔案有重疊，不能把兩者當成獨立備份。本輪沒有第二份獨立儲存副本，**不要刪除原目錄**。後續若安排備份，須另確認目的地、容量與完整性檢查；所有大型處理依 PSC 規則走 Slurm。

## 如何理解「可重現」

目前保存的是可追溯基線，不是 clone 後一條命令即可重跑的正式套件。完整重跑仍依賴：

1. R7 凍結來源、runner／case 與實際表面輸入。
2. PETSc／MPI／編譯器環境，以及 VTK 驗證用 Python 環境。
3. 原 runner 要求的歷史參考目錄、source archive 及 checksum。
4. Slurm allocation；不可在登入節點直接編譯大型目標或執行模擬。

原入口：[run_r7_final.sh](../../../benchmarks/fsi-time-opt/run_r7_final.sh)、[cases-r7-final.json](../../../benchmarks/fsi-time-opt/cases-r7-final.json)。這些檔案與原數值來源在本輪整理前仍有未提交項目；後續 B／C 批需要一併保存。

不要直接執行本目錄的原 `*.sha256` 並期待歷史絕對路徑都存在。本目錄的本地驗證入口是上方 `SHA256SUMS`；原 hashes 用來驗證取回並定位好的原始 artifacts。

## 下一階段，而非本輪已完成的工作

- 保存所有必要數值來源與 regression／runner，再建立 Git 基線 tag。
- 標註歷史 runner 的環境假設與舊監測腳本；本目錄不會啟動監測。
- 規劃獨立備份與大型 artifact 的完整性稽核。
- 整合 main 後重新驗證；本次副本檢查不能替代整合測試。
- 獨立處理正式物理守恆限制，不能透過改變文件狀態使其「通過」。
