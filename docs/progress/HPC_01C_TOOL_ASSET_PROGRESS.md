# HPC-01C：MPI 工具的資料庫一致性與 owner 漏檢

日期：2026-09-08。狀態：本批指定驗收完成；HPC-01C 及整份清單尚未完成。
接續 [Legacy transport 報告](HPC_01C_LEGACY_TRANSPORT_PROGRESS.md)。
本批只修改兩個 CPU 工具、測試、Makefile dependencies 與文件，沒有改變
Jacobian／組裝公式、packed database 格式或其他 solver 的實作。

## 問題與修正

`iga_mesh_check` 與 `iga_assembly_smoke` 原先不比對 rank-local `.ntiga`
內容。兩者現在在開檔解析前讀取 fingerprint，並在分散式檢查／組裝前要求
內容一致；相同內容的副本可使用不同路徑。缺檔、FIFO、目錄與 malformed
header 會共同拒絕，不讓其他 ranks 留在後續 collective 等待。
完整 fingerprint 增加每 rank 的輸入掃描；輸入需在執行期間保持不變，
這不是檔案鎖，也不允許在各副本中省略「當前 rank 沒用到」的不同資料。

Mesh checker 額外修正 owner 漏檢。舊版在以下三種情況都會退出 0 並印出
`elements=1 minimum_detJ=inf bad_elements=0 bad_samples=0`：

1. 兩份資料庫各自合法，但 owner 指派互相矛盾，結果兩個 ranks 都沒載入元素。
2. 所有副本一致，但 ownership index 的 owner 超過 communicator 範圍。
3. 所有副本一致，但 ownership index 與 element record 的 owner 不一致。

第一種情況由 fingerprint agreement 拒絕；後兩者新增 owner 範圍及載入記錄
的 owner 檢查。正常的空 rank 仍被接受：同一份合法資料庫在 2 ranks 中
只有一個 owner 也能通過。這不是完整 `.ntiga` 格式驗證器，未宣稱涵蓋所有
offset、ID、connectivity 或 extraction 的可能損毀。

Assembly smoke 在本體中安裝 scoped PETSc return-error handler，並比對
有效 PETSc options，讓已檢查的 API 錯誤能返回共同錯誤協議。
正常矩陣清理移到結果摘要之前；MatrixOwner 不可複製。兩個工具的最後
logging 也放入 local stage。沒有在 local callback 內包住 PETSc collective。

## 驗收

Ignored 證據位於 `outputs/hpc01/tool-assets/`。

| 驗收 | 結果與範圍 |
|---|---|
| 新增 native regression | 38 筆案例、67 份 rank reports，全部符合預定退出碼與診斷 |
| 健康相容性 | 兩個工具修改前／後各執行 1、2 ranks，結果摘要精確相同；使用不同 rank-local 路徑 |
| 資產錯誤 | 每個工具均拒絕單 rank 的內容差異、不同幾何、缺檔、FIFO、目錄及 malformed header |
| Owner 漏檢 | 上述三種錯誤在舊版重現空 coverage 的成功誤報；新版各自在指定階段共同退出 1 |
| 真實 PETSc 回報 | 對同一矩陣先 ADD 再 INSERT，實際 PETSc 回傳 error 73；新版 1、2 ranks 皆共同回報 element insertion 錯誤，無結果摘要 |
| Returned-error 注入 | assembly end、MatGetInfo、MatDestroy 在真實操作返回後對 rank 0 注入錯誤；新版 1、2 ranks 皆共同退出 1，沒有先印結果 |
| 舊清理順序 | 舊版 1、2 ranks 的 MatDestroy 錯誤雖退出 1，卻已印出 `global_rows=...`；兩份重現保留 |
| 既有工具回歸 | `hpc_tools_regression.py` 的 39 項全數通過，含 legacy transport、非法 fields、rank／resource 差異與退化幾何 |
| Mesh 核心回歸 | `make mesh-test` 退出 0，`mesh_core_test: PASS` |

健康 mesh 的 minimum determinant 為 `0.125`；assembly 的 global rows 為
`128`、nz_used／nz_allocated 均為 `16384`、mallocs 與 missing diagonal 為 0。
原本的退化幾何仍退出 2，minimum determinant=0、bad elements=1、bad samples=64，
沒有把幾何品質不合格混同為輸入協議失敗。
同版 legacy transport 場檔仍與封存參考逐位元組相同；其 1／2 ranks relative
L2 為 `3.9648210454771405e-13`，符合既有 `1e-6` 門檻。

新 native 測試使用 Linux／Open MPI 專用的 `LD_PRELOAD` shim；controller
明確指定故障模式，production code 沒有故障開關。混合 insertion mode 使用
真正 PETSc 錯誤；其餘三類在 API 真實操作完成後才更改返回碼，避免把已卡在
collective 內的狀態當作可恢復案例。每個預期失敗均有限時間退出，沒有以
timeout 作為通過。這不提供程序失聯恢復，也不證明任意 PETSc 內部部分失敗
或析構失敗都可恢復。

## 重現命令與來源

環境：WSL 工作站、GCC 11.4、OpenMPI 4.1.2、PETSc 3.15.5 real64／Int32；
OMP／BLAS 各 1。建置啟用 warnings，沒有新增 compiler warning。
在原生 MPI 授權環境執行，未更改系統 MPI／PETSc 安裝。

```bash
make -C solvers/cpu iga_mesh_check iga_assembly_smoke \
  PETSC_DIR=/usr/lib/petscdir/petsc3.15/x86_64-linux-gnu-real -j2
mpicxx -std=c++17 -O2 -Wall -Wextra -Wpedantic -fPIC -shared \
  -I/usr/lib/petscdir/petsc3.15/x86_64-linux-gnu-real/include \
  solvers/cpu/tests/smoke_failure_preload.cpp -ldl \
  -o outputs/hpc01/tool-assets/smoke-failure-preload.so
python3 scripts/hpc_tool_asset_regression.py \
  --fixture outputs/hpc01/tools/accepted/fixture \
  --baseline-dir outputs/hpc01/tool-assets/before-binaries \
  --preload outputs/hpc01/tool-assets/smoke-failure-preload.so \
  --output-dir /tmp/tool-assets-new
python3 scripts/hpc_tools_regression.py \
  --fixture-source outputs/hpc01/tools/fixture \
  --baseline-dir outputs/hpc01/tool-assets/before-binaries \
  --output-dir /tmp/tools-compatibility-new
```

結果目錄必須是新目錄。Fixture 是已驗證的小型 v5 database，1 element、
64 nodes，包含 serial／2-rank 版本；測試在副本中改動內容，原始 fixture
保持不變。Native 作業有外層 90 秒／每 rank 60 秒 timeout，保存 rank 身分、
退出碼、資源與 log 雜湊。

最終結果為 `native/summary.json`、`tools/summary.json`，來源／binary／input
與日誌稽核見 `acceptance.json`。`before/`、`before-binaries/` 保留舊版，
`source-final.json`／`source-final.tar.gz` 與 `after-binaries/` 保留測試版本。
Git HEAD 加未提交工作樹構成本批版本；完整來源依雜湊辨識。

## 接續範圍

接續檢查 flow runtime、其餘 adapter／executor 的配置與錯誤邊界，並完成
支援入口的覆蓋稽核，再判斷 HPC-01C 是否可完成。HPC-01D 的型別／能力矩陣
仍另行追蹤。沒有新增 GPU、ParaView、跨節點、效能或完整 binary format
驗收；本批小案例的 MPI 成功不能當作擴展性證據。
