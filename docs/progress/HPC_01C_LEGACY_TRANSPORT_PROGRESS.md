# HPC-01C：Legacy transport 的輸入、時間步與清理

日期：2026-09-08。狀態：本批指定驗收完成；HPC-01C 及整份清單保持未完成。
接續 [其餘 MPI 工具報告](HPC_01CD_TOOLS_PROGRESS.md)，本次修改限定於
`iga_transport` 及其測試／文件，沒有更改共用物理元素核心或其他 solver。

## 修正

原有 CLI 會共同處理部分讀取錯誤，但沒有比對各 rank 實際讀到的內容。
在同一個雙 rank 問題中，把 rank 1 的 `D` 從 `0.1` 改為 `0.2`，舊版仍
退出 0 並印成功摘要；與一致輸入的場相比，relative L2 為
`0.2378188687660649`。這份錯誤結果保留在驗收證據內。

新版在開啟資料庫前讀取 fingerprint，並在解析前比對 parameters、mesh、
實際使用的 velocity 及 optional `case_config.json`。邏輯資產名稱相同且內容
相同的 rank-local 副本可以使用不同路徑；optional config 是否存在也需一致。
必要資產的 FIFO／目錄會先被拒絕；顯式 velocity override 的未使用預設檔案
不會開啟。輸入仍須在執行期間保持不變，fingerprint 不提供檔案鎖。

同時加入有效 PETSc options 的精確比對；return-error handler 涵蓋 solver
本體，而非只有 setup／solve。元素積分與插入、boundary rows／values 配置
放在 local stage；PETSc collective 在 callback 外執行，返回後協調狀態。
矩陣組裝、施加邊界、vector 操作、KSP 設定／求解、gather 與正常清理皆檢查
回傳碼。這不能解救已卡在 PETSc／MPI collective 內的程序。

頂層 Mat／Vec／KSP 改由不可複製的 owner 管理，正常路徑明確清理後才印
成功摘要；錯誤退出使用相同順序的 best-effort 清理。輸出仍 gather 到 root，
檢查 vector 尺寸並拒絕非 regular target，保留原有 close 後的寫入檢查。
本體抽成檔案內的 `RunLegacyTransport`，使用傳入的 borrowed communicator；
原生 main 仍以 `PETSC_COMM_WORLD` 呼叫，未新增公開 runtime API。

位置參數、預設 GMRES／block Jacobi、第一步零 initial guess／後續步非零
guess、時間步公式、兩個場的順序與輸出格式保持原有語義。
固定的邊界值只建立一次，初始化改以 owned-row `VecSetValues` 加 assembly
填入相同數值；每一步重用該陣列。總迭代數使用 checked 64-bit 無號累計，
避免以 PetscInt 累加時溢位。組裝與求解時間區間保留，新增的資產掃描在其外。

## 驗收

Ignored 證據根目錄：`outputs/hpc01/legacy-transport/`。

| 測試 | 結果與範圍 |
|---|---|
| 原生函式故障與重試 | `legacy_transport_failure_test` 在 world 2 ranks 及兩個獨立單 rank 程序群共驗證 46 個故障、46 次健康重試 |
| 故障類型 | local insertion 的 allocation exception／PETSc 錯誤、matrix assembly、MatZeroRows、initial／boundary values、MatMult、VecCopy、KSP options／solve、VecNorm、scatter、array read／restore、KSP destroy，以及有效選項差異 |
| 物件清理 | 對 native 顯式建立的 Mat、Vec、KSP、scatter 保留測試引用；退出後先釋放 parent，再檢查每個資料物件僅剩測試引用，沒有遺留 native ownership |
| 獨立程序群 | 使用不同入口值；其中一群多執行一次健康求解，所有故障與重試仍各自完成；沒有用相同呼叫次數掩蓋 world collective |
| Native CLI | 43 筆執行、82 份 rank reports；含修改前／後 1、2 ranks、optional config、velocity override、未使用的 FIFO、零步、無輸出及輸入／選項／輸出錯誤 |
| 舊行為重現 | 保留舊版不一致參數仍成功、場比較不通過的結果；新版相同不一致輸入在 `legacy transport asset parameters` 階段共同拒絕 |
| 既有工具回歸 | 39 項通過；包含三個 MPI 工具的輸入與資源錯誤、geometry、matrix 摘要、legacy field、preonly 不相容案例 |

Native CLI 的 11 份相容性場比較通過 relative `1e-6`／零參考 absolute
`1e-12`，最大 relative L2 為 `4.1385674506546456e-13`；其中包含 baseline
自比，並非 11 份獨立參考。另有一份刻意不一致的舊版場比較失敗，不能計為
數值相容通過。既有工具測試的 legacy 場同 rank 修改前後逐位元組相同，
1／2 ranks relative L2 為 `3.9648210454771405e-13`。沒有放寬任何數值門檻。

函式測試用 GNU linker `--wrap`：local 操作可在呼叫前回報錯誤；collective
則在所有 rank 的真實操作正常返回後，改變指定 rank 的返回碼。
這驗證錯誤返回後的協調與清理，不宣稱 PETSc 內部部分失敗能恢復。
每次重新求解都建立新的 runtime objects，與同群組健康場逐位元組比較；
這不是從失敗的時間步原地續算。retained-reference 檢查針對觀測的頂層
ownership，不能代替對 PETSc 內部所有配置的 memory leak 稽核。

## 測試修正與保留紀錄

初次封存 binary 使用 `copyfile` 而未保留執行權限，導致 `tools/` 與 `cli/`
在 solver 執行前失敗，legacy runner 退出 126。補回原執行權限後使用新結果
目錄重跑；沒有覆蓋那些失敗紀錄。

`cli-verified/` 曾預期 rank 1 環境中的 `PETSC_OPTIONS` 差異會被拒絕，實際
卻退出 0。在此 PETSc 3.15.5 安裝上，初始化之後捕捉的有效 options 一致，
雖然 rank reports 仍記錄不同環境字串。因此該情境改記為
`petsc-env-different` 相容性測試；另在 MPI 初始化後直接修改 rank 1 的實際
options database，驗證共同拒絕、清除差異後同群組重試及場一致性。
前者不能當作有效選項差異的拒絕證據。

最終驗收位於 `unit-final.log`、`cli-accepted/summary.json`、
`tools-accepted/summary.json`，使用 `build-final.log` 對應的同版 binary。
`acceptance.json` 保存稽核結果；`source-final.json`／`source-final.tar.gz`
保存來源與雜湊，`after-binaries/` 保存最終執行檔。Git HEAD 加未提交工作樹
才是本批版本，不能僅用 HEAD 當作完整來源識別。

## 重現命令

環境：WSL 工作站、GCC 11.4、OpenMPI 4.1.2、PETSc 3.15.5 real64／Int32、
MUMPS；OMP／BLAS 各 1。建置開啟 warnings，沒有新增 compiler warning。
這些小案例可重疊執行作正確性驗證，不能當作隔離效能量測。

```bash
make -C solvers/cpu iga_transport legacy_transport_failure_test \
  PETSC_DIR=/usr/lib/petscdir/petsc3.15/x86_64-linux-gnu-real -j2
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 \
  timeout --kill-after=5s 180s mpiexec --map-by core --bind-to core -np 2 \
  solvers/cpu/legacy_transport_failure_test \
  outputs/hpc01/tools/accepted/fixture /tmp/legacy-unit-new
python3 scripts/hpc_legacy_transport_regression.py \
  --fixture outputs/hpc01/tools/accepted/fixture \
  --reference-binary outputs/hpc01/legacy-transport/iga_transport-before \
  --output-dir /tmp/legacy-cli-new
python3 scripts/hpc_tools_regression.py \
  --fixture-source outputs/hpc01/tools/fixture \
  --baseline-dir outputs/hpc01/legacy-transport/tools-before \
  --output-dir /tmp/legacy-tools-new
```

所有結果目錄都必須是新目錄。第一個 fixture 已將 legacy 不支援的 label 3
映射為 2；最後一個 tool harness 自行在 raw fixture 副本上做同樣適配。
fixture 包含 1 個元素、64 個 nodes、serial／2-rank 的 packed databases，
因此只能作小型正確性驗證。所有 native CLI 有外層 90 秒及每 rank 60 秒
timeout，預期失敗皆為共同非零退出而非 timeout。

## 下一步與限制

接續檢查 `iga_mesh_check`／`iga_assembly_smoke` 的完整資產一致性，以及
flow runtime、其餘 adapters／executor 的錯誤邊界。HPC-01C 仍未完成。
本批沒有更改 legacy 共用 parser 的格式政策，也沒有驗證所有極端浮點／
整數輸入；型別與能力矩陣繼續由 HPC-01D 追蹤。
PETSc/MPI 初始化本身、程序失聯、collective 內部停住、析構重試失敗、完整
checkpoint 發布、分散式輸出與跨節點擴展性不屬於本批的完成聲明。
