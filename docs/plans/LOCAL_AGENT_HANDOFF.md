# 本地開發接手：保存基線 → 整合 main → 分階段新架構

更新：2026-09-18。這是目前接手入口；舊的 FSI P0 執行交接文件是歷史紀錄，不是現在的待辦。

## 1. 目前目標與授權邊界

使用者決定轉到本地繼續開發，並要求保存、推送現有工作及提供交接 Prompt。
沿用同一個 repository，不另外拆新 repo。先保存既有成果、核對 main 的差異，
再透過隔離的 integration branch/worktree 整合；不要在研究基線分支直接進行大規模重寫。

本次交接不代表已完成 main 整合、建立 release/tag、通過整合回歸或授權一次完成整份 TODO。
本地接手第一輪先完成環境／來源核對和整合方案，回報具體下一步，再按使用者授權執行。
保留本地所有既有修改；不 reset --hard、clean、強推或覆寫遠端歷史。

## 2. Repository 與進度

- Remote：`git@github.com:EngineerEricXie/TubularFlowIGA.git`。
- 接手分支：`nextgen-phase-0-runtime-foundation`。
- `c53c2db`：已保存多尺度 [todo.md](../../todo.md) 與 [技術矩陣](../../TECHNOLOGY_MATRIX.md)。
- `5430071`：A 批，`.gitignore` 與 27 個精選 evidence 檔案。
- `00b78a5`：B 批（核心／測試／demo）；`87edb85`：C 批（benchmark 工具與歷史狀態）。
- 本交接文件所在版本為 D 批（歷史計畫／進度與本地交接）。
  用 `git log -6 --oneline` 核對實際 commit；以推送後的分支 tip 為接手版本，不停在 A 批。
- 2026-09-18 fetch 時，`origin/main` 是 `dfaede86692a70949565d254861d203f8b2f4012`。
  當時相對 A 批 HEAD，main-only 15 commits、研究分支-only 214 commits；之後保存批次會增加後者。
  這是分歧分支，不能假設 fast-forward。接手時重新 fetch 和計數。
- 沒有 merge main、沒有建立 tag、沒有改寫數值程式以配合本地環境。
- 本輪驗證是來源／證據 SHA-256、Git 差異及輕量靜態檢查，不是新一次模擬或整合回歸。

## 3. 已完成的 PSC 加速，不要重做 P0

PSC FSI 時間改善已完成 R7 最終交付，來源狀態以最新紀錄為準：

1. [精選證據說明](../evidence/bifurcation-fsi-r7/README.md)。
2. [加速結果與歷史紀錄](../progress/BIFURCATION_FSI_TIME_OPTIMIZATION_RESULTS.md)（按最新終態閱讀，不能停在早期 RUNNING 記錄）。
3. [原 run state](../../benchmarks/fsi-time-opt/20260915T035654Z/state.json)。
4. [原 run 敘述紀錄](../../benchmarks/fsi-time-opt/20260915T035654Z/STATE.md)。

最終 Job `46243889`：20 steps，solver 約 6968.953 s；19 項原等價性／輸出／封裝檢查通過。
歷史 13603 s 是跨 job／節點參考，不是配對效能證據。
**17/20 步守恆門檻仍失敗，formal_physics_acceptance=false，僅 visualization delivery。**
不得降低物理容差、改 dt、網格或積分規格來宣稱修復／加速。

`BIFURCATION_FSI_GPT_SOL_HANDOFF.md` 與 PSC 優化計畫保留當時「尚未實作」等歷史描述；
不能用它們覆蓋最新 state，也不要重新建立 P0 goal 或啟動舊 job 監測。
`monitor_job.sh` 只是保存的歷史工具，不是使用者目前同意的監測方式。

## 4. 新架構方向與未完成能力

完整閱讀 [todo.md](../../todo.md) 及 [TECHNOLOGY_MATRIX.md](../../TECHNOLOGY_MATRIX.md)。
新計畫是 **T0–T10**，與舊效能計畫 P0–P7 分開。

- 保留現有 0D、1D、貼體 IGA、immersed IGA、耦合與狀態管理。
- 新增 surface→貼體 FEM 體積網格／流體路線，以及 ALE。
- IGA 薄殼為特色主線；現有 P1 法向預張力膜不是 IGA 殼。
- 器官依物理需求選腔室流體、心肌固體、肝臟降階交換或 Darcy／多孔彈性，不以器官名称暗中推斷方程。
- ALE、通用 FEM 體積流程、IGA 殼、通用三維固體／多孔介質仍需實作與驗證。
- 首先 T0 現況／物理基準，再 T1 輸入契約、T2 固定壁對照；未決 backend／QoI 要先與使用者確認。

## 5. 本地第一輪工作清單

1. 讀取本地環境與 repo 的適用 `AGENTS.md`；核對作業系統、CPU/RAM、Git、編譯器、MPI、PETSc、Python、HDF5 等實際可用工具。
2. 記錄 `git status --short`、分支、HEAD、upstream；若工作樹不乾淨，不直接 pull、switch 或覆寫。
3. Fetch origin。已有 checkout 且分支可 fast-forward 時才 `git pull --ff-only`；若分歧先回報。
4. 在尚未改程式的基線上驗證下節兩組 hashes；缺檔／不匹配要說明，不能重寫清單消除失敗。
5. 閱讀 CPU build 文件與 Makefile，辨識 PSC 硬編碼路徑／module／Slurm 依賴；先提出最小本地建置方案。
   不直接執行舊 `run_r7_final.sh`，它依賴 Git 外的 source archives／歷史結果和 PSC 配置。
6. 重新核對 `origin/main...HEAD` 的 commit 與檔案差異，特別保留 main 的幾何修正、版本化輸出與 reduced-order 工作。
7. 回報隔離 worktree/integration branch 的起點、預期衝突、最小 regression gates 與成本。
   未獲使用者確認前不要合併到 main，也不要開始大規模重寫。

新 clone 可使用（本地已存在同名目錄時不要覆寫）：

```bash
git clone --branch nextgen-phase-0-runtime-foundation git@github.com:EngineerEricXie/TubularFlowIGA.git
cd TubularFlowIGA
git status --short
git log -6 --oneline
```

## 6. 可攜完整性檢查與建置限制

Repository 根目錄：

```bash
(cd docs/evidence/bifurcation-fsi-r7 && sha256sum --check SHA256SUMS)
sha256sum --check --quiet docs/evidence/bifurcation-fsi-r7/provenance/R7Final-source-files.sha256
```

前者核對精選證據，後者核對 442 個已測来源項目。macOS 若缺 `sha256sum` 可用 `shasum -a 256 -c`；
跨平台 checkout 換行轉換也可能影響 hashes，先診斷 Git EOL，不改寫原始 evidence。
正常開始重寫後來源 hashes 不匹配是預期，不能再宣稱來源與 R7 相同。

建置入口：[CPU README](../../solvers/cpu/README.md)、[CPU Makefile](../../solvers/cpu/Makefile)、
[benchmark runner](../../benchmarks/fsi-time-opt/runner.mk)、[MPI test wrapper](../../benchmarks/fsi-time-opt/mpi_test_wrapper.cpp)。
PETSc／MPI 與 C++ 編譯鏈必須一致；多執行緒 MPI 測試保留 `MPI_THREAD_FUNNELED` gate。
不要假定一個 `make all` 包含全部優化 regression，也不要一開始全量編譯／重跑 20 步。

使用者轉到本地主要是接續開發與整合。先前「所有計算與大型編譯走 PSC Slurm」的約束
不可因遷移自動當作取消：本地先做讀寫／Git／輕量靜態檢查；開始本地數值測試或大型建置前，
確認使用者同意的本地資源預算與例外。PSC 登入節點仍只能協調，PSC 計算走 Slurm（通常 `mch260002p`）。
若需 PSC Python，依實際環境使用 `module load anaconda3`，不要在本地機器盲目套用 module 指令。

## 7. Git 沒有搬走的資料

PSC 原 repo：`/ocean/projects/mch260002p/thsieh1/TubularFlowIGA`。
Git 只帶走來源、測試、工具、文件及精選 evidence，不包含完整模擬輸出、執行檔或 source tar archives。

定位索引：[external-artifacts.json](../evidence/bifurcation-fsi-r7/external-artifacts.json)。
它列出 49 個外部檔案與 2 個 evidence 目錄；其中有 45 份 source archives、最終 ParaView package、
執行檔等。目錄和檔案有重疊，不能計為獨立備份。**索引不是備份，獨立備份仍未驗證。**
歷史 `/local/...` 暫存路徑可能失效；PSC `/ocean/...` 路徑不是本地路徑。
不要刪除 PSC 原資料。若確需重跑／取回完整場，先選定最小必要 artifacts、確認 PSC SSH host／目的地／容量，
再傳輸並核對 hashes；本輪沒有下載這些大型資料到使用者電腦。

## 8. 等待與回報偏好

若後續提交 Slurm 工作且無其他工作可做，用**對話端 adaptive sleep**，依估計剩餘時間調整；
不要在 SSH connection sleep，不啟動背景監測腳本，不每分鐘查狀態，不重讀相同日誌。
遵守當前工具允許的等待上限；無法低成本等待時記錄 Job ID／狀態並交接，不能假稱持續監測。
回報時區分歷史證據、本輪新驗證、尚未驗證和物理限制。
