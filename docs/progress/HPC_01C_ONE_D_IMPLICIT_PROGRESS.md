# HPC-01C：1D implicit PETSc 呼叫與 callback 錯誤協調

日期：2026-09-08。狀態：**本報告範圍通過；HPC-01C 整體仍部分完成。**

## 變更

四種 1D 隱式方法的 graph preparation、owned-row assembly、convergence 與
state update 現在協調本地例外。`AdvanceImplicitOneD` 在分支前比較 formulation，
矩陣建立前比較 unknown count；graph builder 拒絕無效端點、cell 數及超出
既有 int 索引範圍的布局。步長必須有限且為正，inlet flow 必須有限。
這些檢查不替代完整外部幾何與配置身分的一致性驗證。

原有未檢查的 Mat／Vec 插入、KSP／SNES 設定、求解與收斂查詢均增加回傳碼
處理。本地插入在 local stage 中檢查；collective 在所有成員返回後才協調
回傳碼，不將 collective 包進本地 work callback。成功路徑的操作名稱使用
`string_view`，避免為每次檢查建立字串。

兩個 SNES callback 現在為 `noexcept`。它們先協調局部 residual／Jacobian
組裝錯誤，保存共同的 `exception_ptr`，再回傳 `PETSC_ERR_USER`；C++ 例外
不跨越 PETSc 的 C 呼叫堆疊。`SNESSolve` 返回後，入口優先回報保存的 callback
診斷，並檢查其他 PETSc 錯誤與 SNES convergence reason。

`OneDPetscSupport.hpp` 提供借用 communicator 的檢查與 stack-owned PETSc
資源。已成功建立的 solver、matrix、vector、scatter 在共同失敗時依序釋放。
vector gather 的本地讀取使用 RAII array view，配置／複製失敗時也先 restore
再協調；成功後才 swap 到輸出。初始 vector 的長度錯誤在 VecAssembly 前拒絕。

各求解器在暫存 flow state 中產生更新，成功協調並關閉 PETSc 資源後才發布。
`AdvanceImplicitOneD` 的 inlet 更新也納入此規則。受控錯誤不改動呼叫者的
原有狀態，可解除故障後重新求解。沒有修改離散公式、KSP／PC 預設、
時間積分、耦合方法或 checkpoint／場檔案格式。

## 驗收與重現

環境：WSL、Intel i9-14900KF、GCC 11.4、Open MPI 4.1.2、PETSc 3.15.5
real/double、32-bit PetscInt、MUMPS；HEAD
`ee8da2ab528849814b6d12c8186476b5f7f59124` 加未提交工作樹，精確來源以 inventory 為準。
MPI 一般使用每 rank 1 OpenMP／OpenBLAS thread；trial unit 另測 4 OpenMP
threads 且關閉 dynamic teams。使用既有 CPU relative L2 `1e-6`、零參考
absolute L2 `1e-12`；restart clock 沿用既有 `1e-12` 規則，沒有放寬門檻。

```bash
make -C solvers/one_d core-test iga_1d one_d_petsc_test \
  one_d_implicit_failure_test one_d_subcommunicator_test one_d_trial_failure_test \
  PETSC_DIR=/usr/lib/petscdir/petsc3.15/x86_64-linux-gnu-real -j2
make -C solvers/coupling petsc multidomain_subcommunicator_test \
  PETSC_DIR=/usr/lib/petscdir/petsc3.15/x86_64-linux-gnu-real -j2
make -C solvers/one_d implicit-failure-test MPIEXEC='mpiexec --oversubscribe' \
  PETSC_DIR=/usr/lib/petscdir/petsc3.15/x86_64-linux-gnu-real
python3 scripts/hpc_one_d_cli_regression.py --output-dir /tmp/hpc-implicit-cli-new
python3 scripts/hpc_one_d_checkpoint_regression.py --output-dir /tmp/hpc-implicit-checkpoint-new
```

`implicit-failure-test` 的實際執行命令由 `run_units.py` 記錄，包含 180 秒
timeout、3 ranks 與原生 straight SWC。其測試以 world 及獨立 1+2 groups
執行，單 rank group 一次、兩 rank group 兩次，檢查不同 collective 次數無串擾。

| 測試 | 實際結果 |
|---|---|
| 新 implicit solve 故障與 retry | 四種方法，各測 step、graph、assembly、publication、layout、method、PETSc options、KSP convergence；world 32、單 rank group 24、兩 rank group 64 次通過 |
| SNES callback | residual 與 Jacobian 分別透過實際 `SNESSolve` C stack 觸發單 rank 短缺 old-flow；共 8 次 group 執行，均回傳錯誤並保存共同診斷 |
| vector 初始化 | 共 4 次 group 執行；單 rank 長度錯誤共同拒絕，同一 Vec 修正輸入後可重用，gather 精確一致 |
| 失敗不發布狀態 | packed fields、inlet、outlet node、step／time／internal-substeps 保持原值；retry 逐物理欄位比較 COMM_SELF 參考通過 |
| 四種 implicit subcommunicators | 1+3 groups 的 area／flow／pressure／node pressure／segment flow 均通過，最大 relative L2 `1.16694e-15` |
| 既有 combined／staged trial | 3 ranks、每 rank 1／4 threads 均通過，含 rollback／retry 與 NaN 物種故障 |
| PETSc／checkpoint unit | 1 與 3 ranks 通過 |
| C++ core／coupling／runtime | 全部通過，核心 runtime 仍不需要 MPI 相依 |
| 原生 CLI | 21 項全部通過，含正常單／三 rank 與 13 項一致失敗案例 |
| 原生 checkpoint | 25 項全部通過，含三個物理案例續跑、14 項故障及 legacy PETSc binary 相容性 |
| 最終 flow／species graph | 獨立 1+2 groups 通過，最大 relative L2 分別 `7.45082e-15`／`3.35727e-15` |

單 rank 沒有跨 rank method／layout mismatch，因此該 group 每種方法少兩項。
以上故障數以程序群執行計算，不將每個 rank 算成獨立案例。KSP options 錯誤
與不收斂案例由整個群組設定；其餘局部故障只施加於指定 rank。unit 預期捕捉
故障，因此成功的 unit 作業退出 0；原生 CLI 故障另外驗證所有 rank 與 launcher
退出非零，且不發布成功標記。

CLI 逐欄位比較誤差為 0。checkpoint 續跑的非時間欄位誤差為 0；最大 clock
relative difference `6.78149e-17`，仍由既有 clock gate 驗收。

CLI batch 最大單 rank wall `0.414230 s`／RSS `40,902,656 bytes`；checkpoint
batch 分別 `0.414231 s`／`40,890,368 bytes`。這些包含啟動／輸出，且多個驗證
批次同時執行，不是獨立效能量測，不能據此宣稱加速。新增 state copy、gather
暫存與協調成本仍需由 HPC-00C／06A／09D 量測。

## 證據與剩餘工作

後續 1D adapter 本地操作與準備階段的補強，見
[adapter 進度](HPC_01C_ONE_D_ADAPTER_PROGRESS.md)；以下保留本報告版本當時的剩餘範圍。

ignored `outputs/hpc01/one-d-implicit/` 保存 `unit-summary.json`、八個 unit／graph
logs、`run_units.py`、`cli/summary.json`、`checkpoint/summary.json`、各 rank
輸入／log hashes、build logs、`evidence-summary.json` 及 `inventory.json`。
初次 test 使用錯誤 checkpoint pack 參數的編譯失敗 log 保留；權威建置是
`final-source-build` 與 `coupling-final-source-build`，無 compiler warning／error。

仍未承諾修復 PETSc backend 內部已卡住的 collective、MPI process loss、
OOM-killer，或在部分 rank 只建立一半 PETSc 物件時原地恢復。這些與已建立
物件之間可捕捉的本地錯誤不同。呼叫者須在同一 communicator 進入相同入口，
並提供一致的物理配置；本次只新增 method／dimension 的防護，不驗證所有
外部資產的內容身分。HPC-01C 的 adapter port／BeginStep／bookkeeping、
其他 CPU／coupling CLI 及 runtime 邊界仍需補齊，因此保持未勾選。

沒有新的 GPU、跨節點或 scaling 驗收。全域 checkpoint、分散式 FSI 與其他
階段仍依 [完整待辦清單](../WORKSTATION_HPC_TODO.md) 推進。
