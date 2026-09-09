# HPC-01C：建構失敗與 PETSc 資源清理

狀態：已併入正式來源，正式路徑通過 229 個建構故障／重試案例與耦合 CLI
回歸。HPC-01C 的其他邊界仍未完成，保持未勾選。
日期：2026-09-08。

## 已確認的缺口

貼體 flow／transport 的建構子取得 Mat、Vec、KSP、IS、VecScatter 後，
若後續步驟拋出例外，C++ 不會呼叫該物件的解構子。既有建構流程因此缺少
這些 handle 的清理，也有未檢查的 PETSc 回傳碼。flow 的 boundary label
目錄與兩個 runtime 的 halo 配置還有本地配置失敗後，其他 rank 繼續進入
collective 的風險。

## 候選修正

候選來源保存在 `outputs/hpc01/construction/candidate/`；
`outputs/hpc01/construction/base.json` 記錄複製當時 249 份原始檔的 SHA-256。
正式 FSI 效能矩陣、無 OpenMP FSI 回歸與貼體同 rank 單 thread 記憶體基準
完成後，已逐檔核對 249 份 base 與正式來源，再整合七個來源／建置檔案。
整合後 251 份來源與驗收候選一致；舊來源另存回復用 archive。

- flow／transport 建構子加入失敗清理，解構子共用清理程序；禁止複製
  擁有 PETSc handle 的 runtime，避免重複釋放。
- OwnedRowAssembler 的矩陣／向量工廠在錯誤時清理尚未交付的 handle；
  檢查共同回傳狀態與 nonzero-pattern 選項一致性。
- boundary label 目錄分成準備、gather layout、發布三個本地階段，
  在進入 MPI 前協調錯誤並檢查整數溢位。
- halo 的本地準備先共同決策，再依序建立並檢查 IS／Vec／scatter。
- transport 初始化陣列在本地例外後歸還借用，再進行共同錯誤決策。
- KSP 配置前比對 communicator 內的 PETSc options，並檢查 setup 回傳碼。
- 大案例 fieldsplit 預設選項先記錄這次新增的 key；建構失敗時只清除
  這些 key。各 rank 即使在不同插入位置失敗，也各自復原到原有選項。
  成功建構保留既有預設值行為，不改使用者原先指定的 override。

`RuntimeConstruction.hpp` 提供僅在測試編譯旗標下啟用的階段故障注入與
物件觀察器。測試額外保留一個 PETSc reference，建構失敗後依相反建立
順序檢查 reference count 必須只剩測試持有的一個，再釋放該 reference。
此方法要驗證清理責任，不能只以「程式未崩潰」宣告沒有遺漏。

## 候選原生驗收

新測試 `test_runtime_construction_failure.cpp` 已在 world 三 ranks 與
獨立一／二-rank groups 執行：

1. 對 flow／transport 的輸入、目錄、矩陣、不同序號的向量、halo、KSP
   及最後發布階段，在群組最後一個 rank 注入失敗。
2. 使用真正無效的 `-ksp_type`，以及同群不同的 solver options，確認共同拒絕。
3. 每次清理後重新建構；正常 64-basis cube 執行 flow／transport 求解，
   與 COMM_SELF 參考比較速度、壓力和物種場。CPU relative L2 沿用 `1e-6`，
   零參考 absolute L2 沿用 `1e-12`；常數 source 另驗證 `2 + dt` 的解析值。
4. 另一個單一常數 basis 的 cube 專門產生空 owned rows、空 required-element
   與空 halo ranks，重複故障／清理／建構及初始場檢查。這不是流場數值驗收案例。
5. 使用既有 C0 square-duct fixture 建立 27 elements／1,000 nodes，實際進入
   fieldsplit 預設分支。分別在第一／後續預設值插入、預設值準備完成、runtime
   ready 注入失敗，另用真正無效的 KSP 選項。每次比對完整 PETSc options
   快照、保留使用者 `-fieldsplit_0_pc_type jacobi`，再驗證重新建構及預設值。
   此案例檢查建構、所有權與初始場，不宣稱 fieldsplit 求解效能已驗收。

另驗證同群不同的 matrix pattern policy 在進入條件式 collective 前共同拒絕，
之後可正常求解。

| communicator ranks | 正常 fixture：flow／transport | 空列 fixture：flow／transport | fieldsplit 回復 | pattern policy | 合計 |
|---:|---:|---:|---:|---:|---:|
| 3 | 19／17 | 19／17 | 5 | 1 | 78 |
| 1 | 18／16 | 18／16 | 5 | 0 | 73 |
| 2 | 19／17 | 19／17 | 5 | 1 | 78 |

合計 **229 個故障／重試案例**通過，原生 MPI 程序退出 0，stderr 為空。
候選以 GCC 11.4／OpenMPI 4.1.2／PETSc 3.15.5 real64／Int32 建置，
啟用編譯警告與斷言，沒有 compiler warning。測試使用單 thread 與 MUMPS；
實際 fieldsplit 預設分支的建構另由 1,000-node fixture 覆蓋。

第一次執行 `native-candidate-1/` 正常退出失敗：空列 fixture 只含 label 0，
卻沿用 label 0–3 的傳輸邊界，健康重試被輸入驗證拒絕。已修正 fixture，
沒有放寬驗收門檻。失敗版來源、binary 與日誌已保留。

修正後證據在 `outputs/hpc01/construction/native-candidate-2/`：
`evidence.json`、`native-acceptance.json`、build／launcher logs、
`candidate-source.tar.gz` 及當次 binary。另行重新核對全部候選來源、binary
和日誌雜湊，並保留執行程式 `run_candidate.py`。

既有回歸也已重建並通過：trial 的 93 個故障案例；獨立 flow graph 群組
最大 relative L2 `7.45082e-15`、species graph `3.35727e-15`；貼體
subcommunicator 的 velocity／pressure／transport／mass 最大誤差分別為
`3.48315e-15`、`9.94533e-16`、`1.70503e-14`、`1.62206e-14`。
證據在 `existing-regressions/` 與 `body-regressions/`。兩／四-rank 的
1／2-thread MPI/OpenMP 回復測試也已通過；兩-rank 場差異為 0，
四-rank 最大 velocity relative L2 `2.98772e-15`、pressure `9.57441e-16`。
同 rank 單 thread 記憶體基準結束後，已核對來源並併入正式 headers。

四-rank 的第一次執行因啟動設定使各 rank 主執行緒都綁在 CPU 0，
已在核實 `/proc` affinity 後主動停止，證據保存在 `parallel-4/placement-failure.json`。
這不是數值失敗或觀察 timeout。重跑 `parallel-4-bound/` 使用
`--map-by slot:PE=2 --bind-to core` 與 `OMP_PLACES=cores`，
四份 rank profile、八個互不重複的配置核心、實際執行緒 masks 及原生回復
gates 已核對通過。這些回歸與無 OpenMP FSI 有重疊，不作效能統計。

## 正式路徑確認

`outputs/hpc01/construction/native-root/` 保留正式來源、binary 與 229 案例
通過證據。CPU `petsc`／`iga_navier_stokes_openmp` 和 coupling `petsc`
目標已重建，確認再次執行 make 為最新狀態；建置日誌沒有 compiler warning。
`root-cli-smoke/` 的正式 bifurcation CLI 通過單 rank explicit、單 rank Aitken、
雙 rank explicit、單／雙 rank port 數值比較，以及失敗 trial 不發布完成標記。
外層 timeout 180 秒，退出 0；`integration.json` 記錄來源及測試／binary 雜湊。
這是正確性回歸，不重新替代 HPC-02 已封存的效能測量。

## 尚未涵蓋的邊界

- 建構呼叫端的參數複製與物件儲存空間配置、其他 runtime／adapter 邊界。
- 成功建構仍沿用全域 PETSc options；不同 domain 的 solver options prefix／
  database 隔離不由此次失敗回復修正宣告完成。
- 外部資產內容身分及其餘 HPC-01C 工作。
- MPI／PETSc collective 內部無法返回的失敗、程序遺失或 OOM killer，
  不能由本地 exception 協調證明可恢復。注入點位於已返回的本地階段或
  全群已成功建立物件之後，不能據此宣稱涵蓋所有底層配置失敗。

呼叫端配置與輸入準備的後續補強見
[建構輸入進度](HPC_01C_CONSTRUCTION_INPUT_PROGRESS.md)。本紀錄保留該補強前
229 案例的版本化證據，不宣告分散式 FSI 已完成。
