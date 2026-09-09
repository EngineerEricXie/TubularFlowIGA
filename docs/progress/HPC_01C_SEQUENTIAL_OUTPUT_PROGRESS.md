# HPC-01C：Sequential CSV 與 manifest 關檔

日期：2026-09-08（本機 EDT）。本批接續 [registry 驗收](HPC_01C_REGISTRY_PROGRESS.md)，
HPC-01C 與整份 goal 仍未完成。

## 問題與修正

`iga_1d_3d_explicit` 的 explicit CSV、strong history／iterations CSV 與兩種
coupling manifest 原先只在串流仍開啟時檢查狀態。若最後關檔失敗，析構不會
回報錯誤，程式仍可能寫出 graph binding manifest 並印出 completed。

現在五個 writer 都明確 close，再檢查串流；CSV 必須成功關檔，才能開始寫
coupling manifest。Graph binding manifest 保持原來的 tmp／close／rename
協議，仍最後發布。Root 輸出改由 `CollectiveLocalStage` 協調，所有 ranks
取得同一 stage／root 診斷；也不再於錯誤協調前配置 `output_error` 字串。
此 callback 僅做本地檔案 I/O，不含 MPI 或 PETSc 呼叫。

修改不改變 CSV／JSON 內容、物理計算、收斂門檻或 CLI 參數。
關檔失敗會留下可診斷的部分輸出，沒有增加整個輸出目錄的交易回復，也沒有
承諾 fsync／斷電持久性；完整 checkpoint 發布另由 HPC-05 追蹤。

## 驗收

`scripts/hpc_sequential_output_regression.py` 使用原有 smoke fixture，覆蓋
explicit、fixed、Aitken，以及 1／2 ranks。Linux 測試 library
`text_close_preload.so` 只攔截指定完整路徑的 `fclose`，真正關閉檔案後回傳
一次 `EOF`／`EIO`。每個故障案例必須觀察到注入訊息，避免攔截失效而誤判。
該 library 不連入 production binary，也不新增 production 故障環境選項。

| 驗收 | 結果 |
|---|---|
| 原生矩陣 | 51 個案例、82 份 rank reports，全部依預期返回且沒有 timeout |
| 新版故障 | 22 次關檔錯誤共同退出 1，沒有 completed 或 graph binding 完成標記 |
| 舊版重現 | 2 ranks 的 8 次 CSV／coupling manifest 關檔錯誤被忽略，仍返回 0 並發布 graph binding；原有 graph binding close 的 3 次錯誤則正確拒絕 |
| 健康執行／重試 | 6 次舊版、6 次新版及故障序列後的 6 次新版健康執行，全部成功 |
| 格式與數值相容 | 10 份 CSV 加 12 份 JSON，總計 22 份逐位元組相同；健康重試也與相應新版結果相同 |
| 原有 sequential smoke | 修改後全套通過，包含 positional／graph、1／2 ranks、subcycling、fixed／Aitken 等價與既有拒絕案例 |

每份 rank report 記錄 exit status、timeout、輸出與資源使用量；多 rank
controller 等待全部 reports 落盤後才退出，避免 launcher 提前殺掉其他
rank 的驗收記錄。圖形完成標記的故障目標是實際 `.json.tmp` 路徑。
CSV 關檔錯誤的結果亦確認未產生後續 coupling manifest。

本機環境沿用 WSL／GCC 11.4／Open MPI 4.1.2／PETSc 3.15.5 real64 Int32；
OMP／BLAS 各 1，使用原有 smoke 的 `-ksp_type preonly -pc_type lu`。
Compilation 啟用 `-Wall -Wextra -Wpedantic`，沒有新 compiler warning。
沒有執行新的 GPU、多節點或效能驗收。

```bash
make -C solvers/coupling iga_1d_3d_explicit text_close_preload.so \
  PETSC_DIR=/usr/lib/petscdir/petsc3.15/x86_64-linux-gnu-real
python3 scripts/hpc_sequential_output_regression.py \
  --case-dir /path/to/retained-explicit-smoke-fixture \
  --reference-binary /path/to/before/iga_1d_3d_explicit \
  --output-dir /path/to/new-results
```

Fixture 來自以 `TUBULARFLOWIGA_KEEP_TEST_OUTPUT=1` 執行的原有
`explicit_coupling_smoke_test`，需要其 `one.ntiga`、`two.ntiga`、三個 domain
資料夾及 schema-v5 configuration。Controller 複製輸入再調整 scheme／rank
database，不更改來源 fixture。

Ignored 證據位於 `outputs/hpc01/sequential-output/`；`accepted/summary.json`
保存完整案例與比較結果，`before/` 保存原始 writer 與 binary。
`smoke.log` 保存原有 smoke 的成功結果；`source-final.json`／tar 與
`after-binaries/` 保存最後受測版本，`acceptance.json` 與 `audit.py` 可核對證據。

## 接續入口

下一批應補齊 sequential 強耦合及其共用初始化，而非只檢查收斂分支：

- 前置 graph／port catalog、waveform／boundary、1D open-loop 初始化與
  history 配置，目前仍有本地例外直接進入外層 catch 的路徑。
- Strong-fixed／Aitken 的 Begin／input／port observation／local result、
  全群收斂、prepare 與 postcommit bookkeeping，必須在 ranks 間協調。
- Abort 應嘗試清理所有 runtime 後再組合診斷；單一 1D abort 失敗不能讓
  其他 ranks 提前進入 3D collective。Runtime 的內部協調仍由各 runtime 負責。
- 完成上述項目後，繼續支援入口覆蓋稽核；不能以這次輸出回歸宣告 HPC-01C 完成。
