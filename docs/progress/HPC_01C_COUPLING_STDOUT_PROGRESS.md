# HPC-01C：1D、耦合與序列工具的 stdout 邊界

日期：2026-09-08。基準工作樹 `ebd924e`。HPC-01C 仍未完成。
本批完成 required stdout 的失敗檢查，不改數值模型或輸出檔案格式。

## 實作與驗收範圍

原生 1D 的 `--check` 摘要改在共同本地階段輸出；1D 最終摘要、
multidomain／bifurcation 與 sequential 的完成訊息均 flush 並檢查串流狀態。
失敗時所有 rank 返回 1，再完成既有清理／PETSc finalize。
不改全域 ostream exception mask，亦不在含 PETSc collective 的區塊外套本地 stage。

四個序列工具 `iga_inspect`、`iga_case_check`、`iga_config_check`、
`iga_flow_validate` 在每條成功返回路徑檢查 stdout；寫入失敗返回 1。
新增的兩個控制器封存 binary／input SHA-256、argv、退出碼與比較結果。
修改前 binaries 在重建前封存，實際 binary 身分以 summary hashes 為準；
不把基準工作樹 revision 當作所有封存 binary 的完整建置來源證明。

## 結果

本機 GCC 11.4、Open MPI 4.1.2、PETSc 3.15.5 real64／Int32／MUMPS。
MPI 測試使用 1／2／3 ranks、OMP／OpenBLAS／MKL 各 1 thread，
`-ksp_type preonly -pc_type lu -pc_factor_mat_solver_type mumps`。

`outputs/hpc01/coupling-stdout/native-final/summary.json` 全部通過：

- 216 個原生 MPI 作業：18 個修改前參考、18 個健康作業、90 次故障及 90 次新作業重試。
- 348 份 rank 報告均核對退出碼、timeout 與資源紀錄；全部故障各 rank 返回 1，無 timeout。
- 五種故障為 short-write、flush failure、runtime_error、bad_alloc、非標準例外。
  每次均確認 rank 0 注入成功與對應共同階段的診斷。
- 660 個檢查摘要／輸出比較通過。1D 的 summary 排除實測 timing／RSS 後完全相同，
  同時驗證這四個量測有限、非負且 RSS 大於零；其餘輸出逐位元相同。
- 覆蓋 1D check／solve、schema-v5 flow／0D graph、schema-v6 species graph、
  獨立 bifurcation executable，以及 sequential explicit／fixed／Aitken 三種方法。

第一次 `native/` 控制器將實測 RSS 也當作固定摘要欄位比較，因此正確中止；
已保留 failed summary。修正比較範圍後以全新目錄重跑完整矩陣，沒有放寬物理量標準。

`outputs/hpc01/utility-stdout/native/summary.json`：48 個原生序列工具作業通過。
12 條成功路徑的健康摘要前後逐位元相同；12 個舊版 `/dev/full` 案例皆靜默返回 0，
12 個新版案例皆返回 1 並報告 `cannot write complete text stream`。
涵蓋 legacy／configured case、3D／1D／flow graph／species graph config，以及
flow validator 的 field、compare、manifest、compare-manifests、Womersley 模式。
Womersley 僅驗證 CLI 摘要相容，不將解析點樣本誤當成已驗證的 IGA 係數近似誤差。

`make -C solvers/one_d core-test` 三個既有測試通過；
`make -C solvers/coupling petsc-test` 既有耦合 smoke test 返回 0。
本批建置無新增 compiler warnings，`git diff --check` 通過。

## 重現

```bash
make -C solvers/one_d iga_1d one_d_stdout_failure_test
make -C solvers/coupling petsc graph_stdout_failure_test \
  bifurcation_stdout_failure_test sequential_stdout_failure_test
python3 scripts/hpc_coupling_stdout_regression.py \
  --baseline-dir outputs/hpc01/coupling-stdout/before \
  --graph-fixtures outputs/hpc01/registry/native-final \
  --sequential-fixtures outputs/hpc01/sequential-output/accepted \
  --output-dir /path/to/new-coupling-output
make -C solvers/cpu iga_inspect iga_case_check iga_config_check iga_flow_validate
python3 scripts/hpc_utility_stdout_regression.py \
  --baseline-dir outputs/hpc01/utility-stdout/before --fixture-root outputs \
  --output-dir /path/to/new-utility-output
```

原始案例由既有 registry、sequential-output、solver-stdout、serial-tools-audit
驗收工具生成；控制器需要保留的 fixture 目錄。各次作業 logs、逐 rank RSS／
wall time 與建置／整合 logs 在上述 evidence roots。本批屬失敗處理驗收，
沒有宣稱效能提升、跨節點驗收或輸出 bundle 原子發布。
最終摘要失敗時，數值檔案可能已經寫出；重試使用新作業與新輸出目錄。

剩餘：CUDA stdout 邊界、F05 診斷、F06 完整入口／清理核對與 HPC-01D；
所有尚未驗收的主清單項目繼續保持未勾選。
