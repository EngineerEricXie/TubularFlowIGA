# HPC-08A／B／C 完成與驗收報告

- 狀態：工作站實作與驗收完成；跨節點驗收由 HPC-09 追蹤。
- 基準 revision：`d4a148b` 加本報告所列工作區變更。
- 日期／主機：2026-09-11，`TsungYehLab`，本機 Open MPI／PETSc real64。
- 完整機器可讀證據：
  `outputs/hpc08/domain-group-scaling-v3/summary.json`。

## 實作範圍

`simulation_config.json` 新增可選的 `resources`。省略時仍使用原有 shared
communicator；`domain_groups` 以 manifest 順序配置互斥的連續 ranks，要求每個
domain 恰好被一個 group 擁有，並允許未配置的 ranks 只參與全域協調。

原生 multidomain runner 只在 owner group 建立 0D／1D／body-fitted 3D PETSc
runtime。全域 registry 改持有 collective proxy；port state、species accounting
與 group-local failure 由 owner root 傳回全域。子域內所有 PETSc、MPI reduction
及 solver context 都收到 group communicator，未新增 `PETSC_COMM_WORLD` 的子域
呼叫。

Hydraulic executor 依 flow dependency level 建立 batch；同層且位於不同 group
的 domains 同時求解，同 group 則維持確定順序。Species executor 在每一步解析
donor 後，以實際 transport direction 重建 batch，因此 direction reversal 與
staged hydraulic／transport trial、rollback、commit 語義均保留。完整使用方式及
重建規則見 [MULTIDOMAIN_RESOURCES.md](../MULTIDOMAIN_RESOURCES.md)。

## 正確性驗收

執行：

```bash
make -C solvers/coupling multidomain_config_test pressure_flow_executor_test \
  species_pressure_flow_executor_test
make -C solvers/coupling domain_group_runtime_proxy_test \
  multidomain_domain_groups_test \
  PETSC_DIR=/usr/lib/petscdir/petsc3.15/x86_64-linux-gnu-real
mpiexec -np 3 solvers/coupling/domain_group_runtime_proxy_test
mpiexec -np 5 solvers/coupling/multidomain_domain_groups_test /tmp/hpc08-groups
```

結果：

- parser／planner 拒絕重複、遺漏、超額及未知 group 配置；shared schema 相容測試通過。
- 3-rank proxy 測試通過 2-rank group、single-rank owner、port/species 傳輸、
  group-local failure 的全域拒絕與健康重試。
- schema-v5 的兩個 3D groups 對 serial 參考最大 normalized difference
  `9.76074e-8`。
- schema-v6 staged transport 最大 normalized difference `3.45439e-7`。
- 0D owner 的八步狀態與 conservation accounting 最大 normalized difference
  `6.81272e-9`。
- executor 單元測試驗證 branch batch、相依違規拒絕、transport donor 反轉及
  staged callback 路徑；所有結果低於既有 `1e-6` graph 門檻。

## 資源配置實測

執行：

```bash
python3 scripts/hpc_domain_group_scaling.py \
  --output-dir outputs/hpc08/domain-group-scaling-v3 --repeats 3
```

案例使用同樣五個 ranks 與兩個各 128 elements／539 nodes 的 3D ducts。shared
模式讓每個 3D runtime 使用全部五 ranks 並依序求解；grouped 模式以一 rank
執行所有小型 0D／1D domains，兩個 3D domains 各用兩 ranks，且在同一 graph
dependency level 同時求解。每種配置執行三個新程序 repetitions。

| 指標（三次中位數） | shared | grouped | grouped/shared |
|---|---:|---:|---:|
| max-rank process wall | 25.7044 s | 13.0003 s | 0.5058 |
| max-rank application elapsed | 25.3742 s | 12.6856 s | 0.5000 |
| individual peak RSS 加總 | 320,487,424 B | 294,842,368 B | 0.9200 |
| max-rank peak RSS | 65,380,352 B | 64,913,408 B | 0.9929 |
| max-rank explicit communication／等待 | 1.3437 s | 12.6788 s | 9.435 |

三次 grouped/shared 數值比較皆為 `5.48528163302372e-7`。程序分組使兩個 3D
solve 重疊，wall time 約減半，並減少約 8.0% 的 individual peak RSS 加總；
小型 owner rank 幾乎整段等待，communication phase 明確呈現這項成本。因此
shared 保持預設，grouped 只在 domain concurrency 或記憶體收益足以抵銷協調時啟用。

另以相同大小的前後相依雙 3D 案例量測到 grouped/shared wall ratio `1.0062`、
RSS ratio `0.9195`；此負面結果保留在
`outputs/hpc08/domain-group-scaling-v2/summary.json`，證明 scheduler 沒有把具相依
關係的 domains 錯誤重疊。

## 限制與持久化策略

Completion manifest 記錄 `resource_mode` 與
`contiguous_manifest_group_order`；原始 graph manifest 是重建 group mapping 的
權威資料。Grouped native checkpoint/restart 在 runtime 建立前明確拒絕，需續跑的
作業目前使用 shared 模式。Grouped immersed 3D 尚未支援。跨節點 rank placement、
網路成本與 scheduler policy 尚未在工作站冒充驗收，統一由 HPC-09 完成。
