# HPC-04A legacy transport options

本報告保留各開發階段當時的狀態；目前 HPC-04A 已完成，整合來源與驗收範圍見
[完成稽核](HPC_04A_COMPLETION_AUDIT.md)。HPC-04B／C 仍進行中。

- 狀態：legacy CLI 限定範圍通過；HPC-04A 整體仍進行中。
- 日期：2026-09-09（美東）；基準 `d00e2fa`，改動限於 `iga_transport`、build
  dependencies、驗證 harness 與文件。PETSc 3.15.5／Open MPI 4.1.2／GCC 11.4，
  工作站 1／2 ranks，OpenMP 未啟用，OMP／OpenBLAS threads=1。

`iga_transport` 現使用 `domain_neuron_transport_transport_`，與 legacy 轉換所得
`neuron_transport` system 身分一致。保留未加 prefix 的共同選項，新 prefix 優先。
私有 snapshot 涵蓋 SetFromOptions、SetUp、Solve，生命週期延續到 KSP 銷毀。
每步新增實際 KSP／PC／factor backend、iterations／reason 記錄；輸出錯誤仍集體返回。
既有位置參數、GMRES restart=50、rtol=1e-8、block-Jacobi 及 warm start 不變；
沒有改變 `.ntiga` 或最終場格式。此舊 CLI 使用 `PETSC_OPTIONS` 配置，沒有擴充位置參數解析。

## 驗收

`scripts/hpc_legacy_solver_options.py` 使用既有 legacy 相容的 1-element、64-node
fixture（`outputs/hpc01/tools/accepted/fixture`），兩步、1／2-rank packed database。
本機 `outputs/hpc04/legacy-options/final-v2/acceptance.json` 保存來源 binary／fixture
hashes、argv、選項、實際配置與場比較。

20 個作業通過：14 正向、6 預期負向，21 份成功 rank report。原預設及原未加 prefix
GMRES／LU／MUMPS 的改動前後場逐 byte 相同（4 次比較）。Scoped FGMRES 加 viewer
的 4 次比較全部通過相對 L2 ≤ `1e-6`（zero-reference absolute ≤ `1e-12`），
最大相對 L2 `1.0777783045643756e-08`。正向 solve 均有正收斂原因。
無效 KSP、無效 factor backend、viewer 無法開檔各在 1／2 ranks 退出 1，
無最終場、無成功摘要、無 MPI_ABORT；隨後正常作業可再次成功。

初版 `final-v1` 中 FGMRES／block-Jacobi 沿用 `rtol=1e-8`，雙 rank 場相對 L2
`9.406299986579029e-06` 未通過，保留失敗紀錄。Final 候選配置使用 `rtol=1e-12`，
未更改 production defaults 或場驗收門檻。GMRES 與 FGMRES 的停止 norm 不同，
不能從相同 tolerance 值推論場精度相同。

原 `legacy_transport_failure_test` 也重新建置執行：world 2 ranks 的 16 次故障／重試，
以及各 singleton 子群的 15 次故障／重試通過，ownership=PASS。證據在 `failure-ranks`。
另以重新建置的 `legacy_stdout_failure_test` 注入新 `solver_configuration` 診斷的
short write、exception、bad_alloc、nonstandard exception、flush failure；5 個雙 rank
作業均退出 1、無場輸出、返回 `legacy transport solver diagnostics: rank 0:`。
證據為 `stdout-acceptance.json` 及 `stdout-*` logs。

主要 20-job suite 的最大單 rank host RSS 42,033,152 bytes。雙 rank 原預設的
assembly `0.0036217 s`、solve `0.000422856 s`、兩步合計 14 iterations；這是並行
功能測試期間的小案例觀測，不能作 speedup／scaling 證據。未涉及 CUDA allocation。
本次不增加物理模型驗證或傳輸 budget 宣稱，僅檢查既有離散場與錯誤／生命週期契約。

```bash
make -C solvers/cpu iga_transport legacy_transport_failure_test legacy_stdout_failure_test \
  PETSC_DIR=/usr/lib/petscdir/petsc3.15/x86_64-linux-gnu-real
python3 scripts/hpc_legacy_solver_options.py \
  --fixture outputs/hpc01/tools/accepted/fixture \
  --reference-binary outputs/hpc04/legacy-options/iga_transport-before \
  --output-dir NEW_OUTPUT
```

下一步：完整 moving 回歸、多層子求解器診斷，以及 HPC-04B／C 的 PDE 預條件器與
網格／硬體擴展評估。HPC-04A 仍維持未勾選。
