# HPC-01C／01D：Standalone FSI exporter 執行入口

日期：2026-09-08（EDT）。多 rank 拒絕與單 rank 完整數值驗收均通過。

`phase8_compliant_channel_fsi_paraview` 使用單程序 FSI runtimes，原本由每個
MPI rank 各跑一份模型並寫入同一路徑。現在在解析輸出路徑、建立目錄、清除
既有 diagnostic 或建構 runtime 前確認 communicator size 為 1。
多 rank 啟動會明確返回錯誤，沒有把重複執行當成分散式 FSI。

`outputs/hpc01/exporter-entry/rejection.json` 已驗證：

- 2 ranks、新輸出路徑：全部返回 1，未建立輸出目錄。
- 3 ranks、既有輸出目錄：全部返回 1，既有 diagnostic 與 sentinel 均未改變。
- 1 rank、缺少參數：仍返回原有 usage 錯誤，沒有誤判為 rank 限制。

共 6 份 rank reports，均未 timeout。重新建置無 compiler warning。
上述拒絕路徑沒有進入物理求解，不能取代單 rank 完整案例的驗收。

完整案例初次設定 300 秒 timeout，已返回 124；既有
[FSI 單執行緒量測](HPC_02_FSI_PERFORMANCE.md)約需 703 秒，因此這個上限不足。
原始 `healthy.log` 與輸出目錄保留；確認作業終止後，已以 1800 秒上限及新
目錄 `healthy-final/` 重跑，日誌為 `healthy-final.log`。該次已返回 0，
4 次耦合迭代收斂、94312 個流體採樣點、9 個膜節點與 8 個三角形通過既有
數值／守恆／epoch gates。Moving mass 為 `0.015561258177167596`，wall leakage
為 `0.015561258176288955`，continuity 為 `1.8222488628470318e-14`，均沿用原有
門檻。產生的 PVD、VTM 及兩份 VTU 均已解析，且沒有 export failure marker。

該次 wall time 為 668.793 秒，peak RSS 為 176432 KiB；執行期間有其他開發／
建置活動，這不是隔離效能量測。`healthy-files.json` 保存五個產物的 hashes，
`executions.json` 保留初次 timeout 與最終成功結果。驗收 binary 與當時來源
對應 `../transport-visualization-init/after-binaries/` 及該批 `source-final.tar.gz`；
後續共用文字讀取器變更另行驗收，不把此舊 binary 的結果當成重新測試全部新來源。

有效命令如下，重跑須使用新的輸出路徑：

```bash
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 \
  timeout --kill-after=5s 1800s mpiexec --map-by core --bind-to core -np 1 \
  solvers/cpu/phase8_compliant_channel_fsi_paraview /absolute/path/to/new-output
```

分散式 FSI 仍依 HPC-07 的介面、所有權、全域收斂與跨節點驗收要求開發。
2026-09-11 的 [HPC-01D 最終簽核](HPC_01D_COMPLETION_REPORT.md)另在此入口建立
輸出前加入共同 PETSc ABI／thread resource preflight；兩種非法 thread 設定的
原生拒絕案例通過。既有單 rank 完整物理結果未被改寫。
