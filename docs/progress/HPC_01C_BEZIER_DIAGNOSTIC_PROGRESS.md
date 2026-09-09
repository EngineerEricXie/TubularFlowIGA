# HPC-01C：Bezier 診斷格式化

日期：2026-09-08。基準 `5ce448e`。F05 完成，HPC-01C 保持未勾選。

`BuildBezierVisualizationMesh` 在兩個 extraction signature 相同但座標不相容時，
原本即會拒絕幾何；惟診斷 formatter 的 allocation／num_put 例外可能被
ostringstream 吞掉，最後僅拋出截斷的幾何錯誤訊息。
現在該局部 formatter 啟用 bad/fail exceptions，保留原始格式化例外。
正常 diagnostic 文字、正確幾何、數值與資料格式不變。

延伸既有 `bezier_visualization_test`，生成真實共享節點但座標相差 0.01 的
兩元素資料庫。正常拒絕訊息必須含兩個 element ID、local point 與 distance；
三種 formatter 故障必須保留 runtime_error、bad_alloc 或非標準 int 例外，
每次恢復 locale 後再要求診斷文字完全一致。

先以舊實作執行新測試，確實在原始例外保存斷言失敗（返回 -6，關閉 core dump）；
加入一行 exception mask 後相同測試返回 0，三次故障／三次恢復均通過，
原有 conforming／coincident／overlap／point-array 測試亦一起通過。
本機 GCC 11.4、純 C++、無 compiler warnings：

```bash
make -C solvers/cpu bezier_visualization_test
timeout 30 solvers/cpu/bezier_visualization_test
```

本批證據為 `outputs/hpc01/cuda-stdout/bezier-{before,after}*`（與前一批 evidence
共用目錄，但獨立 logs／退出紀錄）。此改動僅在幾何已被拒絕的診斷路徑，
組裝／求解 timing、RSS、CPU/GPU 場誤差及 mesh 平滑驗收為 N/A；
沒有把前一批 CUDA binaries 宣稱為此 header 修改後的重建結果。
完整入口／constructor／destructor／早退邊界簽核仍由 F06 追蹤。
