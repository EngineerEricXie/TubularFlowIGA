# HPC-02：有界元素批次執行器準備

後續更新（2026-09-08）：已接入實際浸入式暫態體積積分，並通過 27-cell
元素與失敗回復測試；完整驗收見 [整合進度](HPC_02_VOLUME_PROGRESS.md)。
以下保留 helper 準備階段的原始範圍與證據。

本報告記錄早期執行器與機制測試。後續已完成求解器整合與 HPC-02A／B／C
驗收，見 [體積積分整合進度](HPC_02_VOLUME_PROGRESS.md) 與
[FSI 效能驗收](HPC_02_FSI_PERFORMANCE.md)；以下保留當時的機制證據。
這是接續已量測組裝瓶頸的準備，沒有以 helper 測試宣稱流體或 FSI 已加速。

[ParallelElementBatch.hpp](../../solvers/cpu/include/ParallelElementBatch.hpp)
將一批元素的工作分為三段：呼叫執行緒準備輸入、OpenMP workers 計算、
呼叫執行緒按原順序消費結果。`capacity` 限制同時保留的元素輸入／結果
數量；呼叫端仍須計算每個元素及 worker scratch 的 byte 預算。

worker 接收 const 輸入，每個元素保留獨立結果與 exception slot。
所有 workers 結束後先檢查錯誤；任一元素失敗時不消費該批次，以最小索引
的例外向呼叫端回報。先前批次若已插入，仍由 runtime 的 rollback／重新
組裝處理。動態排程不改變消費順序。

函式本身沒有 PETSc 或 MPI 依賴。整合時 compute callback 也必須沒有
PETSc／MPI 呼叫或共享可變 scratch；泛用模板不能替呼叫端保證這件事。
原生無 OpenMP 建置只接受一個 thread，明確拒絕多 thread 要求。

```bash
make -C solvers/cpu parallel-element-batch-test
```

同一測試來源另以直接 GCC 命令編譯無 OpenMP 及 `-fopenmp` 版本，
使用 `-O3 -std=c++17 -Wall -Wextra -Wpedantic`。兩者通過，無新增警告。
無 OpenMP 版本 406 項檢查；OpenMP 版本對 1／2／4 threads 共 1,215 項
檢查，驗證實際 team size、主執行緒 prepare／consume、私有 scratch、
批次上限、固定消費順序，以及 worker 例外不逃出 OpenMP 區域。
這些是機制檢查數量，不是物理案例數量。
直接編譯／執行日誌保存在 `outputs/hpc02/batch/build-test.log`。

下一步優先整合 `ImmersedTransientFlowRuntime::Assemble` 的體積積分。
其 body force 已凍結，適合 worker 只讀；現有 `Gather` 仍呼叫 PETSc，
須留在 prepare，矩陣／向量 Scatter 及壁面／port 工作留在 consume。
保留原元素順序及 volume system 供壁面 delta 使用，避免重新積分或重算。
已有完整 FSI 參考，可比較速度、壓力、膜與材料位移、速度、牽引力及力。

整合後仍須執行真實 1／2／4 thread 全場及原生 gate、錯誤回復、
RSS 與分階段時間比較；之後再完成純 MPI／混合模式的固定核心數比較。
