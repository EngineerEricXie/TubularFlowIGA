# HPC-04C 固定網格 rank 比較

狀態：量測流程已驗證，正式 128 元素重複量測待完成。

`scripts/hpc_duct_rank_scaling.py` 序列執行同一 C2 duct 網格的 1／2／4 ranks，
每個配置預設兩次，候選為 LU／MUMPS 與 block-Jacobi／ILU。底層仍使用原
Phase 9 validator 與相同場門檻，每個 rank 配置的 LU reference／候選比較之外，
還將每次候選的兩步 checkpoint 速度／壓力與該候選第一次單 rank 比較。

每次執行保存全部 rank 的 RSS、CPU affinity、iterations、wall，以及 assembly、
setup、linear solve、communication、output 的最大 rank exclusive seconds。
摘要提供 median 與 min/max wall、觀察 speedup／efficiency；它不是 benchmark
信賴區間。各 phase 最大值可能來自不同 rank，不能相加作單 rank wall。

候選 harness 新增 `--launcher`，以 argv 執行，既有預設仍是 `mpiexec`。
Rank sweep 預設 `mpiexec --map-by core --bind-to core`，OMP／OpenBLAS=1。
這個 launcher 介面會追加 `-np`，不是通用 Slurm `srun` 介面。

```bash
python3 scripts/hpc_duct_rank_scaling.py \
  --binary outputs/hpc04/immersed/native-v4-binary \
  --output-dir NEW_OUTPUT --transverse 4 --axial 8 \
  --ranks 1 2 4 --repeats 2 --timeout 2400
```

Native frozen binary 為 `01aa4c8`；新 harness 記錄自身、依賴與 binary hashes，
逐次保存 machine-readable acceptance，失敗不列入完成的 scaling 摘要。
原大 duct 的並行作業已結束，啟動前以 process list 確認無殘留 MPI 測試；
但本機工作站並非獨占 allocation，其他使用者／系統活動仍可能影響數字。

## 流程驗證

`outputs/hpc04/rank-scaling/smoke-v1/acceptance.json`：2 元素、1／2 ranks、
一次 repetition、兩候選，共四次求解，原 validator、候選場比較及跨 rank
checkpoint 比較均通過。CPU affinity 確認 binding 生效。這只驗證 harness，
不能作為 strong scaling 證據。Python 語法與 diff whitespace 檢查亦通過。

正式數據仍須完成、核對 rank 工作量及 source，再解讀收益與負收益。
本流程只有 fixed-mesh strong scaling 的局部工作站配置；其他網格、weak scaling、
跨節點與不同硬體仍依 HPC-04C／HPC-09 原要求驗收。
