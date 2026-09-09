# HPC-03B：浸入式積分工作量分區

日期：2026-09-09。基準 `64497fd` 加本批修改。

## 分區與量測

[WeightedWorkPartition.hpp](../../include/WeightedWorkPartition.hpp) 提供不依賴 MPI／
PETSc 的連續分區。輸入是按穩定 cell ID 排序的正整數工作量，以二分搜尋最小
可行的最大 group 工作量，再貪婪建立不超過 rank 數的連續 groups。保持 cell
空間順序；允許空 ranks；拒絕零權重、無效 rank 數及工作量加法／乘法溢位。
此最適性只針對給定權重與連續 groups，不表示實際 wall time 或通訊量最適。

正式 static operator／runtime 的 constructor 新增選用參數：

```cpp
iga::ImmersedStaticDistributedRuntime runtime(comm, domain, volume, surface,
    ghost, options, iga::ImmersedWorkPartition::WeightedContiguous);
```

預設 `CellCount` 保留既有分配。兩種模式共用相同 active nodes／row offsets、
物理控制、gauge／port scalar IDs、數值 kernels 與資料格式。分區方式納入
collective configuration agreement；單 rank 無效 enum 或同群有效但不同模式
都在建立 distributed PETSc 資源前共同拒絕。

初始模型估計每次組裝的 matrix-entry 訪問量：

- 每個 stencil 的 dense scratch 清零成本為列數平方。
- Volume 以 catalog 的 logical quadrature point 數乘四場 dense block 大小。
  Compact cut rule 的實際積分點數會反映細分工作；不把已複製建立的 geometry
  catalog／切割分類成本當作 owner 的數值積分工作。
- Conservative trace 與選定 wall 以 surface point 數及 dense block 大小計費。
  Wall builder 即使不發布 wall delta 仍會計算，因此估計沿用此實際行為。
- Ghost face 使用 16 點、same-field block 大小；成本附在 minus cell owner。
- Gauge 計入 scratch；ports 計入其 surface 點與 block，measurement-only stencil
  使用量測向量長度。廉價 controller target 保持在末 rank，不納入分區權重。

Volume、gauge、ports、ghost 的成本先依 cell 聚合；沒有體積但仍需要 port 量測的
cell 亦可納入。Weighted 模式重新指定 stencil owner，row owner 維持原有分配。
這是可解釋的初始算術量模型，尚未以硬體時序擬合各 kernel 的成本係數。

Backend 記錄每次成功組裝的 halo、局部積分／PETSc 插入、Mat／Vec stash exchange
時間。局部時間在共同錯誤協調前結束，能呈現各 rank 工作不平衡；halo 包含
scatter、local copy 與共同檢查，stash 包含 PETSc assembly 及共同檢查。
Operator 總 assembly time 另含其他協調、全域診斷等成本，不能只加上述三項
當成總時間。失敗時不發布新 timing。另輸出 estimated work、required／remote
halo rows，配合實際時間評估工作量與通訊的取捨。

## 驗收方法

[獨立分區測試](../../solvers/cpu/tests/test_weighted_work_partition.cpp) 將長度
0–7、每項權重 1／3／11、1–5 ranks 的 16,400 組案例與窮舉 contiguous cuts 的
最小最大工作量比較，另檢查 uint64 上界與拒絕行為。

正式物理測試使用既有 27 active cells／54 ghost faces、compact depth 2、兩個
flow controllers／gauge、非零四場及 body force。每 rank 的序列 reference 僅
用於測試；比對 residual／Jacobian action、全域 port／wall／pressure diagnostics，
要求 relative L2 `1e-6` 及 exact declared preallocation、零動態配置。

加權 Newton 另跑四 ranks 的 flow 與 empty-work，保留 HPC-03A 的逐場誤差、
全域守恆、controller、candidate／prepare 失敗及 commit／rollback gates。
Weighted empty-work 的兩個 cells 由 rank 0／1 積分，rank 2／3 沒有 stencil／halo，
但仍持有各自 rows 並參與共同求解。

[比較控制器](../../scripts/hpc_immersed_partition_benchmark.py) 在相同核心綁定與
OMP／BLAS 1 下順序執行：單 rank 跑 A／B；2／4 ranks 跑 A／B／B／A，其中
A 是 CellCount、B 是 WeightedContiguous。每個作業組裝四次，排除第一次暖機，
保留其餘三次的各 rank 原始時間；彙整每次最大 rank 時間的中位數。兩次不同
模式總 estimated work 必須相同，weighted 最大 work 不得超過 count partition。
Controller 檢查每 rank CPU affinity 非空且互不重疊；數值與退出狀態仍由原有
物理回歸 gate 驗證。時間改善不是正確性 gate，不會隱藏加權模式較慢的結果。

```bash
make -C solvers/cpu weighted-work-partition-test
make -C solvers/cpu immersed_distributed_physics_test \
  immersed_distributed_static_flow_test immersed_distributed_assembly_test \
  PETSC_DIR=/usr/lib/petscdir/petsc3.15/x86_64-linux-gnu-real
python3 scripts/hpc_immersed_partition_benchmark.py \
  --output-dir /path/to/new-partition-comparison
python3 scripts/hpc_immersed_static_regression.py --partition weighted \
  --mode flow --ranks 4 --output-dir /path/to/new-weighted-newton
python3 scripts/hpc_immersed_static_regression.py --partition weighted \
  --mode empty-work --ranks 4 --output-dir /path/to/new-weighted-empty
python3 scripts/hpc_immersed_assembly_regression.py --kind physics \
  --physics-mode faults --partition weighted --ranks 2 \
  --output-dir /path/to/new-weighted-faults
```

工作站 GCC 11.4、Open MPI 4.1.2、PETSc 3.15.5 real64/int32／MUMPS。
Evidence 留於 ignored `outputs/hpc03/work-partition/`；輸出目錄須使用新路徑。

## 單機結果

Benchmark 的 10 個作業、26 份 rank reports 全部退出 0、無 timeout，
所有物理與稀疏配置 gates 通過。以下為暖機後每次最大 rank 時間的中位數：

| Ranks | 分區 | Assembly s | 局部積分／插入 s | Halo s | Stash s | 局部工作 max／mean |
|---:|---|---:|---:|---:|---:|---:|
| 1 | cell-count | 8.073803 | 8.073474 | 7.368e-06 | 1.8729e-05 | 1 |
| 1 | weighted | 8.044446 | 8.044092 | 7.326e-06 | 1.7739e-05 | 1 |
| 2 | cell-count | 4.097081 | 4.070312 | 5.834e-06 | 0.02616772 | 1.004242 |
| 2 | weighted | 4.085158 | 4.058063 | 5.4205e-06 | 0.02674936 | 1.004233 |
| 4 | cell-count | 2.21904 | 2.193324 | 7.6325e-06 | 0.02591951 | 1.077075 |
| 4 | weighted | 2.224692 | 2.195763 | 8.2935e-06 | 0.02829316 | 1.079207 |

| Ranks | 分區 | Estimated work max／mean | Σ required rows | Σ remote rows | 最大 rank peak RSS bytes |
|---:|---|---:|---:|---:|---:|
| 1 | cell-count | 1 | 867 | 0 | 66064384 |
| 1 | weighted | 1 | 867 | 0 | 65859584 |
| 2 | cell-count | 1.000071 | 1529 | 662 | 86204416 |
| 2 | weighted | 1.000071 | 1529 | 662 | 86085632 |
| 4 | cell-count | 1.075318 | 2661 | 1818 | 80023552 |
| 4 | weighted | 1.075318 | 2649 | 1806 | 86429696 |

最大 residual／Jacobian-action relative L2 分別為 `2.844863e-16`, `7.279741e-16`。

這個幾何的原始 cell-count 分配已達到連續分區的最小最大估計工作量；
weighted 在 2 ranks 得到相同分配，在 4 ranks 改變中間的 cell 分界，但最大
預估負載不變。4-rank assembly ratio（weighted／count）約 `1.00255`，沒有
明確加速；總 remote halo rows 從 1818 降到 1806，但 stash 中位時間略增。
不以約千分之幾的時間差宣稱效能改善，亦不更換既有預設。

這些是固定小案例的單機量測。RSS 包含複製的 geometry catalogs 與每 rank
序列 reference；remote rows 是一次 required-state scatter 的索引數，不是
PETSc 內部實際傳輸位元組。Benchmark 僅組裝、不執行 KSP；求解時間不適用，
不是零秒求解。Newton correctness 作業另保存 assembly／linear timings。
沒有 CUDA allocation 或跨節點結果；大型、非均勻幾何與其他硬體仍須量測，
初始模型不能保證任意案例加速。

Weighted flow／empty-work 的四 rank Newton 回歸通過，完整場最大 relative L2
分別 `3.33575e-11`／`3.38179e-15`；逐場、全域守恆、controller 與交易回復 gates
均通過。Backend unit 的 1／2／4／split 1+2 通過。最終 physical executable 的
weighted faults 2-rank 驗收涵蓋 physics／partition 無效與不一致，以及積分後
失敗與同物件重試。所有 affected targets 以 warnings 啟用建置，無新 compiler
warnings；Python syntax 與 whitespace checks 通過。

`acceptance.json` 彙整正式 benchmark、Newton、unit、final-faults 的 hashes；
較早 `faults/` 是加入 remote-row 輸出前的 executable 結果，不混入最終 binary
驗收。HPC-03B 的可選工作量分區與比較工具完成；後續以實際案例調整模型，
case／graph／暫態入口依 HPC-03C 接入。整份清單尚有 27 項未勾選。
