# HPC-09A–E 本機部署進度

- 狀態：本機歷史驗收已完成；後續 HPC-09A–D 真實跨節點驗收亦通過，見 [Bridges-2 完成報告](HPC_09_BRIDGES2_CROSS_NODE_ACCEPTANCE.md)。
- 基準 revision：`66fd930` 加本報告所列 scaling case 產生器變更。
- 日期／主機：2026-09-11，`TsungYehLab`。
- 機器可讀摘要：`outputs/hpc09/local-deployment-v3/summary.json`。

## 已完成實作

- `check_dependencies.sh COMPONENT --report FILE` 與 `hpc_build_manifest.py`
  記錄 source state、C++／MPI wrapper、PETSc scalar/index 與 MUMPS／HYPRE／HDF5
  能力、HDF5 link configuration、CUDA compiler／目標架構、Slurm allocation 和
  可選 binary linkage。本機 CUDA probe 能找到未加入 `PATH` 的
  `tubularflow-cuda` CUDA 12.6 compiler。
- `hpc_test_tiers.py` 分離 unit、1／2／4-rank MPI、GPU hardware smoke 與
  scheduled two-node tier。每項命令都有外層 timeout、獨立 stdout／stderr、非零
  exit 傳遞和 `result.json`；缺少必要資源時只產生具理由的 `skipped`。
- `multinode_graph.sbatch` 對齊 rank／thread／BLAS 綁定，將相同輸入 staging 到
  每個節點的 NVMe，結果與 checkpoint 保留在 shared filesystem。Slurm 的 batch
  shell 預警會明確轉送至 Open MPI ranks，完成下一個 accepted step checkpoint 後
  才重提交流程。
  它可在 attempt zero 自動建立 16,384-element、8-step graph，並以可控延遲將
  SIGUSR1 送至與 Slurm `B:` 預警相同的 shell trap；第一個作業 checkpoint 後
  requeue，下一個 attempt 自動 discovery／restart。
- `hpc_cross_node_scaling.py` 與 `cross_node_scaling.sbatch` 提供 exclusive Slurm
  strong／weak scaling。Collector 要求至少三次 repetition、strong 一 rank 數值
  參考、每次 native physical validation，以及 weak elements/rank 容許範圍；輸出
  wall、speedup、efficiency、iterations、phase／communication 與 RSS。
  `hpc_prepare_scaling_cases.py` 可直接產生 16,384-element fixed strong case 與
  每 rank 256 elements 的 weak cases，免除提交前手工建立未稽核案例。
- `cross_node_fsi.sbatch` 將 1-rank reference、跨兩節點 4-rank nonzero strong
  FSI writer 與新 2-rank read-only pair restore 串成單一驗收；既有 checkers 驗證
  完整 field／surface／history／ports／守恆、4→2 重分區及 checkpoint bytes 不變。
  同一作業先以每節點一 rank 執行 scheduled tier，避免只驗大型 fixture 而缺少
  小型跨節點測試層。
- 根 README、Bridges-2 guide、文件索引及新的
  [HPC deployment guide](../HPC_DEPLOYMENT.md) 已更新能力、限制與可重現命令。
  Grouped checkpoint／immersed、native moving／FSI graph restart 與跨節點證據均明列。
- `hpc_finalize_cross_node.py` 統一核對三個作業的 clean、相同 revision、兩節點
  hostname、graph signal/requeue/restart、scheduled tier、FSI 4→2 restart、完整
  repetitions、物理解與數值比較，以及 wall／iteration／RSS／communication 指標。
  只有它產生 `status: passed` 才能關閉尚待 allocation 的七項工作。

## 本機驗收

```bash
PETSC_DIR=/usr/lib/petscdir/petsc3.15/x86_64-linux-gnu-real \
  ./scripts/check_dependencies.sh cpu --report /tmp/cpu-build.json
./scripts/check_dependencies.sh cuda --report /tmp/cuda-build.json
python3 -m unittest discover -s scripts/tests
python3 scripts/hpc_test_tiers.py --tier unit --timeout 300 --output-dir /tmp/unit
PETSC_DIR=/usr/lib/petscdir/petsc3.15/x86_64-linux-gnu-real \
  python3 scripts/hpc_test_tiers.py --tier mpi --timeout 300 --output-dir /tmp/mpi
LD_LIBRARY_PATH=/home/tsungyeh/anaconda3/envs/tubularflow-cuda/targets/x86_64-linux/lib \
  python3 scripts/hpc_test_tiers.py --tier gpu --timeout 120 --output-dir /tmp/gpu
python3 scripts/hpc_test_tiers.py --tier scheduled --output-dir /tmp/scheduled
```

結果：CPU／CUDA dependency reports 通過；57 個 Python tests 通過；unit tier
四組命令、MPI build 加 1／2／4-rank 數值案例、CUDA execution contract 加真實
RTX 4080 SUPER `device-info` 都通過。非 Slurm 工作站的 scheduled tier 正確寫出
`skipped`，理由為缺少至少兩節點 allocation，沒有誤報為 pass。Slurm scripts 經
`bash -n`，Python scripts 經 `py_compile`，全工作區 `git diff --check` 通過。

既有 shared multidomain smoke suite 亦在相同本機程式狀態完整通過；species 最大
edge／domain／global residual 分別為 `6.7763e-20`、`9.5563e-15`、`8.8902e-15`。

Scaling case generator 與 collector 另以 128-element fixed case、每 rank 32
elements 的工具自測執行 strong／weak 各三次。六次實際 PETSc process、geometry、
physical validation、strong fields、profile／RSS 與統計全部通過；strong 三次各
30 次、weak 三次各 25 次線性迭代。證據在
`outputs/hpc09/scaling-collector-self-test-v3/results/summary.json`。此案例只驗證
collector 全路徑，不作效能或 scaling 宣稱；正式 scheduler 預設使用 16,384-element
strong 與每 rank 256 elements 的 weak cases。

## 原本待實際硬體驗收的項目（現已完成）

以下記錄原工作站驗收當時缺少 Slurm／跨節點 allocation 的項目；後續已在 Bridges-2 完成，數值與 job IDs 見上述完成報告：

1. 同 revision 的 Bridges-2 build manifest、two-node scheduled tier 與 rank binding。
2. 大型 native graph 的跨節點 checkpoint、`SIGUSR1` 轉送、完整 generation 及
   requeue 後續跑；這同時是 HPC-05C 的剩餘驗收。
3. distributed moving flow／FSI 的跨節點場、守恆、強耦合、surface owner 與
   restart；這同時是 HPC-07C／D 的剩餘驗收。
4. 非微型 strong／weak cases 的三次以上獨占量測及完整效能／記憶體報告。

上述四組跨節點證據經 finalizer passed 後，HPC-09A–D 與相關 HPC-05C／07C／07D 已勾選；全部 38 項完成。原工作站的 GPU 與本機結果仍屬各自歷史 revision，沒有改寫為本批實測。
