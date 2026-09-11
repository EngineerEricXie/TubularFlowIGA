# HPC-09A–E 本機部署進度

- 狀態：HPC-09E 完成；HPC-09A–D 的本機工具與回歸通過，等待實際跨節點 allocation 驗收。
- 基準 revision：`ed4333f` 加本報告所列工作區變更。
- 日期／主機：2026-09-11，`TsungYehLab`。
- 機器可讀摘要：`outputs/hpc09/local-deployment-v2/summary.json`。

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
- `hpc_cross_node_scaling.py` 與 `cross_node_scaling.sbatch` 提供 exclusive Slurm
  strong／weak scaling。Collector 要求至少三次 repetition、strong 一 rank 數值
  參考、每次 native physical validation，以及 weak elements/rank 容許範圍；輸出
  wall、speedup、efficiency、iterations、phase／communication 與 RSS。
- `cross_node_fsi.sbatch` 將 1-rank reference、跨兩節點 4-rank nonzero strong
  FSI writer 與新 2-rank read-only pair restore 串成單一驗收；既有 checkers 驗證
  完整 field／surface／history／ports／守恆、4→2 重分區及 checkpoint bytes 不變。
- 根 README、Bridges-2 guide、文件索引及新的
  [HPC deployment guide](../HPC_DEPLOYMENT.md) 已更新能力、限制與可重現命令。
  Grouped checkpoint／immersed、native moving／FSI graph restart 與跨節點證據均明列。

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

結果：CPU／CUDA dependency reports 通過；56 個 Python tests 通過；unit tier
四組命令、MPI build 加 1／2／4-rank 數值案例、CUDA execution contract 加真實
RTX 4080 SUPER `device-info` 都通過。非 Slurm 工作站的 scheduled tier 正確寫出
`skipped`，理由為缺少至少兩節點 allocation，沒有誤報為 pass。Slurm scripts 經
`bash -n`，Python scripts 經 `py_compile`，全工作區 `git diff --check` 通過。

既有 shared multidomain smoke suite 亦在相同本機程式狀態完整通過；species 最大
edge／domain／global residual 分別為 `6.7763e-20`、`9.5563e-15`、`8.8902e-15`。

## 尚待實際硬體驗收

目前主機沒有 `sbatch`／`srun`、`SLURM_JOB_ID` 或跨節點 allocation，因此沒有執行
或宣稱下列結果：

1. 同 revision 的 Bridges-2 build manifest、two-node scheduled tier 與 rank binding。
2. 大型 native graph 的跨節點 checkpoint、`SIGUSR1` 轉送、完整 generation 及
   requeue 後續跑；這同時是 HPC-05C 的剩餘驗收。
3. distributed moving flow／FSI 的跨節點場、守恆、強耦合、surface owner 與
   restart；這同時是 HPC-07C／D 的剩餘驗收。
4. 非微型 strong／weak cases 的三次以上獨占量測及完整效能／記憶體報告。

因此 HPC-09A–D 保持未勾選。HPC-09E 已完成，因為 README 能力矩陣、限制、命令
及所有硬體缺口均已收尾且沒有將整份清單誤標完成。
