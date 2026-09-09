# HPC-01D：CUDA 啟動資源與容量檢查

日期：2026-09-09。基準 HEAD `8a7e549` 加本批工作樹。
HPC-01D 部分完成，整份清單仍維持原範圍。

## 實作

[ExecutionEnvironment.hpp](../../include/ExecutionEnvironment.hpp) 共用既有 CPU
thread-setting parser，保留 namespace、數值及錯誤語義，新增不依賴 MPI 的
單程序啟動檢查。[CudaExecution.hpp](../../solvers/cuda/include/CudaExecution.hpp)
在 CUDA 五個 CLI 入口讀檔及 GPU 查詢前驗證資源設定並檢查 stdout flush。

- Open MPI 或 PMI 的 rank／size 必須完整、合法且 size 為 1；拒絕不完整、
  負值、溢位及 rank 越界。MPI metadata 優先於外層 Slurm allocation。
- Slurm step 使用 `SLURM_STEP_NUM_TASKS`，缺少時使用 `SLURM_NTASKS`；
  batch／extern shell 及只有 allocation task count 的環境不當成多程序啟動。
  依據 [Open MPI 環境變數](https://docs.open-mpi.org/en/v5.0.8/tuning-apps/environment-var.html)、
  [srun](https://slurm.schedmd.com/srun.html) 與
  [sbatch](https://slurm.schedmd.com/sbatch.html) 的區分。
- OMP thread 清單／limit 與 OpenBLAS、MKL、BLIS requests 使用 CPU 相同 parser。
  摘要記錄 requests、單程序／單 GPU、real64／int index 位寬及 host OpenMP
  編譯旗標；未設定為 -1，BLAS 0 保留 library default，沒有改環境或綁定政策。
- 裝置摘要新增 CUDA build／linked runtime、driver API 版本、HDF5 build／linked
  版本，以及 serial VTKHDF、all-element database loading 能力。
- Mesh counts 必須為正且不超過 INT_MAX；solver node count 額外受
  `INT_MAX / max(fields, 3)` 限制，涵蓋 cuBLAS vector length、field 與 XYZ velocity
  的 int 索引。Configured fields 仍為 1–8。檢查先於大型標籤／場配置。
- Block IDs 維持 int，block-value offsets 維持 size_t；配置前檢查 bytes 乘法溢位，
  不把原本 size_t 可表示的 scalar-value count 錯限縮為 INT_MAX。

沒有更改弱式、solver tolerance、資料格式或 GPU 讀取 CPU 多分區資料的行為。
此檢查不涵蓋未知 launcher 或沒有受支援 metadata 的 PMIx-only launcher，
也不是防止兩個獨立作業寫入相同目錄的鎖定機制。

## 驗收

| 範圍 | 結果 |
|---|---|
| 無 CUDA／PETSc 的 host checks | 84 項通過：launcher、thread、stdout、索引邊界；無大型配置 |
| CPU 共用 parser 回歸 | 150 項通過：OpenMP 1／4 threads、無 OpenMP，world3 與 split1+2 |
| 原生 CUDA 啟動 | 23 作業通過，含實際 Open MPI 多程序拒絕及單程序接受 |
| CUDA 數值／stdout 回歸 | 38 作業、27 場比較通過；包含 17 次故障與新作業重試，無 timeout |
| CPU 兩分區資料庫由單 GPU 讀取 | 1,005 nodes、720 elements；bad samples 0，min detJ `9.0103244285259191e-5` |

[原生控制器](../../scripts/hpc_cuda_resources_regression.py) 先用封存 binary
重現 `mpiexec -np 2 device-info` 兩程序皆返回 0，再確認修改後五個入口兩程序
皆返回 1。Rank wrapper 等待各程序紀錄完成，核對逐程序退出碼及 diagnostic。
thread／PMI／Slurm 拒絕案例同時隱藏 GPU，證明設定檢查早於 GPU 查詢。
`/dev/full` 測試資源摘要失敗；四種超限 node header 在小型 database 副本中
修改，證明容量拒絕不需巨大記憶體。所有拒絕皆在 project device buffer
peak/live 為 0 時發生，未產生指定輸出。合法原生 Slurm／PMI 環境由測試合成，
沒有宣稱在真實 Slurm allocation 或 PMI launcher 執行。

數值回歸沿用 [CUDA stdout 控制器](../../scripts/hpc_cuda_stdout_regression.py)。
CPU 場使用已接受 HPC-00 矩陣，重新核對原輸入及場 hash，沒有重跑 CPU solver。
保留 GPU 前後 relative `1e-6`、CPU/GPU `1e-5`、零參考 absolute `1e-12` 門檻，
實测最大 relative L2 分別為 `1.09172e-12` 與 `4.06867e-6`。
速度／壓力分開比較，transport 保留 node IDs。獨立 `iga_flow_validate` 檢查
before／healthy／retry 三份流場，最大 relative mass imbalance `3.43214e-8`，
符合 `1e-6`；命令與結果見 `mass-validation.json`。所有 solver 作業都有正值
project allocation peak，釋放後 live bytes 為 0。

以下為本批單次健康作業的 exclusive 時間，僅供驗收記錄，不作加速結論。
Host RSS 包含執行期間峰值；project CUDA peak 不包含 driver/library 自行配置。

| 路徑 | 組裝 s | 線性求解 s | Host peak bytes | CUDA project peak bytes |
|---|---:|---:|---:|---:|
| Flow | 6.335321 | 36.686502 | 305872896 | 52207740 |
| Configured transport | 0.085441 | 1.162254 | 296132608 | 24744428 |
| Legacy transport | 0.000792 | 0.759427 | 283435008 | 976084 |

## 重現與證據

環境：本機 RTX 4080 SUPER、SM89、CUDA 12.6；實測 runtime build/linked
`12060`、driver API `13010`、HDF5 build/linked `1.12.1`。Driver API 數字不是
NVIDIA driver package version。CPU parser 使用 GCC 11.4、Open MPI 4.1.2、
PETSc 3.15.5 real64/int32。建置 logs 無新 compiler warning。

```bash
make -C solvers/cuda execution-test
conda run -n tubularflow-cuda make -C solvers/cuda CUDA_ARCHS=89 \
  iga_cuda cuda_stdout_failure_test
conda run --no-capture-output -n tubularflow-cuda env \
  LD_LIBRARY_PATH=/home/tsungyeh/anaconda3/envs/tubularflow-cuda/targets/x86_64-linux/lib \
  python3 scripts/hpc_cuda_resources_regression.py \
  --baseline outputs/hpc01/cuda-resources/before/iga_cuda \
  --fixture-root outputs --output-dir /path/to/new-startup-results
conda run --no-capture-output -n tubularflow-cuda env \
  LD_LIBRARY_PATH=/home/tsungyeh/anaconda3/envs/tubularflow-cuda/targets/x86_64-linux/lib \
  python3 scripts/hpc_cuda_stdout_regression.py \
  --baseline-dir outputs/hpc01/cuda-resources/before \
  --fixture-root outputs --output-dir /path/to/new-numerical-results
```

需先保留 HPC-00 flow／transport 矩陣及 HPC-01 legacy fixture。
`outputs/hpc01/cuda-resources/` 保存 host／CPU tests、build-final log、
`native/summary.json`、逐程序 logs、binary／input hashes 及數值比較。
封存 binary 以實際 SHA256 識別，不聲稱由此 HEAD 重新編譯。
首次 `numerical/` 程序中斷、缺少完整報告，保留作中斷紀錄；重新執行於
`numerical-final/`，不把中斷作業納入通過案例。

HPC-01D 尚需其餘 embedding 入口、PETSc 64-bit／complex／backend 與輸出能力矩陣、
真實 scheduler 資源／綁定核對；本批不是效能、跨節點或多 GPU 驗收。
