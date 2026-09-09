# HPC-01C：CUDA stdout 邊界

日期：2026-09-08。基準工作樹 `7e31f86`；HPC-01C 仍未完成。
證據位於 `outputs/hpc01/cuda-stdout/`，含建置 log、封存 binary、輸入／
binary SHA-256、逐作業 argv、timeout、RSS、allocation 及場比較。

## 修改與測試

CUDA 的裝置摘要、mesh-check 提早返回、solver 最後 allocation／profile
摘要均明確 flush 並檢查串流狀態。裝置摘要失敗在配置 GPU buffer 前返回 1；
其餘 stdout 失敗最晚在對應 CLI 成功返回前被檢出。
保留單程序／單 GPU 模式、既有 exception mask 與數值／檔案格式。
中途 logging 失敗可能在完成求解後才回報；不宣稱每行都立即 fail-fast。

原生 wrapper 在指定訊息片段替換 stdout buffer；CUDA 自己的 main 捕捉
受檢查串流產生的例外並返回 1，wrapper 另外證明注入實際命中。
`device-info` 測五種 short／flush／runtime_error／bad_alloc／非標準例外；
mesh-check、flow、configured／legacy transport、allocation 與 profile
各測 short／flush。每個入口都有修改前、健康與新作業重試。

`native/summary.json` 的 38 個作業全部通過，含 17 次故障及 21 次健康作業；
沒有 timeout。27 個場比較全部通過：修改前後 GPU 最大 relative L2
`4.81687e-13`（門檻 `1e-6`），相對既有 CPU 基準最大 `4.06867e-6`
（門檻 `1e-5`）；零參考 absolute 門檻均為 `1e-12`。
CPU 使用 HPC-00 已接受場與原輸入，控制器重新核對 hash；沒有重跑 CPU binary。
傳輸比較保留節點 ID，流場速度／壓力分開比较，不調整 pressure gauge。

另對修改前、健康與重試三份 flow 執行獨立 `iga_flow_validate`，均返回 0，
最大 relative mass imbalance `3.43214e-8`，通過原 `1e-6` 門檻；
記錄在 `mass-validation.json`。所有 solver 作業（包括最後摘要故障）
都有正值 project-buffer peak 且 live bytes 為 0。

## 環境與量測

本機 RTX 4080 SUPER、compute capability 8.9、CUDA 12.6、SM89；
單程序、OMP／BLAS 1 thread，`IGA_PROFILE=1`。建置無 compiler warnings。
下列為單次健康作業的 exclusive 階段時間與實測記憶體，僅供正確性驗收記錄，
不作為獨立效能或加速基準；CUDA driver/library allocation 不在 project peak 內。

| 路徑 | 組裝 s | 線性求解 s | Host peak bytes | CUDA project peak bytes |
|---|---:|---:|---:|---:|
| Flow | 4.016868 | 11.950671 | 302583808 | 52207740 |
| Configured transport | 0.052607 | 0.392991 | 295821312 | 24744428 |
| Legacy transport | 0.000646 | 0.271425 | 280166400 | 976084 |

## 重現與限制

```bash
conda run -n tubularflow-cuda make -C solvers/cuda CUDA_ARCHS=89 \
  iga_cuda cuda_stdout_failure_test
conda run --no-capture-output -n tubularflow-cuda env \
  LD_LIBRARY_PATH=/home/tsungyeh/anaconda3/envs/tubularflow-cuda/targets/x86_64-linux/lib \
  python3 scripts/hpc_cuda_stdout_regression.py \
  --baseline-dir outputs/hpc01/cuda-stdout/before --fixture-root outputs \
  --output-dir /path/to/new-cuda-output
```

需要保留的 HPC-00 flow／transport 矩陣與 HPC-01 legacy fixture；
flow 使用既有較嚴格 nonlinear／mass controls，不放寬誤差門檻。
GPU runtime 需要本機 GPU 執行權限，並非在 cluster login node 運行。
健康重試是新程序／新目錄，不代表損壞的 writer 可恢復，也不提供輸出原子發布。
尚未驗證的 F05／F06、HPC-01D 與其他清單項目繼續保留。
