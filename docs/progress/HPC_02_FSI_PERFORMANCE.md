# HPC-02：FSI 單機 OpenMP 隔離量測

日期：2026-09-08。1／4-thread 隔離矩陣及預定效能／記憶體門檻已通過；
修正後無 OpenMP binary 的完整場回歸也已通過。貼體同 rank 單 thread
記憶體比較亦已完成，HPC-02A／B／C 的必要驗收已通過。

## 結果

同一個 FGMRES 版本、完整 compliant-channel FSI fixture、相同輸出範圍，
每種配置執行首次一次及三次正式重複。首次不納入統計，各次以新程序執行，
1／4 threads 交錯排列；沒有其他編譯或數值測試與矩陣並行。

| threads | 端到端中位數 s | assembly 中位數 s | solver setup 中位數 s | linear solve 中位數 s | 正式重複的最大 peak RSS bytes |
|---:|---:|---:|---:|---:|---:|
| 1 | 703.224653 | 692.263471 | 0.759159 | 0.005827 | 71,868,416 |
| 4 | 305.748723 | 294.870696 | 0.754370 | 0.005804 | 73,424,896 |

4-thread 端到端速度比為 **2.300008**，RSS ratio 為 **1.021657**。
兩者通過預先指定的時間至少改善 10%、RSS ratio 不超過 1.25 門檻。
RSS 採每次程序峰值，再取三次正式重複的最大值；不是把不同時間峰值相加。
phase 表列互斥時間的各自中位數，不應將各 phase 的中位數相加當成某一次總時間。

八次原生 FSI 收斂、force／moment、moving mass、wall leakage 與 continuity
gates 均通過。與 HPC-00D 封存參考的 **72 個完整場比較**全部通過，
最大 relative L2 為壓力的 `6.961917358512401e-12`，沿用 CPU `1e-6`
與零參考 absolute `1e-12` 門檻。沒有改變 fixture、壓力 gauge 或物理參數。

## 環境與核對

WSL 工作站、GCC 11.4、OpenMPI 4.1.2、PETSc 3.15.5 real64／Int32。
CPU 集合為 1-thread `{0}`、4-thread `{0,2,4,6}`，經系統 topology 核對
為不同回報核心；OpenMP `PROC_BIND=close`、`PLACES=threads`，BLAS threads 為 1。
實際 worker team 與有界 batch 由原生日誌確認。

```bash
python3 scripts/hpc_openmp_matrix.py \
  --reference outputs/hpc00/reference-states/fsi-reference \
  --output-dir outputs/hpc02/volume/isolated-matrix-fgmres \
  --binary solvers/cpu/compliant_channel_fsi_openmp_test \
  --threads 1 4 --cpus 0 2 4 6 --repetitions 3 --timeout 1800
```

矩陣 controller 已退出 0。另以
`outputs/hpc02/volume/verify_fgmres_matrix.py` 重新讀取每份 collection、
原生與 launcher 日誌、GNU time 資源記錄、實際 argv／CPU affinity／環境、
場檔及來源雜湊，重算摘要後與矩陣結果一致。

完整證據位於 `outputs/hpc02/volume/isolated-matrix-fgmres/`：
`matrix.json`、`accepted-evidence.json` 與 `source-binary-and-evidence.tar.gz`。
封存包含來源、binary、參考場、各次輸出及核對程式；accepted evidence 記錄
封存及各檔雜湊。`partial-verification.json` 保留量測進行時的部分核對，
最終結果以八次執行的 accepted evidence 為準。

這是單一工作站、單 MPI rank 之 FSI 體積組裝加速證據，不代表已完成
分散式 FSI。固定總核心數的貼體純 MPI／混合比較見
[HPC_02_HYBRID_PROGRESS.md](HPC_02_HYBRID_PROGRESS.md)。

## 無 OpenMP 回歸

重新建置的 `compliant_channel_fsi_test` 已退出 0，原生 gates 與九個完整場
比較全部通過，最大 relative L2 同為 `6.961917358512401e-12`。建置命令
沒有 `-fopenmp`，binary 的 undefined symbols 也沒有 OpenMP 呼叫。
來源、binary、build log 與場檔已核對並封存在
`outputs/hpc02/volume/fsi-fgmres-no-openmp/`。該次與其他正確性測試重疊，
不納入上述隔離效能結果。
