# HPC-06D：大型前處理完成報告

狀態：完成。日期：2026-09-11。驗收機器：單一 workstation，16 個可見 CPUs；
量測使用 GNU time wrapper，每個 stage 為 fresh process。機器可讀證據位於
`outputs/hpc06/preprocess-scaling-v2/summary.json`。

## 案例與方法

`hpc_preprocess_scaling.py` 由同一 schema-v4 straight-tube 設定產生三個幾何長度，
依序執行 C++ mesh pipeline、cache-only spline（1／4 OpenMP threads）、4-way
METIS partition、sparse-cache packer 與 database inspector。每個 stage 記錄 wall、
user/system CPU、peak RSS、stdout/stderr hash；輸入、工具及主要產物亦記錄 SHA256。
所有案例通過 mesh scaled-Jacobian、bad-element、surface-intersection、spline
element count、packer 及 inspector gates，且 peak RSS 均低於明定的 2 GiB 上限。

| 案例 | 控制點 | 元素 | control mesh | cache | `.ntiga` |
|---|---:|---:|---:|---:|---:|
| small | 3,417 | 2,880 | 263,130 B | 29,435,068 B | 29,588,512 B |
| medium | 13,065 | 11,520 | 1,067,218 B | 119,632,828 B | 120,200,992 B |
| large | 38,793 | 34,560 | 3,340,626 B | 360,160,188 B | 361,834,272 B |

產物隨元素數近似線性成長。Large cache 與 database 都超過 343 MiB，但 spline
peak RSS 為 168.0 MiB、packer peak RSS 為 19.7 MiB，證明既有 256-element chunks、
ordered streaming 與 packer sequential record read 沒有把完整輸出保留在記憶體。

## Stage 成本

| large stage | wall | peak RSS |
|---|---:|---:|
| mesh pipeline | 0.315 s | 20.4 MiB |
| spline，1 thread | 6.677 s | 167.1 MiB |
| spline，4 threads | 2.718 s | 168.0 MiB |
| METIS partition | 0.364 s | 140.5 MiB |
| sparse-cache packer | 1.916 s | 19.7 MiB |
| database inspector | 0.665 s | 4.2 MiB |

Spline 是量測中最長的 stage，既有 OpenMP 在 large case 提供 2.456× wall speedup，
且 1／4-thread 的 `bzmeshinfo.txt`、`spline_cache.igacache`、`bzmesh.vtk` 與
`geometry_transform.json` SHA256 完全相同。保留目前 chunked OpenMP 實作。
Mesh 與 packer 分別只占 0.315 與 1.916 秒，本機證據不支持增加其平行化複雜度。
METIS 的 RSS 高於 mesh／packer，但 wall 很短且為外部工具，因此本批不更換 partitioner。

Small case 額外執行 legacy text spline 與 `iga_pack --legacy-text`。`cmat.txt` 與
`bzpt.txt` 合計 47,340,745 bytes，cache 為 29,435,068 bytes；legacy pack wall
0.665 秒，cache pack wall 0.164 秒。兩條路徑輸出的 `.ntiga` SHA256 都是
`5cd8264e5bd19ec75373825fcae7b001af3de6d2b9aaf91e8a9069a5a7630ab6`，所以保留
legacy 介面不會犧牲 native cache 路徑的預設效率。

## 重現

```bash
make mesh spline cpu EIGEN_DIR=/usr/include/eigen3
python3 scripts/hpc_preprocess_scaling.py \
  --output-dir NEW_OUTPUT --threads 1 4 --memory-limit-gib 2
```

這是單工作站、直管拓撲的工作量與記憶體驗收，不代表複雜分岔的幾何常數或共享
檔案系統效能。既有 NMO 31,680-element 報告提供分岔案例的相近規模交叉證據；
跨節點與 scheduler filesystem 行為仍由 HPC-09 驗收。
