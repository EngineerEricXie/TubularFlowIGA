# HPC-04C 固定網格 rank 比較

狀態：128 元素、1／2／4 ranks、兩次重複的本機量測與證據核對完成；
HPC-04C 已於 2026-09-11 由 [完成報告](HPC_04BC_COMPLETION_REPORT.md)簽核。

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

正式本機數據與 source 核對見下一節，rank 工作量與 assembly 瓶頸仍須進一步定位。
本流程只有 fixed-mesh strong scaling 的局部工作站配置；其他網格、weak scaling、
跨節點與不同硬體仍依 HPC-04C／HPC-09 原要求驗收。


## 128 元素正式重複量測

`outputs/hpc04/rank-scaling/fixed128-v1/acceptance.json` 的六個 rank／repeat
配置、每個 LU 與 block-Jacobi 共 12 次求解全數通過。固定 4×4×8 元素、
539 control nodes、2156 flow DOFs；仍是原 C2 duct graph 與兩個 accepted steps。
同 rank 候選對 LU 場比較、跨 rank／repeat checkpoint 場比較及原 Phase 9
物理 validator 均保持原門檻。跨 rank 的速度／壓力最大相對 L2 為
`2.7998249980599e-11`。

後驗 audit 重新核對 672 個 hashes、28 份 rank reports／logs、六次候選對 LU
場比較與 12 次跨 rank checkpoint 比較，並重新執行 12 次物理 validator，全數
通過。Audit 與重現腳本分別為同目錄的 `audit.json`、`audit_scaling.py`。
Frozen native binary SHA256 為
`9771d5cebd123b95ab4329d8d5d3a4efc2b2e64d34ae9435966de49027fad3a1`；
它代表 `01aa4c8` 的歷史求解器，不代表後續 FSI 元件的 current HEAD binary。

下表時間均為秒；wall 是每次最大 rank wall 的兩次 median，括號為 min–max；
各 phase 也是每次最大 rank exclusive seconds 的 median，可能來自不同 rank。

| 候選 | ranks | wall median（範圍） | assembly | setup | linear solve | communication | iterations／次 | speedup |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| lu | 1 | 278.164（277.131–279.196） | 258.426 | 15.614 | 0.160 | 0.005 | 42 | 1.000× |
| lu | 2 | 267.539（266.678–268.399） | 251.714 | 12.055 | 0.121 | 2.045 | 42 | 1.040× |
| lu | 4 | 268.249（268.179–268.319） | 252.858 | 11.252 | 0.114 | 6.706 | 42 | 1.037× |
| bjacobi | 1 | 264.674（262.907–266.441） | 256.466 | 3.180 | 1.088 | 0.006 | 588 | 1.000× |
| bjacobi | 2 | 255.073（253.811–256.335） | 248.888 | 0.816 | 2.129 | 1.763 | 3060 | 1.038× |
| bjacobi | 4 | 261.596（260.633–262.559） | 255.467 | 0.121 | 1.331 | 6.782 | 3724 | 1.012× |

RSS 以下取兩次執行、所有 ranks 的最大單 rank peak，單位 bytes；不是同時的
aggregate peak，也不代表 root 以外的每個 rank 都有同樣用量。

| 候選 | 1 rank | 2 ranks | 4 ranks |
|---|---:|---:|---:|
| LU | 132,435,968 | 120,258,560 | 107,106,304 |
| block-Jacobi | 86,769,664 | 64,651,264 | 55,422,976 |

Assembly 仍約 249–258 秒，增加 ranks 未出現接近線性的加速。4 ranks 的
observed efficiency 為 LU 25.9%、block-Jacobi 25.3%；2 ranks 的 wall 甚至略
低於 4 ranks。Block-Jacobi 在此網格 wall 較 LU 低，但 iterations 由 588 增至
3060／3724。這是需要保留的迭代數與 scaling 負收益，不足以變更預設策略。
數字指出 assembly 值得進一步定位，尚不能單靠 aggregate phase 時間判定
是重複工作、負載不均、同步等待或其他特定程式原因。

## 執行條件與解讀限制

每個 rank 的 affinity 核對為 `[0,1]`、`[2,3]`、`[4,5]`、`[6,7]`（依
rank 數取前幾組），OMP／OpenBLAS=1。測量期間另有綁定 CPU 14、15 的
FSI 編譯與正確性測試，已逐批記錄在 `background*.json`；它們仍共享記憶體、
頻寬與系統資源。因此這批資料不是獨占 allocation 的效能認證，兩次重複也
不能提供穩健的信賴區間。沒有把本機 MPI 結果延伸為跨節點或其他硬體結論。

完成報告另彙整 16／128／1024-element mesh growth 與目前 HEAD 小案例；
HPC-09 的 scheduler／跨節點正式 scaling 驗收仍未完成。沒有更改任何預設求解策略。
