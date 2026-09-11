# HPC-04B／C：預條件器與擴展性完成報告

日期：2026-09-11。完成 audit 基準 `dc63fc8`；機器可讀結果位於
`outputs/hpc04/completion-bc-v1/acceptance.json`。

## 候選範圍

Body-fitted C2 duct 以 LU／MUMPS、block-Jacobi／ILU、A11 full Schur／LU、
velocity GAMG＋pressure Jacobi，以及三種 selfp Schur 變體完成評估。網格從
16、128 到 1024 elements；中、大案例使用具 SHA256 身分的 frozen binary，
目前 HEAD 另重跑 16-element、2-rank 四候選，確認現行求解路徑仍符合相同
物理、場值與 solver diagnostics 門檻。

Immersed 路徑涵蓋 closed pressure gauge、兩個 flow-controller ports、pressure
reference，以及小切割元素。Gauge／controller 額外 rows 沒有被誤當成均勻
四欄位 saddle-point block；AMG 只用於速度子區塊。

| 幾何／案例 | 評估結果 |
|---|---|
| 16-element duct，目前 HEAD | 四候選皆通過；最大場 relative L2 `1.08564e-10` |
| 128-element duct，4 ranks | 四候選皆通過；最大場 relative L2 `1.67480e-10` |
| 1024-element duct，4 ranks | LU／block-Jacobi 通過；最大場 relative L2 `3.68215e-11` |
| selfp Schur 小 duct | LU、selfp-LU、selfp-GAMG-LU 通過；selfp-GAMG 不收斂 |
| immersed regular cut | closed／flow／pressure 各三候選，共九項通過 |
| immersed small cut | LU 通過；block-ILU／block-LU-shift 均以 `DIVERGED_ITS` 2000 次停止 |

所有成功案例保留既有 pressure、flow、3D mass、全系統體積平衡與 checkpoint
field gates。三個不收斂組合完整保留為候選評估結果，沒有把它們改列成功。

## 網格與 rank 擴展

小／中／大 duct 均記錄 assembly、solver setup、linear solve、communication、
output、iterations、wall 與每 rank peak RSS。1024-element、4-rank 作業中，LU
peak rank RSS 為 412,999,680 bytes，block-Jacobi 為 131,657,728 bytes；兩者
wrapper wall 分別 2318.256／2171.397 秒。這批執行有背景活動，只作大型功能、
記憶體與成本觀察。

固定 128-element strong-scaling sweep 以 1／2／4 ranks、LU 與 block-Jacobi、
各兩次重複，共 12 次求解。28 份 rank reports、672 個檔案 hash、12 次重新執行
的物理 validator 與跨 rank checkpoint 場比較通過；最大 relative L2
`2.7998249980599e-11`。

| 候選 | 1→4 rank wall speedup | 4-rank efficiency | 1→4 rank iterations | 1→4 peak rank RSS |
|---|---:|---:|---:|---:|
| LU | 1.037× | 25.9% | 42→42 | 132,435,968→107,106,304 bytes |
| block-Jacobi | 1.012× | 25.3% | 588→3724 | 86,769,664→55,422,976 bytes |

Assembly 在各 rank 約 249–258 秒，communication 隨 rank 增加；本機沒有呈現
接近線性的加速。Block-Jacobi 降低 factor setup 與單 rank RSS，但 Krylov 次數
隨 ranks 增加。Schur-LU 在小／中網格分別累計 62,388／84,362 次線性迭代；
速度 GAMG／壓力 Jacobi 雖通過場門檻，仍有明顯成本。這些負收益支持維持現有
預設策略，候選只透過已文件化 prefix 選用。

## 證據稽核與界線

完成稽核共核對 38 個候選／scaling evaluations：35 個成功及三個預期不收斂；
同時驗證 13 份 acceptance／audit manifest。1024-element 的歷史 harness 由
`launch_commit` 的 Git blob 重建 hash，避免錯把後續 harness 更新判成舊證據損毀；
其 frozen native runtime 核心與目前 `TransientFlowRuntime.hpp` 相同。目前 HEAD
的小案例則使用現行 binary，四候選全部通過。

本機工作站數字不延伸為跨節點或其他硬體效能結論；正式 scheduler allocation、
跨節點 binding 與 scaling 由 HPC-09 驗收。HPC-04C 要求的網格加密、rank 增加、
建立／求解成本、記憶體、迭代數及負收益已具可追溯證據。
