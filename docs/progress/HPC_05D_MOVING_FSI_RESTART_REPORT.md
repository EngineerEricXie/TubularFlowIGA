# HPC-05D：moving／FSI restart 與重分區完成報告

日期：2026-09-10。此報告完成 CPU 浸入式／moving-flow 與 bounded
moving-FSI runtime 的 checkpoint／restart；native graph CLI 與跨節點排程
仍分別由 HPC-05C、HPC-07D／09 追蹤。

## 交付內容

- Moving-flow bundle 保存 accepted clock、owned PETSc field、前一／目前材料幾何、
  geometry predecessor、history extension provenance、ports、可變 controls 與
  conservation state。Fresh runtime 先驗證完整 epoch，再發布候選。
- Field restore 依 stable row ID 串流到目標 ownership；4→2、4→1 與 changing-layout
  4-rank continuation 均已通過，不以改啟動參數假裝重分區。
- Paired FSI catalog 在同一 epoch 加入每個來源 rank 的 owned traction publication，
  以及單一 owner 的 membrane metadata／numerical state。任何缺片、checksum 或
  context 錯誤都在回傳 fresh pair 前清理候選。
- Surface restore 驗證來源 reference identity、共同 layout／composition／context、
  不重複 partition identity，以及 stable node ID 的精確全域覆蓋；接著只保留目標
  rank 擁有的 traction／nodal force，並建立新的 partition、producer 與 projection
  provenance。Membrane owner 可隨目標 communicator 改變。
- 相同 rank 與 ownership 的路徑仍使用嚴格 parser，保留 publication bytes 的精確
  往返；跨 rank 路徑以既定數值門檻比較，不要求不同 reduction 順序逐位元相同。

## 驗收結果

GCC 11.4、Open MPI 4.1.2、PETSc 3.15.5 real64／Int32／MUMPS；本機
OMP／OpenBLAS 均為 1 thread。跨 rank 作業使用已終止 writer 所發布的 nonzero
bundle，reader 只在 checksum fault injection 期間暫時改動 membrane shard，並在
繼續前還原；checker 再要求所有輸入檔案 hash 與執行前完全相同。

| 來源→目標 ranks | accepted field scaled L2 | 下一步 field scaled L2 | 下一步迭代 | wall time／每 rank peak RSS | 稽核 |
|---|---:|---:|---:|---|---|
| 4→2 | `5.370410638499683e-13` | `1.4317113742544146e-12` | 7 | 1701.46 s；121.1／149.9 MiB | `paired-repartition-fixed-4-to-2-audit.json` |
| 4→1 | `3.303457280884551e-13` | `1.679465531835937e-12` | 7 | 3332.96 s；135.5 MiB | `paired-repartition-fixed-4-to-1-audit.json` |
| 1→4 | `5.958747452066621e-13` | `1.159705383816217e-12` | 7 | 1011.63–1011.68 s；98.3–136.0 MiB | `paired-repartition-layout-validation-1-to-4-audit.json` |

三組都低於預先固定的 `1e-8` field 門檻，並通過 accepted surface、下一步
displacement／velocity／traction／consistent force、完整 Aitken history、ports、
五項 conservation diagnostics、rank reports、timeout、stdout／stderr digest 與
來源 bundle 不變檢查。先前相同 rank 的 1／2／4-rank pair、新作業 zero／nonzero
reader、moving flow 4→2／4→1 與可變 port control 續跑證據，詳見
[HPC-07 進度](HPC_07A_MATERIAL_COMPOSITION_PROGRESS.md)。

第一輪 1→4 正確拒絕了以目標局部 ownership 誤建的單 rank 來源 layout；四份
reports 均 exit 1，沒有 continuation 數值宣稱。失敗及修正來源綁定於
`paired-repartition-fixed-1-to-4-failure-audit.json`。修正後改為比對已認證來源
shard 的共同 layout identity，再以 reference identity 與完整 stable-node coverage
綁定目標；上述 1→4 結果使用修正後 frozen binary。舊 4→2／4→1 binary 的正常
路徑已使用同一 coverage／provenance 實作；最後的修正只移除無效的來源-layout
重建，兩批精確來源分別由 `paired-repartition-fixed-source-manifest.json` 與
`paired-repartition-layout-validation-source-manifest.json` 保存。

重現 checker：

```bash
python3 scripts/hpc_check_paired_fsi_restart.py \
  --run outputs/hpc03/moving-graph-v1/paired-repartition-fixed-4-to-2 \
  --source-ranks 4 \
  --input-hashes outputs/hpc03/moving-graph-v1/paired-repartition-source-4-input-hashes.json
```

4→1 使用相同來源與 `--source-ranks 4`；1→4 改用
`paired-repartition-layout-validation-1-to-4`、來源 rank 1 的 input-hash catalog。
驗收器另已用健康 fixture 與把 accepted error 改成 `1` 的負向 fixture 自測。
本批來源、forced rebuild、codec／sanitizer／zero smoke、三組 matrix、負向案例
與限制的最終索引為 `paired-repartition-completion-audit.json`。

## 限制

這是本機功能驗證，並行執行時另有工作負載，wall time 不作 scaling 結論。
Moving／FSI bundle 目前是 library／test integration；正式 native graph 使用者 CLI、
大型 scheduler signal forwarding 及跨節點結果尚未完成。Membrane 數值矩陣仍在
單一 bounded owner 求解，沒有宣稱結構矩陣已分散。
