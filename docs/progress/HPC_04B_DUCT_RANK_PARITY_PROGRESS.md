# HPC-04B 方管跨 rank 場一致性

- 狀態：同一 16-element 案例的 1／2／4-rank 數值驗證通過；整體 HPC-04 尚未完成。
- 日期：2026-09-09（美東）；本次工具基準 `eedda3e`。所有 native 作業仍使用
  frozen `01aa4c8` binary，沒有改 production headers、solver 預設或 checkpoint 格式。

新增 `hpc_compare_checkpoint_fields.py`，依 native graph 的 domain hash 選擇已驗證
manifest 中的 body-fitted owned field shards；檢查 global rows、每段 begin/end、
payload 長度、完整覆蓋、無重疊／缺口，再按全域 coefficient 順序比較。
Flow 額外分開 interleaved velocity 與 pressure 的 L2；沿用逐值及整體
`1e-12 + 1e-6*reference` 門檻。工具不替 caller 判定兩案例是否有相同幾何、basis、
物理單位或 pressure reference；本次另比對九個非 partition case 檔案逐 byte 一致，
使用同一 fixture builder、相同 2×2×4 元素及 175 nodes，僅 packed rank 數不同。

兩步 0D／3D／1D／RCR 強耦合的 1／4-rank 作業已重新執行原 Phase 9 validator，
2-rank 使用先前同一 frozen binary 的 LU 結果。三者皆通過原 physical gates。
比較各步完整 700-row 3D field，結果如下：

| 比較 | 最大 velocity 相對 L2 | 最大 pressure 相對 L2 |
|---|---:|---:|
| 1 vs 2 ranks | 3.78398e-14 | 1.19317e-14 |
| 1 vs 4 ranks | 4.05959e-14 | 1.60775e-14 |

這是 accepted checkpoint 場的比較，不宣稱支援不同 rank 間 checkpoint restart，
也沒有將 0D／1D shard 合併冒充 3D 場；它們的物理 gate 由原 validator 檢查。

## 測試與重現

三個 Python unit tests 涵蓋重排／空 owned range、重疊／缺口／global size 不同、
越界／截斷／不完整資料，以及非有限值／錯誤數值的 gate。另複製真實 bundle、
修改一個 shard byte，CLI 返回 machine-readable `status=failed`、退出 1。
初次 integration probe 使用可讀 domain 當檔名前綴而拒絕空 shard 集合；已依
`GraphCheckpointDomainPrefix` 的既有 hash 規則修正，沒有改資料格式。

```bash
python3 -m unittest discover -s scripts -p test_hpc_compare_checkpoint_fields.py
python3 scripts/hpc_duct_solver_candidates.py \
  --binary outputs/hpc04/immersed/native-v4-binary \
  --output-dir NEW_OUTPUT --ranks 4 --candidates
python3 scripts/hpc_compare_checkpoint_fields.py \
  ONE_RANK_BUNDLE FOUR_RANK_BUNDLE --domain root3d
```

`--candidates` 空清單新增為只執行 LU 參考；省略該參數仍跑原四候選。
比較工具目前會載入完整欄位，且沿用既有 manifest verifier；不宣稱它是大型資料的
串流或分散式比較器。成功／失敗均輸出 JSON，失敗退出 1。

本機證據位於 `outputs/hpc04/pde-candidates/one-rank-v1`、`four-rank-small-v1`、
原 `duct-final-v1/lu`；`crossrank-acceptance.json`、`rank-one-two.json`、
`rank-one-four.json` 保存結果，`audit_crossrank.py` 核對原物理 validator、case hashes、
shard 數／總 rows、故障拒絕與比較工具 hash。

## 資源與限制

工作站可見 16 logical CPUs（Intel Core i9-14900KF）、PETSc 3.15.5 real64/int32、
Open MPI 4.1.2、GCC 11.4，OMP／OpenBLAS=1。表中為最大 rank exclusive phase 秒：

| ranks | assembly | setup | linear solve | communication | output | peak rank RSS bytes |
|---:|---:|---:|---:|---:|---:|---:|
| 1 | 38.7478 | 1.27663 | 0.0318074 | 0.00383974 | 0.180167 | 61,480,960 |
| 2 | 30.1305 | 0.913038 | 0.0266137 | 0.531543 | 0.127022 | 59,211,776 |
| 4 | 37.9054 | 1.11640 | 0.0333100 | 3.94875 | 0.225025 | 55,476,224 |

這些數值來自有其他背景回歸的功能作業，2-rank 亦非同時重跑；不可據此算平行效率。
各 phase 的最大值不一定來自同一 rank，不得相加冒充 wall time。未涉及 CUDA allocation。

128-element、539-node、4-rank 原候選 sweep 的 LU 與 block-Jacobi 已通過，
Schur 候選仍執行中；完整 moving 回歸也尚未終結。後續須收集它們的正式結果，
繼續其他幾何／pressure reference 條件及無干擾、具排程資源的大小案例評估。
