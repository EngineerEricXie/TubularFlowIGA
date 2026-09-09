# HPC-05A：Checkpoint 狀態契約稽核

日期：2026-09-09。基準 revision：`70b1715992912871eebb4a6458634bedb1b9a542`。
狀態：**HPC-05A 已完成**；HPC-05B／C／D 尚未完成，整份 HPC TODO 繼續執行。

## 交付與驗收範圍

交付 [耦合 checkpoint 狀態契約](../architecture/COUPLED_CHECKPOINT_CONTRACT.md)。
本項是 TODO 明列的狀態盤點／介面稽核，沒有修改 numerical runtime、CLI 或檔案格式。
程式現況和未來必須履行的契約在文件中分開標示。

本次驗收標準是：列齊時間／步數、流場／history、species、0D／outlet／VCA、
structure／FSI、moving geometry 與 donor hysteresis；每項有實際 state owner／方法
依據及保存或重建的決策；明定 all-domain accepted macro-step 邊界；列出既有
serializer／restore 的缺口和 05B／C／D 可執行的後續 gates。

| 已核對項目 | 結論 |
|---|---|
| graph clock／commit 順序 | 在 domain finalize、donor 發布、accepted result／pressure 更新與 accepted clock 更新全部完成後擷取；使用實際累加時間，保存 domain 自有 counters |
| pressure-flow／species Aitken | residual／omega 在每次 macro-step 重建；外部下一步 pressure map 必須保存 |
| donor hysteresis | 保存所有 `(edge,species)` donor；近零 flow 不能由空 map 重算 |
| 1D | 完整 hydraulic state、substeps、species concentrations／flux／accounting、LastInlet 及 dynamic radius；area0／resistance 可重建，configuration 本身未被 vasodilation 修改 |
| body-fitted 3D | BE 下一步會從 accepted field 複製 previous velocity；transport 的 `previous_` 是每步重建的 matrix，真正持久化欄位為 concentration 與 step counter |
| 0D／outlet／VCA | 保存 pressure storage、RCR capacitor、accepted ports／accounting、reservoir 與 last arterial species；輸出帳目不能用當前狀態猜測 |
| immersed／moving | owned fields、controllers／gauge／clock、material source topology／kinematics、predecessor publication provenance、conservation primitives；trial extension／force 按下一步重建 |
| membrane／FSI | committed displacement／velocity、reference／clamp model、兩側完整 committed publications；Strong FSI 每步 reset Aitken，predictor 使用已恢復的結構狀態 |
| 可觀察歷史 | accepted records prefix 有 epoch、checksum、範圍與筆數；新 job 接續輸出不得重複或只留尾段；執行時間另計 |

## 找到的實作缺口

- `CoupledDomainRuntime` 尚無 checkpoint capability；native graph runner 沒有完整
  bundle checkpoint／restart 入口。不能把各 domain 的 standalone serializers
  依序呼叫，便宣稱具備多域原子發布。
- `OneDCheckpoint v2` 與現有 `RestoreCommittedState` 未保存／恢復完整 LastInlet
  和 species boundary flux／step accounting。3D transport `ReadState()` 把 counter
  設為 1，必須新增可驗證的真實 counter restore。
- 0D constructor 未恢復 committed index／count／publication；species executor
  目前只有 donor getter，沒有完整 restore。浸入式的場 setter 也未同時恢復
  adapter clock／publication／commit count，transient adapter 仍要求 fresh index 0。
- moving geometry 的 publication hash 依賴 predecessor geometry digest 與 cell
  transitions；一般 constructor 建出的 genesis publication 不相同。
  `SetCommittedGlobalState()` 還會清除 conservation，不能冒充完整恢復。
- membrane、FSI lifecycle、fluid material／traction publication 需要配對 candidate
  restore。現有 visualization snapshot 不包含上述完整持久化狀態。
- VCA 格式缺完整 circuit cache／accepted history；1D CLI 明確拒絕 closed-loop
  checkpoint。此稽核沒有替尚未存在的 VCA／CUDA graph provider 宣稱支援。

## 檢查與證據

使用本工作站的唯讀原始碼稽核，單一 Python 3 程序及 Git；未啟動 MPI 或 GPU。
input identity 為上述基準 revision 及契約所連結原始碼的 SHA-256 清單，保存在
`outputs/hpc05/contract/audit.json`。同目錄的 `audit.py` 是本次執行的文件檢查程式，
檢查契約／本報告的本機連結、38 個任務的核取方塊與來源檔案身分；不屬於 solver 測試。

```bash
python3 outputs/hpc05/contract/audit.py
git diff --check
```

結果：文件連結、來源清單與 TODO 計數檢查通過，退出碼 0；`git diff --check`
退出碼 0。清單為 13／38 項完成，25 項未完成，其中 2 項已有部分進度。

數值 relative L2、守恆、收斂原因、MPI／thread scaling、assembly／solve 時間、
communication／I/O 吞吐、host peak RSS、CUDA peak allocation：**N/A**。
原因是本次只完成狀態與介面契約，沒有新的求解或 checkpoint I/O 實作，依 TODO
「僅文件、配置或介面稽核任務」規則記錄；沒有將舊數值測試視為本次 restart 驗收。

## 剩餘工作與下一步

HPC-05B 先實作版本化 manifest、有限大小的嚴格解析、temporary shards、checksum、
sync、完成標記最後發布，以及混 epoch／截斷／缺片／錯配置拒絕。
接著 05C 補齊既有 graph 的 accepted capture／candidate restore，執行獨立 MPI job
的 uninterrupted 對 restart 全歷史比較及三個中斷點測試。05D 在相應 runtime 可用後
加入 immersed／moving／FSI 的 publication 與 repartition 恢復。

契約第 6 節列出的所有 restart／故障／效能 gates 仍待上述實作後執行，未在本次勾選。
