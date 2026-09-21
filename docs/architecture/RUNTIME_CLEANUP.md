# 貼體 runtime 的終止清理契約

`TransientFlowRuntime` 與 `TransientTransportRuntime` 借用 communicator；
呼叫端必須讓該 communicator 與 PETSc 存活到 runtime 完成清理。
runtime 不擁有或釋放 communicator。

## 正常結束

呼叫端先完成所需的 field gather、診斷與 checkpoint，再由同群所有 rank
以相同順序呼叫 `Close()`。此為終止操作，之後只能再次 `Close()` 或析構，
不能重用 runtime 求解或查詢 PETSc 狀態。需要重試模擬時建立新 runtime。

Close 依固定順序嘗試全部 destroy，以不配置記憶體的 result 保存第一個
回傳錯誤及物件名稱；完成所有 destroy 後，才進入 `CollectiveLocalStage`
協調失敗 rank／錯誤訊息。各 destroy 不放入僅允許本地工作的 callback。
多個 rank 失敗時沿用最小 rank 的診斷，該 rank 保留第一個失敗物件。

| Runtime | 釋放順序 |
|---|---|
| Flow | solver、scatter、destination IS、ghost previous、ghost state、source IS、rhs、update、previous、committed state、state、Jacobian |
| Transport | solver、scatter、destination IS、ghost state、source IS、rhs、next、committed、current、forcing、previous matrix、left matrix |

每個 runtime 最多執行一次這組釋放。重複 Close 重新協調已保存的結果；
首次失敗不會在第二次變成成功，也不在某些 rank 上重試不確定的 destroy。
析構保留 noexcept 的資源釋放後備路徑，但無法向呼叫端回傳成功／失敗；
因此 embedding caller 必須明確 Close，才能把清理列入作業成功條件。

## 呼叫端與發布順序

| 入口 | 最後使用與成功回報之間的清理 |
|---|---|
| CPU flow／VCA | 完成最終 field、index、VTKHDF close 後，先 Close 可選 transport，再 Close flow，最後列印 completion／profile |
| Multidomain／bifurcation | 全部 accepted step bookkeeping 後，按 domain map 順序 Close transport／flow；然後寫 graph CSV／completion manifest 與 stdout |
| Sequential explicit／fixed／Aitken | 最後 accepted history bookkeeping 後 Close 3D runtime；然後寫 CSV／completion manifest 與 stdout |

Graph／sequential 最後的 writer 只使用已保存的 accepted histories、配置、
資產與容器數量，不再讀取已關閉的 PETSc state。
Registry 擁有的貼體 adapters 借用這些 runtime，析構不呼叫其求解或 abort。

CPU flow 的既有 accepted checkpoint／field 可能早於清理失敗寫出；
它們表示已接受的步狀態，不是整個作業成功標記。原子 checkpoint bundle、
最後完整 epoch 與 job restart 行為遵循
[coupled checkpoint contract](COUPLED_CHECKPOINT_CONTRACT.md)。

## 失敗範圍

建構或求解已失敗時，析構嘗試同一固定清理順序，保留原本失敗退出語義。
這不替代 runtime 內部的 collective ordering：各 rank 必須先完成共同失敗
協議再一起退棧；單純在最外層 catch 中處理無法救回已卡在 collective 的 peer。

驗收注入在真實 destroy 返回後覆寫回傳碼，證明共同退出、剩餘釋放、
引用計數與新 runtime 重試。它不模擬 destroy 內部永久等待、MPI 失聯、
程序死亡或已毀損 PETSc 內部物件；不宣稱能恢復這些故障或保證其無洩漏。

其他 CLI 的成功路徑已有個別明確清理，例如 configured／legacy transport、
assembly smoke、四種 1D implicit 方法。各入口仍須由對應 regression 驗證。
