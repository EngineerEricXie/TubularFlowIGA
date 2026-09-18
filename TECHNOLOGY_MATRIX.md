# 血管／器官 3D–1D–0D 技術、網格與現況矩陣

更新日期：2026-09-18。

配套文件：[詳細 TODO 與驗證路線](todo.md)。

範圍：本文件整理本專案與目前推薦架構相關的主要技術，不是所有 CFD／生物力學方法的完整百科。**推薦／可選方法不代表已在專案中實作。**本輪為程式與文件的唯讀核對，未重新執行數值測試。

## 1. 狀態與術語

| 標記／術語 | 含義 |
|---|---|
| 已有 | 已核對到程式、文件或既有案例；仍受其已驗證範圍限制 |
| 基礎／受限 | 有 vertical slice、專用路線或部分功能，不能推論通用化已完成 |
| 待新增 | 目前盤點沒有足以宣稱完整功能的實作／驗證依據 |
| 外部 option | 可考慮採用的外部後端；不是已安裝、已整合或已驗證的專案能力 |
| FEM／IGA／FVM | 空間離散方法；IGA 常採用與有限元素相同的弱式框架，但以樣條等基底離散 |
| 貼體 | 計算網格邊界貼合物理幾何；不限定 FEM、IGA 或 FVM |
| ALE | 移動流體網格與相應方程描述；通常仍需要流體體積網格 |
| Immersed | 物理邊界不必沿背景網格；不同 immersed 技術有不同積分與耦合成本 |
| Shell／membrane | 三維空間中的結構模型，但可用二維中面網格；殼可含彎曲，簡化膜不能自動等同完整殼 |
| 1D／0D | 1D 沿路徑／網路有空間離散；0D 為無空間網格的集中參數狀態 |

推薦總體定位：**保留貼體 IGA、低維網路與耦合基礎，發展 IGA 薄殼特色；新增通用 surface→貼體 FEM、ALE 與按需選用的組織模型。**

## 2. 3D 物理區域、方法與 mesh

同一個器官通常對應多列。心臟可包含血液、心肌、瓣膜；肝臟可包含大血管、灌流區室與組織骨架。器官名稱本身不決定方程。

| 物體／物理區域 | 推薦方法 | Mesh／輸入資料 | 其他 options | 目前能力與缺口 |
|---|---|---|---|---|
| 血管樹內血流，固定管壁 | 保留貼體 IGA Navier–Stokes；新增通用 FEM 路線 | IGA：六面體控制網格、樣條體積與 Bézier extraction；FEM：四面體／六面體流體體積網格；出入口標籤 | FVM；固定背景 immersed IGA | **已有** CPU／CUDA 支援子集的貼體 IGA，穩態／Backward Euler 暫態與 VMS；通用四面體流體 FEM 待新增 |
| 血管樹內血流，管壁變形 | 貼體 FEM 或 IGA＋ALE，耦合結構 | 隨壁面變形的流體體積網格；另有壁面殼／實體網格 | Moving immersed IGA；其他 immersed 方法 | **已有受限** moving immersed＋膜 FSI；ALE／remeshing 待新增 |
| 血管薄壁 | IGA 薄殼，作為特色主線 | 樣條中面、厚度、材料、纖維／初始應力等按需資料 | FEM 殼；簡化膜；厚壁則用實體 FEM／IGA | 現有為 **P1 三角形預張力膜**，法向位移未知數；不是 IGA 殼或通用血管材料 |
| 血管厚壁／局部壁內應力 | 三維固體 FEM，必要時混合近不可壓縮形式 | 壁厚內部體積網格、材料／纖維與約束 | 體積 IGA；薄壁近似但需驗證其可接受性 | 通用三維固體後端 **待新增** |
| 心室等器官腔室血流 | 貼體 FEM／IGA＋ALE | 腔室流體體積網格；壁面運動或心肌結構；出入口標籤 | 指定運動流體；immersed 流體 | **已有理想化指定運動 LV immersed 案例**；不是完整心肌驅動心臟 FSI |
| 心肌／厚壁器官組織 | 非線性固體 FEM；按需纖維各向異性與主動收縮 | 組織體積網格、材料、纖維、初始應力、外部支撐／約束 | 體積 IGA；僅研究血流時使用指定運動 | 通用固體、主動心肌模型 **待新增** |
| 瓣膜、薄膜與器官包膜 | 薄殼假設適用時採 IGA 殼 | 樣條中面＋厚度；接觸／夾持等條件 | FEM 殼、實體元素；ALE／immersed 流體依運動選擇 | 有簡化膜基礎；彎曲殼、接觸／閉合及完整瓣膜模型 **待新增** |
| 肝臟等器官內大血管 | 與一般血管相同：固定壁或 ALE 流體＋可選結構 | 血管腔體積網格；FSI 另需管壁／周圍支撐資料 | 部分區域保留既有 1D；3D 血管＋降階組織負載 | 管狀血流核心已有；不等於完整肝臟血管—組織模型 |
| 肝臟等組織三維灌流 | Darcy／多區室多孔介質，FEM 或 FVM | 器官組織體積網格、滲透率、區室、交換／來源／匯 | 降階區域模型；需要較多流動細節時再評估其他方程 | 三維組織灌流 **待新增**；现有代謝／壁面交換不等同此模組 |
| 肝臟等組織變形與液體互動 | 多孔彈性／雙相 FEM；不關心液體時用固體 FEM | 組織體積網格、固體材料、孔隙與滲透參數 | 體積 IGA；固定域灌流；降階力學 | 通用組織力學與多孔彈性 **待新增** |
| 3D 氧、藥物、營養物傳輸 | 對流—擴散—反應，在相應物理域求解 | 該域體積網格、濃度／源項／反應／介面資料 | 穩定化 FEM／IGA；FVM | **已有貼體 IGA** 多物種、反應、壁面交換與部分生理輔助量；所有 moving／FSI／組織組合並未完成 |

物理選項依據：[ALE FSI 官方理論](https://doc.comsol.com/6.3/doc/com.comsol.help.sme/sme_ug_theory.06.079.html)、[IGA 生物薄殼原始研究](https://biomechanics.stanford.edu/paper/CMAME15.pdf)、[雙相組織官方說明](https://febiosoftware.github.io/febio-docs/features/modules/module_biphasic/)、[肝臟多區室灌流原始研究](https://arxiv.org/abs/1605.09162)。這些來源支持方法選項，不證明本專案已完成那些功能。

## 3. 1D 方法與 mesh

| 物理用途 | 方法／可選形式 | Mesh／資料 | 目前狀態與限制 |
|---|---|---|---|
| 剛性血管網路 | `steady_poiseuille`：段阻力与壓力／流量分配 | 帶半徑 SWC／line-OBJ、節點與線段拓樸 | **已有**；隨時間入口可準靜態重算，但不表示有脈波傳播 |
| 剛性血管慣性 | `rigid_inertance`：阻力＋慣性，Backward Euler | 同一個固定截面網路 | **已有**；有慣性／相位，不具彈性管的有限速度壓力波 |
| 彈性血管脈波 | 顯式有限體積 A/Q＋Rusanov flux | 線段內的多個 1D cells、線性／Olufsen 壁法則 | **已有** `explicit_rusanov`，含 CFL substeps |
| 隱式血流 | `pressure_network`、`linearized_aq`、`nonlinear_aq`、`implicit_1d_pde` | 依形式使用集總節點／分支或多 cell 網路 | **已有** PETSc 路線；各形式保留的物理不同，不能視為同一精度模型 |
| 物種傳輸 | 守恆 `A*C`；迎風對流、擴散、反應、來源、Robin／指定壁通量 | 與血流一致的 1D cell network | **已有**多物種、代謝、氧容量、部分血氣衍生量與血管擴張回饋 |
| 組織降階描述 | 沿路徑分布的 1D 模型，或區域交換近似 | 必須定義物理路徑／區域與血管端點映射 | 可利用既有框架，但**尚無足夠依據宣稱獨立通用組織 1D 已完成**；純 well-mixed 區室應歸類 0D |

`iga_1d` 是可執行檔名稱，並不表示上述所有方法使用 IGA 基底。現有原生輸入要求帶半徑、受支援的有根網路拓樸；不將任意 centerline 自動當成合法模型。

依據：[1D 指南](docs/ONE_D.md)、[1D 範例](examples/one_d/README.md)、[傳輸核心](solvers/one_d/include/OneDTransport.hpp)。指南中的部分歷史總述已落後於新 graph 功能；方法表與較新實作需交叉核對。

## 4. 0D 方法與資料

| 物理用途 | 方法／options | Mesh／資料 | 目前狀態與限制 |
|---|---|---|---|
| 末端血管床 | RCR／Windkessel，Backward Euler 更新 capacitor pressure | 無空間網格；近端／遠端阻力、順應性、參考壓力、port | **已有 terminal RCR**；1D 邊界另支援固定壓力／阻力 |
| 流量來源／儲槽 | Source reservoir，儲存＋泵入＋阻力關係 | 無空間網格；元件與壓力狀態 | **已有 source reservoir** |
| 體外循環／VCA | Well-mixed 儲槽、泵、可選氧合器／透析／注入等專用降階部件 | 無空間網格；迴路、物種、灌流液及裝置參數 | **已有專用路徑**；native 1D 與部分 CPU 3D 支援，模式／物種／restart 各有契約 |
| 完整閉環心臟／循環 | 時變彈性心腔、瓣膜、循環元件 ODE／代數系統等 | 無空間網格；完整元件圖、初值與參數 | **待新增／整合**；VCA closed-loop 不能當成此功能 |
| 通用區室物種交换 | 儲存、交換、反應、來源／匯的物質量 ODE | 無空間網格；區室容積、濃度與交換參數 | VCA 有專用功能；**generic 0D graph species 尚未完成** |
| 器官區域平均負載／代謝 | 经校準的阻力、順應性、交換／代謝區室 | 無空間網格；區域與血管端點對應 | 可由現有基礎擴充，但不可自動宣稱有肝臟等器官特異模型 |

依據：[0D 核心](include/ZeroDFlowDomain.hpp)、[0D 多尺度驗證](docs/progress/PHASE_9_REPORT.md)、[VCA 範例](examples/vascular_flow/vca_bifurcation/README.md)。

## 5. 輸入、前處理與 mesh 路由

| 輸入 → 目標 | 主要技術 | 目前狀態 | 重要限制 |
|---|---|---|---|
| Centerline → 1D | 半徑／拓樸驗證、線段、cells | **已有** SWC／半徑 line-OBJ | 不是一般 surface OBJ；無半徑不足以定義截面 |
| Centerline → 貼體 3D IGA | 平滑／重採樣、六面體控制網格、樣條、Bézier extraction、METIS、`.ntiga` | **已有** | 管狀、有根樹、受支援二分岔與幾何品質限制；非通用任意器官 |
| Surface → immersed IGA | 封閉三角表面、Cartesian cubic B-spline、cut volume／surface quadrature、Nitsche、ghost | **已有** | 幾何需滿足封閉、方向／標籤等契約；背景也有體積 cell，不是只算表面 |
| Surface → 貼體 FEM | 修復／cap／標記、四面體或混合體積網格、資料轉換 | **完整流程待新增** | Gmsh／TetGen 等為候選，不代表已整合；需另查授權與環境 |
| Surface → 貼體體積 IGA | 樣條擬合、patch 分解、體積參數化、正 Jacobian | **通用流程待新增** | 光滑表面不等於有效體積樣條；擬合不是精確還原影像 |
| Surface／樣條曲面 → IGA 殼 | 分析用中面、patch 接合、厚度與材料 | **待新增** | 現有體積 extraction 不等於殼求解器；高階連續性需處理 |
| 器官外表面 → 組織體積模型 | 組織 meshing、內部區域、血管定位、材料與交換映射 | **待新增** | 外表面不能唯一決定內部血管、纖維、孔隙、滲透或初始應力 |

依據：[控制網格生成器](preprocessing/mesh/README.md)、[Bézier 前處理](preprocessing/spline/README.md)、[樣條體積參數化研究](https://arxiv.org/abs/1902.00650)。

## 6. 耦合與共同基礎設施

| 技術層 | 推薦設計／options | 現有能力 | 待新增／需保留限制 |
|---|---|---|---|
| 0D↔1D↔3D 壓力／流量 | 顯式或分區強耦合、SI、outward-positive ports | 異質 graph、transactional runtime、收支檢查 | 受支援 acyclic topology；不可稱任意循環 graph 全部可用 |
| 1D↔3D 物種 | 依實際流向交換濃度／守恆通量，處理逆流 | Native 1D 與貼體 3D 的 schema-v6 路線 | 不涵蓋通用 0D、immersed、所有 moving／組織 |
| 流體↔結構 | 強 Dirichlet–Neumann＋Aitken；可選 IQN／monolithic | 已有 bounded FSI 與分散式流體／單 owner 膜 | 通用殼／實體、非匹配介面、IQN／monolithic 待新增 |
| 血管↔組織 | Port／分布來源與區域映射、質量与物種交換 | 可借用既有 port／狀態基礎 | 三維灌流／組織的完整耦合未完成；需防止區室重複計數 |
| 幾何運動 | 固定域、ALE、immersed 明確分路 | 已有固定域與 moving immersed | ALE、mesh quality controller、remeshing／守恆場轉移待新增 |
| 時間積分 | 先沿用已驗證 Backward Euler；高階／multirate 獨立立項 | 3D／0D 的 BE、1D 多種顯隱式形式 | 全框架通用 multirate 或高階 FSI 不應視為已有 |
| CPU | C++、MPI／PETSc、OpenMP；配置適合的 solver／PC | 貼體與 immersed 分散式路線、局部平行組裝 | 結構膜矩陣仍單 owner；不能將小測 parity 當大型 scaling |
| GPU | 適用數值核心可用 CUDA | 單 GPU 貼體 standalone 支援子集 | CUDA FSI／ALE／通用多 domain／組織尚不能宣稱可用 |
| 狀態／重啟 | Trial、prepare、finalize、rollback、checkpoint | 多條路線已有；moving／FSI 有跨 rank restart 證據 | Library／test 與使用者 CLI 完整度不同；需逐路線核對 |
| I/O 與診斷 | VTK／ParaView、QoI、守恆、收斂、source/input hashes | 已有基礎與多個驗證報告 | 新物理仍要新增測試；動畫不是驗收證據 |

較新分散式狀態依據：[moving／FSI restart](docs/progress/HPC_05D_MOVING_FSI_RESTART_REPORT.md)。不要沿用早期 Phase 文件將現在所有 immersed／FSI 都說成 `COMM_SELF`；也不要反過來把單一 heartbeat driver 的單 rank 結果說成已測量多節點加速。

## 7. 代表性應用組合

| 使用情境 | 可組合的模型 | 目前差距 |
|---|---|---|
| 剛壁血管樹＋末端床 | 貼體 3D IGA＋1D 網路＋0D RCR | 核心已有；實際幾何／graph／誤差需求仍須驗證 |
| 可變形血管樹，保留 IGA 特色 | 貼體 FEM 或 IGA＋ALE，耦合 IGA 殼，再接 1D／0D | ALE、IGA 殼及新介面待新增 |
| Surface 輸入的複雜固定血管 | Surface→FEM volume→流體；或現有 immersed IGA | 前者待建；後者有幾何與成本限制 |
| 指定運動心室血流 | ALE 流體＋指定運動；或現有 prescribed immersed | 既有理想化 immersed 基礎；ALE 與通用輸入待新增 |
| 主動心肌與血流互動 | ALE 流體＋非線性心肌固體＋適用的 0D／1D 邊界 | 心肌材料、主動收縮、ALE 與完整整合待新增 |
| 肝臟血管流場，組織只作負載 | 3D 血管＋經定義／校準的降階組織交換模型 | 血流基础可用；組織模型不能憑既有代謝功能直接宣稱完成 |
| 肝臟三維灌流 | 3D／1D 血管＋3D 多區室 Darcy＋0D 邊界 | 多孔介質、來源／匯映射及器官參數待新增 |
| 肝臟變形影響灌流 | 上述模型＋多孔彈性／雙相組織 | 通用組織力學與雙向交換待新增 |

## 8. IGA 應保留在哪裡？

| 角色 | 建議定位 | 證明價值的方式 |
|---|---|---|
| 現有貼體 IGA 流體 | 保留，作為已有後端與比較基準 | 相同 QoI 誤差下比較總成本，而非僅比 DOF |
| IGA 薄殼 | 優先特色主線，但需新實作 | 彎曲／膜變形、材料、patch 接合、與基準解的收斂 |
| 樣條血管幾何 | 保留中心線生成與平滑曲面能力 | 幾何誤差、可靠性、前處理與參數研究成本 |
| 體積 IGA 器官 | 選配，不要求所有器官都使用 | 體積參數化品質、精度與成本是否優於替代路線 |
| Immersed IGA | 保留研究分支與適合案例 | 免除貼體 meshing 的收益，是否抵銷 cut 積分／幾何更新成本 |
| Divergence-conforming IGA 等高階研究 | 獨立 option | 相容空間、映射與守恆證據；不能視為現有 B-spline 自動具備 |

## 9. 驗證聲明與來源索引

最近分岔 FSI 加速交付完成 20 步，solver 約 6968.953 s，19 項既定等價性／輸出等檢查通過；但仍保留 **17/20 步原守恆失敗**，屬 visualization delivery，非正式物理驗收。此限制不應被擴大成其他所有案例都失敗，也不能被省略後宣稱本案例高擬真已驗證。

- [專案總覽與支援限制](README.md)
- [CPU 求解器](solvers/cpu/README.md)
- [CPU 數值架構](solvers/cpu/ARCHITECTURE.md)
- [CUDA 支援](solvers/cuda/README.md)
- [1D 方法](docs/ONE_D.md)
- [0D 核心](include/ZeroDFlowDomain.hpp)
- [實際 P1 膜核心](solvers/cpu/include/PretensionedMembrane.hpp)
- [理想化 prescribed LV 驗證](docs/progress/PHASE_7_REPORT.md)
- [較新的 moving／FSI restart 與限制](docs/progress/HPC_05D_MOVING_FSI_RESTART_REPORT.md)
- [分岔 FSI 最終加速與物理限制](docs/progress/BIFURCATION_FSI_TIME_OPTIMIZATION_RESULTS.md)

本文件的現況是上述核對日期的 snapshot。新增後端、修正守恆或改變可執行拓樸後，應同時更新本矩陣與 [TODO](todo.md)，不可只更改「已有」標籤而沒有對應 evidence。
