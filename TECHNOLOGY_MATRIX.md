# 血管／器官 3D–1D–0D 技術、網格與現況矩陣

更新日期：2026-09-19。

範圍：本文件整理本專案與目前推薦架構相關的主要技術，不是所有 CFD／生物力學方法的完整百科。**推薦／可選方法不代表已在專案中實作。**各項「已有」仍以其連結契約與本地 audit 的明示範圍為限。

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

人工正方體器官案例包含互不重疊的動脈／靜脈樹、固定組織 Darcy、FEM 血管壁 FSI，以及明示來源／匯與守恆再現。它只屬功能驗證，不是肝臟解剖或生理成果。真實完整肝臟、組織變形、瓣膜與心肌材料／纖維／主動收縮尚未實作。原生 ALE／固體／分區 FSI 的受限元件可用於血管壁路線，但尚未構成雙樹加 Darcy 的完整案例。既有 [肝臟幾何審查](docs/LIVER_GEOMETRY_CANDIDATES.md)與局部 ROI 功能證據保留，不把局部結果稱為全肝灌流。

## 2. 3D 物理區域、方法與 mesh

同一個器官通常對應多列。心臟可包含血液、心肌、瓣膜；肝臟可包含大血管、灌流區室與組織骨架。器官名稱本身不決定方程。

| 物體／物理區域 | 推薦方法 | Mesh／輸入資料 | 其他 options | 目前能力與缺口 |
|---|---|---|---|---|
| 血管樹內血流，固定管壁 | 保留貼體 IGA Navier–Stokes；新增通用 FEM 路線 | IGA：六面體控制網格、樣條體積與 Bézier extraction；FEM：四面體流體體積網格；出入口標籤 | FVM；固定背景 immersed IGA | **已有** CPU／CUDA 支援子集的貼體 IGA；另有原生 C++/PETSc tetra P2/P1 FEM 功能 vertical slice、固定域 Backward-Euler 一階時間收斂、manufactured P2 三階 L2 空間收斂及理想管／Y 數值證據。FEM 尚缺三層 strict Poiseuille 空間 gate、大網格可擴展 Schur solver 與生理驗證 |
| 血管樹內血流，管壁變形 | 貼體 FEM 或 IGA＋ALE，耦合結構 | 隨壁面變形的流體體積網格；另有壁面殼／實體網格 | Moving immersed IGA；其他 immersed 方法 | **已有受限** moving immersed＋膜 FSI；原生 tetra ALE 已完成運動學、harmonic mesh motion、GCL、rollback、暫態與小型 MPI runtime，matching native ALE-flow/solid 全鏈亦通過零狀態 invariant 與交易 gate；原生 matching ALE-flow/solid 已有非零牽引驅動與理想化 compliant-channel 多步功能 smoke，但限小型 dense 單 partition；生產 PETSc/distributed FSI、彈性管物理基準、雙樹－Darcy 共同案例與 remeshing 仍未完成 |
| 血管薄壁 | IGA 薄殼，作為特色主線 | 樣條中面、厚度、材料、纖維／初始應力等按需資料 | FEM 殼；簡化膜；厚壁則用實體 FEM／IGA | **已有受限**單 patch NURBS Kirchhoff–Love 殼的原生 element／小型 static runtime 與膜、彎曲、tangent、收斂、traction／constraint 測試；另保留 P1 預張力膜。多 patch、血管材料與 FSI 尚未完成 |
| 血管厚壁／局部壁內應力 | 三維固體 FEM，必要時混合近不可壓縮形式 | 壁厚內部體積網格、材料／纖維與約束 | 體積 IGA；薄壁近似但需驗證其可接受性 | **受限原生路線已有** total-Lagrangian P1 tetra compressible neo-Hookean 與 stabilized mixed P1/P1 `u-p` formulations、analytic tangents、dense global Newton、locking verification、nodal load／約束／反力、inverse-elastostatics prestress 及 matching FSI structure adapter；surface load 的 production integration、空間收斂與 PETSc scalability 尚未完成 |
| 心室等器官腔室血流 | 貼體 FEM／IGA＋ALE | 腔室流體體積網格；壁面運動或心肌結構；出入口標籤 | 指定運動流體；immersed 流體 | **已有理想化指定運動 LV 功能週期**：原生 ALE（96／384／672 tetra）與單程序 immersed 共用 surface／motion／端口 `β=0.5`，各完成 16 步與逐步共點場量對照；三層的速度相對 L2 峰值為 `71.37%／58.15%／54.22%`，改善但未達空間收斂。三層原生 ALE 已量測平均壓力、壓力—容積積分與動能。96／384-tetra ALE 再跑第 2 個週期後，相同 ED 幾何的速度場回返相對 L2 仍為 `27.40%／37.22%`，尚非週期流場。immersed cut-cell ED／ES 體積積分較共同表面體積高 `1.36%`／`1.43%`。空間／時間收斂、完整能量平衡、分散式 immersed 回流、生理驗證及心肌 FSI 仍未完成 |
| 心肌／厚壁器官組織 | 非線性固體 FEM；按需纖維各向異性與主動收縮 | 組織體積網格、材料、纖維、初始應力、外部支撐／約束 | 體積 IGA；僅研究血流時使用指定運動 | 有各向同性 compressible tetra、stabilized mixed P1/P1 `u-p` near-incompressible locking verification，以及明示 load/material/support 的小型 inverse-elastostatics prestress 初始化；patient-specific 參數識別、纖維、主動收縮與 distributed mixed solver **待新增**，因此不可稱心肌模型完成 |
| 瓣膜、薄膜與器官包膜 | 薄殼假設適用時採 IGA 殼 | 樣條中面＋厚度；接觸／夾持等條件 | FEM 殼、實體元素；ALE／immersed 流體依運動選擇 | 有簡化膜與受限原生 KL 彎曲殼基礎；多 patch、接觸／閉合及完整瓣膜模型 **待新增** |
| 肝臟等器官內大血管 | 與一般血管相同：固定壁或 ALE 流體＋可選結構 | 血管腔體積網格；FSI 另需管壁／周圍支撐資料 | 部分區域保留既有 1D；3D 血管＋降階組織負載 | 管狀血流核心已有；不等於完整肝臟血管—組織模型 |
| 肝臟等組織三維灌流 | Darcy／多區室多孔介質，FEM 或 FVM | 器官組織體積網格、滲透率、區室、交換／來源／匯 | 降階區域模型；需要較多流動細節時再評估其他方程 | **已有受限自主 P1 tetra 單區室穩態 Darcy 計算層、JSON CLI 與 surface／volume 工作流**，含解析／manufactured、1／2／4 MPI、rank-owned 壓力／通量 VTU、fTetWild／Gmsh 標籤網格、守恆面通量及 RT0 H(div) 後處理功能測試；通用 port／matching-facet 流量→cell source 映射通過立方體／Y 功能測試；同病例局部 ROI 亦有明示人工 BC／材料的自主流體→Darcy、逐面守恆與可重跑證據，但尚無完整肝臟開口／材料契約、雙樹共同案例或多區室，故不可宣稱器官灌流完成 |
| 肝臟等組織變形與液體互動 | 多孔彈性／雙相 FEM；不關心液體時用固體 FEM | 組織體積網格、固體材料、孔隙與滲透參數 | 體積 IGA；固定域灌流；降階力學 | 通用組織力學與多孔彈性 **待新增** |
| 3D 氧、藥物、營養物傳輸 | 對流—擴散—反應，在相應物理域求解 | 該域體積網格、濃度／源項／反應／介面資料 | 穩定化 FEM／IGA；FVM | **已有貼體 IGA** 多物種、反應、壁面交換與部分生理輔助量；另有自主 tetra ALE P1 物種保守 weak form、PETSc/MPI 稀疏代數、外部標籤網格＋指定速度的版本化 CLI／P1 VTU、可選隱式一階衰減、指定外部濃度的逐標籤壁面交換，及單／多標籤 native FEM↔固定體積 0D 組織儲槽同一步保守交換；後者已有 CLI／surface／volume 工作流及成對跨行程重啟。多物種反應、任意案例生產 graph restart、FSI／器官週期及一般高 Péclet positivity 驗證仍缺 |

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
| 通用區室物種交换 | 儲存、交換、反應、來源／匯的物質量 ODE | 無空間網格；區室容積、濃度與交換參數 | 已有獨立 variable-volume well-mixed 多物種更新、source／terminal RCR staged graph adapters；受控正／逆流、回滾與 source→native FEM→RCR 移動兩步功能鏈及固定測例五 shard 檔案重啟通過。**反應、任意案例生產 checkpoint 及生理混合體積校準未完成**，VCA 專用功能不得替代 |
| 器官區域平均負載／代謝 | 经校準的阻力、順應性、交換／代謝區室 | 無空間網格；區域與血管端點對應 | 可由現有基礎擴充，但不可自動宣稱有肝臟等器官特異模型 |

依據：[0D 核心](include/ZeroDFlowDomain.hpp)、[通用物種契約](docs/T7_GENERIC_ZERO_D_SPECIES.md)、[0D 指南](docs/ZERO_D.md)、[VCA 範例](examples/vascular_flow/vca_bifurcation/README.md)。

## 5. 輸入、前處理與 mesh 路由

| 輸入 → 目標 | 主要技術 | 目前狀態 | 重要限制 |
|---|---|---|---|
| Centerline → 1D | 半徑／拓樸驗證、線段、cells | **已有** SWC／半徑 line-OBJ | 不是一般 surface OBJ；無半徑不足以定義截面 |
| Centerline → 貼體 3D IGA | 平滑／重採樣、六面體控制網格、樣條、Bézier extraction、METIS、`.ntiga` | **已有** | 管狀、有根樹、受支援二分岔與幾何品質限制；非通用任意器官 |
| Surface → immersed IGA | 封閉三角表面、Cartesian cubic B-spline、cut volume／surface quadrature、Nitsche、ghost | **已有** | 幾何需滿足封閉、方向／標籤等契約；背景也有體積 cell，不是只算表面 |
| Surface → 貼體 FEM | 共用 closed-surface preflight、標記、四面體網格、資料轉換 | **已有受限** Gmsh exact-boundary 與 fTetWild envelope-remesh adapter，皆輸出原生 reader 使用的 labelled ASCII Gmsh 4.1 | 單一 fluid region、P1 geometry；不自動 cap／repair；PSC 安裝與大型案例仍需另驗 |
| Surface → 貼體體積 IGA | 樣條擬合、patch 分解、體積參數化、正 Jacobian | **通用流程待新增** | 光滑表面不等於有效體積樣條；擬合不是精確還原影像 |
| Surface／樣條曲面 → IGA 殼 | 分析用 NURBS 中面、patch 接合、厚度與材料 | **已有受限**原生單 patch KL shell data/element/runtime | 一般三角 surface 尚無自動 NURBS fitting；非 C1 seams 目前 fail closed，多 patch continuity 待新增 |
| 器官外表面 → 組織體積模型 | 組織 meshing、內部區域、血管定位、材料與交換映射 | **待新增** | 外表面不能唯一決定內部血管、纖維、孔隙、滲透或初始應力 |

依據：[控制網格生成器](preprocessing/mesh/README.md)、[Bézier 前處理](preprocessing/spline/README.md)、[樣條體積參數化研究](https://arxiv.org/abs/1902.00650)。

## 6. 耦合與共同基礎設施

| 技術層 | 推薦設計／options | 現有能力 | 待新增／需保留限制 |
|---|---|---|---|
| 0D↔1D↔3D 壓力／流量 | 顯式或分區強耦合、SI、outward-positive ports | 異質 graph、transactional runtime、收支檢查；原生 tetra ALE 與 1D，以及 source→tetra→RCR 的 1／2／4-rank 功能圖測試 | 受支援 acyclic topology；原生 0D 強耦合在低慣性及調整阻抗後的高密度功能組通過，均非生理驗證；不可稱任意循環 graph 全部可用 |
| 1D↔3D 物種 | 依實際流向交換濃度／守恆通量，處理逆流 | Native 1D 與貼體 3D 的 schema-v6 路線；另有原生 tetra ALE 兩步移動 1D／0D source→tetra→RCR 功能圖及固定測例五 shard restart | immersed、混合方向單一 port、所有 moving／組織與任意案例生產 graph/file restart 尚未完成 |
| 流體↔結構 | 強 Dirichlet–Neumann＋Aitken；可選 IQN／monolithic | 已有 bounded immersed FSI；另有原生 matching tetra ALE-flow/solid dense reference vertical slice，通過總力／力矩／離散功率、zero-state invariant、非零 traction-driven response、idealized compliant-channel 單步功能 smoke、獨立殘差、reject、prepare-abort 與 atomic commit gate | 尚缺 native compliant-channel／elastic-tube 的收斂與 added-mass 物理驗證、production PETSc adapter、distributed ownership；IGA 殼 map、非匹配介面、IQN／monolithic 待新增 |
| 血管↔組織 | Port／分布來源與區域映射、質量与物種交換 | 可借用既有 port／狀態基礎；單／多 wall label 的 native tetra species↔固定體積 0D 儲槽已驗證同一步正反向物質收支、使用者 CLI 與成對重啟；雙區相反方向交換與跨區響應亦通過 | 幾何到多區域的自動映射、三維灌流／組織的完整耦合未完成；需防止區室重複計數 |
| 幾何運動 | 固定域、ALE、immersed 明確分路 | 已有固定域、moving immersed，以及原生 tetra ALE/harmonic mesh-motion/quality rollback vertical slice；idealized LV 的共享幾何／壁運動、流場 QoI 與共點場量已對照 | deforming-chamber 場量收斂、remeshing／守恆歷史場轉移待新增 |
| 時間積分 | 先沿用已驗證 Backward Euler；高階／multirate 獨立立項 | 3D／0D 的 BE、1D 多種顯隱式形式 | 全框架通用 multirate 或高階 FSI 不應視為已有 |
| CPU | C++、MPI／PETSc、OpenMP；配置適合的 solver／PC | 貼體與 immersed 分散式路線、局部平行組裝 | 結構膜矩陣仍單 owner；不能將小測 parity 當大型 scaling |
| GPU | 適用數值核心可用 CUDA | 單 GPU 貼體 standalone 支援子集 | CUDA FSI／ALE／通用多 domain／組織尚不能宣稱可用 |
| 狀態／重啟 | Trial、prepare、finalize、rollback、checkpoint | 多條路線已有；moving immersed FSI 有跨 rank restart 證據；native tetra ALE-FSI 有 checksum/identity 綁定的 single-partition paired restart parity | Native tetra ALE-FSI 尚缺 distributed shard ownership／repartition restart；library／test 與使用者 CLI 完整度不同 |
| I/O 與診斷 | VTK／ParaView、QoI、守恆、收斂、source/input hashes | 已有基礎與多個驗證報告 | 新物理仍要新增測試；動畫不是驗收證據 |

Moving／FSI runtime 使用分散式 PETSc rows 與 MPI collective；單一 heartbeat driver 的結果不代表多節點加速。重啟介面與限制見 [coupled restart](docs/COUPLED_RESTART.md)。

## 7. 代表性應用組合

| 使用情境 | 可組合的模型 | 目前差距 |
|---|---|---|
| 剛壁血管樹＋末端床 | 貼體 3D IGA＋1D 網路＋0D RCR | 核心已有；實際幾何／graph／誤差需求仍須驗證 |
| 可變形血管樹，保留 IGA 特色 | 貼體 FEM 或 IGA＋ALE，耦合 IGA 殼，再接 1D／0D | ALE、IGA 殼、matching tetra transfer 與 native solid lifecycle adapter 已有分離的原生 vertical slices；IGA 殼 map、共同案例與完整雙向 gate 尚待新增 |
| Surface 輸入的複雜固定血管 | Surface→FEM volume→流體；或現有 immersed IGA | 前者已有受限原生 P2/P1 水力工作流與雙出口 Y 功能例（Gmsh／fTetWild），但尚非任意器官、廣泛阻抗穩健或生理驗證；後者有幾何與成本限制 |
| 指定運動心室血流 | ALE 流體＋指定運動；或現有 prescribed immersed | 兩端共用 idealized LV surface／motion／可選回流端口，各完成 16-step 功能週期；端口流量、體積 RMS 及共點速度 L2 已量化，但空間／時間收斂、壁面差異分解與通用輸入仍待新增 |
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
- [理想化 prescribed LV 驗證](docs/T6_NATIVE_LV_QOI.md)
- [Moving／FSI restart 與限制](docs/COUPLED_RESTART.md)
- [數值基準與限制](docs/BENCHMARKS.md)

本文件的現況是上述核對日期的 snapshot。新增後端、修正守恆或改變可執行拓樸後，應同步更新本矩陣與對應的技術及驗證文件。
