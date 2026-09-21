# 肝臟首版幾何來源審查與單病例取得紀錄

初查日期：2026-09-19；單病例取得日期：2026-09-20 UTC；其他來源複核：2026-09-20。此文件記錄候選來源、原始影像／標籤及受限的單區域與多區域 ROI 試網格，**不代表已產生合格的全肝血管—組織多區域網格或灌流案例**。不可把理想化圓管、立方體或 Y 分岔稱為肝臟成果。原始醫學影像與衍生病人網格不納入程式庫。

| 來源 | 授權與取得 | 已明示幾何／標籤 | 首版判斷 |
|---|---|---|---|
| [TCIA Colorectal-Liver-Metastases](https://www.cancerimagingarchive.net/collection/colorectal-liver-metastases/) | 公開 [NBIA REST API](https://wiki.cancerimagingarchive.net/display/Public/NBIA+Search+REST+API+Guide) 可取得原始 DICOM/SEG；單病例的 series metadata 和 CT ZIP 內 `LICENSE` 均明示 CC BY 4.0 | [資料論文](https://www.nature.com/articles/s41597-024-02981-2) 說明同病例 liver、vessels、tumours、future liver remnant masks。已核實 `CRLM-CT-1046` 與 `CRLM-CT-1072` 各自同一 Study 的 CT 與 SEG，見下節 | 可作首版真實標籤幾何候選；L0 與局部互斥共面 tetra 的功能試驗完成，**全肝連通性／開口／BC／材料仍未驗收**，尚不得稱為灌流案例。 |
| [3D-IRCADb-01](https://www.ircad.fr/research-and-development/data-sets/liver-segmentation-3d-ircadb-01/) | 官方頁標示 CC BY-NC-ND 4.0；可下載個別病人資料，但衍生網格的分享權利需另行確認 | 官網明列 CT、labelled DICOM、mask DICOM、分區 VTK surface；頁面未列各病例完整血管 label 名稱 | 暫不採作可再散布的首版輸入；不能由「有肝表面」推論有合格灌流邊界。 |
| [Medical Segmentation Decathlon](https://registry.opendata.aws/msd/) | 官方 AWS 登錄標示 CC BY-SA 4.0、免 AWS 帳號的公開 S3 | Task03 Liver 與 Task08 Hepatic Vasculature 為不同任務；不能把兩者未配對的病例湊成單一個案 | 可用於各自幾何／標籤工具測試；尚不能證明同病例的肝體積＋血管灌流幾何。 |
| [LIRCAD／Inria Zenodo](https://zenodo.org/records/13897086) | 2026-09-20 複核：官方 Zenodo 記錄可公開瀏覽，但頁面 `Rights → License` 未顯示具體授權文字；13.6 GB ZIP 未下載、未核對原始檔內 LICENSE 或衍生權利 | 官方說明含 77 個病例的 CT、Portal／Hepatic vessel dual labels 及部分分支名稱（main／left／right 等），但未明列同病例肝實質 mask、血管 cap 或流體 BC；主分支名稱**不是**可直接求解的 inlet／outlet surface label | 分支語義可能有研究價值，但授權與肝組織／端口契約尚未核實；不取代已驗證的 TCIA 病例，也不下載或混接不同病例資料。 |
| [HVA-CT／Zenodo](https://zenodo.org/records/19850108)（[資料論文](https://www.nature.com/articles/s41597-026-08346-1)） | 2026-09-20 複核：官方 Zenodo 頁列出可下載檔名，但 `Rights → License` 欄未顯示授權條款；未下載或使用 | 論文及 Zenodo 說明為 MSD Task08 衍生的同病例 Portal／Hepatic vein 改良 mask、另有自動產生的 liver mask；未列出可求解的 capped inlet／outlet surface、病人 BC 或 Darcy 材料 | 仍是**分割研究資料候選**，不能以 open-access 論文代替資料授權，也不能僅憑血管／肝 mask 升級成灌流案例。 |
| [VSNet 重標資料／GitHub](https://github.com/XXYZB/VSNet) | 2026-09-20 複核：README 稱專案採 MIT 並提供重標 ZIP，但原始影像另指向 MSD；未確認原始影像及衍生 mask 的逐項權利，未下載 | README 明列 303 例 CT 的 Portal／Hepatic vein 標註，並未明列同病例肝實質 mask、flow caps、BC 或 Darcy 材料 | 可作日後血管分割候選，尚不足以替代目前已核實同病例授權／標籤的 TCIA 路線。 |

既有 `cases/liver_vessels_simple_implicit_pde/` 是 **1D OBJ network** 案例，設定檔內的入口 Fourier 流量、`0 Pa` 出口、黏度、密度及壁材料數字未附對應 TCIA 病例來源／量測依據；也沒有同病例肝組織 tetra、血管—組織匹配界面或 Darcy 來源映射。保留其程式與輸出作歷史 1D 功能測試，但不把名稱中的 `liver` 當作本首版 L1–L4 或生理 BC 證據，不將數字移植到 `CRLM-CT-1072`。

## 已取得的 TCIA 單病例：CRLM-CT-1046

取得時間：2026-09-20 UTC；來源：[TCIA collection](https://www.cancerimagingarchive.net/collection/colorectal-liver-metastases/) 的公開 [NBIA API](https://wiki.cancerimagingarchive.net/display/Public/NBIA+Search+REST+API+Guide)。授權：**CC BY 4.0**，由 CT 與 SEG 的 `getSeries` metadata 各自確認，CT ZIP `LICENSE` 亦明示；使用／衍生結果應註明來源及授權，且遵守檔內禁止再識別的使用條款。參照 [資料論文](https://www.nature.com/articles/s41597-024-02981-2)。本地原件只存於 `/tmp/tcia-crlm-1046-AIxqT5/`，未納入 git；`/tmp` 可清空，重跑須重新下載並比對 hash。

| 欄位 | CT | SEG |
|---|---|---|
| SeriesInstanceUID | `1.3.6.1.4.1.14519.5.2.1.9203.8273.876266207167921740530708709916` | `1.3.6.1.4.1.14519.5.2.1.9203.8273.198419865300741306933067795121` |
| 取得 API | `getImage?SeriesInstanceUID=<CT UID>`（ZIP） | `getSingleImage?SeriesInstanceUID=<SEG UID>&SOPInstanceUID=1.3.6.1.4.1.14519.5.2.1.9203.8273.207996444071145243031318450804` |
| 本次下載 SHA-256 | `7c6ed88e5a3014b3916c637f511010bd1491d8ff5561b9e431dafa86d83f31a8`（ZIP） | `e3db49bae8841517c3cbede6afdff21594ef1369e1837aaaee67db64dce2dc76`（DICOM） |
| 實際內容 | 47 張 512×512 CT DICOM、ZIP 10,916,177 bytes | 1 張 DICOM SEG、107 frames、3,574,118 bytes |

兩者 StudyInstanceUID 都是 `1.3.6.1.4.1.14519.5.2.1.9203.8273.333449295436763574579232065427`。SEG 明確參照上述 CT series 的 **31 個 SOPInstanceUID**，且 31 個均在下載的 47 張 CT 內；SEG 的 31 個切片位置全部與 CT 一致，CT 另有 16 張在分割範圍外。兩者都是 512×512、像素間距 `0.933594×0.933594 mm`；SEG z 座標為 `−400` 至 `−250 mm`，間距 `5 mm`，CT z 範圍為 `−455` 至 `−225 mm`。CT orientation 為 DICOM patient coordinates 的 row `(1,0,0)`、column `(0,1,0)`；CT 與 SEG in-plane origin 約 `(-247.399994,-239) mm`（SEG Y 欄位格式四捨五入）。轉換時必須保留 DICOM LPS 與 mm，不可默認為 m 或 RAS。

SEG 內的原始 SegmentLabel（**非**自動推論的邊界／生理語義）：

| number | SegmentLabel | 非零 voxel | 6-neighbour component 數 | 最大 component voxel |
|---:|---|---:|---:|---:|
| 1 | `Liver` | 373238 | 2 | 373219 |
| 2 | `Liver Remnant` | 223010 | 1 | 223010 |
| 3 | `Hepatic` | 5346 | 12 | 5235 |
| 4 | `Portal` | 6965 | 13 | 6790 |
| 5 | `Tumor_1` | 410 | 2 | 375 |
| 6 | `Tumor_2` | 3215 | 1 | 3215 |

Component 數以原始 SEG 的 31 個 z 平面按 DICOM 位置排列，使用 `scipy.ndimage.label` 預設 6-neighbour connectivity 計算；小孤島不能無記錄地刪除。`Hepatic`/`Portal` 的 mask 有資料論文中的血管語境，但 DICOM SegmentSequence 的 `SegmentedPropertyCategoryCodeSequence` 類別為 `Tissue`；只靠這些名稱尚**不能**定義流向、入口／出口、肝動脈、截面、壁面與血管—組織交換位置。5 mm 層厚也使細血管與開口品質有待檢查。不得自行添補任何缺失解剖或材料參數。

### WSL 幾何轉換試驗（L1 未通過）

新增 `scripts/dicom_seg_to_surface.py`，只接受規則、未翻轉的 axial DICOM LPS SEG；以原始 mask 的 `0.5` voxel isosurface 導出單一標籤的 STL。需本地 `pydicom`、NumPy、SciPy、VTK；測試時 `pydicom 3.0.1` 僅安裝在 `/tmp/tcia-crlm-1046-AIxqT5/pylib`。STL 不攜帶標籤，必須將 manifest 的 `stl_default_boundary_id` 傳給下游 preflight；這也**不是**完整多區域幾何。沒有**自動**平滑或補洞；連通分量選擇及平滑均須明確指定和受量化上限約束，不作解剖推斷。

```bash
PYTHONPATH=/tmp/tcia-crlm-1046-AIxqT5/pylib python3 scripts/dicom_seg_to_surface.py \
  /tmp/tcia-crlm-1046-AIxqT5/seg.dcm 1 /tmp/tcia-crlm-1046-AIxqT5/liver-repro.stl \
  --manifest /tmp/tcia-crlm-1046-AIxqT5/liver-repro.json \
  --keep-largest-component --smooth-iterations 10 --smooth-passband 0.1 \
  --max-displacement-mm 1.0 --max-volume-relative-error 0.005
solvers/cpu/surface_fem_preflight /tmp/tcia-crlm-1046-AIxqT5/liver-repro.stl \
  /tmp/tcia-crlm-1046-AIxqT5/liver-surface.msh \
  --default-boundary-id 1 --length-scale-to-m 0.001 --max-triangles 200000
```

原樣 `Liver` mask 產生 167,952 三角形；preflight **安全拒絕** `surface face adjacency is disconnected`，對應原始 mask 的 373,219＋19 voxel 兩個連通分量。只有明確加上 `--keep-largest-component` 才會排除 19 voxel（約原始 Liver mask 的 0.0051%；manifest 記錄原始分量與排除量）；未平滑的 167,840 三角形仍因自交被拒絕，診斷位置約 `(-0.1288, 0.0555, -0.3800) m`。`vtkFlyingEdges3D`、`vtkDiscreteFlyingEdges3D`、`vtkMarchingCubes` 的未平滑試驗均被拒絕；不能靠更換演算法免除檢查。

上方**明示的修復**使用 VTK windowed-sinc 10 iterations、passband `0.1`，對最大分量的表面位移中位數 `0.117 mm`、95 百分位 `0.326 mm`、最大 `0.973 mm`，皆由腳本寫入 manifest 並受 `1.0 mm` 上限約束；表面體積 `1,625,590.83 mm³` 與保留的 voxel 體積 `1,626,484.22 mm³` 相差 `0.0549%`，低於明定 `0.5%` 上限。產出的 STL SHA-256 為 `4041f6d71597426a749df1c10e7daeea720751b0fca960fd8deee8e66fa3fd7f`，與初次通過 preflight 的試驗輸出逐位元相同；preflight 通過，得 83,922 頂點、167,840 三角形，canonical SHA-256 `7de7f20e522d9b604a40a3c77795e061df64d0a14fea366a35ea7492c0dfa67c`。這個修復版**不是原始分割不變**，使用時必須保留兩者差異與來源。

用此表面跑 Gmsh `--target-size-m 0.005`、minimum scaled Jacobian `0.001`、maximum volume relative error `0.01` 時，Gmsh 仍有 127 個 ill-shaped tetra，實際最小 scaled Jacobian `2.05785e-05`，adapter **安全拒絕**，沒有發布肝 `.msh`。因此 L1 的完整標籤血管／組織 tetra 及其界面仍未完成；不可把單一肝表面 preflight 通過延伸為灌流幾何已通過。

同一表面用本機 fTetWild `d7d99bb4387a07895b9adce058dc7305f6b6e5ab`、`--target-size-m 0.005 --envelope-m 0.0005 --max-optimization-passes 20 --max-threads 4` 試跑兩次：**第一次** adapter 在輸出邊界重新投影標籤時因 `fTetWild boundary orientation disagrees with source surface` 安全拒絕，**第二次**相同設定通過，輸出 198,435 tetra、42,014 nodes、30,610 個 `boundary_label_1` triangles、最小 scaled Jacobian `0.0943184`、體積相對誤差 `0.0004618`、最大來源表面距離 `0.0003059 m`、輸出網格 SHA-256 `b12a232e8560db51bb26fd3dac24c58176179ac4ea95bee2dd8280736a093108`。第二次的 mesher peak RSS 約 `519 MB`、mesher wall `62.5 s`；兩次都另有約數分鐘的表面 preflight。結果顯示 **4-thread 路線不具穩定重跑證據**，不能把單次成功當作 L1 完成。另以 `--max-threads 1`、其餘設定相同重試，仍被方向 gate 拒絕；具體 facet `(10142,10111,10141)`、centroid 約 `(-0.003468,-0.097089,-0.339199) m`、最近來源距離 `0.0002728 m`、normal dot `−0.484`。故問題不只來自多執行緒排程；還需檢查局部幾何與投影規則，並在不放鬆方向、距離、標籤與 Jacobian gate 的前提下重試。adapter 的拒絕訊息已增加 facet／位置／距離／normal-dot 診斷；`scripts/tests/test_ftetwild_to_fem_volume.py` 通過。

另外，fTetWild 會在 envelope 內重網格，不能將其輸出 manifest 的 `current_geometry` 誤記為原始輸入表面。已修正後續輸出的 `geometry_contract`：來源表面 identity 存於 `source_geometry`，完成的 tetra 網格 hash 作固定 reference/current identity，stable IDs 指向輸出 MSH tags。**上列第二次試跑的舊 manifest 在此修正前產生，雖有正確 mesh hash 與 `geometry_change`，其 `geometry_contract` identity 欄位仍錯，不應直接作最終可重跑證據。**

該單表面 mesher adapter 將唯一體積標為 `fluid`，只是現有通用網格格式的預設，**不是**此病人的肝組織物理角色，更沒有分出血管與 tissue region。輸出 `.msh` 可讓原生 tetra reader 做幾何功能測試（本次 42,014 nodes／198,435 cells 讀取通過），不得直接當作肝 Darcy 或血管流體的案例輸入。L1 真正完成仍需同病例多區域 meshing、界面標籤、正 Jacobian 與穩定重跑。

### 更窄 envelope 的肝表面試驗（仍非 L1 完成）

對**同一**修復版 `Liver` 表面縮小 fTetWild envelope 至 `0.0001 m`，target size `0.005 m`、20 optimization passes、1 thread，保留所有原有的方向、距離、正 Jacobian、標籤與體積 gate。先直接用同一個已通過 preflight 的 canonical surface 啟動 fTetWild，結果經 adapter 的 `oriented_boundary`／`label_boundary` 檢查通過：482,506 tetra，最小 scaled Jacobian `0.07129`，最大來源距離 `0.00005967 m`、最小 normal dot `0.5806`。再用正式 `scripts/ftetwild_to_fem_volume.py` 從 STL 重跑，通過並發布 481,398 tetra、110,839 nodes、108,038 個 `boundary_label_1` triangles，最小 scaled Jacobian `0.07445`、體積相對誤差 `2.20×10⁻⁵`、最大來源距離 `0.00006031 m`、最小 normal dot `0.8411`；正式輸出 SHA-256 `891a9c2dd16b32161d6405213c76b6e6f9c3b9ec176012131de846036f7ed47d`，原生 reader 讀取通過。兩次體積相對差約 `8.53×10⁻⁷`，但 mesh hash／元素數不同，所以**只有幾何 QoI 的重跑一致性，沒有 bitwise 決定性**。正式 manifest 的 source、reference/current identity 已按重網格語義修正。這仍是單一肝外形體積，沒有組織／血管 region 和界面。

### 同病例血管標籤的表面可行性

原始 SEG 的 6-neighbour 最大分量：`Portal` 6,790／6,965 voxels（明示排除 175，約 2.51%），`Hepatic` 5,235／5,346 voxels（排除 111，約 2.08%）。未平滑的 `0.5` voxel isosurface 分別有 23,336／16,484 三角形；兩者各自通過專案的 closed/manifold/self-intersection preflight（canonical SHA-256 分別 `064c1145e9813cd0e2d1adc0dfe30abd9dc9785e2f741e395c063adf304184b9`、`73ab8bdd10ffcaff6c6549accb8be83f02a53e5ef9fc1e062e9fedec43afd1ae`）。但相對保留分量 voxel 體積的 isosurface 體積偏差達 **5.01%／4.09%**；這是 5 mm 層距下細血管幾何的重要限制，不能以表面 preflight 通過當作形態精度已驗證。上述數值只用於界定功能試驗，尚未訂出可接受的病人幾何誤差契約。

在共同 31-slice SEG 格點逐 voxel 比對，`Hepatic` 5,346 voxels 中 3,617 位於 `Liver` mask 內、1,729 在外；`Portal` 6,965 中 4,094 在內、2,871 在外，兩者彼此零重疊。最大分量也各有肝內及肝外部分。這證明這些標籤**不是**純肝內的封閉組織 region；肝表面與血管表面若分別 tetrahedralize，必須明確處理交叉、開口、組織挖除或非符合界面映射，不能直接把三個獨立網格拼接稱為守恆多區域網格，也不能由交叉位置推斷流向或生理 BC。

上述重疊數字現在可由專案 `scripts/audit_dicom_seg_region_overlap.py` 從**原始同一 SEG** 重跑；它先檢查 LPS 方向、共同平面格點、規則 z 間距、binary frame、來源 hash，再輸出每個原始標籤在 `Liver` mask 內／外的 voxel 數。`CRLM-CT-1046` 與 `CRLM-CT-1072` 均重跑吻合；缺標籤會拒絕（exit 2）。兩血管 mask 彼此均不重疊；純集合差 `Liver \ (Hepatic ∪ Portal)` 分別有 365,527／2,745,757 voxels，但**只是未網格化、未驗證的數學集合差**，不得宣稱已得到組織 Darcy region。新增逐 voxel-face 接觸稽核，只對互斥的 `Liver \ vessel_union` 與各血管 label 計算共面鄰接（不計邊／角接觸），並以原始 spacing 給出格點面積。`1046` 的 Hepatic／Portal 各有 6,788／9,440 接觸面（14,765.75／19,579.04 mm²）；`1072` 各有 1,876／33,862 面（1,034.50／19,045.27 mm²）。兩病例不同 spacing 與分割內容，這些**格點接觸面積不是平滑血管壁面積，也不能互相比作生理面積**。它只證明原始互斥 voxel 集合具有面鄰接，尚無經驗收的共同 tetra 界面、開口或流向；因此 L1/L3 仍未完成。

原始 SEG 的**六鄰接分量與 tissue 接觸覆蓋**亦已逐分量量化，並交叉檢查逐分量接觸面總數等於整體 mask 的接觸面數。`1046` 的 Hepatic／Portal 分別 12／13 個分量，均有 tissue-only 面接觸；扣除血管後的 tissue-only 有 2 個分量，最大 365,508 voxels、另一個 19 voxels。`1072` 的 Hepatic 為 1 個分量且有接觸；Portal **116 個分量**，其中 **108 個**有 tissue-only 面接觸、**8 個**無接觸（共 41 voxels）；最大 Portal 分量的 32,248 個接觸面佔整體 33,862 個，其餘分量不可默默併入它。`1072` 的 tissue-only 也有 2 個分量，最大 2,745,756 voxels、另一個 1 voxel。這是原始標籤的拓樸限制；局部 ROI 通過不代表全域血流域連通。全域求解前必須明示每個分量保留／排除的依據及相應 BC，或在缺乏依據時安全停用該模型，不能用一個壓力 gauge 或單一端口替全部不相連血管分量定義流向。

直接用現有多區域網格命令對兩份原始 SEG 執行 `--roi full` 預檢（`1046 --max-voxels 1000000`、`1072 --max-voxels 4000000`）時，兩者均以 exit 2 在 `region 'tissue_candidate' has 2 disconnected 6-neighbour ROI components` 停止，且未寫出 `.msh`。這是**預期的安全拒絕**，不是可透過放寬求解 tolerance 解決的問題；後續血管分量及 `1046` 的跨血管相鄰面也仍需獨立政策。沒有經來源支持的分量處置，不繞過此 gate 來宣稱 L1 全肝網格完成。

原生 C++ 求解器另加入獨立的**共三角面** tetra 分量 gate：流體與組織網格各自必須恰有一個 face-connected 分量，僅共點／共邊不算相連；重複節點、越界節點及超過兩個 owner 的面也安全拒絕。此檢查在設定 BC、組裝與 PETSc 求解前執行，避免直接呼叫 binary 時繞過 Python 前處理。合成 C++ 單元測試通過；新 binary SHA-256 `7ccde60436d6ced487ea0baa4758e52df6e015061130d2ac3f2f88b35ce4ce80` 在 `/tmp/tcia-crlm-1072-native-component-gate-1` 與 `-2` 各完成 1／2／4-rank 原始 SEG→局部 FEM 重跑，跨重跑數值與 Global IDs 解場比較通過。舊封存 field 證據的 binary hash 仍保留，不回寫為新 binary。這仍非全肝分量政策。

上述原生 gate 再補上**同一四節點集合的重複 tetra 與重複 cell ID** 拒絕，避免體積重複積分及 GlobalCellIds 撞號；合成負向測試通過。更新 binary SHA-256 為 `14ef3e9e97aeb7c956f916538cf84b99d145c6a9f2e0539bf2c4c8384e524041`，在 `/tmp/tcia-crlm-1072-native-duplicate-gate-1` 重跑 1／2／4-rank 局部案例通過；與前一 binary 的 1／2-rank QoI 完全相同，4-rank 差異只在約 `10⁻²³ m³/s` 的歸約浮點末位及約 `10⁻²⁵ m³/s` 的 cell defect，遠小於守恆 gate。這不表示 bitwise 跨 binary 相同，也不升級為全肝驗證。

現在的原始 SEG 稽核另列**每個**六鄰接分量的 `zyx` 半開 voxel bounding box、voxel 數、血管—tissue 共面數與是否碰到影像格點邊緣（**影像邊緣不是血管 port**）。`1072` 最大 Portal 分量 ID 1 為 47,283 voxels，範圍 `[[14,187],[163,300],[109,313]]`；第二個 tissue-only 分量 ID 2 是 `[[205,206],[220,221],[316,317]]` 的 1 voxel。這使小島可被精確審查，但未授權排除它們。指定局部 ROI `[[10,20],[210,226],[101,117]]` 的全域分量交集也已稽核：Portal 只含全域 ID 1 的 **32 voxels**、tissue-only 只含全域 ID 1 的 **2,188 voxels**，Hepatic 無 voxel。故這個 2,220-voxel 功能案例**只截取主分量極小部分**，不能由 ROI 求解外推全肝流量、開口或分量處置。來源 SEG hash 與 ROI audit JSON hash `866efff60a36aaafb9e1fc730a8706f3ffa0c0021c7b9bbba1b180edd2718dfa` 已進功能證據卡；稽核輸出拒絕覆寫。

新增 `scripts/audit_seg_exposed_vessel_faces.py` 將每個血管六鄰接分量的 voxel 表面按 `tissue-only`、影像內背景、另一血管標籤、SEG 格點邊緣四類分帳；每分量以 `6N−2E` 獨立核對總表面數，輸出不可覆寫且綁定 SEG／工具 SHA-256。`1072` Hepatic 的四類面數依序為 `1,876/3,048/0/0`，Portal 為 `33,862/6,376/0/0`；`1046` Hepatic 為 `6,788/1,625/50/171`，Portal 為 `9,440/2,746/50/0`。這些是**原始 voxel 拓樸分類**，不是 cap、可用入口／出口、血流方向或壓力邊界。特別是 `1046` 兩血管標籤有 50 個相鄰面，需避免在合併模型中重複或漏算。輸出僅留在 repo 外：`/tmp/tcia-crlm-1072-lXgLlB/exposed-faces-audit-final.json`（SHA-256 `bb0c2e9ead4e5c461a422374104e378d038db5348c91df0f595a7494da5615be`）及 `/tmp/tcia-crlm-1046-AIxqT5/exposed-faces-audit-final.json`（`73878a9fd6e0d7ad6bfb901fd1dfd037944351367f3e61cb832203fbcf096793`）。

現在 `scripts/run_liver_roi_from_raw.py` 在生成局部 tetra 前，也會從**同一份複製的原始 SEG** 執行全域血管暴露面稽核，將 `source-face-audit.json` 的 SHA-256 納入 `summary.json`，並在新的 fresh-series 驗證器中核對來源、標籤、分量數與四類面總數。新摘要為 schema v2，缺此面稽核會拒絕；舊 v1 證據仍可按原契約驗證。`/tmp/tcia-crlm-1072-source-face-bound-v2` 已完成 v2 工作流的 1／2／4-rank 重跑，與先前結果按數值及全域 ID 解場比較通過；新工作流不把全域暴露面誤稱為已定義的生理端口。

v2 對外 fresh-series 驗證另要求**完整摘要檔**、來源面稽核 hash、各階段時間／資源表，以及每個 native MPI rank 的資源紀錄同時存在；原始 SEG→結果腳本在寫摘要前先做內部數值驗證，寫入後再做一次完整 v2 驗證。缺稽核或 stage 表的負向單元測試通過，新的 `/tmp/tcia-crlm-1072-summary-gate-v2` 端到端重跑及與前版跨重跑比較通過。沒有摘要的中途失敗目錄不再能透過對外 fresh-series 完成驗證。

```bash
PYTHONPATH=/tmp/tcia-crlm-1046-AIxqT5/pylib python3 \
  scripts/audit_seg_exposed_vessel_faces.py \
  /tmp/tcia-crlm-1072-lXgLlB/seg.dcm 1 3 4 \
  --manifest /tmp/new-exposed-faces-audit.json
```

```bash
PYTHONPATH=/tmp/tcia-crlm-1046-AIxqT5/pylib python3 \
  scripts/audit_dicom_seg_region_overlap.py \
  /tmp/tcia-crlm-1072-lXgLlB/seg.dcm 1 3 4 \
  --roi 10:20,210:226,101:117 \
  --manifest /tmp/tcia-crlm-1072-lXgLlB/component-roi-membership.json
python3 scripts/validate_liver_roi_functional_evidence.py \
  --component-roi-audit /tmp/tcia-crlm-1072-lXgLlB/component-roi-membership.json
```

```bash
PYTHONPATH=/tmp/tcia-crlm-1046-AIxqT5/pylib python3 scripts/audit_dicom_seg_region_overlap.py \
  /tmp/tcia-crlm-1046-AIxqT5/seg.dcm 1 3 4 \
  --manifest /tmp/tcia-crlm-1046-AIxqT5/overlap-audit.json
```

`Portal` 最大分量以 fTetWild target `0.002 m`、envelope `0.0001 m`、1 thread 試網格時，adapter 因輸出 boundary 超出設定 envelope **安全拒絕**，未發布網格。這與單一肝表面較窄 envelope 成功不矛盾：血管形態及源格點解析度不同。後續嘗試更寬 envelope 時仍須獨立驗證邊界方向、來源距離、體積偏差和正 Jacobian，不以工具能產出 tetra 視為合格。

後續 `Portal` 試驗仍未過 gate：fTetWild target `0.002 m`、envelope `0.00025`／`0.0005 m` 的兩次輸出皆在邊界法向 gate 被拒絕；前者拒絕 facet 的 source distance 約 `0.0000751 m`、normal dot `−0.511`，後者約 `0.000146 m`、normal dot `−0.0553`。Gmsh target `0.002`／`0.001 m` 的最小 scaled Jacobian 分別 `0.0001548`／`0.0003121`，均低於不變的 `0.001` gate，因此沒有發布 `Portal` tetra。不能只放寬品質數字以稱為可用網格。

`Hepatic` 最大分量以 Gmsh target `0.002 m` 試網格則通過：31,360 tetra、9,575 nodes、`boundary_label_3` 的 16,484 原始三角 facets 全數保留，最小 determinant `6.57×10⁻¹¹ m³`、最小 scaled Jacobian `0.001352`、對已導出表面體積的相對誤差約 `5.53×10⁻¹⁴`；輸出 SHA-256 `24e6c1d7934f9ea33e88441b1605707b56649b6d793bf18cf5779d79e0e5e9ce`，原生 tetra reader 通過。Gmsh 仍報告 2 個 ill-shaped tets，但專案量化門檻確實通過；這是單一標籤血管分量的**幾何功能**證據，沒有開口分類、血流 BC、與肝組織的共同界面或生理驗證。它不使 L1 整體完成。

重取時使用上表 UID 與 API base `https://services.cancerimagingarchive.net/nbia-api/services/v1/`；解析 ZIP 前先檢查成員路徑，核對 `LICENSE`、下載 hash、Study/Series/SOP UID、影像座標與非零 SEG label。ZIP 的壓縮位元組若因服務端重打包而改變，應逐一核對 DICOM 內容和 UID，不把 ZIP hash 差異直接解讀成影像差異。

不將原始病人影像納入 repo；重取的最小命令（先建立 case 專屬暫存目錄，URL 完全指定同一病例的 series／SOP）：

```bash
case_dir=$(mktemp -d -p /tmp tcia-crlm-1046-XXXXXX)
curl -fsSL --output "$case_dir/seg.dcm" 'https://services.cancerimagingarchive.net/nbia-api/services/v1/getSingleImage?SeriesInstanceUID=1.3.6.1.4.1.14519.5.2.1.9203.8273.198419865300741306933067795121&SOPInstanceUID=1.3.6.1.4.1.14519.5.2.1.9203.8273.207996444071145243031318450804'
curl -fsSL --output "$case_dir/ct.zip" 'https://services.cancerimagingarchive.net/nbia-api/services/v1/getImage?SeriesInstanceUID=1.3.6.1.4.1.14519.5.2.1.9203.8273.876266207167921740530708709916'
sha256sum "$case_dir/seg.dcm" "$case_dir/ct.zip"
unzip -Z1 "$case_dir/ct.zip"  # 解壓前確認只有 LICENSE 與 00000001.dcm–00000047.dcm
```

解壓與處理前仍須比對上表 hash、metadata 與 `LICENSE`；不要對未經檢查的 ZIP 直接解壓。

## 較細切片的同病例替代候選：CRLM-CT-1072

於 2026-09-20 UTC 從同一 TCIA collection 的公開 API 取得去識別病例 `CRLM-CT-1072` 的 CT ZIP 與 SEG（只存 `/tmp/tcia-crlm-1072-lXgLlB/`）。CT series UID `1.3.6.1.4.1.14519.5.2.1.9203.8273.152783424215299005378585248107`，SEG series UID `1.3.6.1.4.1.14519.5.2.1.9203.8273.160655374625974558520205684891`，SEG SOP UID `1.3.6.1.4.1.14519.5.2.1.9203.8273.178516068093830030868820769150`，共同 Study UID `1.3.6.1.4.1.14519.5.2.1.9203.8273.186150960588321456990090403242`。CT 與 SEG 的 `getSeries` metadata 均明示 **CC BY 4.0**，CT ZIP `LICENSE` 再次確認；下載資訊應與上方 [TCIA collection](https://www.cancerimagingarchive.net/collection/colorectal-liver-metastases/) 和 [NBIA API](https://wiki.cancerimagingarchive.net/display/Public/NBIA+Search+REST+API+Guide) 一起保留。CT ZIP 55,543,830 bytes、SHA-256 `4ad0142f2eae3320d31a14dfe9270bd2c7e8f9d201f944c41f4ef332018cfaba`；SEG DICOM 21,928,022 bytes、SHA-256 `afe7a2bdf63cb65a5f4cea53cde8f2271bac51f5f53369b2c47a028b083c1654`。原件不納入 repo。

CT ZIP 含 240 張 CT 與 `LICENSE`；SEG 參照其中 217 張 SOP UID，**全部匹配**且同一 Study。CT 與 SEG 都是 512×512、平面 spacing `0.724609 mm`、切片位置步長約 `0.8 mm`，DICOM LPS axial orientation `(1,0,0;0,1,0)`，in-plane origin 約 `(-195.800003,-195.5) mm`；CT SliceThickness metadata 為 `1.25 mm`，與 `0.8 mm` 的相鄰 ImagePositionPatient 步長不同，轉換以**實際 frame 位置**建立座標，不把 thickness 當 z 步長。SEG 共 657 frames、217 個 z 位置，原始標籤 `Liver`、`Liver Remnant`、`Hepatic`、`Portal`、`Tumor_1`。這是較細格點，不等於血管解剖／開口精度已臨床驗證。

重取的最小命令；先核對 hash、ZIP 成員與 `LICENSE`，不將原件提交到 repo。若 ZIP 重打包導致 hash 改變，仍須逐個驗證 DICOM UID、影像座標與內容，不能僅依 ZIP hash 判斷影像變動。

```bash
case_dir=$(mktemp -d -p /tmp tcia-crlm-1072-XXXXXX)
curl -fsSL --output "$case_dir/seg.dcm" 'https://services.cancerimagingarchive.net/nbia-api/services/v1/getSingleImage?SeriesInstanceUID=1.3.6.1.4.1.14519.5.2.1.9203.8273.160655374625974558520205684891&SOPInstanceUID=1.3.6.1.4.1.14519.5.2.1.9203.8273.178516068093830030868820769150'
curl -fsSL --output "$case_dir/ct.zip" 'https://services.cancerimagingarchive.net/nbia-api/services/v1/getImage?SeriesInstanceUID=1.3.6.1.4.1.14519.5.2.1.9203.8273.152783424215299005378585248107'
sha256sum "$case_dir/seg.dcm" "$case_dir/ct.zip"
unzip -Z1 "$case_dir/ct.zip"  # 解壓前核對 LICENSE 和 240 個 CT DICOM 成員路徑
```

| label | 原始非零 voxel | 6-neighbour 分量 | 最大分量 voxel | 本地單標籤幾何試驗 |
|---|---:|---:|---:|---|
| `Liver` | 2,776,613 | 2 | 2,776,612 | 明示排除 1 voxel；409,188 三角形，未平滑版因自交拒絕；受限平滑版通過完整表面 preflight，Gmsh 單一 `Liver` mask 體積試網格 6,694,418 tetra、min scaled Jacobian `0.002007649`，原生 reader 通過；未扣除血管或建立組織 region。 |
| `Hepatic` | 7,972 | 1 | 7,972 | 原樣 9,852 三角形表面通過 preflight；Gmsh 24,769 tetra、min scaled Jacobian `0.02777`、正 determinant，原生 reader 通過。 |
| `Portal` | 47,823 | 116 | 47,283 | 明示排除 540 voxels；未平滑版因自交拒絕，10-iteration 受限平滑版 76,628 三角形通過 preflight；Gmsh 181,427 tetra、min scaled Jacobian `0.007410`、正 determinant，原生 reader 通過。 |

`Liver` 修復版最大點位移 `0.1542 mm`、對保留分量 voxel 體積差 `0.00896%`，表面 SHA-256 `49d4d35d0821f8a4ba6448ad657062a70c8d9bbefca35d4197e1ecf24e805b2d`，專案 preflight canonical SHA-256 `775f21661d54b9b85761021b356cb453b7883ebe6e7f81ab29baf23e28ed75e7`；此為**明示修復的幾何功能結果**，不是原始解剖不變或生理驗證。`Portal` 修復版最大點位移 `0.1686 mm`、對保留分量 voxel 體積差 `1.817%`；`Hepatic` 原樣等值面差 `0.743%`。正式 `Portal` Gmsh mesh 的 boundary label 4 保留 76,628 個來源三角 facet，mesh SHA-256 `99035c52abd901435c0d2546f463f4f80cc433bdc90428d8f5c4a9cb8d31084a`；`Hepatic` label 3 保留 9,852 facets，mesh SHA-256 `f79cdd04677201a8d2c16c7af70788736a9fb3a898cecfb77a4faa9755872252`。兩者皆只是**分開的單標籤血管分量**；沒有血管口分類、壁面—組織界面或灌流 BC。`Hepatic` 7,972 voxels 中 2,085 在 `Liver` mask 內、5,887 在外；`Portal` 47,823 中 28,771 在內、19,052 在外；不能不處理交叉就合併網格。

```bash
PYTHONPATH=/tmp/tcia-crlm-1046-AIxqT5/pylib python3 scripts/audit_dicom_seg_region_overlap.py \
  /tmp/tcia-crlm-1072-lXgLlB/seg.dcm 1 3 4 \
  --manifest /tmp/tcia-crlm-1072-lXgLlB/overlap-audit.json
```

對 `Liver` 修復表面最初嘗試 Gmsh target `0.005 m`、來源單位轉換 `0.001 m/mm`、boundary label 1、原有最小 scaled Jacobian `0.001` gate；整個 adapter 以 `900 s` 限時。adapter 重新執行上述 409,188-triangle 完整 preflight，至時限仍未開始產出可驗收 tetra，程序以 `124` 超時結束，未發布 `.msh`／manifest。**該次**網格結果未定（inconclusive），不是幾何不合格；後續已用保留安全檢查的 hash-verified preflight 重用路徑解決，見下段。

此重複 preflight 成本已由帶 hash 驗證的重用路徑解決：更新後的專案 C++ preflight 重新執行通過，其 canonical MSH SHA-256 `3dc0d73146885ed595e98b893194f4b2b55afa0ebf3ba6e56bd05d2ee6c68c7f` 與上次獨立輸出**完全相同**。Gmsh adapter 核對原始 STL SHA、canonical MSH SHA、選項和所有幾何通過旗標後，以相同 target `0.005 m` 成功產生 `6,694,418` tetra、`1,170,873` nodes；最小 determinant `1.243079282536043×10⁻¹² m³`、最小 scaled Jacobian `0.002007649`（門檻 `0.001`）、409,188 個 boundary label 1 facet 全數保留，tetra 邊界集合精確吻合來源面，對封閉表面體積相對差 `2.01×10⁻¹⁴`。輸出 MSH SHA-256 `9c7ece1efc0efebd2fc43fd6bdb3d96100c5d7a4356853dcfeb07dc4df57ada9`，單次重用／meshing／專案品質稽核約 `329 s`、峰值 RSS 約 `6.99 GB`。專案原生 tetra reader／P2 拓樸建構在 WSL 用此網格讀取通過（`8,069,884` edges、約 `59 s`、峰值 RSS 約 `1.89 GB`），**僅為讀入功能測試，沒有求解 PDE**。Gmsh 曾警告 2 個 ill-shaped tets，但實際所有元素仍通過專案正 Jacobian 與 `0.001` scaled-Jacobian gate；不把警告當成驗收，也沒有調低門檻。這**只**是 `Liver` mask 的單區域幾何網格；Gmsh adapter 的唯一 volume physical name 預設為 `fluid`，在此**不是**組織物理角色，且原始 mask 與血管重疊。沒有互斥分區、共同界面、開口或 Darcy 案例。

同一 STL、同一帶 hash 的 preflight pair、同一 Gmsh 設定與本機環境再次完整重跑，第二次也是 `6,694,418` tetra、`1,170,873` nodes、相同最小 determinant／scaled Jacobian／體積與 boundary facet 數，輸出 MSH SHA-256 **完全相同** `9c7ece1efc0efebd2fc43fd6bdb3d96100c5d7a4356853dcfeb07dc4df57ada9`（第二次約 `329.3 s`）。這證明本機此環境下的**單區域幾何重跑一致**，不證明跨平台 bitwise 決定性或肝臟灌流模型可重現。

```bash
solvers/cpu/surface_fem_preflight /tmp/tcia-crlm-1072-lXgLlB/liver-smoothed.stl \
  /tmp/tcia-crlm-1072-lXgLlB/liver-validated.msh \
  --manifest /tmp/tcia-crlm-1072-lXgLlB/liver-validated.json \
  --default-boundary-id 1 --length-scale-to-m 0.001 --max-triangles 450000
python3 scripts/surface_to_fem_volume.py \
  /tmp/tcia-crlm-1072-lXgLlB/liver-smoothed.stl \
  /tmp/tcia-crlm-1072-lXgLlB/liver-cached-gmsh.msh \
  --manifest /tmp/tcia-crlm-1072-lXgLlB/liver-cached-gmsh.json \
  --target-size-m 0.005 --length-scale-to-m 0.001 \
  --default-boundary-id 1 --max-triangles 450000 \
  --validated-surface-msh /tmp/tcia-crlm-1072-lXgLlB/liver-validated.msh \
  --validated-surface-manifest /tmp/tcia-crlm-1072-lXgLlB/liver-validated.json
```

### 原始 SEG 的互斥多區域 tetra ROI（幾何功能驗證）

新增 `scripts/dicom_seg_to_multiregion_tet.py`，直接從同一原始 SEG 的 `Liver`、`Hepatic`、`Portal` masks 建立 `Liver \ (Hepatic ∪ Portal)` 的 `tissue_candidate` 與原始血管標籤互斥分區；逐 voxel 用固定六 tetra 分割、共享同一格點節點，並明列 `roi_cut_not_a_flow_port`、`unclassified_exterior_not_a_flow_port` 和 tissue–vessel 精確共面界面。這是**原始 voxel-faithful 階梯幾何**，不使用獨立平滑表面，不得宣稱血管細節已達解剖精度或切面是入口／出口。血管 masks 若重疊、ROI 任一存在區域非六鄰接單分量、血管彼此共面而無明確契約、區域表面有非 manifold edge／vertex、無 tissue–vessel 共面、超過指定體素上限或已存在輸出會安全拒絕。發布 manifest 前由多區域 tetra audit 驗證正 Jacobian、min scaled Jacobian `0.001`、face-connected region、所有內外面分區、來源 voxel 體積和每個 voxel 接觸面對應的兩個 interface 三角形。

在 `CRLM-CT-1072` 的明示 half-open ROI `z=10:20, y=210:226, x=101:117`，由原始 SEG 的 2,220 個標記體素產出 2,811 節點、13,320 tetra（tissue-candidate 13,128，Portal 192）、158 個 Portal–tissue 共面 interface 三角形；最小 determinant `4.200465623×10⁻¹⁰ m³`、最小 scaled Jacobian `0.37417548`，其餘外表面按 ROI 裁切與未標籤 SEG 外部分類。各 region tetra 體積與原始 voxel 數及 spacing 在浮點誤差內吻合；兩次相同來源／ROI／工具在本機產出相同 mesh SHA-256 `75201bd15c126d5fd8ab71162e4f7f8b242972cc65f502dfb431951ee02acd94`。來源 SEG SHA-256 為 `afe7a2bdf63cb65a5f4cea53cde8f2271bac51f5f53369b2c47a028b083c1654`。`--max-voxels 100` 對此 ROI 以 exit 2 拒絕，合成兩 voxel 案例通過確定性、品質、分區及 nonmanifold／重疊拒絕測試。

```bash
PYTHONPATH=/tmp/tcia-crlm-1046-AIxqT5/pylib python3 \
  scripts/dicom_seg_to_multiregion_tet.py \
  /tmp/tcia-crlm-1072-lXgLlB/seg.dcm 1 3 4 \
  /tmp/tcia-crlm-1072-lXgLlB/roi-portal-connected.msh \
  --roi 10:20,210:226,101:117 --max-voxels 50000
python3 scripts/audit_multiregion_tet_mesh.py \
  /tmp/tcia-crlm-1072-lXgLlB/roi-portal-connected.msh \
  /tmp/tcia-crlm-1072-lXgLlB/roi-portal-connected.contract.json
```

同一個已稽核多區域 `.msh` 現可由 `scripts/split_multiregion_tet_for_native.py` 抽成兩份原生 reader 可讀的單區域網格；它重做原網格正 Jacobian／完整面分區稽核，為 interface 兩側指定相同明確數字邊界 label，並驗證兩側三角面實體座標集合完全一致，輸出前不推斷入口／出口。此 ROI 的 tissue-candidate 子網格有 13,128 tetra、2,807 節點、label 4 上 158 面，SHA-256 `9af3cc058de3de3df0fa0ea2f0d848146765c25dc2fbc3da1698e0bc39fcf1f5`；Portal 子網格有 192 tetra、88 節點、同 label 4 上 158 面，SHA-256 `6cec8db10f064c01b13d02e2478da2f655e0f86417b33b640ae663740a573a59`。兩次分割得到相同 hash；兩份均由專案 C++ `ReadNativeTetMeshGmsh41`、P2 topology 與正 Jacobian 檢查通過。子網格唯一 volume physical name 為原生 reader 所需的 `fluid`，**這不是 tissue 的物理角色**；角色仍以來源 contract／split manifest 為準。

用專案 `MapNativeTetMatchingInterfaceToTissueSource` 的幾何／守恆路徑，給 Portal 子網格**未求解的製造仿射場** `u_x=(x-x_min)/(0.01 m) m/s`（其餘分量為零），實際 158 個界面面核對精確座標、唯一 cell owner、相反法向後，vessel-outward 及沉積 tissue-cell 的流量均為 `1.34414899937535×10⁻⁶ m³/s`；錯誤 interface label 以 exit 2 拒絕。這個場帶非零散度、不是血流 PDE 解，數字完全不能解讀為病人流量。

同一 manufactured source 現亦已送入專案自主 P1 Darcy／RT0 求解器，在 tissue 子網格 `boundary_label_2`（ROI cut）和 `boundary_label_3`（未分類 SEG exterior）**僅為功能測試**指定 `0 Pa`，interface label 4 保持自然零面通量，明示功能材料 `K=1 m²/(Pa s)`；不是肝組織實測 mobility 或生理壓力。WSL 1／2／4 MPI ranks 均通過：界面來源 `1.3441489993753505×10⁻⁶ m³/s`，RT0 保守外流約 `1.34414899937515×10⁻⁶ m³/s`，最大逐 cell 收支缺陷約 `1.064×10⁻²⁰ m³/s`，Darcy KSP converged reason `3`。此功能路線已到「真實標籤局部幾何→製造速度場→逐面來源→自主 Darcy→保守外流」，但**沒有解血管流體方程，不能冒稱患者灌流案例或 L4 端到端完成**。錯誤 interface label 也安全拒絕。

進一步的**實際求解雙物理功能案例**使用相同 ROI 子網格，從命令列明示 label 1（此 ROI 的人工截面）給定 `u=(0,0,-0.01) m/s`、密度 `1000 kg/m³`、黏度 `0.004 Pa·s`，在 label 4（血管—組織界面）給定人工 `0 Pa` 作血管壓力 outlet；以專案自主 P2/P1 固定域穩態流體 FEM 解速度／壓力後，逐面轉移該**已求解**速度的通量至組織 Darcy cell。組織同樣使用只供功能驗收的 `K=1 m²/(Pa s)` 和 label 2/3 的 `0 Pa` 人工外邊界，RT0 還原保守通量。1／2／4 ranks 均通過：血管入口約 `−1.5751746086430047×10⁻⁸ m³/s`，界面來源 `+1.5751746086430044×10⁻⁸ m³/s`，Darcy 保守外流約 `1.5751746086512×10⁻⁸ m³/s`，最大逐 cell defect 約 `5.020×10⁻²² m³/s`；流體 Newton 4 次更新、最終線性收斂原因 `4`，Darcy 線性收斂原因 `3`。錯 interface label 以 exit 2 拒絕。此案例**確有自主血流 PDE 與 Darcy PDE 的數值解及守恆交換**，但入口與壓力／材料值純為明示的人工功能輸入，ROI 不是全肝，亦非量測灌流／生理驗證。

`scripts/run_native_matching_roi_functional.py` 在啟動 MPI 前逐一比對原始 SEG、互斥多區域 mesh／contract、兩份子網格與求解器 binary 的 SHA-256，重新執行多區域 Jacobian／面分區 audit，再以顯式 SI 輸入生成每個 rank count 的 JSON（含數值 QoI、收斂原因、wall time 與 MPI launcher 的子程序 RSS）。來源 contract 不匹配時實測在求解前拒絕；既有結果檔也拒絕覆寫。可明示 `--field-output-dir`，以專案 rank-owned VTU/PVTU 輸出 fluid P2 速度／壓力與 tissue P1 壓力、Darcy RT0 通量、逐面保守流量和體源，wrapper 逐檔 hash 並拒絕重複 cell ID／不完整覆蓋。1／2／4-rank JSON 和 VTU 已留於病例專屬 `/tmp`；摘要證據卡記錄來源與 field bundle hash，驗證器按 Global IDs 重組比較跨 rank 解場。wall time／RSS 和浮點輸出不宣稱 bitwise 相同；此 rank-owned **局部功能輸出**仍不使全肝 L4 完成。

相同輸出目錄另有 `interface_facets.json`：專案 mapping kernel 對每個匹配面輸出兩側 cell ID、三角形座標和 `vessel_outward_m3_s`。wrapper 從兩份原始子網格獨立核對 158 個幾何面完整且不重複、每面兩側 owner 存在、逐面有號流量總和等於入口／界面守恆來源，並逐 tissue cell 檢查 `ΣQ_f = source_s_inv × volume_m3`；證據驗證器比較 1／2／4 ranks 的逐面幾何、owner 與有號流量。在這個**人工條件**案例有 64 面正向、94 面反向，淨流量仍約 `+1.5751746086430×10⁻⁸ m³/s`。這是功能測試中的局部回流數值，不可解讀為病人血流或組織灌流方向。

```bash
python3 scripts/split_multiregion_tet_for_native.py \
  /tmp/tcia-crlm-1072-lXgLlB/roi-portal-connected.msh \
  /tmp/tcia-crlm-1072-lXgLlB/roi-portal-connected.contract.json \
  /tmp/tcia-crlm-1072-lXgLlB/roi-native-connected-split
make -C solvers/cpu native_tet_mesh_inspect native_tet_matching_mesh_inspect
solvers/cpu/native_tet_mesh_inspect \
  /tmp/tcia-crlm-1072-lXgLlB/roi-native-connected-split/seg_4_portal.msh
solvers/cpu/native_tet_matching_mesh_inspect \
  /tmp/tcia-crlm-1072-lXgLlB/roi-native-connected-split/seg_4_portal.msh \
  /tmp/tcia-crlm-1072-lXgLlB/roi-native-connected-split/tissue_candidate.msh 4
make -C solvers/cpu native_tet_matching_darcy_smoke PETSC_DIR=/usr/lib/petsc
mpiexec -np 2 solvers/cpu/native_tet_matching_darcy_smoke \
  /tmp/tcia-crlm-1072-lXgLlB/roi-native-connected-split/seg_4_portal.msh \
  /tmp/tcia-crlm-1072-lXgLlB/roi-native-connected-split/tissue_candidate.msh \
  4 2 3 1
make -C solvers/cpu native_tet_matching_solved_roi_smoke PETSC_DIR=/usr/lib/petsc
python3 scripts/run_native_matching_roi_functional.py \
  --source-seg /tmp/tcia-crlm-1072-lXgLlB/seg.dcm \
  --geometry-manifest /tmp/tcia-crlm-1072-lXgLlB/roi-portal-connected.manifest.json \
  --multiregion-mesh /tmp/tcia-crlm-1072-lXgLlB/roi-portal-connected.msh \
  --multiregion-contract /tmp/tcia-crlm-1072-lXgLlB/roi-portal-connected.contract.json \
  --split-manifest /tmp/tcia-crlm-1072-lXgLlB/roi-native-connected-split/manifest.json \
  --vessel-region seg_4_portal --tissue-region tissue_candidate \
  --interface-name interface_seg_4_portal_tissue_candidate \
  --inlet-label 1 --tissue-exit-labels 2 3 \
  --inlet-velocity-m-s 0 0 -0.01 --density-kg-m3 1000 \
  --viscosity-pa-s 0.004 --mobility-m2-pa-s 1 \
  --interface-pressure-pa 0 --tissue-exit-pressure-pa 0 \
  --ranks 2 \
  --field-output-dir /tmp/tcia-crlm-1072-lXgLlB/solved-fields-direct-r2 \
  --output /tmp/tcia-crlm-1072-lXgLlB/solved-functional-direct-r2.json
python3 scripts/validate_liver_roi_functional_evidence.py \
  --results-dir /tmp/tcia-crlm-1072-lXgLlB \
  --field-results-dir /tmp/tcia-crlm-1072-lXgLlB
```

最後一行核對**這次留存的**原始碼版本、1／2／4-rank JSON 與解場 hash／ownership；完整重跑時須以各自 `--ranks` 產生三份輸出，使用新的病例專屬目錄，並以新證據卡記錄當次 hash，不能期望不同 rank 的浮點 field 位元相同。上述路徑示範已存在輸出時，生成器與 wrapper 會拒絕覆寫。可攜的摘要證據另存於 `benchmarks/liver_roi_functional_evidence.json`，不包含可取代原始病例輸入的解場。

同一人工功能輸入另凍結為 `cases/liver_roi_functional.json`。它只含已稽核來源的相對路徑／hash、ROI／label 與**明示人工** SI 邊界／材料數值，不附病人 SEG，也不宣稱病人 BC。`scripts/run_liver_roi_functional_case.py` 先逐檔核對，再呼叫上述求解與守恆 wrapper；2-rank 實測通過，案例 SHA-256 `0027687c0d8f952f459a37496b3e1d3878a99da9d6af0f3349355eab877c939a` 寫入結果。錯資料目錄在 MPI 前以 exit 2 拒絕。案例檔或任何輸入變更均需新證據，不能套用本次通過卡。

`scripts/run_liver_roi_from_raw.py` 現可從**原始同病例 SEG** 在全新、repo 外目錄直接重做 ROI 多區域 tetra、原生子網格、案例 hash gate、1／2／4-rank 兩個自主 FEM 解、逐面來源與 rank-owned VTU。`--output-dir` 已存在或原始 SEG hash 錯時在寫出前 exit 2 拒絕；中途失敗的目錄會保留供診斷，不假裝已完成。兩次完整成功重跑 `/tmp/tcia-crlm-1072-full-rerun-2` 與 `-3` 均由 2,220 標記 voxels 得 13,320 tetra，兩次 `summary.json` SHA-256 同為 `64a7121d0b71114cb4d7f4b73e53669d0167d65fc553bf8dbf5600e507a87ba1`。非封存式 `--fresh-series-dir` 驗證兩次各自的來源／案例／求解器、1／2／4-rank 解場與面流量，再跨重跑按 Global IDs 和數值容差比較；通過但**只限局部人工功能案例**。需本地 pydicom／SciPy／Gmsh 和 PETSc／MPI，並遵守不把病人衍生網格放進 git 的要求。

後續同案例新增每階段 GNU time 紀錄（`summary.json` 的 `stage_metrics`），本機一次受限 ROI 重跑的 mesh／split／1／2／4-rank solve wall time 分別約 `2.98/0.41/0.77/0.68/0.67 s`；mesh／split 最大 RSS 約 `1.37 GB/97.8 MB`。solve 階段約 `47–49 MB` 是**外層 Python／MPI launcher 命令的 GNU time RSS**，不是所有 MPI rank 的 RSS 總和，更不是 GPU allocation 或大型全肝效能測量。耗時與記憶體不可做 bitwise 重跑 gate；網格 hash、守恆量和全域 ID 解場才是再現性 gate。

另在 `/tmp/tcia-crlm-1072-rank-resource-2` 的新重跑中，`rankN/result.json` 以 `rank_solver_resources` 記錄每個**原生求解器 MPI rank 程序**的 GNU time wall 與 peak RSS，並以 SHA-256 核對不可覆寫的 `fields.rank-metrics/rankN.time.txt`。1-rank 為 `47,132 KB`；2-rank 各為 `48,620/48,156 KB`；4-rank 各為 `46,144/45,940/46,172/46,368 KB`。這是個別 rank 的最高 RSS，**不能相加當同時總量**。新舊兩種結果的跨 rank／跨重跑數值與場量比較通過；此量測只補局部人工功能案例的資源證據，不代表全肝效能或生理驗證。

```bash
PYTHONPATH=/tmp/tcia-crlm-1046-AIxqT5/pylib python3 \
  scripts/run_liver_roi_from_raw.py \
  --source-seg /tmp/tcia-crlm-1072-lXgLlB/seg.dcm \
  --output-dir /tmp/my-new-crlm-1072-run-A
PYTHONPATH=/tmp/tcia-crlm-1046-AIxqT5/pylib python3 \
  scripts/run_liver_roi_from_raw.py \
  --source-seg /tmp/tcia-crlm-1072-lXgLlB/seg.dcm \
  --output-dir /tmp/my-new-crlm-1072-run-B
python3 scripts/validate_liver_roi_functional_evidence.py \
  --fresh-series-dir /tmp/my-new-crlm-1072-run-A \
  --compare-fresh-series-dir /tmp/my-new-crlm-1072-run-B
```

```bash
python3 scripts/run_liver_roi_functional_case.py \
  --data-dir /tmp/tcia-crlm-1072-lXgLlB --check-only
python3 scripts/run_liver_roi_functional_case.py \
  --data-dir /tmp/tcia-crlm-1072-lXgLlB --ranks 2 \
  --field-output-dir /tmp/tcia-crlm-1072-lXgLlB/solved-fields-ledger-final-r2 \
  --output /tmp/tcia-crlm-1072-lXgLlB/solved-functional-ledger-final-r2.json
python3 scripts/validate_liver_roi_functional_evidence.py \
  --case-result /tmp/tcia-crlm-1072-lXgLlB/solved-functional-ledger-final-r2.json
```

這只建立**局部真實標籤上的功能耦合證據**。ROI 的兩個區域現在均先經六鄰接單分量 gate，輸出 contract 設 `require_connected=true` 並由 tetra audit 再驗證 face-connected cell graph；重跑 mesh hash 不變。這不代表**全肝**互斥多區域網格已連通：原始全域 Portal 有 116 個分量，未明示分量選擇前應被此 gate 拒絕。血管切面沒有被判定為真正入口／出口，且沒有病人 BC 或材料；上段雖實際求解兩個 PDE，所有外部數值均為人工功能測試。不可將上述 ROI 解讀成 L1 全域完成、L3 病人生理耦合或肝臟灌流驗證。全解析度 `1072` 的肝臟相關約 280 萬 voxels 會生成約 1,680 萬 tetra，需另設可承受的記憶體／輸出與全域 manifold gate；目前 CLI 預設 50,000-voxel 上限，不能默默放大。

實測 `--roi full --max-voxels 3000000`（明示提高上限、只驗證 fail-closed）在**產生任何 mesh 前**以 exit 2 拒絕：`region 'tissue_candidate' has 2 disconnected 6-neighbour ROI components`，沒有寫出 `.msh`。因此當前 CLI 不會把全域 116 個 Portal 分量或 tissue-only 小島默默包成一個看似可求解的域；需先有記錄在案的分量處置契約，然後才可重新評估全肝記憶體、manifold 與 BC。

## 後續 fail-closed gate

1. 已記錄本次 collection、下載日期、原始檔 hash、授權與引用、病人／series 去識別 ID、DICOM SEG 每個 segment 的原始名稱及參考影像 UID；重跑時再驗證來源版本和 hash。
2. 使用同病例的肝臟與血管 mask；已核對座標、體素單位、CT/SEG 配準、原始分量數與局部 ROI 的正 Jacobian／共面標籤。**全肝**互斥分區仍須先決定不連通小分量的保留／排除政策，再驗收全域 manifold、正 Jacobian、標籤保留與可用開口；不能把局部 ROI gate 外推為全域通過。
3. 僅在可從資料與明示操作辨識入口、出口、血管—組織區域對應時建立多尺度案例。若 vessel mask 只有未分類的合併標籤，就只驗證幾何／體積網格，不推斷門靜脈、肝動脈或肝靜脈。
4. 流體黏度、Darcy mobility、壓力／流量、交換係數均由案例檔明示來源與量綱；缺值就停止，不用器官名稱補預設值。理想化參數可用於功能測試，但不可稱為生理驗證。

現狀：**L0 的授權、同病例影像／標籤與座標核實已完成；L1–L4 未完成。** 在可用開口、界面、分區網格與明示 BC／材料參數取得前，肝臟專屬灌流案例仍暫停。可先用這一個案測試真實標籤 surface／tetra 幾何流程，同時繼續通用血管—組織耦合及解析測例；均不可冒稱肝臟灌流或生理驗證。PSC 大型計算不以 WSL 測試代替。
