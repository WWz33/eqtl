# Changelog

## [Unreleased]

### Fixed
- **`nb_irls` 收敛标志** — 外层 φ 轮询实际检测到收敛失败时正确报告 `glm_converged = 0`。
- **协变量共线性**：`prep_lm` 前新增重构一致性检查（`‖XtX·XtX⁻¹ − I‖∞ < 1e-6`)，共线时将该基因标记为不可测（`n_tested=0, p=NaN`)；若整个数据集共线则 fail-fast。
- **主输出缓冲**:`pubsetbuf` 在 `open()` 之后调用对 libstdc++ 无效。改为 set buf → open,pairs/top/region 实际生效 4MB 缓冲（此前仅后两个生效，pairs 走了默认 8KB)。
- **命令行解析严格**:`--window`、`--perm` 等参数用 `strtol` 全程校验；`--window 1e6` 将不再静默截断为 1。
- **scan_cis**：加载共享 GRM 基谱后，置换循环内的 `basis_ref`/`K_ref` 正确透传（此前空指针）；post-error 清理 tmp 文件。
- **`fission`**:`write_tsv` 精度提升为 `setprecision(17)`,Y1/Y2 可直接回读为 eqtl `-e` 输入。

### 性能（LMM）
- **cis 置换**：谱空间缓存 `g_til`（每 SNP 只投影一次 `Qᵀg`，不再每个置换重投影）、复用 `chi0` 与 `X_til`、置换内只对 |t| 最大的 SNP 求一次不完全 beta、黄金分割每步取一个内点。`--thread 4`、200×1000 测试面板实测：`--perm 100` 117.9 s → 12.5 s，`--perm 20` 26.3 s → 5.3 s，缺失面板 `--perm 20` 30.0 s → 6.3 s；除黄金分割一项外均逐字节相同。
- **trans 缺失面板**：mixed-keep 分支恢复并行、按 keep 分组共享 `Qᵀg`、GRM 基谱按 keep 记忆化（逐字节相同）。实测 `--perm 0`：missH 230 s → 22 s，变长面板 230 s → 48 s。

### 运行接口变更
- **trans 置换 p 停用**：`--perm` 在 trans/gw 下不再触发 stage-2，`p_emp`/`p_beta`/β 形参输出 NA（运行时 warn 一次）。原因：其 top-K 置换 null 按观测 p 选取，观测端带选择偏倚，`p_emp` 对 ~97% 基因恒等于 1/(B+1)（简并而非保守）。cis 置换不受影响。分析见 `eqtl-gene-level-p-plan.md`。
- **删除 `--perm-trans-thr` / `--perm-trans-top`**：stage-2 停用后已无任何代码读取，属于死参数；传入会报未知选项。README 表格同步移除。
- **`--disable-beta-approx` → `--perm-mode {beta,exact}`**：同一个概念（基因级置换 p 的算法）一个旋钮；`exact` 等价于原开关（不算 β 拟合，`p_beta`/β 形参为 NA），`beta` 为默认，行为与改前完全一致。
- **删除 `--fast`**：一个开关装了三种不同语义，且没有一种站得住——LMM 半边是对 GRM 做稀疏近似（实测零加速、结果改变），GLM 半边只写了一个从未被读取的字段（无效果），GLMM 半边冻结 σ² 属可省步骤。默认路径（不带该开关）行为完全不变；稀疏 GRM 与 basis 缓存里为此设的键均一并移除。
- `--perm-freeze-delta`（默认关闭）：LMM 置换链路复用观测 δ、跳过每次抽样的 REML 搜索。cis 置换快 14–30%、trans stage-2 快 23%；`p_emp`/`p_beta` 有轻微漂移（nominal 结果不变），既有输出默认不受影响。
- `EQTL_PROF=1`：退出时打印 REML 求值次数与耗时（排查用，默认关闭，无行为影响）。
- `--model glm` / `glmm` 现在自动拒绝 `--mode trans/gw`(count 模型只有 cis 路径，否则会按数百基因 × 500k SNP PQL 重拟合）。
- `make test` 变绿色。

### 数值准确性
- LMM δ 的黄金分割单点求值改变了内点复用的低位：δ 有 ULP 级漂移，`acat_p` 相对漂移 ≤ 4.2e-5（cis）/ 9.4e-4（trans），`p_emp`/`n_sig`/`n_tested` 不变。
- `p_from_t`、`beta_cdf`、`pnorm_two_sided` 与 scipy 一致（相对能差 < 1e-12)。
- bed 2-bit 编码解码与 PLINK2 / tensorQTL / GCTA 完全一致（位操作 + 正读双跑）。
- GRM 计算与 `gcta64 --make-grm` 匹配（sub-divider 正确包含所有 non-missing markers)。

## [0.1.0] — 首次提交
