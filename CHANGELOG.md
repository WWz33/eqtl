# Changelog

## [Unreleased]

### Fixed
- **`nb_irls` 收敛标志** — 外层 φ 轮询实际检测到收敛失败时正确报告 `glm_converged = 0`。
- **协变量共线性**：`prep_lm` 前新增重构一致性检查（`‖XtX·XtX⁻¹ − I‖∞ < 1e-6`)，共线时将该基因标记为不可测（`n_tested=0, p=NaN`)；若整个数据集共线则 fail-fast。
- **主输出缓冲**:`pubsetbuf` 在 `open()` 之后调用对 libstdc++ 无效。改为 set buf → open,pairs/top/region 实际生效 4MB 缓冲（此前仅后两个生效，pairs 走了默认 8KB)。
- **命令行解析严格**:`--window`、`--perm` 等参数用 `strtol` 全程校验；`--window 1e6` 将不再静默截断为 1。
- **scan_cis**：加载共享 GRM 基谱后，置换循环内的 `basis_ref`/`K_ref` 正确透传（此前空指针）；post-error 清理 tmp 文件。
- **`fission`**:`write_tsv` 精度提升为 `setprecision(17)`,Y1/Y2 可直接回读为 eqtl `-e` 输入。

### 运行接口变更
- `--model glm` / `glmm` 现在自动拒绝 `--mode trans/gw`(count 模型只有 cis 路径，否则会按数百基因 × 500k SNP PQL 重拟合）。
- `make test` 变绿色。

### 数值准确性
- `p_from_t`、`beta_cdf`、`pnorm_two_sided` 与 scipy 一致（相对能差 < 1e-12)。
- bed 2-bit 编码解码与 PLINK2 / tensorQTL / GCTA 完全一致（位操作 + 正读双跑）。
- GRM 计算与 `gcta64 --make-grm` 匹配（sub-divider 正确包含所有 non-missing markers)。

## [0.1.0] — 首次提交
