# eDNA Fixed-Tree Placement Backend

环境 DNA 短序列放置后端：在**固定参考树**（拓扑与原枝长不变）上，用
JC69 模型和 Felsenstein 递推对每条待测序列独立枚举所有枝，找到整树
对数似然最大的插入位置，并输出 jplace v3 结果。

## 范围

- 参考序列 3–20 条，待测序列 1–5 条，全部等长且 <= 500 列。
- 允许 ACGT、IUPAC 歧义码（RYSWKMBDHVN）与缺口 -（视为未知）。
- 空序列、非法字符、重名 ID、非正/非有限枝长、树叶与参考 ID 不一致、
  非无根二叉树等错误会**定位并整单拒绝**（HTTP 422）。
- 有根二叉 Newick（根度为 2）会被自动压缩为无根树，两段子枝长度相加，
  无根拓扑与总枝长保持不变。
- 似然计算严格按 Felsenstein 递推对祖先状态求和（含逐列缩放避免下溢），
  不使用汉明距离或最近叶替代。
- 歧义码对所有兼容状态求和；缺口不提供碱基证据。全部由缺口/完全未知
  码组成的待测序列不强制归属，返回无法放置原因。
- 对每条枝：插入点把原枝拆成两段（非负、总长不变），新叶枝长在
  [0, 2] 内优化（粗网格 + 黄金分割），最大化整树对数似然。
- 返回全部枝的结果（枝号、对数似然、相对似然权重、distal/pendant
  长度），权重在全枝集合上归一化，按似然降序、同值按枝号升序排列。
  权重是**相对支持度**，不是物种归属正确的概率。

## 运行

    .venv/bin/python -m uvicorn app.main:app --port 8000

## 接口

- POST /api/place — JSON 放置结果。
- POST /api/place/jplace — 相同计算，返回可下载的 jplace v3 文件。
- GET /api/health — 健康检查。

请求体：

    {
      "reference_fasta": ">A\nACGT...\n>B\n...",
      "query_fasta": ">q1\nACGT...",
      "newick": "((A:0.08,B:0.12):0.15,(C:0.10,(D:0.06,E:0.09):0.12):0.18);"
    }

jplace 输出含五个标准放置字段
edge_num, likelihood, like_weight_ratio, distal_length, pendant_length，
树中枝以 {edge_num} 标注，与结果中的枝号一致。

## 示例与自测

examples/ 下有 5 条参考、3 条待测（含一条全缺口不可放置序列）和参考树。

    PYTHONPATH=. .venv/bin/python tests/selftest.py

## 代码结构

- app/schemas.py — 字母表、IUPAC、输入校验
- app/io_parsers.py — FASTA / Newick 解析（Biopython）
- app/tree.py — 无根二叉树结构、稳定枝号、带枝号 Newick 序列化
- app/likelihood.py — JC69 + Felsenstein 递推（缩放、方向性缓存）
- app/optimize.py — 全枝枚举与插入点/新枝长优化、权重归一化
- app/jplace.py — jplace v3 组装
- app/main.py — FastAPI 路由与整单校验流程
