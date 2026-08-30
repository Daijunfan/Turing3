# GermSynth-R：多表型、可证明、可再生的递归算法胚种原型

GermSynth-R 从精确双线性张量规范出发，完成以下闭环：

1. 从朴素分解出发搜索低秩局部结构；
2. 将 `GF(2)` 支撑精确提升为整数系数恒等式；
3. 输出独立可检查的 JSON 证明证书；
4. 将固定大小局部恒等式递归生成为任意 `2^k` 规模算法；
5. 构造多个功能等价的算法表型，并允许递归树中逐节点切换；
6. 精确求解最小再生覆盖，在抽象局部资源失效时选择仍可执行的表型。

当前原型精确发现并验证：

- 二项式乘法的秩 3 胚种（Karatsuba 结构）；
- `2×2` 矩阵乘法的秩 7 胚种；
- 16 个不同的整数系数表型；
- 一个基数最小的三表型单故障全覆盖；
- 任意节点表型混合与逐节点单资源故障下的精确递归计算。

> 边界：本项目验证的是一个新的“多表型算法胚种与局部再生”计算内核。秩 7 的 `2×2` 公式本身不是新的矩阵乘法指数；所有结论均严格限定在证书定义的精确整数/环语义与抽象资源故障模型内。

## 环境

- Python 3.11+
- NumPy
- pytest
- 可选：支持 C++17 的 `g++` 或 `clang++`

## 一条命令复现

```bash
cd germsynth
./reproduce.sh
```

完整复现包括：搜索、整数提升、证书生成、单元测试、负向证书测试、两个独立 Python 检查器、两个独立 C++17 检查器。

## 分步运行

```bash
# 1. 搜索与全量实验
PYTHONPATH=. python3 run_all.py

# 2. 单元测试与故意破坏证书的负向测试
PYTHONPATH=. pytest -q

# 3. 不导入 germsynth 包的独立 Python 检查器
python3 independent_verify.py certificates/karatsuba_germ.json
python3 independent_verify.py certificates/matrix_rank7_germ.json
python3 independent_verify_pool.py certificates/matrix_phenotype_pool.json

# 4. 从证书生成并编译独立 C++ 检查器
python3 emit_cpp.py certificates/matrix_rank7_germ.json generated/matrix_rank7_germ.cpp
g++ -std=c++17 -O2 -Wall -Wextra -pedantic \
    generated/matrix_rank7_germ.cpp -o build/matrix_rank7_germ_verify
./build/matrix_rank7_germ_verify

python3 emit_cpp_pool.py certificates/matrix_phenotype_pool.json generated/matrix_pool_verify.cpp
g++ -std=c++17 -O2 -Wall -Wextra -pedantic \
    generated/matrix_pool_verify.cpp -o build/matrix_pool_verify
./build/matrix_pool_verify
```

## 核心目录

```text
germsynth/
├── germsynth/                     # 搜索、提升、递归执行、谱证书、覆盖优化
├── certificates/                  # 独立可检查的 JSON 证书
├── generated/                     # 由证书生成的独立 C++17 实现/检查器
├── tests/                         # 正向与负向测试
├── results/                       # 原始结果和验证日志
├── run_all.py                     # 完整研究流水线
├── independent_verify.py          # 不依赖包实现的单胚种检查器
├── independent_verify_pool.py     # 不依赖包实现的表型池检查器
├── emit_cpp.py                    # 单胚种 C++ 生成器
├── emit_cpp_pool.py               # 多表型/故障覆盖 C++ 生成器
├── RESEARCH_REPORT.md             # 研究报告、定理、结果与边界
└── reproduce.sh                   # 全量复现脚本
```

## 核心形式化对象

一个秩为 `R` 的双线性胚种是：

```text
T = Σₛ uₛ ⊗ vₛ ⊗ wₛ
```

局部张量恒等式一旦成立，将标量替换为大小减半的块并结构归纳，即得到全部 `n=2^k` 的正确性。标量乘法次数为：

```text
M(2^k) = R^k = n^(log₂ R)
```

若多个表型都实现同一个 `T`，递归树的每个内部节点可独立选择任意表型；局部恒等式在每个节点成立，因此全局正确性不依赖表型选择序列。

## 已复现的关键数字

- Karatsuba：625 个整数基例全穷举通过；`n=1024` 时严格执行 `3^10=59049` 次标量乘法。
- 矩阵乘法：6561 个整数基例全穷举通过；`128×128` 时严格执行 `7^7=823543` 次标量乘法。
- 20 个确定性搜索种子全部找到秩 7 支撑并成功整数提升，去重后得到 16 个系数表型。
- 64 阶矩阵递归树逐节点随机切换表型，结果精确通过。
- 44 种抽象局部资源的全部单点失效，由一个经穷举证明最小的三表型集合完全覆盖。
- 16 表型池对 946 个双故障集合覆盖 871 个，对 13244 个三故障集合覆盖 9695 个。

详细数据见 `results/results.json` 与 `RESEARCH_REPORT.md`。
