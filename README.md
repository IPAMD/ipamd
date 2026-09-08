# IPAMD

**IPAMD**是一个基于插件架构的 Python 包，用于生物大分子的**粗粒化（CG）分子动力学**模拟与分析。它提供简洁高效的 API，支持从序列/结构构建体系、GPU 加速模拟、轨迹分析到批量组学数据处理的全流程工作流。

- **版本**：0.0.31
- **Python**：>= 3.12
- **许可证**：[LGPL-3.0](LICENSE)
- **详细文档**：[doc/IPAMD Documation.md](doc/IPAMD%20Documation.md)
- **示例**：[demo/](demo/)
- **论文引用**：doi:[10.1021/acs.jctc.5c00147](https://doi.org/10.1021/acs.jctc.5c00147)

## 特性

- **序列分析模块（Sakuanna）**：蛋白质/DNA/RNA 序列读取、统计与格式化输出
- **多种 CG 力场**：内置 Calvados、HPS-Urry、HPST 等；可通过插件扩展 Mpipi 等力场
- **GPU 加速**：基于 Numba CUDA 与底层cuda引擎构建，支持gpu加速的分子动力学模拟
- **插件化架构**：构建器、分析算子、数据处理等功能均可通过插件扩展

## 安装

### 环境要求

推荐使用 [Miniforge](https://github.com/conda-forge/miniforge) 创建 Python 3.12 环境：

```bash
mamba create -n ipamd python==3.12
mamba activate ipamd
```

### 从 PyPI 安装

```bash
pip install ipamd
```

### 从源码安装

```bash
git clone https://github.com/secretqsan/ipamd.git
cd ipamd
pip install .
```

### 依赖

核心依赖：`numpy`、`numba`、`rich`、`periodictable`。部分插件（如绘图、Pandas 导出）需要额外安装 `matplotlib`、`pandas` 等，安装插件时会自动提示并安装。

## 快速开始

```python
from ipamd import App
from ipamd.public.models.md import Unit

app = App(name='example')
box = app.builder.new_simulation_box(20, 20, 20).in_solvent('pure water')

app.builder.load_example_molecule('fus')
mol = app.builder.protein_from_pdb('fus.pdb', rigid_from_plddt=True)
box.place_molecule_periodically(mol, 1, 1, 1)

simulation = app.simulation.new_simulation(
    box,
    'simulation',
    total_time=1 * Unit.TimeScale.ns,
    snap_shot=10,
)
simulation.run()
```

首次运行时会自动创建配置文件 `~/.config/ipamd/config.json`，模拟结果默认写入 `~/ipamd_output/<项目名>/`。

## 使用说明

`doc`目录下提供了详细的文档，详细内容请参考[IPAMD Documation](doc/IPAMD%20Documation.md)。

## 示例

`demo/` 目录包含完整示例，可通过 `ipamd -s` 查看链接：

| 示例 | 说明 |
|------|------|
| [single_chain_simulation](demo/single_chain_simulation/) | 单链蛋白模拟 |
| [basic_example](demo/basic_example/) | 液滴模拟与接触、RDF、MSD、密度分析（需 `plotting` 插件） |
| [sequence_analysis](demo/sequence_analysis/) | 蛋白质序列分析 |
| [slab_simulation](demo/slab_simulation/) | Slab 模拟（需 `slab_simulation`、`plotting` 插件） |
| [RNA+Protein](demo/RNA+Protein/) | RNA + 蛋白复合体系（需 `mpipi_forcefield` 插件） |
| [custom_system](demo/custom_system/) | 自定义分子体系构建与模拟 |

## 支持与反馈

- **作者**：Xiaoyang Liu
- **邮箱**：[liuxiaoyang_Q@outlook.com](mailto:liuxiaoyang_Q@outlook.com)
- **Issues**：[GitHub Issues](https://github.com/secretqsan/ipamd/issues)
