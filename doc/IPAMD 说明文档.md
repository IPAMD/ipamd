# IPAMD 说明文档
*v0.0.32*\
作者：Xiaoyang Liu

---

## 目录
- [简介](#简介)
- [安装](#安装)
- [使用方法](#使用方法)
- [插件包](#插件包)
- [二次开发](#二次开发)
- [支持](#支持)

## 简介
IPAMD 是一个基于插件架构的 Python 包，用于使用粗粒化（CG）方法模拟生物大分子。它提供了简洁高效的方式来构建、模拟和分析生物分子体系。目前，IPAMD 支持多种一珠一氨基酸 CG 模型，例如 HPS 模型和 Calvados 模型，以及双珠 RNA 模型（Calvados RNA）。更多 CG 模型可以通过插件包添加（例如 Mpipi 模型）。IPAMD 易于使用和扩展。用户可以通过编写插件方便地为该软件包添加新功能。编写插件的详细说明见[二次开发](#二次开发)一节。

## 安装
### 环境要求
IPAMD 是一个 Python 包。使用本软件包需要先在系统中安装 Python。可以从官方网站下载 Python，但更推荐通过 Miniforge 安装。Miniforge 可通过下面的链接下载。
[conda-forge/miniforge: A conda-forge distribution.](https://github.com/conda-forge/miniforge)
本软件包**必须**在 Python 3.12 下运行。可以用 mamba 创建虚拟环境。
```bash
mamba create -n ipamd python==3.12
mamba activate ipamd
```

### 从源码安装
从源码安装时，需要先将工作目录切换到软件包的根目录，然后运行下面的命令。
```bash
python setup.py install
#或
pip install .  #推荐
```

### 从 PyPI 安装
目前本软件包也已发布到 PyPI。可以用 pip 安装，这样可以得到最新的稳定版本。
```
pip install ipamd
```

核心依赖包括 `numpy`、`numba`、`rich`、`pybioseq`、`periodictable` 和 `pypdbio`。部分可选插件可能还需要 `matplotlib`、`pandas`、`scipy` 或 `biopython` 等额外软件包。安装对应插件包时会自动安装这些依赖。

## 使用方法
在 Python 脚本中导入 ipamd 后，即可使用本软件包提供的功能。各函数的详细用法见后续章节。示例代码也放在 demo 目录中。可以通过 `ipamd -s` 列出官方示例。

### 数据存储
默认情况下，输出数据由 ipamd 管理在另一个目录中，以保持工作目录整洁。数据目录位于 `~/ipamd_output`。如果不希望这样，可以修改位于 `~/.config/ipamd/config.json` 的配置文件。修改配置文件的详细说明见[配置](#配置)一节。

### App
App 是本软件包的主模块。首先需要从本软件包导入 App 类并创建实例。创建 App 实例时可以设置项目名称。项目名称会用作输出目录的名称。如果未设置名称，默认会使用一个随机字符串。
```python
from ipamd import App

app = App(name='protein')
```
创建 app 对象后，默认会自动加载全部插件。随后即可使用本软件包提供的功能。

App 类还提供了若干用于初始化模拟的函数，如下所示。

- `__init__`：初始化 app 对象。该函数有三个参数。
- `name`：项目名称。默认是一个随机字符串。

- `init_working_dir`：是否初始化工作目录。默认是 True。通常情况下不应修改该参数。

- `gpu_id`：要使用的 GPU 编号。默认是 0。该编号与 `nvidia-smi` 中显示的值一致。

- `gen_link`：该函数可以生成指向数据目录的符号链接。当需要使用本软件包提供的数据文件时，这个函数会很有用。

- `available_ff`：打印系统中可用的力场。

- `use`：加载力场。第一个参数可以是力场名称（`app.available_ff()` 打印出的名称之一），也可以是与力场文件格式相同的字典。如果省略名称或设为 `'default'`，将使用配置文件中 `default_ff` 指定的力场。可选参数 `override`（默认 True）控制新力场是否替换当前力场。如果 `override` 为 False，新的原子定义和参数会合并到当前力场中。这在需要把蛋白质力场与额外力场（例如 Calvados RNA 或弹性网络）组合使用时很有用。

```python
app.use('Calvados3')
app.use('Calvados_rna', override=False)
```

- `switch`：切换当前项目。工作目录会切换到与新项目名称对应的目录。该函数主要用于开发。

- `load_file`：把当前工作目录中的数据文件复制到本软件包的数据目录。如果选择不由本软件包管理数据文件，可以忽略该函数。如果在配置文件中设置了 `$decentralized`，则不应使用该函数。
- `export_file`：把本软件包数据目录中的数据文件复制到当前工作目录。如果选择不由本软件包管理数据文件，可以忽略该函数。第一个参数是要导出的文件名，第二个参数是目标目录。如果未设置目标目录，文件会被复制到当前工作目录。

App 对象还提供若干模块：`builder`、`simulation`、`analysis`、`sakuanna`、`mdanalysis` 和 `data_process`。其中 `sakuanna`、`mdanalysis` 和 `data_process` 会在首次访问时延迟加载。

### 模型
MD 模拟中的一些常见概念定义在 `ipamd.public.models.*` 中。大多数情况下不需要手动创建这些对象，因为本软件包提供的函数会为你创建它们。因此这里只做简要介绍。
#### Atom
Atom 对象定义原子的速度、电荷和质量。如果设置了这些值，它们会覆盖力场中定义的默认值。

#### Molecule
Molecule 对象定义所有原子的位置以及原子之间的成键关系。分子还可以包含键角项和刚性基团。可以通过 `molecule.transform()`、`molecule.rotate()` 和 `molecule.move()` 对分子进行变换、旋转或平移。

#### Frame
Frame 是包含若干分子的对象。读取帧中分子数据有两种方式。第一种是直接读取 frame 对象的分子数组。第二种是运行 `frame.properties()`。返回值是包含体系性质的 Dist。

#### Environment
Environment 对象定义模拟环境，例如温度、介电常数等。
可以通过下面的方式创建环境对象：
```python
from ipamd.public.models.md import Environment
env = Environment()
```
也可以在初始化时提供环境名称参数，此时会从对应模板创建环境对象。目前可用的名称包括 `'pure water'`、`'normal saline'` 和 `'cytoplasm'`。

还可以通过以下方法进一步修改环境：
- `set_temperature(temperature, unit="K")`：设置温度。可用单位为 `K` 和 `C`。介电常数会自动更新。
- `set_ionic_strength(ionic_strength)`：设置离子强度，单位为 mol/L。
- `set_ph(ph)`：设置 pH。
- `set_pressure(pressure, unit="atm")`：设置压强。可用单位为 `atm` 和 `Pa`。

盒子可以通过 `box.in_solvent(env)` 或 `box.in_solvent('pure water')` 应用环境。

#### Box
Box 对象定义模拟盒子。盒子包含若干帧以及环境。如果要读取某一帧的数据，应先运行 `box.frame(n)` 设置当前帧，再运行 `box.current_frame()`。也可以通过 `box.compute(method, title='', target_frame='', **kwargs)` 在选定帧上计算分析算子。`target_frame` 参数接受诸如 `'50-100'` 这样的范围字符串。

*注意*：通过 box 进行计算是旧功能，将逐步移除。

#### Unit
IPAMD 中默认的长度和时间单位是纳米、皮秒。如果要使用其他单位，可以用 Unit 类进行换算。下面是在 IPAMD 中使用埃和飞秒的示例。
```python
from ipamd.public.models.md import Unit
time = 100 * Unit.TimeScale.fs  
length = 100 * Unit.LengthScale.A
```

#### Sequence
DNA、RNA 和蛋白质序列以类似的方式表示。可以直接从 `ipamd.public.models.sequence` 创建序列对象，也可以通过 sakuanna 模块从 fasta 文件加载。序列对象以类似 List 的方式编写。可以像列表一样修改、切片和遍历序列对象。
```python
from ipamd.public.models.sequence import ProteinSequence, RNASequence

seq = ProteinSequence(name='ELN', sequence='VPGAGVPGAGVPGAG')
rna = RNASequence(name='polyU', sequence='U' * 30)
```

### 运行模拟
需要先创建模拟盒子，然后通过 `simulation = app.simulation.new_simulation(parameters)` 创建模拟对象。最后可以通过 `simulation.run(timesteps)` 运行模拟。创建模拟对象时可用的参数如下。
- `simulation_box`：模拟盒子对象。
- `job_name`：任务名称。默认是 `simulation`。所有输出文件都会使用该前缀。
- `dt`：模拟的时间步长。默认是 0.01 ps。
- `total_time`：模拟总时长。默认是 0，表示不运行模拟。
- `snap_shot`：写入快照文件的间隔。默认是 0，表示不写入快照文件。
- `run_step`：要运行的步数。默认是 0，表示不运行模拟。
- `period`：写入轨迹的间隔。如果设置了 `total_time` 和 `snap_shot`，则 `run_step` 和 `period` 会被忽略。
- `thermo_bath`：要使用的积分器。默认是 `langevin_nvt`。
- `res_auto_read`：是否自动读取重启文件。默认是 True。如果希望在模拟结束后做分析，这个选项会很有帮助。
- `minimize_energy`：是否在运行模拟前进行能量最小化。默认是 True。
- `fixed_particle`：应从积分器中排除的粒子索引列表。默认是空列表。

`simulation.run()` 可以接受一个可选参数以跳过交互式提示。例如，`simulation.run('n')` 不会重新运行已完成的模拟，`simulation.run('y')` 则会重新运行。

典型工作流如下：
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

### 构象转换
IPAMD 还提供了一些构象转换工具。可以通过执行 `box.converter_function()` 调用这些工具。这些工具包括：

`align_center`：将整个体系的质心对齐到盒子中心。该函数不需要参数。

`merge_nearest`：有时液滴中的分子会被模拟盒子分开。该函数可以把这些分子合并回来。该函数不需要参数。

`read_xml`：从 .xml 文件读取构象。唯一的必需参数是 `filename`，即 .xml 文件的名称。注意所提供的文件名不应带有 .xml 扩展名。

`unwrap`：在 MD 模拟中通常使用周期性边界条件（pbc）。所有跑出盒子的粒子都会被折回盒子内。该函数可以把粒子展开到原始位置。该函数不需要参数。在分析蛋白质扩散时这个函数会很有用。

`to_pdb`：将构象保存为 .pdb 文件。必需参数是 `filename`，即不含扩展名的 .pdb 文件名。可选参数包括 `ignoring_pbc`（默认 True）和 `atom_type_override`。

`to_xml`：与 `read_xml` 类似，但将构象保存为 .xml 文件。XML 是 IPAMD 的原生格式，因此可以用该函数保存模拟盒子并在之后继续模拟。可选参数包括 `ignoring_pbc`（默认 False）和 `indent`（默认 True）。

### 力场
IPAMD 中的所有力场都以 .xml 格式提供。力场文件位于本软件包的数据目录中。每个力场文件包含两部分，第一部分是原子定义，第二部分是相互作用势。原子定义部分包含每个 CG 珠的电荷和质量。相互作用势部分包含每种相互作用的形式和参数。力场文件的格式如下。
```xml
<?xml version="1.0" encoding="UTF-8" ?>
<ff>
    <atom_definition>
        <ATOM charge="0" mass="0"/>
    </atom_definition>
    <ff_param>
        <FORCE_NAME global_parameter="value">
            <ATOM parameter="value"/>
        </FORCE_NAME>
    </ff_param>
</ff>
```
也可以定义随模拟条件变化的参数，方法是将值设为以 `compute:` 开头的字符串。前缀之后可以写一个 python 表达式来计算该值。环境中的运行时物理量可以通过 `{physical_value}` 使用。

在运行时，也可以创建一个与力场文件格式相同的字典，并通过 `app.use(dict)` 加载到 app。

内置力场文件包括 `Calvados1`、`Calvados2`、`Calvados3`、`Calvados_lc`、`Calvados_rna`、`HPS-Urry`、`HPST`、`droplet_prepare` 和 `slab_prepare`。Mpipi 等额外力场可以通过插件包安装。

### 模块
本软件包的主要功能组织在 app 对象的若干模块中。各模块及其包含的函数如下。

#### builder
builder 模块包含一些用于构建体系的函数。各函数的详细信息如下。

`app.builder.new_simulation_box`：该函数创建一个空的模拟盒子。必需参数是 `x`、`y` 和 `z`，即盒子边长，单位为纳米。返回的盒子已经包含一个空帧。可以在该函数后链式调用 `.in_solvent(...)`。

```python
box = app.builder.new_simulation_box(20, 20, 20).in_solvent('pure water')
```

`app.builder.box_from_xml`：该函数可以从 .xml 文件生成模拟盒子。唯一的必需参数是 `filename`，即 .xml 文件的名称。当需要构建包含多个分子的复杂体系，或从先前的模拟继续时，这个函数会很有用。

`app.builder.molecule_from_xml`：该函数从 .xml 文件读取一个分子。如果文件包含多个分子，只返回第一个。该分子会对齐到临时盒子的中心。

`app.builder.cg_molecule_from_pdb`：该函数可以从 pdb 文件生成粗粒化分子。CG 珠的名称与 pdb 文件中的残基名相同。必需参数是 `file_name`，即 pdb 文件的名称。其他可选参数包括 `rigid_range`，该参数指定哪些残基应视为刚体。该参数的值应类似 `'5-20,30-50'`，表示残基 5 到 20 以及 30 到 50 是刚体。默认是 None，表示所有残基都是柔性的。

`app.builder.download_pdb`：该函数可以从 RCSB PDB 数据库下载 pdb 文件。唯一的必需参数是 `pdb_id`，即 RCSB PDB 数据库中的编号。pdb 文件会写入工作目录。

`app.builder.load_example_molecule`：该函数把内置示例分子复制到工作目录。目前可用的名称包括 `'fus'` 和 `'polyU'`。调用该函数后，可以通过 `protein_from_pdb` 或 `calvados_rna_from_pdb` 构建分子。

```python
app.builder.load_example_molecule('fus')
mol = app.builder.protein_from_pdb('fus.pdb', rigid_from_plddt=True)
```

`app.builder.gen_elastic_network`：该函数可以为蛋白质生成弹性网络。该函数有两个参数：`molecule` 是目标分子，`max_gap` 是来自两个不同结构区域的残基之间的最大间隔。`max_gap` 的默认值是 4。该函数会返回两个值，第一个是目标分子，第二个是该函数生成的额外力场。需要在模拟中使用这个额外力场，弹性网络才会生效。

示例：
```python
mol, extra_ff = app.builder.gen_elastic_network(
    molecule,
    max_gap=4
)
app.use(extra_ff, override=False)
```

`app.builder.linear_protein`：该函数可以按指定序列生成线性蛋白质。必需参数是 `protein`，即蛋白质序列对象。

`app.builder.spiral_protein`：该函数可以生成阿基米德螺旋形状的蛋白质。必需参数与 `linear_protein` 相同。另一个可选参数是 `d`。该值表示径向距离。默认值是 0.38 nm。

`app.builder.stacked_protein`：该函数可以从序列生成紧凑的、填充空间的蛋白质构象。残基按组沿类似三维 Hilbert 曲线的路径放置。必需参数是 `protein`。可选参数 `span`（默认是 2）是每组中的残基数。

`app.builder.molecule_from_curve`：该函数沿用户定义的曲线放置序列。必需参数是 `mol`（序列对象）和 `curve`（把残基索引映射到 `(x, y, z)` 的函数）。可选参数包括：
- `start_add`、`end_add`：添加到首尾珠的额外质量。默认是 0。
- `cg`：粗粒化标签。默认是 `MC`。
- `rename_map`：把残基名重映射为原子类型的字典。
- `new_bond`：相邻珠之间的键类型。默认是 `B-B`。

`app.builder.protein_from_pdb`
该函数可以从 pdb 文件自动生成蛋白质对象。pdb 文件应位于工作目录中。该函数的参数包括
- `file_name`：pdb 文件的名称。
- `ignoring_h`：是否忽略 pdb 文件中的氢原子。默认是 True。
- `cg`：粗粒化珠应位于何处。可用值为 `alpha` 和 `mass_center`。默认是 `mass_center`。
- `rigid_range`：该参数指定哪些残基应视为刚体。可参考 `cg_molecule_from_pdb`。多个刚性基团可以用 `;` 分隔。
- `rigid_from_plddt`：是否根据 alphafold2 的 plddt 分数判断残基是否为刚体。默认是 `False`。
- `threshold`：当 `rigid_from_plddt` 设为 true 时，该参数决定 plddt 分数大于阈值时残基可视为刚体，即具有结构的蛋白质。默认是 70。
- `max_gap`：可参考 `gen_elastic_network`。

示例：
```python
protein = app.builder.protein_from_pdb(
    'protein',
    rigid_from_plddt=True
)
```

`app.builder.update`：该函数更新分子中每个原子所附着的力场。在更改 app 的力场之后，这个函数会很有用。

`app.builder.calvados_rna_from_pdb`：该函数从 pdb 文件生成双珠 Calvados RNA 分子。骨架珠放在磷原子上（`NMPB`），侧链珠放在 N1（U/A/C）或 N9（G）上。同时会添加沿骨架的键角项。

`app.builder.calvados_rna_from_sequence`：该函数从 RNA 序列和曲线函数生成双珠 Calvados RNA 分子。必需参数是 `mol`（`RNASequence` 对象）和 `curve_func`。可选参数包括 `reverse` 和 `fallback_ref_vec`。

示例：
```python
from ipamd.public.models.sequence import RNASequence
seq = RNASequence('polyU', 'U' * 30)
mol_rna = app.builder.calvados_rna_from_sequence(
    seq,
    lambda index: (0, 0, index * 0.59)
)
```

#### genbox
这部分包含一些用于在模拟盒子中放置分子的函数。基本用法是运行 `box.gen_function(parameters)`。可用函数如下。

`place_molecule_randomly`：在模拟盒子中随机放置分子。必需参数如下。
- `molecule`：要添加的分子对象；n，分子数量。
- `threshold`：两个粒子之间的最小距离，默认是 1。
- `max_tries`：函数尝试放置一个分子的最大次数，默认是 5。
- `strict`：如果该参数设为 true，程序会持续尝试直到达到目标 `n`，默认是 false。
- `allow_out_of_box`：是否允许分子放置在盒子边界之外，默认是 true。如果要通过 `sub_box` 生成盒子，这个选项会很有用。

`place_molecule_periodically`：在模拟盒子中周期性放置分子。必需参数如下。
- `molecule`：要添加的分子对象。
- `nx`、`ny`、`nz`：分别沿 x、y、z 方向的分子数量。

`place_molecule_at`：在指定位置放置一个分子。必需参数如下。
- `molecule`：要添加的分子对象。
- `x`、`y`、`z`：放置分子的坐标。

`enneutral`：放置离子直到体系电中性。仅在使用 Calvados-lc 力场运行模拟时才需要该函数。该函数不需要参数。

`sub_box`：可以用该函数通过组合已有盒子来生成模拟体系。参数如下。
- `sub_box`：要放入的已有模拟盒子。
- `x`、`y`、`z`：放置 `sub_box` 的坐标。

#### analysis
analysis 模块包含一些用于分析模拟结果的算子。可以用这些算子获取体系信息，例如温度、势能、回转半径等。这些算子可以通过模拟盒子的 `compute` 函数调用。各算子的详细信息如下。

`app.analysis.contact_map`：该函数可以计算两类分子之间的接触图。该函数有三个参数。
- `type1`：第一类分子的名称。
- `threshold`：如果两个 CG 珠之间的距离小于该阈值，则这两个 CG 珠视为接触。默认是 4。
- `type2`：第二类分子的名称。默认是 ''，表示与 type1 相同。

`app.analysis.momentum`
该函数是计算体系动量的算子。不需要其他参数。函数会返回体系沿 x、y、z 轴的动量。
通常体系动量应为零或非常接近零。非零动量表明体系不正确，应重新运行模拟。

`app.analysis.temperature`
该函数是计算体系温度的算子。不需要其他参数。函数会返回体系温度。

`app.analysis.rg`
该函数是计算体系回转半径的算子。可选参数 `type_` 可用于选择分子类型。函数会返回体系的回转半径。

`app.analysis.potential`
该函数是计算体系势能的算子。该函数需要把模拟对象作为参数。函数会返回体系的势能。
*注意：* 使用该函数前必须先完成模拟。

`app.analysis.rmsd`
该函数是计算体系均方根偏差（RMSD）的算子。函数会返回体系的 RMSD。

`app.analysis.rmsf`
该函数是计算每个珠的均方根涨落（RMSF）的算子。它需要多帧。可选参数 `n` 用于选择分子索引。默认是第一个分子。

`app.analysis.flory_monomer`
该函数从一系列帧估计单链聚合物的 Flory 标度指数。可选参数包括 `type_`（分子类型）和 `sub`（采样密度，默认是 1）。分子应至少包含 30 个珠。帧数过少可能导致结果不准确。

`app.analysis.flory_multimer`
该函数从不同长度链的混合物估计 Flory 标度指数。它需要多帧。

`app.analysis.ripley`
该函数计算粒子分布的 Ripley K 函数。参数包括 `start_d`、`step`、`end_d`（距离范围，默认 1、1、10）、`l`（是否返回 L 函数）以及 `ref`（是否返回随机分布的理论参考值）。

示例：
```python
flory = box.compute(
    app.analysis.flory_monomer,
    title='flory',
    target_frame='50-100'
)
flory.print()
```

#### mdanalysis
`app.mdanalysis` 是推荐用于轨迹分析的模块。与 `app.analysis` 相比，这些算子返回带类型的数据对象（`Scalar`、`Vector`、`Matrix`），可以由 `app.data_process` 处理。大多数算子接受 `box` 和 `target_frame`。可以直接计算一帧，也可以通过 `batch_compute` 对一段帧范围取平均。

`app.mdanalysis.batch_compute`：在一段帧范围上计算某个 mdanalysis 算子并返回平均值。第一个参数是算子名称。其他参数会传递给该算子。

```python
rg = app.mdanalysis.batch_compute(
    'rg_v1',
    box=box,
    target_frame='51-200'
)
app.data_process.print(rg)
```

`app.mdanalysis.rg_v1`：计算回转半径。可选参数 `target_molecule` 可以是分子名称，或带残基范围的名称，例如 `'fus:1-100'`。如果未设置，将包含所有分子。结果单位是 nm。

`app.mdanalysis.density_align`：计算沿某一轴的质量密度。参数包括 `target_molecule`、`direction`（`X`、`Y` 或 `Z`，默认是 `Z`）以及 `d`（分箱大小，单位 nm，默认是 1）。结果单位是 g/mL。

`app.mdanalysis.density_radius`：从原点计算径向质量密度。必需参数是 `cutoff`。可选参数包括 `origin`（默认是 `(0, 0, 0)`）、`d` 和 `target_molecule`。

`app.mdanalysis.density_box`：计算子盒子内的质量密度。参数包括 `x0`、`y0`、`z0`、`lx`、`ly` 和 `lz`。如果长度为 0，将使用整个盒子的尺寸。结果单位是 g/cm3。

`app.mdanalysis.contact_map_v1`：计算接触图。参数包括 `threshold`（默认是 4）、`type1`、`type2` 和 `mode`。如果 `type1` 和 `type2` 相同，`mode` 可以是 `inter`（仅分子间接触）、`intra`（仅分子内接触）或空（仅排除自身接触）。

`app.mdanalysis.contact_number`：从接触图计算每个残基的接触数。参数与 `contact_map_v1` 相同。也可以通过 `cm` 传入已有的接触图。

`app.mdanalysis.net_charge`：计算体系净电荷。结果单位是 e。

`app.mdanalysis.slab_free_energy`：从 slab 模拟估计转移自由能以及浓相/稀相密度。参数包括 `direction`（默认是 `Z`）和 `d`（分箱大小）。结果是一个向量，包含自由能（kJ/mol）、浓相密度和稀相密度。

`app.mdanalysis.momentum_v1`：计算体系总动量。结果是沿 x、y、z 的向量，单位为 kg m/s。无其他参数。

`app.mdanalysis.temperature_v1`：计算体系瞬时温度。刚性基团按 6 个自由度、柔性珠按 3 个自由度计数。结果单位是 K。

`app.mdanalysis.potential_v1`：从模拟日志读取势能。必需参数是 `simulation`。结果单位是 kJ/mol。使用前必须先完成模拟。

`app.mdanalysis.ripley_v1`：计算粒子分布的 Ripley K 函数。参数包括 `start_d`、`step`、`end_d`（距离范围，默认 1、1、10）、`l`（是否返回 L 函数）、`ref`（是否返回随机分布的理论参考值）以及 `target_molecule`。

`app.mdanalysis.rmsd_v1`：计算轨迹相对第一帧的 RMSD。必需参数是 `box` 和 `target_frame`（帧范围）。可选参数 `n` 选择分子索引，`target_molecule` 按分子名或 `'name@1-100'` 选择。结果单位是 nm。

`app.mdanalysis.rmsf_v1`：计算每个珠相对第一帧的 RMSF。参数与 `rmsd_v1` 相同。需要至少两帧。结果单位是 nm。

`app.mdanalysis.flory_monomer_v1`：从一段轨迹估计单链聚合物的 Flory 标度指数。参数包括 `box`、`target_frame`、`target_molecule`（或旧参数 `type_`）以及 `sub`（采样密度，默认 1）。分子应至少包含 30 个珠。

`app.mdanalysis.flory_multimer_v1`：从不同长度链的混合物估计 Flory 标度指数。参数包括 `box`、`target_frame` 和 `target_molecule`。

`app.mdanalysis.rdf`：计算径向分布函数 g(r)。参数包括 `type1`、`type2`、`dr`（分箱大小，默认 0.1 nm）、`r_max`（默认盒子最短边的一半）以及 `mode`（`inter` / `intra`，仅当两类粒子相同时生效）。可用 `batch_compute` 对多帧取平均。

`app.mdanalysis.msd`：计算均方位移。必需参数是 `box` 和 `target_frame`。默认跟踪分子质心；`per_particle=True` 时按珠计算。可选参数 `dt` 为相邻所选帧之间的时间间隔，`n` / `target_molecule` 用于选择分子。结果单位是 nm²。

示例：
```python
app.data_process.plot(
    app.mdanalysis.batch_compute(
        'contact_map_v1',
        target_frame='5-10',
        box=box,
        threshold=4,
        mode='inter'
    ),
    save_figure=True
)
```

#### simulation
simulation 模块包含相互作用势和积分器。
该模块的第一部分是相互作用势。相互作用势用于描述 CG 珠之间的相互作用。一般来说，参数会根据所用力场自动设置，不需要手动设置。可用的相互作用势如下。
- `app.simulation.ah`：Ashbaugh-Hatch 势。
- `app.simulation.bond_harmonic`：谐振子键势。
- `app.simulation.angle_harmonic`：谐振子键角势。
- `app.simulation.debye`：Debye 势。
- `app.simulation.pppm`：pppm 势。
- `app.simulation.stacking`：Calvados RNA 使用的堆积势（侧链堆积）。
- `app.simulation.centripetal`：用于制备液滴或 slab 构型的向心力。
- `app.simulation.wf`：Wang-Frenkel 势。该势由 `mpipi_forcefield` 插件包提供。

该模块的第二部分是积分器。创建模拟时可以选择使用哪个积分器。积分器的参数由 app 的 `env` 生成，不需要手动设置。可用的积分器如下。
- `app.simulation.langevin_nvt`：Langevin NVT 积分器。
- `app.simulation.em`：NVE 积分器，可用于能量最小化。
- `app.simulation.noose_hover_nvt`：Noose Hover NVT 积分器。
- `app.simulation.npt_z`：NPTMTK 积分器。仅在 z 方向有压强耦合。
- `app.simulation.andersen_npt`：Andersen NPT 积分器，在所有方向都有压强耦合。

#### sakuanna
Sakuanna 是序列分析模块。它包含一些用于分析蛋白质序列的函数。各函数的详细信息如下。

可以直接创建序列对象：
```python
from ipamd.public.models.sequence import ProteinSequence
sequence = ProteinSequence(
    name='ELN',
    sequence='VPGAGVPGAGVPGAG'
)
```

`app.sakuanna.sequence_from_fasta`：从工作目录中的 fasta 文件加载序列对象。必需参数是 `fasta_path`。可选参数 `mol` 可以是 `'protein'`、`'dna'` 或 `'rna'`。默认是 `'protein'`。如果 fasta 文件包含一条序列，将返回单个序列对象。如果包含多条序列，将返回序列对象组成的元组。

示例：
```python
seq = app.sakuanna.sequence_from_fasta('proteins.fasta', mol='protein')
```

`app.sakuanna.to_fasta`：把序列对象写入 fasta 文件。可选参数包括 `filename` 和 `comment`。如果未设置 `filename`，文件将以序列名称命名。

`app.sakuanna.pretier`：该函数可以更美观地打印蛋白质序列。该函数有三个参数。
- `protein`：要打印的蛋白质序列对象。
- `word3`：是否以三字母代码打印序列。默认是 True。
- `ter`：是否打印末端基团。默认是 True。
示例：
```python
sequence = ProteinSequence(
    name='test',
    sequence='AAAAAAAAAAAAAAAAAAAAAAAA'
)
app.sakuanna.pretier(
    sequence, 
    word3=True, 
    ter=True
)
```

`app.sakuanna.statistic`：该函数可以计算蛋白质序列中特定氨基酸的比例。该函数有三个参数。
- `protein`：要分析的蛋白质序列对象。
- `targets`：要分析的目标氨基酸或氨基酸类别。该值可以是氨基酸和氨基酸类别组成的列表。可用的氨基酸类别包括 `aromatic`、`positive`、`negative`、`charged`、`polar`、`nonpolar`。
- `format_`：返回值的格式。可用值为 `ratio` 和 `count`。如果值为 `ratio`，函数会返回目标氨基酸在蛋白质序列中的比例。如果值为 `count`，函数会返回目标氨基酸在蛋白质序列中的计数。默认是 `ratio`。

示例：
```python
res = app.sakuanna.statistic(
    sequence, 
    targets=['A', 'charged'],
    format_='ratio'
)
app.data_process.plot(res)
```

`app.sakuanna.tag`：该函数通过匹配子序列为蛋白质序列中的残基打标签。参数 `tags` 是一个字典。键是标签名，值是子序列列表。可以把某个标签设为 `'rest'`，用于标记未匹配的残基。返回值是一个 `Vector`，可以通过 `app.data_process.plot` 绘图。

```python
tags = app.sakuanna.tag(
    sequence,
    tags={
        'aromatic': ['F', 'W', 'Y'],
        'other': 'rest'
    }
)
```

#### data_process
mdanalysis 和 sakuanna 的输出数据以带类型的对象提供（`Scalar`、`Vector`、`Matrix`、`Ratio`、`Distribution`、`PointSet`、`String`）。可以通过 `app.data_process` 处理这些对象。

- `print`：更美观地打印数据。可选参数 `precision`（默认是 3）。
- `plot`：显示数据的图。可选参数包括 `style`、`save_figure` 以及 `appearance`（对向量可以是 `line`、`bar`、`heatmap` 或 `discrete_heatmap`）。
- `flatten`：把向量降为标量，或把矩阵降为向量。参数 `by` 可以是 `average`（默认）或 `sum`。对于矩阵，需要提供 `axis`。
- `average`：对若干相同类型的数据对象求平均。
- `sum`：对若干相同类型的数据对象求和。
- `normalize`：归一化数据。对于比例数据，数值会被缩放使总和为 1。对于向量或矩阵，使用最小-最大归一化。对于点集，可以把 `target` 设为 `'x'`、`'y'` 或两者。
- `to_csv`：把数据保存为 csv 文件。必需参数是 `filename`。

较旧的 `box.compute` 路径仍会返回 `AnalysisResult` 对象。可以对该对象调用 `print`、`distribution`、`flatten`、`merge`、`normalize`、`plot` 和 `save`。对于新代码，推荐使用 `app.mdanalysis` 加上 `app.data_process`。

### OmicsLoader
OmicsLoader 是用于加载组学数据的模块。该模块可以让批量蛋白质分析更轻松。它可以从 fasta、pdb 等多种格式加载数据。可以使用该模块提供的函数自动加载数据。
要使用该模块，需要先导入 OmicsLoader 和 batch_run。第一个类用于读取你提供的数据，第二个函数则把不同蛋白质的任务分发到不同 GPU。
```python
from ipamd import OmicsLoader, batch_run
```
首先需要用数据目录的路径创建 OmicsLoader 实例。然后需要创建一个处理数据的函数。该函数应有两个参数，第一个参数是 app，第二个参数是 data。这两个参数会由 batch_run 函数传入。data 参数是一个字典，包含数据的名称、序列、类型和路径。`batch_run` 还接受可选参数 `gpus`，它是 GPU 编号的范围字符串，例如 `'0-3'`。如果未设置，将使用所有可用 GPU。详细示例如下。
```python
loader = OmicsLoader(path='input_data/transcription')
def process_sequence(app, data):
    app.use('Calvados3')
    box = app.builder.new_simulation_box(150, 150, 150)
    mol = app.builder.protein_from_pdb(
        data['name'] + '.pdb', 
        rigid_from_plddt=True
    )
    box.place_molecule_periodically(mol, 1, 1, 1)
    simulation = app.simulation.new_simulation(
        box,
        f'{data['name']}', 
        total_time=20 * Unit.TimeScale.ns, 
        snap_shot=100
    )
    simulation.run('n')
    flory = box.compute(
        app.analysis.flory_monomer, 
        title='flory', 
        target_frame='50-100'
    )
    flory.print()
    with open('results.csv', 'a') as f:
        f.write(f'{data["name"]},{flory.data}\n')

batch_run(loader, process_sequence)
```

### 配置
本软件包的配置文件位于 `~/.config/ipamd/config.json`。首次运行本软件包时会创建该配置文件。可以修改配置文件来改变软件包的行为，例如数据文件路径、额外插件路径以及输出文件路径等。详细配置项如下。

| 配置项 | 说明 | 可用值 |
|--------|------|--------|
| result_dir | 模拟结果的写入位置 | 目标目录，或设为 `$decentralized` 如果你不希望输出由 IPAMD 管理 |
| output_level | 打印多少日志 | 0：verbose，1：info，2：warning，3：error |
| auto_load | 是否自动加载插件 | True 或 False |
| default_ff | 默认力场 | 力场名称，应为 `app.available_ff()` 打印出的可用力场之一 |
| external_plugin_dir | 额外插件目录 | 目标目录列表 |
| sakuanna_plugin_dir | 额外 sakuanna 插件的路径 | 目标目录列表 |
| simulation_plugin_dir | 额外 simulation 插件的路径 | 目标目录列表 |
| builder_plugin_dir | 额外 builder 插件的路径 | 目标目录列表 |
| analyse_plugin_dir | 额外 analyse 插件的路径 | 目标目录列表 |

### 示例
demo 目录中提供了一些示例代码。可以运行这些示例来熟悉本软件包。demo 目录中有若干子文件夹，每个子文件夹包含示例代码以及所需的输入数据。也可以通过 `ipamd -s` 列出官方示例。示例的含义可以参考本软件包的论文（doi: 10.1021/acs.jctc.5c00147）。

| 示例 | 说明 |
|------|------|
| single_chain_simulation | 运行单链模拟 |
| sequence_analysis | 对蛋白质序列进行分析 |
| basic_example | 液滴模拟与接触、RDF、MSD、密度分析 |
| slab_simulation | Slab 模拟，含密度与自由能分析 |
| RNA+Protein | 模拟 RNA 与蛋白质体系 |
| custom_system | 创建自定义分子体系 |
| omics | 使用 OmicsLoader 进行批量分析 |
| calvados_rna | 双珠 RNA 模拟 |

### 命令行工具
IPAMD 提供命令行工具 `ipamd`，用于管理插件包和列出示例。

```bash
ipamd -s                  # 显示官方示例
ipamd -l                  # 列出已安装的插件包
ipamd -i mpipi_forcefield # 安装插件包（从本地路径或从 GitHub）
ipamd -r mpipi_forcefield # 移除插件包
ipamd -p plugin_dir       # 将插件目录打包为 zip 文件
ipamd -v                  # 显示版本
```

如果 `-i` 的参数不是本地文件或目录，IPAMD 会尝试从官方插件包仓库下载 `{name}.zip`。

## 插件包
部分功能以可选插件包的形式提供。可以通过 `ipamd -i <pack_name>` 安装它们。安装后，相应的插件和数据文件会被复制到 IPAMD 的安装目录。

### mpipi_forcefield
该插件包提供 Mpipi 力场。安装后可以通过 `app.use('mpipi')` 加载该力场。该插件包还会添加 Wang-Frenkel 势（`app.simulation.wf`）。模拟前必须先用 `mpipi_chtype` 处理蛋白质。

`app.builder.mpipi_chtype`：该函数会转换有结构残基的原子类型，使其能被 Mpipi 力场识别。对于属于刚性基团的每个珠，如果原子类型尚未以 `0` 结尾，就会加上后缀 `0`。唯一的必需参数是 `molecule`。函数返回修改后的分子。

示例：
```python
app.use('mpipi')
mol = app.builder.protein_from_pdb('protein.pdb', rigid_from_plddt=True)
mol = app.builder.mpipi_chtype(mol)
```

### calvados_rna
该插件包仍在开发中。它为 Calvados RNA 模拟提供支持文件，并与 Calvados2 力场兼容。安装后可以通过 `app.use('Calvados_rna', override=False)` 把 RNA 力场合并到当前力场中。

`app.builder.calvados_rna_from_pdb`：该函数从 pdb 文件生成双珠 Calvados RNA 分子。骨架珠放在磷原子上（`NMPB`），侧链珠放在 N1（U/A/C）或 N9（G）上。同时会添加沿骨架的键角项。唯一的必需参数是 `file_name`。

`app.builder.calvados_rna_from_sequence`：该函数从 RNA 序列和曲线函数生成双珠 Calvados RNA 分子。必需参数是 `mol`（`RNASequence` 对象）和 `curve_func`。可选参数包括 `reverse` 和 `fallback_ref_vec`。

该插件包还提供 `app.simulation.stacking` 和 `app.simulation.angle_harmonic`。

示例：
```python
from ipamd.public.models.sequence import RNASequence
app.use('Calvados2')
app.use('Calvados_rna', override=False)
seq = RNASequence('polyU', 'U' * 30)
mol_rna = app.builder.calvados_rna_from_sequence(
    seq,
    lambda index: (0, 0, index * 0.59)
)
```

### plotting
该插件包添加基于 matplotlib 的绘图后端。安装后可以通过 `app.data_process.plot` 绘制分析结果。

`app.data_process.plot`：该函数显示数据的图。必需参数是 `data`。可选参数包括：
- `style`：matplotlib rc 样式字典。
- `save_figure`：是否将图保存为 `{title}.png`。默认是 False。
- `appearance`：向量图的外观。可用值为 `line`、`bar`、`heatmap` 和 `discrete_heatmap`。默认是 `line`。

标量和字符串数据无法绘图。此时会改为打印数据。

### to_pandas
该插件包将分析结果转换为 pandas DataFrame。它依赖 `pandas`。

`app.data_process.to_df`：该函数将数据对象转换为 `pandas.DataFrame`。唯一的必需参数是 `data`。支持标量、向量、矩阵、比例、分布和点集数据。

示例：
```python
df = app.data_process.to_df(rg)
df.to_csv('rg.csv')
```

### cif_format_support
该插件包将构象保存为 CIF 格式。它依赖 `biopython`。

`box.to_cif`：该函数将当前构象保存为 .cif 文件。必需参数是 `filename`。可选参数 `ignoring_pbc`（默认 True）控制写入前是否展开周期性映像。

### cg2all_integrate
该插件包将 CG 蛋白质转换为全原子蛋白质。安装后请重启 shell。

`box.to_aa`：该函数将 CG 蛋白质转换为全原子蛋白质。唯一参数是全原子蛋白质的输出文件名。
*注意：* 该函数使用 cg2all 作为计算后端，因此必须确保计算机上已安装 cg2all，并且可以在终端中调用 `convert_cg2all`。

### alphafpld2_integrate
该插件包集成 AlphaFold2/ColabFold。安装后请重启 shell。

`app.builder.af2`：该函数可以调用计算机上已安装的 alphafold2（colabfold）来生成蛋白质的初始结构。运行该函数后，蛋白质的 pdb 文件会写入工作目录。唯一的必需参数是 `protein_sequence`。该参数是一个蛋白质序列对象。
*注意：* 需要确保系统中可以使用 `colabfold_batch` 命令。

### slab_simulation
该插件包用于支持 slab 模拟。它依赖 `scipy`。安装后还会提供 `slab_prepare` 力场。

`app.mdanalysis.slab_free_energy`：该函数从 slab 模拟估计转移自由能以及浓相/稀相密度。必需参数是 `box` 和 `target_frame`。可选参数包括 `direction`（默认是 `Z`）和 `d`（分箱大小）。结果是一个向量，包含自由能（kJ/mol）、浓相密度和稀相密度。

### sequence_alignment
该插件包提供序列比对工具。

`app.sakuanna.sequence_align`：该函数对两条序列进行比对。必需参数是 `seq1` 和 `seq2`。可选参数包括：
- `algorithm`：比对算法。可用值为 `needleman-wunsch`、`smith-waterman` 和 `diff`。默认是 `needleman-wunsch`。
- `match_score`：匹配得分。默认是 1。
- `mismatch_score`：错配得分。默认是 -1。
- `gap_score`：缺口得分。默认是 -1。

函数返回三个数据对象组成的列表：比对后的参考序列、比对后的目标序列，以及比对得分。

`app.data_process.print_diff`：该函数打印两条已比对序列之间的差异。必需参数是 `ref` 和 `target`，它们应是 `sequence_align` 返回的字符串对象。可选参数包括 `print_ref`（默认 True）和 `seperate_with_bracket`（默认 True）。

示例：
```python
ref, target, score = app.sakuanna.sequence_align(seq1, seq2)
app.data_process.print_diff(ref, target)
```

### adv_data_processing
该插件包提供高级数据处理工具。它依赖 `scipy`。

`app.data_process.gaussian`：该函数对数据应用高斯滤波。唯一的必需参数是 `data`。函数返回带平滑值的新数据对象。

## 二次开发
IPAMD 是基于插件的软件，开发者可以通过编写插件方便地为本软件包添加新功能。编写插件的详细说明见本节。

### 插件
插件是包含一个或多个类或函数的 Python 模块。其中应有一个名为 `func` 的函数，该函数会被 IPAMD 自动检测并加载。每个插件都有一个配置变量。该变量类似：
```python
configure = {
    'type': 'function',
    "schema": 'schema_name',
    "apply": ['applied_value']
}
```
`type` 变量已弃用，将在未来版本中移除。唯一可用的值是 `function`。`apply` 变量指示应将哪个值传给插件。你的函数也应有一个与所应用值同名的参数。`schema` 变量表示一组已应用值的组合。但在大多数情况下 schema 变量并非必需，也不推荐开发者使用。

较新的插件还可以声明 `resource`（注入的资源，例如 `ff` 和 `persistency_dir`）以及 `alias`（调用插件时使用的名称）。

也可以在自己的插件中调用另一个插件。为此，需要从 `ipamd.public.utils` 导入 `PluginBase` 类，然后使用 `PluginBase.call('plugin_name', parameters)` 调用目标插件。`plugin_name` 是要调用的插件名称，`parameters` 是包含目标插件所需参数的字典。*注意：该功能仍为测试版，可能尚不稳定。*

插件包是包含 `meta.json` 文件、插件脚本以及可选数据文件的目录。可以通过 `ipamd -p <dir>` 打包，并通过 `ipamd -i <zip>` 安装。

## 支持
如有任何问题或建议，请通过[邮件](mailto:liuxiaoyang_Q@outlook.com)联系我。也可以在 GitHub 页面提交 issue。我会尽快回复。
