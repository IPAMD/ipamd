# IPAMD Documentation
*v0.0.32*\
Author: Xiaoyang Liu

---

## Contents
- [Introduction](#Introduction)
- [Installation](#Installation)
- [Usage](#Usage)
- [Plugin-packs](#Plugin-packs)
- [Second-time-development](#Second-time-development)
- [Supporting](#Supporting)

## Introduction
IPAMD is a plugin-based Python package designed to simulate biomolecules using CG methods. It provides a simple and efficient way to build, simulate, and analyze biomolecular systems. Currently, IPAMD supports several one-bead-one-amino-acid CG models, such as the HPS model and the Calvados models, as well as two-bead RNA models (Calvados RNA). More CG models can be added through plugin packs (for example, the Mpipi model). IPAMD is designed to be easy to use and extend. The user can easily add new functions to this package by writing a plugin. The detailed instruction on writing a plugin is shown in the [second-time-development](#Second-time-development) section.

## Installation
### prerequisites
IPAMD is a Python package. To use this package, you need to have Python installed on your system. You can download Python from the official website. But installing Python by Miniforge is recommended. Miniforge can be downloaded by the link below.
[conda-forge/miniforge: A conda-forge distribution.](https://github.com/conda-forge/miniforge)
This package *MUST* run in python 3.12. You can create a virtual environment by mamba.
```bash
mamba create -n ipamd python==3.12
mamba activate ipamd
```

### Install from source code
To install this package from the source code, you need firstly to change your working directory to the root directory of the package, and then run the following command.
```bash
python setup.py install
#or
pip install .  #recommended
```

### Install from Pypi
Currently, this package is also available on PyPI. You can also install it by pip. By this way, you can get the latest stable version of this package.
```
pip install ipamd
```

Core dependencies include `numpy`, `numba`, `rich`, `pybioseq`, `periodictable` and `pypdbio`. Some optional plugins may require extra packages such as `matplotlib`, `pandas`, `scipy` or `biopython`. These packages will be installed automatically when you install the corresponding plugin pack.

## Usage
Import ipamd in your python script, and then you can use the functions provided by this package. You can find the detailed usage of each function in the following sections. Example code is also available in the demo folder. You can list the official examples by running `ipamd -s`.

### Data storage
By default, the output data is managed by ipamd in another directory to make your working space clean. The data directory is located at `~/ipamd_output`. If you don't like this behavior, you can change this by modifying the config file located at `~/.config/ipamd/config.json`. The detailed instruction on modifying the config file is shown in the [configuration](#configuration) section.

### App
App is the main module of this package. The first thing you need to do is to import App class from this package and create an instance. When creating an instance of App, you can set the name of the project. The name of the project will be used as the name of the output directory. If the name is not set, the default name will be a random string.
```python
from ipamd import App

app = App(name='protein')
```
After creating an app object, all the plugins will be loaded automatically by default. You can then use the functions provided by this package.

We also provide several functions in the App class to help you init your simulation. The functions are listed below.

- `__init__`: initialize the app object. There are three parameters in this function.
- `name`: the name of the project. Default is a random string.

- `init_working_dir`: whether to initialize the working directory. Default is True. In normal cases, you should not change this parameter.

- `gpu_id`: the id of the GPU to use. Default is 0. The id is the same with the value shown by `nvidia-smi`.

- `gen_link`: this function can generate a symbolic link to the data directory. This function is useful when you want to use the data files provided by this package.

- `available_ff`: print the available force fields in your system.

- `use`: load a force field. The first parameter can be the name of a force field (one of the names printed by `app.available_ff()`), or a dict with the same format as a force field file. If the name is omitted or set to `'default'`, the force field specified by `default_ff` in the config file will be used. The optional parameter `override` (default is True) controls whether the new force field replaces the current one. If `override` is False, the new atom definitions and parameters will be merged into the current force field. This is useful when you want to combine a protein force field with an extra force field, such as Calvados RNA or an elastic network.

```python
app.use('Calvados3')
app.use('Calvados_rna', override=False)
```

- `switch`: switch the current project. The working directory will be changed to the directory corresponding to the new project name. This function is mainly useful in development.

- `load_file`: copy a data file from current working directory to the data directory of this package. If you choose not to manage the data files by this package, i.e. set `$decentralized` in your config file, you can ignore this function and this function shouldn't be used.

- `export_file`: copy a data file from the data directory of this package to the current working directory. If you choose not to manage the data files by this package, you can ignore this function. The first parameter is the name of the file to export, and the second parameter is the target directory. If the target directory is not set, the file will be copied to the current working directory.

The App object also exposes several modules: `builder`, `simulation`, `analysis`, `sakuanna`, `mdanalysis` and `data_process`. `sakuanna`, `mdanalysis` and `data_process` are loaded lazily when they are first accessed.

### Models
Some common concepts in MD simulation and other things are defined in `ipamd.public.models.*`. In most times you don't need to create these objects manually, as the functions provided by this package will create these objects for you. So there would just be a brief of these models.
#### Atom
The Atom object defines the velocity, the charge and the mass of an atom. If these values are set, they will override the default values defined in the force field.

#### Molecule
Molecule object defines the position of all atoms and the bonding relationship between atoms. A molecule can also contain angle terms and rigid groups. You can transform, rotate or translate a molecule by `molecule.transform()`, `molecule.rotate()` and `molecule.move()`.

#### Frame
Frame is the object containing several molecules. There are two ways of reading the molecule data contained in a frame. The first way is to directly read the molecule array of a frame object. The second way is to run `frame.properties()`. The returned value is a Dist containing the properties of the system.

#### Environment
The Environment object defines the simulation environment, such as the temperature, the dielectric constant et al.
You can create an environment object by
```python
from ipamd.public.models.md import Environment
env = Environment()
```
You can also provide a parameter of environment to the initial statement. Then the environment object will be created from the corresponding template. Currently available names include `'pure water'`, `'normal saline'` and `'cytoplasm'`.

You can further modify the environment by the following methods:
- `set_temperature(temperature, unit="K")`: set the temperature. Available units are `K` and `C`. The dielectric constant will be updated automatically.
- `set_ionic_strength(ionic_strength)`: set the ionic strength, in mol/L.
- `set_ph(ph)`: set the pH.
- `set_pressure(pressure, unit="atm")`: set the pressure. Available units are `atm` and `Pa`.

A box can apply an environment by `box.in_solvent(env)` or `box.in_solvent('pure water')`.

#### Box
Box object defines the simulation box. The box contains several frames and the environment. If you want to read the data in a specific frame, you should first run `box.frame(n)` to set the active frame and then run `box.current_frame()`. You can also compute an analysis operator on selected frames by `box.compute(method, title='', target_frame='', **kwargs)`. The `target_frame` parameter accepts a range string such as `'50-100'`.

*Note*: Computing through box is an old feature and will be removed gradually.

#### Unit
The default length and time unit in IPAMD is nanometer, picosecond. If you want to use other units, you can use the Unit class to convert the units. There is an example of using angstrom and femtosecond in IPAMD.
```python
from ipamd.public.models.md import Unit
time = 100 * Unit.TimeScale.fs  
length = 100 * Unit.LengthScale.A
```

#### Sequence
The DNA, RNA and protein sequence are represented by similar way. You can create a sequence object directly from `ipamd.public.models.sequence`, or load it from a fasta file by the functions provided by the sakuanna module. The sequence object was written in a List-like way. You can modify, slice and iterate the sequence object like a list.
```python
from ipamd.public.models.sequence import ProteinSequence, RNASequence

seq = ProteinSequence(name='ELN', sequence='VPGAGVPGAGVPGAG')
rna = RNASequence(name='polyU', sequence='U' * 30)
```

### running a simulation
You need to create a simulation box first, and then create a simulation object by `simulation = app.simulation.new_simulation(parameters)`. Finally, you can run the simulation by `simulation.run(timesteps)`. The parameters available when creating a simulation object are listed below.
- `simulation_box`: the simulation box object.
- `job_name`: the name of the job. Default is `simulation`. All the output files will have this prefix.
- `dt`: the time step of the simulation. Default is 0.01 ps.
- `total_time`: the total time of the simulation. Default is 0, which means no simulation will be run.
- `snap_shot`: the interval of writing the snapshot file. Default is 0, which means no snapshot file will be written.
- `run_step`: the number of steps to run. Default is 0, which means no simulation will be run. 
- `period`: the interval of writing the trajectory. If the `total_time` and `snap_shot` is set, the `run_step` and `period` will be discarded.
- `thermo_bath`: the integrator to use. Default is `langevin_nvt`.
- `res_auto_read`: whether to read the restart file automatically. Default is True. This will be helpful if you want to do some analysis after the simulation.
- `minimize_energy`: whether to minimize the energy before running the simulation. Default is True.
- `fixed_particle`: a list of particle indices that should be excluded from the integrator. Default is an empty list.

`simulation.run()` can take an optional argument to skip the interactive prompt. For example, `simulation.run('n')` will not rerun a finished simulation, and `simulation.run('y')` will rerun it.

A typical workflow looks like:
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

### configuration converter
IPAMD also provides some tools for configuration conversion. You can call these tools by executing `box.converter_function()`. These tools include:

`align_center`: align the mass center of the whole system to the center of the box. No parameter is needed in this function.

`merge_nearest`: sometimes molecules in a droplet would be separated by the simulation box. This function can merge these molecules. No parameter is needed in this function.

`read_xml`: read the configuration from a .xml file. The only required parameter is `filename`, which is the name of the .xml file. Note that the filename you provided shouldn't have a .xml extension.

`unwrap`: in MD simulation, a pbc (periodic boundary condition) is usually used. All the particles out of the box will be wrapped back to the box. This function can unwrap the particles to their original position. No parameter is needed in this function. This function is useful when you analyze the diffusion of the proteins.

`to_pdb`: save the configuration to a .pdb file. The required parameter is `filename`, which is the name of the .pdb file (without extension). Optional parameters include `ignoring_pbc` (default is True) and `atom_type_override`.

`to_xml`: similar to `read_xml`, but save the configuration to a .xml file. XML format is the native format of IPAMD, so you can use this function to save your simulation box and continue your simulation later. Optional parameters include `ignoring_pbc` (default is False) and `indent` (default is True).


### force fields
All the force fields in IPAMD are provided in .xml format. The force field files are located in the data directory of this package. Each force field file contains two parts, the first part is the atom definition and the second part is the interaction potentials. The atom definition part contains the charge and the mass of each CG bead. The interaction potential part contains the interaction form and the parameters of each interaction potential. The format of the force field file is shown below.
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
You can also define parameters that varying with the simulation condition by setting the value with a string starting with `compute:`. After the prefix, you can write a python expression to compute the value. The runtime values in environment can be used by `{physical_value}`.

During the runtime, you can also create a dict with the same format as the force field file, and load it to the app by `app.use(dict)`.

Built-in force field files include `Calvados1`, `Calvados2`, `Calvados3`, `Calvados_lc`, `Calvados_rna`, `HPS-Urry`, `HPST`, `droplet_prepare` and `slab_prepare`. Extra force fields such as Mpipi can be installed from plugin packs.

### modules
The main functions of this package are organized in several modules of app object. The modules and the contained functions are listed as bellow.

#### builder
The builder module contains some functions for building the system. Detailed information of each function is listed as bellow.

`app.builder.new_simulation_box`: this function creates an empty simulation box. The required parameters are `x`, `y` and `z`, which are the box lengths in nanometer. The returned box already contains one empty frame. You can chain `.in_solvent(...)` after this function.

```python
box = app.builder.new_simulation_box(20, 20, 20).in_solvent('pure water')
```

`app.builder.box_from_xml`: this function can generate a simulation box from a .xml file. The only required parameter of this function is `filename`, which is the name of the .xml file. This function is useful when you want to build a complex system with multiple molecules or continue your simulation from a previous simulation. 

`app.builder.molecule_from_xml`: this function reads a molecule from a .xml file. If the file contains more than one molecule, only the first one will be returned. The molecule will be aligned to the center of a temporary box.

`app.builder.cg_molecule_from_pdb`: this function can generate a coarse-grained molecule from a pdb file. The name of the CG beads will be the same as the residue name in the pdb file. The required parameters of this function is `file_name`, which is the name of the pdb file. The other optional parameters include `rigid_range`, this parameter indicates which residues should be regarded as rigid body. The value of this parameter should like `'5-20,30-50'`, which means residues from 5 to 20 and from 30 to 50 are rigid body. Default is None, which means all residues are flexible.

`app.builder.download_pdb`: this function can download a pdb file from RCSB PDB database. The only required parameter of this function is `pdb_id`, which is the id in RCSB PDB database. The pdb file will be written in the working directory.

`app.builder.load_example_molecule`: this function copies a built-in example molecule to the working directory. Currently available names include `'fus'` and `'polyU'`. After calling this function, you can build the molecule by `protein_from_pdb` or `calvados_rna_from_pdb`.

```python
app.builder.load_example_molecule('fus')
mol = app.builder.protein_from_pdb('fus.pdb', rigid_from_plddt=True)
```

`app.builder.gen_elastic_network`: this function can generate an elastic network for a protein. The two parameters of this function are `molecule`, which is the target molecule, and `max_gap`, which is the maximum distance between two residues from two different structural regions. Default value of `max_gap` is 4. This function will return two values, the first one is the target molecule, and the second one is the extra force field generated by this function. You need to use this extra force field in your simulation to make the elastic network work.

Example:
```python
mol, extra_ff = app.builder.gen_elastic_network(
    molecule,
    max_gap=4
)
app.use(extra_ff, override=False)
```

`app.builder.linear_protein`: this function can generate a linear protein with a specific sequence. The required parameters of this function are `protein`, which is a protein sequence object.

`app.builder.spiral_protein`: this function can generate a protein with the shape of Archimedean spiral. The required parameters of this function is the same with `linear_protein`. Another optional parameter is `d`. This value indicates the radial distance. The default value is 0.38 nm.

`app.builder.stacked_protein`: this function can generate a compact, space-filling protein conformation from a sequence. Residues are placed group by group along a 3D Hilbert-like path. The required parameter is `protein`. The optional parameter `span` (default is 2) is the number of residues in each group.

`app.builder.molecule_from_curve`: this function places a sequence along a user-defined curve. The required parameters are `mol` (a sequence object) and `curve` (a function that maps residue index to `(x, y, z)`). Optional parameters include:
- `start_add`, `end_add`: extra mass added to the first and last beads. Default is 0.
- `cg`: coarse-graining label. Default is `MC`.
- `rename_map`: a dict that remaps residue names to atom types.
- `new_bond`: the bond type between neighboring beads. Default is `B-B`.

`app.builder.protein_from_pdb`
This function can automatically generate the protein object from a pdb file. The pdb file should be located in the working directory. The parameters of this function include
- `file_name`: the name of the pdb file.
- `ignoring_h`, whether to ignore the hydrogen atoms in the pdb file. Default is True.
- `cg`: Where should the coarse-grained bead locate. Available values are `alpha` and `mass_center`. Default is `mass_center`.
- `rigid_range`, this parameter indicates which residues should be regarded as rigid body. Reference `cg_molecule_from_pdb` for this parameter. Multiple rigid groups can be separated by `;`.
- `rigid_from_plddt`, determine whether a residue is rigid body by plddt score from alphafold2. Default is `False`.
- `threshold`, when `rigid_from_plddt` is set true, this parameter determines when plddt score is larger than the threshold, the residue can be regarded as rigid body, i.e., protein with structure. Default is 70.
- `max_gap`, reference `gen_elastic_network` for this parameter.

example:
```python
protein = app.builder.protein_from_pdb(
    'protein',
    rigid_from_plddt=True
)
```

`app.builder.update`: this function updates the force field attached to every atom in a molecule. It is useful after you change the force field of the app.

`app.builder.calvados_rna_from_pdb`: this function generates a two-bead Calvados RNA molecule from a pdb file. Backbone beads are placed at phosphorus atoms (`NMPB`), and side-chain beads are placed at N1 (U/A/C) or N9 (G). Angle terms along the backbone will also be added.

`app.builder.calvados_rna_from_sequence`: this function generates a two-bead Calvados RNA molecule from an RNA sequence and a curve function. The required parameters are `mol` (an `RNASequence` object) and `curve_func`. Optional parameters include `reverse` and `fallback_ref_vec`.

example:
```python
from ipamd.public.models.sequence import RNASequence
seq = RNASequence('polyU', 'U' * 30)
mol_rna = app.builder.calvados_rna_from_sequence(
    seq,
    lambda index: (0, 0, index * 0.59)
)
```

#### genbox
This part contains some functions for placing molecules in the simulation box. The basic usage is to run `box.gen_function(parameters)`. Available functions are listed below.

`place_molecule_randomly`: randomly place molecules in the simulation box. The required parameters of this function are listed bellow.
- `molecule`, the molecule object to add; n, the number of molecules.
- `threshold`, the minimum distance between two particles, default is 1.
- `max_tries`, maximum times the function will try to place a molecule, default is 5.
- `strict`, if this parameter is set to true, the program will continue trying until the  target `n` is reached, default is false.
- `allow_out_of_box`, if the molecules are allowed to be placed over the box boundary, default is true. If you want to generate box by `sub_box`, this option will be useful.

`place_molecule_periodically`: periodically place molecules in the simulation box. The required parameters of this function are listed bellow.
- `molecule`, the molecule object to add.
- `nx`, `ny`, `nz`, the number of molecules in x, y and z direction respectively.

`place_molecule_at`: place one molecule at a specific position. The required parameters of this function are listed bellow.
- `molecule`, the molecule object to add.
- `x`, `y`, `z`, the coordinate to place the molecule.

`enneutral`: place ions until the system is neutral. This function is only necessary when you run a simulation through Calvados-lc force field. No parameter is needed in this function.

`sub_box`: you can use this function to generate a simulation by combining existing boxes. The parameters are listed below.
- `sub_box`: the pre-existed simulation box to be placed in.
- `x`, `y`, `z`: the coordinate to place the `sub_box`.

#### analysis
Analysis module contains some operators for analyzing the simulation results. You can use these operators to get the information of the system, such as the temperature, the potential energy, the radius of gyration et al. These operators can be called by the `compute` function of a simulation box. Detailed information of each operator is listed as bellow.

`app.analysis.contact_map`: this function could calculate the contact map between two kinds of molecules. The function has three parameters.
- `type1`: the name of the first kind of molecule.
- `threshold`: if the distance between two CG beads is smaller than the threshold, the two CG beads are regarded as contact. Default is 4.
- `type2`: the name of the second kind of molecule. Default is '', which means the same as type1.

`app.analysis.momentum`
This function is an operator for computing the momentum of the system. No other parameters are needed. The function will return the momentum align x, y and z axe of the system.
Normally, the momentum of the system should be zero or very close to zero, None zero momentum indicates that the system is not correct, and you should rerun the simulation.

`app.analysis.temperature`
This function is an operator for computing the temperature of the system. No other parameters are needed. The function will return the temperature of the system.

`app.analysis.rg`
This function is an operator for computing the radius of gyration of the system. An optional parameter `type_` can be used to select a molecule type. The function will return the radius of gyration of the system.

`app.analysis.potential`
This function is an operator for computing the potential energy of the system. The simulation is needed as a parameter of this function. The function will return the potential energy of the system.
*Note:* before using this function, you must finish the simulation.

`app.analysis.rmsd`
This function is an operator for computing the root-mean-square deviation(RMSD) of the system. The function will return the RMSD of the system.

`app.analysis.rmsf`
This function is an operator for computing the root-mean-square fluctuation (RMSF) of each bead. It needs multiple frames. An optional parameter `n` selects the molecule index. Default is the first molecule.

`app.analysis.flory_monomer`
This function estimates the Flory scaling exponent of a single-chain polymer from a series of frames. Optional parameters include `type_` (molecule type) and `sub` (sampling density, default is 1). The molecule should contain at least 30 beads. A small number of frames may lead to inaccurate results.

`app.analysis.flory_multimer`
This function estimates the Flory scaling exponent from a mixture of chains with different lengths. It needs multiple frames.

`app.analysis.ripley`
This function computes Ripley's K function of the particle distribution. Parameters include `start_d`, `step`, `end_d` (distance range, default 1, 1, 10), `l` (whether to return the L function) and `ref` (whether to return the theoretical reference of a random distribution).

example:
```python
flory = box.compute(
    app.analysis.flory_monomer,
    title='flory',
    target_frame='50-100'
)
flory.print()
```

#### mdanalysis
`app.mdanalysis` is the recommended module for trajectory analysis. Compared with `app.analysis`, these operators return typed data objects (`Scalar`, `Vector`, `Matrix`) that can be processed by `app.data_process`. Most operators take a `box` and a `target_frame`. You can compute one frame directly, or average over a range of frames by `batch_compute`.

`app.mdanalysis.batch_compute`: compute an mdanalysis operator on a range of frames and return the average. The first parameter is the name of the operator. Other parameters will be passed to that operator.

```python
rg = app.mdanalysis.batch_compute(
    'rg_v1',
    box=box,
    target_frame='51-200'
)
app.data_process.print(rg)
```

`app.mdanalysis.rg_v1`: compute the radius of gyration. Optional parameter `target_molecule` can be a molecule name, or a name with a residue range such as `'fus:1-100'`. If it is not set, all molecules will be included. The unit of the result is nm.

`app.mdanalysis.density_align`: compute the mass density along one axis. Parameters include `target_molecule`, `direction` (`X`, `Y` or `Z`, default is `Z`) and `d` (bin size in nm, default is 1). The unit of the result is g/mL.

`app.mdanalysis.density_radius`: compute the radial mass density from an origin. Required parameter is `cutoff`. Optional parameters include `origin` (default is `(0, 0, 0)`), `d` and `target_molecule`.

`app.mdanalysis.density_box`: compute the mass density inside a sub-box. Parameters include `x0`, `y0`, `z0`, `lx`, `ly` and `lz`. If the lengths are 0, the full box size will be used. The unit of the result is g/cm3.

`app.mdanalysis.contact_map_v1`: compute a contact map. Parameters include `threshold` (default is 4), `type1`, `type2` and `mode`. If `type1` and `type2` are the same, `mode` can be `inter` (only intermolecular contacts), `intra` (only intramolecular contacts) or empty (exclude self contacts only).

`app.mdanalysis.contact_number`: compute the contact number of each residue from a contact map. Parameters are the same as `contact_map_v1`. You can also pass an existing contact map by `cm`.

`app.mdanalysis.net_charge`: compute the net charge of the system. The unit of the result is e.

`app.mdanalysis.slab_free_energy`: estimate the transfer free energy and the densities of the dense/dilute phases from a slab simulation. Parameters include `direction` (default is `Z`) and `d` (bin size). The result is a vector containing free energy (kJ/mol), dense-phase density and dilute-phase density.

`app.mdanalysis.momentum_v1`: compute the total momentum of the system. The result is a vector along x, y and z in kg m/s. No other parameters are needed.

`app.mdanalysis.temperature_v1`: compute the instantaneous temperature of the system. Rigid groups contribute 6 degrees of freedom and flexible beads contribute 3. The unit of the result is K.

`app.mdanalysis.potential_v1`: read the potential energy from the simulation log. The required parameter is `simulation`. The unit of the result is kJ/mol. The simulation must be finished first.

`app.mdanalysis.ripley_v1`: compute Ripley's K function of the particle distribution. Parameters include `start_d`, `step`, `end_d` (distance range, default 1, 1, 10), `l` (whether to return the L function), `ref` (whether to return the theoretical reference of a random distribution) and `target_molecule`.

`app.mdanalysis.rmsd_v1`: compute RMSD relative to the first frame of a trajectory. Required parameters are `box` and `target_frame` (a frame range). Optional parameter `n` selects a molecule index; `target_molecule` selects by name or `'name@1-100'`. The unit of the result is nm.

`app.mdanalysis.rmsf_v1`: compute the RMSF of each bead relative to the first frame. Parameters are the same as `rmsd_v1`. At least two frames are required. The unit of the result is nm.

`app.mdanalysis.flory_monomer_v1`: estimate the Flory scaling exponent of a single-chain polymer from a trajectory. Parameters include `box`, `target_frame`, `target_molecule` (or the legacy parameter `type_`) and `sub` (sampling density, default is 1). The molecule should contain at least 30 beads.

`app.mdanalysis.flory_multimer_v1`: estimate the Flory scaling exponent from a mixture of chains with different lengths. Parameters include `box`, `target_frame` and `target_molecule`.

`app.mdanalysis.rdf`: compute the radial distribution function g(r). Parameters include `type1`, `type2`, `dr` (bin size, default 0.1 nm), `r_max` (default is half of the shortest box edge) and `mode` (`inter` / `intra`, used when both selections are the same). You can average over frames with `batch_compute`.

`app.mdanalysis.msd`: compute the mean squared displacement. Required parameters are `box` and `target_frame`. Molecule centers of mass are tracked by default; set `per_particle=True` to use beads. Optional parameter `dt` is the time interval between consecutive selected frames, and `n` / `target_molecule` select molecules. The unit of the result is nm².

example:
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
Simulation module contains the interaction potentials and the integrators. 
The first part of this module is the interaction potentials. The interaction potentials are used to describe the interactions between the CG beads. Generally, the parameters will be set automatically depending on the force field you use, and you don't need to set them manually. The available interaction potentials are listed below.
- `app.simulation.ah`:  Ashbaugh-Hatch potential.
- `app.simulation.bond_harmonic`: harmonic bond potential.
- `app.simulation.angle_harmonic`: harmonic angle potential.
- `app.simulation.debye`: debye potential.
- `app.simulation.pppm`: pppm potential.
- `app.simulation.stacking`: stacking potential used by Calvados RNA (side-chain stacking).
- `app.simulation.centripetal`: centripetal force used to prepare droplet or slab configurations.
- `app.simulation.wf`: Wang-Frenkel potential. This potential is provided by the `mpipi_forcefield` plugin pack.

The second part of this module is the integrators. You can choose which integrator to use in your simulation when creating your simulation. The parameters of the integrator are generated by the `env` of the app, and you don't need to set them manually. The available integrators are listed below.
- `app.simulation.langevin_nvt`: Langevin NVT integrator.
- `app.simulation.em`: NVE integrator, could be used for energy minimization.
- `app.simulation.noose_hover_nvt`: Noose Hover NVT integrator.
- `app.simulation.npt_z`: NPTMTK integrator. There is pressure coupling only in z direction.
- `app.simulation.andersen_npt`: Andersen NPT integrator, with pressure coupling in all directions.

#### sakuanna
Sakuanna is a sequence analysis module. It contains some functions for analyzing the sequence of the protein. Detailed information of each function is listed bellow.

You can create a sequence object directly:
```python
from ipamd.public.models.sequence import ProteinSequence
sequence = ProteinSequence(
    name='ELN',
    sequence='VPGAGVPGAGVPGAG'
)
```

`app.sakuanna.sequence_from_fasta`: load sequence objects from a fasta file in the working directory. The required parameter is `fasta_path`. The optional parameter `mol` can be `'protein'`, `'dna'` or `'rna'`. Default is `'protein'`. If the fasta file contains one sequence, a single sequence object will be returned. If it contains more than one sequence, a tuple of sequence objects will be returned.

example:
```python
seq = app.sakuanna.sequence_from_fasta('proteins.fasta', mol='protein')
```

`app.sakuanna.to_fasta`: write a sequence object to a fasta file. Optional parameters include `filename` and `comment`. If `filename` is not set, the file will be named after the sequence name.

`app.sakuanna.pretier`: this function can print the protein sequence in a prettier way. The function has three parameters.
- `protein`: the protein sequence object to print.
- `word3`: whether to print the sequence in three-letter code. Default is True.
- `ter`: whether to print the terminal group. Default is True.
example:
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

`app.sakuanna.statistic`: this function can calculate the ratio of specific amino acids in the protein sequence. The function has three parameters.
- `protein`: the protein sequence object to analyze.
- `targets`: the target amino acids or amino acid classes to analyze. The value could be a list of amino acids and amino acid classes. Available amino acid classes include `aromatic`, `positive`, `negative`, `charged`, `polar`, `nonpolar`. 
- `format_`: the format of the return value. Available values are `ratio` and `count`. If the value is `ratio`, the function will return the ratio of the target amino acids in the protein sequence. If the value is `count`, the function will return the count of the target amino acids in the protein sequence. Default is `ratio`.

example:
```python
res = app.sakuanna.statistic(
    sequence, 
    targets=['A', 'charged'],
    format_='ratio'
)
app.data_process.plot(res)
```

`app.sakuanna.tag`: this function tags residues in a protein sequence by matching sub-sequences. The parameter `tags` is a dict. The key is the tag name, and the value is a list of sub-sequences. You can set one tag to `'rest'` to label unmatched residues. The returned value is a `Vector` that can be plotted by `app.data_process.plot`.

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
The output data of mdanalysis and sakuanna is provided as typed objects (`Scalar`, `Vector`, `Matrix`, `Ratio`, `Distribution`, `PointSet`, `String`). You can process these objects by `app.data_process`.

- `print`: print the data in a prettier way. Optional parameter `precision` (default is 3).
- `plot`: show the plot of the data. Optional parameters include `style`, `save_figure` and `appearance` (`line`, `bar`, `heatmap` or `discrete_heatmap` for vectors).
- `flatten`: reduce a vector to a scalar, or a matrix to a vector. Parameter `by` can be `average` (default) or `sum`. For a matrix, `axis` is required.
- `average`: average several data objects of the same type.
- `sum`: sum several data objects of the same type.
- `normalize`: normalize the data. For a ratio, the values will be scaled so that the sum is 1. For a vector or a matrix, min-max normalization will be used. For a point set, you can set `target` to `'x'`, `'y'` or both.
- `to_csv`: save the data to a csv file. The required parameter is `filename`.

The older `box.compute` path still returns an `AnalysisResult` object. You can call `print`, `distribution`, `flatten`, `merge`, `normalize`, `plot` and `save` on that object. For new code, `app.mdanalysis` plus `app.data_process` is recommended.

### OmicsLoader
OmicsLoader is a module for loading omics data. This module will give you an easier experience in batch protein analysis. It can load data from various formats, such as fasta, pdb et al. You can use the functions provided by this module to load the data automatically.
To use this module, you need first to import OmicsLoader and batch_run. The first class is for reading the data you provide, while the second function is to dispatch the jobs of different proteins to different GPUs.
```python
from ipamd import OmicsLoader, batch_run
```
You first need to create an instance of OmicsLoader with the path of your data directory. Then you need to create a function to process the data. The function should have two parameters, the first parameter is the app, and the second one is data. These two parameters will be passed by the batch_run function. The data parameter is a dictionary, which contains the name, sequence, type and path of the data. `batch_run` also accepts an optional `gpus` parameter, which is a range string of GPU ids, such as `'0-3'`. If it is not set, all available GPUs will be used. A detailed example is shown below.
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

### configuration
The config file of this package is located at `~/.config/ipamd/config.json`. The configuration file will be created when you first run this package. You can modify the configuration file to change the behavior of the package, such as the path of the data file, the path of extra plugins, and the path of the output file, et al. The detailed configuration items are listed below.

| Item         | Description                           | Available values                                                                                         |
|--------------|---------------------------------------|----------------------------------------------------------------------------------------------------------|
| result_dir   | where to write the simulation results | the target directory, or set to `$decentralized` if you don't want to get your output managed by IPAMD   |
| output_level | how much log is printed               | 0: verbose, 1: info, 2: warning, 3: error                                                                |
| auto_load    | whether to load the plugins automatically | True or False                                                                                            |
| default_ff  | the default force field               | the name of the force field, should be one of the available force fields printed by `app.available_ff()` |
| external_plugin_dir | extra plugin directories          | a list of target directories                                                                             |
| sakuanna_plugin_dir | the path of extra sakuanna plugins | a list of target directories                                                                             |
| simulation_plugin_dir | the path of extra simulation plugins | a list of target directories                                                                             |
| builder_plugin_dir | the path of extra builder plugins | a list of target directories                                                                             |
| analyse_plugin_dir | the path of extra analyse plugins | a list of target directories                                                                             |

### examples
Some example code is provided in the demo folder. You can run these examples to get familiar with this package. There are some subfolders in the demo folder, each subfolder contains an example code and the input data needed. You can also list the official examples by `ipamd -s`. The meaning of the examples could refer to the paper of this package(doi: 10.1021/acs.jctc.5c00147).

| Example | Description |
|---------|-------------|
| single_chain_simulation | Run a single chain simulation |
| sequence_analysis | Perform analysis on protein sequences |
| basic_example | Droplet simulation with contact, RDF, MSD and density analysis |
| slab_simulation | Slab simulation with density and free-energy analysis |
| RNA+Protein | Simulate RNA and protein systems |
| custom_system | Create custom molecular systems |
| omics | Batch analysis with OmicsLoader |
| calvados_rna | Two-bead RNA simulation |

### CLI
IPAMD provides a command-line tool `ipamd` for managing plugin packs and listing examples.

```bash
ipamd -s                  # show official examples
ipamd -l                  # list installed plugin packs
ipamd -i mpipi_forcefield # install a plugin pack (from local path or from GitHub)
ipamd -r mpipi_forcefield # remove a plugin pack
ipamd -p plugin_dir       # pack a plugin directory into a zip file
ipamd -v                  # show version
```

If the argument of `-i` is not a local file or directory, IPAMD will try to download `{name}.zip` from the official plugin pack repository.

## Plugin-packs
Some functions are provided as optional plugin packs. You can install them by `ipamd -i <pack_name>`. After installation, the corresponding plugins and data files will be copied into the IPAMD installation directory.

### mpipi_forcefield
This pack provides the Mpipi force field. After installation, you can load the force field by `app.use('mpipi')`. This pack also adds a Wang-Frenkel potential (`app.simulation.wf`). Proteins must be processed by `mpipi_chtype` before the simulation.

`app.builder.mpipi_chtype`: this function converts the atom types of structured residues so that they can be recognized by the Mpipi force field. For every bead that belongs to a rigid group, a suffix `0` will be added to the atom type if it is not already present. The only required parameter is `molecule`. The function returns the modified molecule.

example:
```python
app.use('mpipi')
mol = app.builder.protein_from_pdb('protein.pdb', rigid_from_plddt=True)
mol = app.builder.mpipi_chtype(mol)
```

### calvados_rna
This pack is in development. It provides supporting files for Calvados RNA simulations, and is compatible with the Calvados2 force field. After installation, you can merge the RNA force field into the current force field by `app.use('Calvados_rna', override=False)`.

`app.builder.calvados_rna_from_pdb`: this function generates a two-bead Calvados RNA molecule from a pdb file. Backbone beads are placed at phosphorus atoms (`NMPB`), and side-chain beads are placed at N1 (U/A/C) or N9 (G). Angle terms along the backbone will also be added. The only required parameter is `file_name`.

`app.builder.calvados_rna_from_sequence`: this function generates a two-bead Calvados RNA molecule from an RNA sequence and a curve function. The required parameters are `mol` (an `RNASequence` object) and `curve_func`. Optional parameters include `reverse` and `fallback_ref_vec`.

This pack also provides `app.simulation.stacking` and `app.simulation.angle_harmonic`.

example:
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
This pack adds a plotting backend based on matplotlib. After installation, you can plot analysis results by `app.data_process.plot`.

`app.data_process.plot`: this function shows a plot of the data. The required parameter is `data`. Optional parameters include:
- `style`: a matplotlib rc style dict.
- `save_figure`: whether to save the figure as `{title}.png`. Default is False.
- `appearance`: the appearance of a vector plot. Available values are `line`, `bar`, `heatmap` and `discrete_heatmap`. Default is `line`.

Scalar and string data cannot be plotted. In that case, the data will be printed instead.

### to_pandas
This pack converts analysis results to a pandas DataFrame. It depends on `pandas`.

`app.data_process.to_df`: this function converts a data object to a `pandas.DataFrame`. The only required parameter is `data`. Scalar, vector, matrix, ratio, distribution and point-set data are supported.

example:
```python
df = app.data_process.to_df(rg)
df.to_csv('rg.csv')
```

### cif_format_support
This pack saves a configuration to CIF format. It depends on `biopython`.

`box.to_cif`: this function saves the current configuration to a .cif file. The required parameter is `filename`. An optional parameter `ignoring_pbc` (default is True) controls whether periodic images are unwrapped before writing.

### cg2all_integrate
This pack converts a CG protein to an all-atom protein. After installation, please restart your shell.

`box.to_aa`: this function converts a CG protein to an all-atom protein. The only parameter is the output filename of the all-atom protein.
*Note:* This function uses cg2all as its computing backend, so you must ensure that cg2all is installed on your computer and that `convert_cg2all` can be called in the terminal.

### alphafpld2_integrate
This pack integrates AlphaFold2/ColabFold. After installation, please restart your shell.

`app.builder.af2`: this function can call the alphafold2 (colabfold) installed on your computer to generate the initial structure of the protein. After running this function, a pdb file of the protein will be written in the working directory. The only required parameter is `protein_sequence`. This parameter is a protein sequence object.
*Note:* You should make sure that the `colabfold_batch` command is available in your system.

### slab_simulation
This pack supports slab simulations. It depends on `scipy`. After installation, a `slab_prepare` force field will also be available.

`app.mdanalysis.slab_free_energy`: this function estimates the transfer free energy and the densities of the dense/dilute phases from a slab simulation. The required parameters are `box` and `target_frame`. Optional parameters include `direction` (default is `Z`) and `d` (bin size). The result is a vector containing free energy (kJ/mol), dense-phase density and dilute-phase density.

### sequence_alignment
This pack provides sequence alignment tools.

`app.sakuanna.sequence_align`: this function aligns two sequences. The required parameters are `seq1` and `seq2`. Optional parameters include:
- `algorithm`: the alignment algorithm. Available values are `needleman-wunsch`, `smith-waterman` and `diff`. Default is `needleman-wunsch`.
- `match_score`: the score for a match. Default is 1.
- `mismatch_score`: the score for a mismatch. Default is -1.
- `gap_score`: the score for a gap. Default is -1.

The function returns a list of three data objects: the aligned reference sequence, the aligned target sequence, and the alignment score.

`app.data_process.print_diff`: this function prints the difference between two aligned sequences. The required parameters are `ref` and `target`, which should be the string objects returned by `sequence_align`. Optional parameters include `print_ref` (default is True) and `seperate_with_bracket` (default is True).

example:
```python
ref, target, score = app.sakuanna.sequence_align(seq1, seq2)
app.data_process.print_diff(ref, target)
```

### adv_data_processing
This pack provides advanced data processing tools. It depends on `scipy`.

`app.data_process.gaussian`: this function applies a Gaussian filter to the data. The only required parameter is `data`. The function returns a new data object with smoothed values.

## Second-time-development
IPAMD is a plugin-based software, and the developer could easily add new functions to this package by writing a plugin. The detailed instruction on writing a plugin is shown in this section.

### Plugin
A plugin is a Python module that contains one or more classes or functions. There should be a function named `func`, and this func will be detected and loaded by IPAMD automatically. There is a configuration variable for each plugin. This variable looks like:
```python
configure = {
    'type': 'function',
    "schema": 'schema_name',
    "apply": ['applied_value']
}
```
The `type` variable is deprecated and will be removed in future versions. The only available value is `function`. The `apply` variable indicates which value should be passed to the plugin. Your function should also have a parameter that has the same name with the applied value. The `schema` variable indicates a combination of applied value. But the schema variable is not necessary in most cases and not recommended for developers.

Newer plugins can also declare `resource` (injected resources such as `ff` and `persistency_dir`) and `alias` (the name used when calling the plugin).

You can also call another plugin in your plugin. To do this, you need to import the `PluginBase` class from `ipamd.public.utils`, and then use `PluginBase.call('plugin_name', parameters)` to call the target plugin. The `plugin_name` is the name of the plugin you want to call, and the `parameters` is a dict that contains the parameters needed by the target plugin. *Note: This feature is beta and may be not stable.*

A plugin pack is a directory containing a `meta.json` file, plugin scripts and optional data files. You can pack it by `ipamd -p <dir>` and install it by `ipamd -i <zip>`.

## Supporting
Contact me at [email](mailto:liuxiaoyang_Q@outlook.com) if you have any questions or suggestions. You can also submit an issue on the GitHub page. I will reply to you as soon as possible.
