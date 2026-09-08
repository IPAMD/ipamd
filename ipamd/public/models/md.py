"""
class definition for molecular dynamics
"""
import copy
import math
from dataclasses import make_dataclass
import sys
from typing import List
import numpy as np
from ipamd.public.utils.parser import range_to_list, value_of
from ipamd.public.utils.output import error
from ipamd.public import shared_data
from ipamd.public.utils.plugin_manager import PluginBase
from ipamd.public.models.common import AnalysisResult

ForceField = make_dataclass(
    'ForceField', 
    [('atom_definition', dict), ('ff_param', dict)]
)

class Atom:
    """
    Atom class
    """
    def __init__(self, velocity, atom_type, ff, mass=None, charge=None):
        """
        Initialize an atom
        
        :param velocity: velocity of the atom
        :param atom_type: type of the atom
        :param ff: force field
        :param mass: mass of the atom
        :param charge: charge of the atom
        :return: None
        """
        self.__ff = ff
        if atom_type not in self.__ff.atom_definition.keys():
            error('Invalid atom type ' + atom_type)
            sys.exit(1)
        self.__type = atom_type
        self.__v = {
            'velocity': velocity,
            'mass': mass,
            'charge': charge,
        }

    def get(self, target):
        """
        Get the property of the atom
        
        :param target: property of the atom
        :return: property of the atom
        """
        if target == 'type':
            return self.__type
        elif target in self.__v and self.__v[target] is not None:
            return self.__v[target]
        else:
            try:
                return self.__ff.atom_definition[self.__type][target]
            except KeyError:
                error("Invalid atom property")
                sys.exit(1)

    def set(self, target, value):
        """
        Set the property of the atom
        
        :param target: property of the atom
        :param value: value of the property
        :return: None
        """
        if target == 'type':
            self.__type = value
        elif target in self.__v:
            self.__v[target] = value
        else:
            error("Invalid atom property")
            sys.exit(1)

    def update(self, ff):
        """
        Update the force field of the atom
        
        :param ff: force field
        :return: None
        """
        self.__ff = ff

class Residue:
    """
    Residue class
    """
    def __init__(self, type_):
        """
        Initialize a residue
        
        :param type_: type of the residue
        :return: None
        """
        self.atom_list = []
        self.type_ = type_

class Molecule:
    """
    Molecule class
    """
    def __init__(self, type_, cg):
        """
        Initialize a molecule
        
        :param type_: type of the molecule
        :param cg: corse graining method
        :return: None
        """
        self.type_ = type_
        self.cg = cg
        self.atoms = []
        self.angles = []
        self.__activated_index = -1
        self.length = 0

    def add_atom(
            self,
            atom,
            coordinate: List[float],
            bond=None,
            rigid_group=-1
    ):
        """
        Add an atom to the molecule
        
        :param atom: atom to be added
        :param coordinate: coordinate of the atom
        :param bond: bond between the atom and the previous atom
        :param rigid_group: rigid group of the atom
        :return: None
        """
        self.atoms.append({
            "prototype": atom,
            "links": [],
            "offset": coordinate,
            'rigid_group': rigid_group,
        })
        self.length += 1
        if self.__activated_index == -1:
            self.__activated_index = 0
        else:
            previous_index = self.__activated_index
            self.__activated_index += 1
            if bond is not None:
                self.link(previous_index, self.__activated_index, bond)

    def link(self, atom1, atom2, type_=''):
        """
        Link two atoms in the molecule
        
        :param atom1: index of the first atom
        :param atom2: index of the second atom
        :param type_: type of the bond
        :return: None
        """
        atom1, atom2 = min(atom1, atom2), max(atom1, atom2)
        type1 = self.atoms[atom1]['prototype'].get('type')
        type2 = self.atoms[atom2]['prototype'].get('type')
        self.atoms[atom1]['links'].append(
            {
                'to': atom2,
                'type': f'{type1}-{type2}' if type_ == '' else type_
            }
        )

    def add_angle(self, atom1, atom2, atom3, type_=''):
        """
        Add an angle to the molecule
        
        :param atom1: index of the first atom
        :param atom2: index of the second atom
        :param atom3: index of the third atom
        :param type_: type of the angle
        :return: None
        """
        atom1, atom2, atom3 = sorted([atom1, atom2, atom3])
        type1 = self.atoms[atom1]['prototype'].get('type')
        type2 = self.atoms[atom2]['prototype'].get('type')
        type3 = self.atoms[atom3]['prototype'].get('type')
        self.angles.append({
            'atom1': atom1,
            'atom2': atom2,
            'atom3': atom3,
            'type': f'{type1}-{type2}-{type3}' if type_ == '' else type_
        })

    def transform(self, m, by='zero'):
        """
        Transform the molecule
        
        :param m: transformation matrix
        :param by: center or zero
        :return: None
        """
        center = np.array([0.0, 0.0, 0.0])
        if by == 'center':
            for atom in self.atoms:
                offset = np.array(atom['offset'])
                center += offset
            center /= len(self.atoms)
        for atom in self.atoms:
            old_pos = atom['offset']
            x = old_pos[0] - center[0]
            y = old_pos[1] - center[1]
            z = old_pos[2] - center[2]
            coord = np.array([x, y, z])

            new_coord = np.dot(coord, m)
            new_coord += center
            atom['offset'] = (new_coord[0], new_coord[1], new_coord[2])

    def rotate(self, theta_x, theta_y, theta_z, unit='degree', by='zero'):
        """
        Rotate the molecule
        
        :param theta_x: angle of rotation around x axis
        :param theta_y: angle of rotation around y axis
        :param theta_z: angle of rotation around z axis
        :param unit: degree of the angle
        :param by: center or zero
        :return: None
        """
        if unit == 'degree':
            theta_x = math.radians(theta_x)
            theta_y = math.radians(theta_y)
            theta_z = math.radians(theta_z)
        c1 = math.cos(theta_x)
        s1 = math.sin(theta_x)
        c2 = math.cos(theta_y)
        s2 = math.sin(theta_y)
        c3 = math.cos(theta_z)
        s3 = math.sin(theta_z)

        center = np.array([0.0, 0.0, 0.0])
        if by == 'center':
            for atom in self.atoms:
                offset = np.array(atom['offset'])
                center += offset
            center /= len(self.atoms)

        for atom in self.atoms:
            old_pos = atom['offset']
            x = old_pos[0] - center[0]
            y = old_pos[1] - center[1]
            z = old_pos[2] - center[2]
            new_x = c2 * c3 * x - c2 * s3 * y + s2 * z
            new_y = (c1 * s3 + c3 * s1 * s2) * x + (c1 * c3 - s1 * s2 * s3) * y - c2 * s1 * z
            new_z = (s1 * s3 - c1 * c3 * s2) * x + (c3 * s1 + c1 * s2 * s3) * y + c1 * c2 * z
            new_x += center[0]
            new_y += center[1]
            new_z += center[2]
            atom['offset'] = (new_x, new_y, new_z)

    def move(self, t):
        """
        Move the molecule
        
        :param t: translation vector
        :return: None
        """
        for atom in self.atoms:
            atom['offset'] = (
                atom['offset'][0] + t[0],
                atom['offset'][1] + t[1],
                atom['offset'][2] + t[2]
            )

    def __bonds(self):
        """
        Get the bonds in the molecule
        
        :return: list of bonds
        """
        bond_list = []
        for i, atom1 in enumerate(self.atoms):
            neighbor_list = atom1['links']
            for j in neighbor_list:
                bond_list.append({
                    'type': j['type'],
                    'atom1': i,
                    'atom2': j['to']
                })
        return bond_list

    def distance(self, index1, index2):
        """
        Get the distance between two atoms
        
        :param index1: index of the first atom
        :param index2: index of the second atom
        :return: distance between the two atoms
        """
        d = math.sqrt(
            (self.atoms[index1]['offset'][0] - self.atoms[index2]['offset'][0]) ** 2 +
            (self.atoms[index1]['offset'][1] - self.atoms[index2]['offset'][1]) ** 2 +
            (self.atoms[index1]['offset'][2] - self.atoms[index2]['offset'][2]) ** 2
        )
        return d

    def properties(self):
        """
        Get the properties of the molecule
        
        :return: dictionary of properties
        """
        result = {
            'type': [],
            'mass': [],
            'charge': [],
            'position': [],
            'bond': [],
            'velocity': [],
            'rigid_group': [],
            'angle': self.angles
        }
        for atom in self.atoms:
            result['type'].append(atom['prototype'].get('type'))
            result['position'].append(atom['offset'])
            result['velocity'].append(atom['prototype'].get('velocity'))
            result['mass'].append(atom['prototype'].get('mass'))
            result['charge'].append(atom['prototype'].get('charge'))
            result['rigid_group'].append(atom['rigid_group'])
        result['bond'] = self.__bonds()
        return result


    @staticmethod
    def mass_center(prop):
        x_mass = 0
        y_mass = 0
        z_mass = 0
        total_mass = 0
        for index in range(len(prop['mass'])):
            mass = prop['mass'][index]
            position = prop['position'][index]
            x_mass += position[0] * mass
            y_mass += position[1] * mass
            z_mass += position[2] * mass
            total_mass += mass
        return np.array([x_mass / total_mass, y_mass / total_mass, z_mass / total_mass])

    @staticmethod
    def align(mol_prop, target_prop):
        P = np.array(mol_prop['position'])
        Q = np.array(target_prop['position'])
        center_P = np.mean(P, axis=0)
        center_Q = np.mean(Q, axis=0)

        P_centered = P - center_P
        Q_centered = Q - center_Q

        H = np.dot(P_centered.T, Q_centered)
        U, _, Vt = np.linalg.svd(H)
        R = np.dot(U, Vt)
        if np.linalg.det(R) < 0:
            Vt[2, :] *= -1
            R = np.dot(U, Vt)

        t = center_Q - center_P

        return t, R

class Environment:
    """
    Environment class
    """
    @staticmethod
    def __epsilon(t):
        """
        Get the dielectric constant of the environment
        
        :param t: temperature
        :return: dielectric constant
        """
        return 5321 / t + 233.76 - 0.9297 * t + 1.417e-3 * t ** 2 - 8.292e-7 * t ** 3

    def __init__(self, name=None):
        if name=='pure water' or name is None:
            self.values = {
                'ionic_strength': 1e-7,
                'temperature': 298,
                'ph': 7.0,
                'pressure': 1,
                'epsilon': Environment.__epsilon(298)
            }
        elif name=='normal saline':
            self.values = {
                'ionic_strength': 0.154,
                'temperature': 298,
                'ph': 7.0,
                'pressure': 1,
                'epsilon': Environment.__epsilon(298)
            }
        elif name == 'cytoplasm':
            self.values = {
                'ionic_strength': 0.154,
                'temperature': 310,
                'ph': 7.2,
                'pressure': 1,
                'epsilon': Environment.__epsilon(310)
            }
        else:
            error('unknown system')
            quit()
        self.values['epsilon'] = Environment.__epsilon(self.values['temperature'])

    def set_ph(self, ph):
        self.values['ph'] = ph

    def set_ionic_strength(self, ionic_strength):
        self.values['ionic_strength'] = ionic_strength

    def set_temperature(self, temperature, unit="K"):
        if unit == "K":
            self.values['temperature'] = temperature
        elif unit == "C":
            self.values['temperature'] = temperature + 273.15
        else:
            error("Invalid temperature unit")
            raise ValueError("Invalid temperature unit")
        self.values['epsilon'] = Environment.__epsilon(self.values['temperature'])

    def set_pressure(self, pressure, unit="atm"):
        if unit == "atm":
            self.values['pressure'] = pressure
        elif unit == "Pa":
            self.values['pressure'] = pressure / 101325
        else:
            error("Invalid pressure unit")
            raise ValueError("Invalid pressure unit")

class Frame:
    def __init__(self, outer, no):
        self.box = outer
        self.x = self.box.x
        self.y = self.box.y
        self.z = self.box.z
        self.no = no
        self.molecules = []

    def add_molecule(self, molecule, offset=(0, 0, 0)):
        """
        Add a molecule to the frame
        
        :param molecule: molecule to add
        :param offset: offset of the molecule
        :return: None
        """
        self.molecules.append({
            'prototype': copy.deepcopy(molecule),
            'offset': offset
        })

    def molecule(self, index):
        return self.molecules[index]['molecule']

    def move_molecule(self, index, offset):
        self.molecules[index]['offset'][0] += offset[0]
        self.molecules[index]['offset'][1] += offset[1]
        self.molecules[index]['offset'][2] += offset[2]

    def set_size(self, x, y, z):
        if self.x == 0:
            self.box.x = x
        self.x = x
        if self.y == 0:
            self.box.y = y
        self.y = y
        if self.z == 0:
            self.box.z = z
        self.z = z

    def properties(self, ignoring_image=False, filter=None):
        x = self.box.x
        y = self.box.y
        z = self.box.z
        env = self.box.env.values
        if filter is None:
            filter = range(len(self.molecules))
        else:
            filter = range_to_list(filter)
        result = {
            'size': (x, y, z),
            'no': self.no,
            'molecules': []
        }
        for index, molecule in enumerate(self.molecules):
            if index not in filter:
                continue
            prop = molecule['prototype'].properties()
            molecule_prop = {
                'CG': molecule['prototype'].cg,
                'molecule_type': molecule['prototype'].type_,
                'n_atoms': molecule['prototype'].length,
                'n_bonds': len(prop['bond']),
                'n_angles': len(prop['angle']),
                'bonds': prop['bond'],
                'type': prop['type'],
                'mass': [value_of(mass, env) for mass in prop['mass']],
                'charge': [value_of(charge, env) for charge in prop['charge']],
                'position': [],
                'image': [],
                'rigid_group': prop['rigid_group'],
                'velocity': prop['velocity'],
                'angles': prop['angle']
            }
            if not ignoring_image:
                for position in prop['position']:
                    image = [0, 0, 0]
                    new_position = [0, 0, 0]
                    image[0], new_position[0] = divmod(position[0] + molecule['offset'][0], x)
                    image[1], new_position[1] = divmod(position[1] + molecule['offset'][1], y)
                    image[2], new_position[2] = divmod(position[2] + molecule['offset'][2], z)
                    new_position[0] -= x / 2
                    new_position[1] -= y / 2
                    new_position[2] -= z / 2
                    molecule_prop['position'].append(new_position)
                    molecule_prop['image'].append(image)
            else:
                for position in prop['position']:
                    molecule_prop['position'].append([
                        position[0] - x / 2 + molecule['offset'][0],
                        position[1] - y / 2 + molecule['offset'][1],
                        position[2] - z / 2 + molecule['offset'][2]
                    ])
                    molecule_prop['image'].append([0, 0, 0])
            result['molecules'].append(molecule_prop)
        return result

class Box(PluginBase):
    """
    Simulation Box Class
    """
    def __init__(self, x, y, z, force_field, persistency_dir):
        super().__init__([
            shared_data.module_installation_dir + '/plugins/converter',
            shared_data.module_installation_dir + '/plugins/genbox'
        ])

        self.persistency_dir = persistency_dir
        self.x = x
        self.y = y
        self.z = z
        self.force_field = force_field
        self.env = Environment()

        self.__current_frame = 0
        self.frames = []
        self.def_schema(
            'frame',
            {
                'frame': lambda: self.frames[self.__current_frame]
            }
        )
        self.add_resource('ff', self.force_field)
        self.add_resource('persistency_dir', self.persistency_dir)
        if shared_data.config.get('auto_load'):
            self.load_all()

    def set_box_size(self, x, y, z):
        self.x = x
        self.y = y
        self.z = z

    def in_solvent(self, environment):
        if isinstance(environment, str):
            self.env = Environment(environment)
        else:
            self.env = environment
        return self

    def new_frame(self):
        self.__current_frame = len(self.frames)
        self.frames.append(Frame(self, self.__current_frame))
        return self

    def compute(self, method, title='', target_frame='', **kwargs):
        target_frame = range_to_list(target_frame)
        result = AnalysisResult(
            title=title,
            type_=AnalysisResult.Type.UNKNOWN,
            data={}
        )
        if len(target_frame) == 0:
            frame = self.frames[self.__current_frame]

            res = method(frame, **kwargs)
            result.parse(res)
            return result
        else:
            if "target" in method.attr and method.attr['target'] == 'multi':
                target_frames = []
                for i in target_frame:
                    target_frames.append(self.frame(i).current_frame())
                res = method(target_frames, **kwargs)
                result.parse(res)
            else:
                result.name_as(str(target_frame[0]))
                res = method(self.frame(target_frame[0]).current_frame(), **kwargs)
                result.parse(res)
                if result.type == AnalysisResult.Type.VECTOR:
                    result = result.flatten(new_name=str(target_frame[0]))
                for i in target_frame[1:]:
                    new_result = AnalysisResult(
                        title=str(i),
                        type_=AnalysisResult.Type.UNKNOWN,
                        data={}
                    )
                    res = method(self.frame(i).current_frame(), **kwargs)
                    new_result.parse(res)
                    if new_result.type == AnalysisResult.Type.VECTOR:
                        new_result = new_result.flatten(new_name=str(i))
                    result = result.merge(other=new_result)
            return result.name_as(title)

    def frame(self, index):
        self.__current_frame = index
        return self

    def current_frame(self):
        return self.frames[self.__current_frame]

    def clean(self):
        """
        clean all loaded frames
        """
        self.__current_frame = -1
        self.frames = []
        return self

class Unit:
    """
    scale factor for units
    """
    class TimeScale:
        """
        time scale factor
        """
        fs = 0.001
        ps = 1
        ns = 1000
        us = 1000000

    class LengthScale:
        """
        length scale factor
        """
        nm = 1
        A = 0.1

    class Unitless:
        """
        unitless scale factor
        """
        rad = 1
        degree = math.pi / 180
