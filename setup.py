"""
IPAMD Setup Script
"""
import os
from setuptools import setup, find_packages
home_dir = os.path.expanduser('~')
config_dir = os.path.join(home_dir, '.config/ipamd')

with open('README.md', 'r', encoding='utf-8') as f:
    long_description = f.read()
setup(
    name='ipamd',
    python_requires='>=3.12',
    license='LGPL-3.0',
    long_description=long_description,
    long_description_content_type="text/markdown",
    entry_points={
        'console_scripts': [
            'ipamd = ipamd.cli:main'
        ],
    },
    description='A Python package for MD simulations and analysis of biomolecules',
    version='0.0.31',
    include_package_data=True,
    package_data={
        'ipamd': [
            "data/**", 
            "app/lib/*",
            "plugins/**",
            "meta/*"
        ]
    },
    install_requires=[
        'numpy',
        'rich',
        'numba',
        'numpy',
        'pybioseq',
        'periodictable'
    ],
    packages=find_packages(),
    author='Xiaoyang Liu',
    author_email='liuxiaoyang_Q@outlook.com'
)
