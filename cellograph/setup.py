# Copyright (C) 2023 Purvis Lab, CompCy Lab, University of North Carolina at Chapel Hill

import os
from setuptools import setup, find_packages

install_requires = [
    "tensorflow>=2.7.0",
    "tensorflow.keras>=2.6.0",
    "phate>=1.0.7",
    "scanpy>=1.9.1",
    "scvelo>=0.2.4",
    "pandas>=1.2.2",
    "anndata>=0.8.0",
    "numpy>=1.20.0",
    "scipy>=1.6.1",
    "scikit-learn",
    "seaborn>=0.11.2"
    "spektral"
]

version = '0.0.1' 

readme = open("README.md").read()

setup(
    name="cellograph",
    version=version,
    description="cellograph",
    author="Jamshaid Shahir, Purvis Lab, University of North Carolina at Chapel Hill",
    author_email="jashahir@live.unc.edu",
    packages=find_packages(),
    license="MIT License - See LICENSE file",
    python_requires=">=3.8",
    long_description=readme,
    long_description_content_type="text/markdown",
    url="https://github.com/jashahir/cellograph",
    keywords=["big-data", "manifold-learning", "computational-biology", "graph neural networks", "single-cell", "genomics"],
    classifiers=[
        "Development Status :: 5 - Production/Stable",
        "Intended Audience :: Developers",
        "Intended Audience :: Science/Research",
        "Natural Language :: English",
        "Operating System :: MacOS :: MacOS X",
        "Operating System :: POSIX :: Linux",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
    ],
)

# get location of setup.py
setup_dir = os.path.dirname(os.path.realpath(__file__))




VERSION = '0.0.1' 
DESCRIPTION = 'A graph neural network model to learn latent information in multi-condition single-cell genomics data'
LONG_DESCRIPTION = '''
This is a two-layer graph neural network that takes as input a neighborhood graph of single-cell data spanning two or more biological conditions/treatments, and learns the difference between these conditions in such a way that''s amenable to visualization and clustering
'''
# Setting up
setup(
       # the name must match the folder name 'verysimplemodule'
        name="cellograph", 
        version=VERSION,
        author="Jamshaid Shahir",
        author_email="<jashahir@live.unc.edu>",
        description=DESCRIPTION,
        long_description=LONG_DESCRIPTION,
        packages=find_packages(),
        install_requires=[], # add any additional packages that 
        # needs to be installed along with your package. Eg: 'caer'
        
        keywords=['python', 'gnn', 'single-cell', 'genomics'],
        classifiers= [
            "Development Status :: 3 - Alpha",
            "Intended Audience :: Researchers",
            "Programming Language :: Python :: 3",
            "Operating System :: MacOS :: MacOS X",
            "Operating System :: Microsoft :: Windows",
        ]
)