
StellarFB-IC
============
This is a python-based initial condition genereting tool for developing 
the stellar feedback implications in the moving-mesh cosmological simulation 
code AREPO. 

Author: Yunwei Deng

Requirement
===========
python3
numpy
scipy
astropy
h5py
matplotlib

Quick start
===========
0. $vim IC_Generator.py
1. choose a type of IC by uncommenting its name
2. Set the stellar parameters
3. Set the gas parameters
4. $:wq
5. $python IC_Generator.py
6. Check the filename.pdf file
