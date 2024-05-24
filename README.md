# PGdet: Determination of symmetry pointgroup from a SMILES organic molecule

[![MIT License](https://img.shields.io/badge/License-MIT-green.svg)](https://choosealicense.com/licenses/mit/)

PGdet is a university project created by students of the Swiss Federal Institute of Technology Lausanne: Aldric Mercier, Alexandre Keller and Miguel Correia. The aim of this code is to take organic molecules written in SMILES and extract their point group.


## 🛠️ Installation

To install our project, you can create a new conda environment with Anaconda with this code:

```bash
  conda create -n pgdet
  conda activate pgdet
```
Then, clone it into your new environment with:
```bash
  git clone https://github.com/alexandreklr/PGdet.git
```
Finally navigate to PGdet and install the dependencies with:
```bash
  cd PGdet
  pip install -e .
```


## 💡 Basic usage

The functions of this project may be imported as follows:

```bash
from PGdet import atom_mapping
```
The main function is atom_mapping(smiles, desc, plot). It takes three arguments. The first one is a smiles as a string or a list of smiles. The second one is called 'desc' and the third one is called 'plot'. The second and third arguments are booleans and are considered False if not precised.
The function return the value of the smiles input. For example, it should return 'D6h' when given the smiles 'c1ccccc1'.
```bash
pg = atom_mapping('c1ccccc1')
print(pg)
```
The function also accepts a list of smiles. In this case it returns a list of pointgroups in the same order. For example, when given ['CC','CCO','C'], it should return ['D3d','Cs','Td'].
```bash
pg = atom_mapping(['CC','CCO','C'])
print(pg)
```
If 'desc' is set to True, the function will also print a description of the molecule. The description contains the molecule's point group, a list of atoms in the molecule, a list of their coordinates, a list of sets indicating which atoms are bonded, a list of their lenghts normalized to approximatively 1 and a list of sets containing atom indexes and the angle formed by them.
```bash
atom_mapping('CC', True)
```
However, the bonds in a ring are not normalised to 1:
```bash
atom_mapping('c1ccccc1', True)
```
If 'plot' is set to True, the function will also plot the molecule in colors. Some atoms, of which the color is not precised in the program, will appear grey.
```bash
atom_mapping('CCO', False, True)
```
The function does not accept inorganic molecules unless they are part of a list of exception. The code below should return an error precising which inorganic molecules are accepted.
```bash
atom_mapping('FS(F)(F)(F)(F)F')
```
The function also does not accept molecules which contain more than 25 non hydrogen atoms.
```bash
atom_mapping('CCCCCCCCCCCCCCCCCCCCCCCCCC')
```


## ⚛️ Authors

- [@alexandreklr](https://github.com/alexandreklr)
- [@aldricmrcr](https://github.com/aldricmrcr)
- [@Migloops](https://github.com/Migloops)

## 📞Feedback

If you have any feedback or suggestions, please reach out to us at alexandre.keller@epfl.ch
