# Protein Powder
Proteins are long strings of aminoacids that are important for several processes in the human body. The function of proteins is determined by the combination of aminoacids as well as its folding. In this case the number of aminoacids is constrained to three: hydrophobic (H), polar (P) and cysteine (C) aminoacids. Hydrophobic aminoacids attract each other, as well as cysteine aminoacids. Cysteine aminoacids can also interact with each other which leads to the biggest increase in stability of the folding. To gain more insight into the folding of proteins, this project contains several algorithms to fold certain proteins in the most optimal way.

## Get everything working
### Requirements
This codebase is written in Python 3.10. Requirements.txt contains all necessary package to run the code succesfully. Requirements.txt can be used true the following instructions:

```
pip install -r requirements.txt
```
or via conda:

```
conda install --file requirements.txt
```

### Use
The main.py document contains an example for all different algorithms on how to run these algorithms. If you run the following command without changing anything in the document, only the random algorithm will be executed for the first protein. You can change which protein you want to explore and which algorithm you want to be executed by commenting or decommenting certain parts.

```
python main.py
```

#### Experiments
**INPUT NEEDED**

### Structure
This list contains the structure of this repository and what different folders or files contain:
- **/code**: contains all code of the project
    - **code/algorithms**: contains the code for the algorithms
    - **code/classes**: contains the two classes that are necessary for this case
    - **code/visualization**: contains the code to visualize the final fold
    - **code/experiments**: contains the code to run different experiments
- **/data**: contains results of experiments

## Authors
- Alaya Storm
- Luna van der Velden
- Sofie Perizonius
