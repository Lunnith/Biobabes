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

Below, a description of the approach of the different algorithms is given. Examples of how to run the algorithms are given in main.py.

### Approach of different algorithms
## Random algorithm

## Greedy

## Depth First
The Depth First algorithm is a constructive algorithm that searches through all the possible foldings of a protein and saves the protein with the lowest score. It creates child states by adding an extended version of the previous state in all possible directions. This algorithm can be used for folding in 2D and 3D. Furthermore, it has two optional forms of pruning that can be used. The first one is P-pruning, which decreases the amount of states that have to be judged by making a rule that several P's in a row can't all fold in the same direction. The second one is directions pruning, which also reduced the amount of states to be judged by making a rule that if only 2 directions have been used after folding the protein for 2/3, these states won't be expanded.

## Important Parts
The Important Parts algorithm is a faster way to use the Depth First Algorithm. For this algorithm, it is possible to either split a protein on the P count or on size. The Depth First algorithm will then be applied to all different parts and the folded parts are then connected. 

## Hill Climber with Simulated Annealing

### Experiments
**INPUT NEEDED**
## Compare all
This experiment compares the outcome (lowest score achieved) and time for all algorithms for both a short and a long protein and creates a combined barchart to visualize the results. This experiment can be run by the following code:

```
python -m code.experiments.compare_all
```

## Depth First vs Greedy (with depth) vs Hill Climber vs Simulated Annealing
This experiment compares the outcome (lowest score achieved) and number of iterations for the Depth First algorithm, Greedy algorithm with different split sizes, Hill Climber algorithm and the Hill Climber algorithm together with Simulated Annealing. This is done for one of the shorter proteins in 2D to compare the scores of the algorithms with the score that is achieved by the Depth First algorithm because this is probably the optimal score. This experiment can be run by the following code:

```
python -m code.experiments.depth_greedy_hill_climber
```

## Greedy combined with Simulated Annealing
This experiment uses a folded protein made by the Greedy with depth algorithm and inputs this in the Simulated Annealing algorithm to see if this will lead to an even better folded protein than just using Greedy with depth. In a plot, the difference in score after the greedy algorithm and after the Simulated Annealing algorithm is shown. This experiment can be run by the following code:
```
python -m code.experiments.greedy_vs_simulated_annealing
```
## Greedy with depth vs Important Parts
This experiment compares the score distribution of running the Greedy with depth algorithm with split size 5 and running the Important Parts algorithm when using splitting on size with a size of 5. This experiment can be run by the following code:

```
python -m code.experiments.greedy_vs_important_parts
```
## Greedy with depth experiment
There are a 2 important input variables for testing greedy. The splits and the before. The script of the experiement consists of a dictionary with the value for split as the key and the number of iterations as values. Because for the larger splits, the time will increase, using the dictionary you can adjust the amount of iterations so it will not run too long. In the script it also states that for every split, the 'before' will take the value for all numbers lower than the split number. 

To do the experiments we have conducted, give the first split a number of iterations (for instance: 2000). Since the split 2 will do the number of iterations twice (for every 'before' possible), give this a number of iterations of 1000. Do this until split 5. Now all splits have the same number of iterations. 

To get the parameters for the best results in the shortest amount of time, we tested different combinations and looked at the results of the experiments. We saw that the best result we could find was -51 and the shortest and most sure way of reaching this was using ...

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
