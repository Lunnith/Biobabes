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
#### Random algorithm
The random algorithm makes a protein using random chosen directions for the bonds at each step. The algorithm stops when it is stuck and will return the shorter protein.

#### Greedy (with depth)
The greedy algorithm makes the protein step by step and choses the bond direction that leads to the best score. When the protein gets stuck, it goes back one step and tries again without being able to go to the stuck position. To use greedy depth, you can determine at how many future bonds you want to look at each step. So if you want to use split = 3, you go trough all possible directions these 3 bonds can go together and see which leads to the best score. If there are multiple directions that lead to the same best score, the algorithm will chose one at random. The larger the split, the slower the agorithm, because it will have to look at all the states for this selection of bonds (similar to depth first algorithm). When there are bonds at the end that can not be devided by the split anymore, these will be seen as one split (but smaller) and the best directions of the bonds will be determined in the same way. It is further possible to set the amount of aminoacids before the splits start, so the aminoacids within one split are not the same each time. When using the input variable 'before' it will use the number in before as the amount of bonds before starting with the normal splits (best direction determined the same way as with the bonds after the splits that do not fit anymore).

#### Depth First
The Depth First algorithm is a constructive algorithm that searches through all the possible foldings of a protein and saves the protein with the lowest score. It creates child states by adding an extended version of the previous state in all possible directions. This algorithm can be used for folding in 2D and 3D. Furthermore, it has two optional forms of pruning that can be used. The first one is P-pruning, which decreases the amount of states that have to be judged by making a rule that several P's in a row can't all fold in the same direction. The second one is directions pruning, which also reduced the amount of states to be judged by making a rule that if only 2 directions have been used after folding the protein for 2/3, these states won't be expanded.

#### Important Parts
The Important Parts algorithm is a faster way to use the Depth First Algorithm. For this algorithm, it is possible to either split a protein on the P count or on size. The Depth First algorithm will then be applied to all different parts and the folded parts are then connected. 

#### Hill Climber with Simulated Annealing
The hill climber algoritm takes a protein (either unfolded, or pre-folded) and changes the direction of a specified amount of bonds randomly. After this change, it checks wether or not the score has improved. If the change had a negative impact on the score, the change gets rejected. In case of a score improvement, this will become the new protein to work with. 
This can result into the algorithm getting stuck into a local optimum. Which is why there is an optional Simulated Annealing function built.
This function makes the hill climber able to sometimes accept a change with a negative impact on the score, in the hopes of finding a better solution later on. The chances of accepting this negative change are relatively high in the beginning of the run, and decrease with each iteration.

### Experiments

#### Randomise
This experiment runs an amount of iterations and creates a random generated protein every time. It will then give a list with all the scores, the best folded protein and the corresponding best score. It will then make a histogram of all the scores and the best score and the used time. This experiment can be run by the following code:

```
python -m code.experiments.randomise_base
```
#### Compare all
This experiment compares the outcome (lowest score achieved) and time for all algorithms for both a short and a long protein and creates a combined barchart to visualize the results. This experiment can be run by the following code:

```
python -m code.experiments.compare_all
```

#### Depth First vs Greedy (with depth) vs Hill Climber vs Simulated Annealing
This experiment compares the outcome (lowest score achieved) and number of iterations for the Depth First algorithm, Greedy algorithm with different split sizes, Hill Climber algorithm and the Hill Climber algorithm together with Simulated Annealing. This is done for one of the shorter proteins in 2D to compare the scores of the algorithms with the score that is achieved by the Depth First algorithm because this is probably the optimal score. This experiment can be run by the following code:

```
python -m code.experiments.depth_greedy_hill_climber
```

#### Greedy combined with Simulated Annealing or Hill Climber
This experiment uses a folded protein made by the Greedy with depth algorithm and inputs this in the Hill Climber or Simulated Annealing algorithm to see if this will lead to an even better folded protein than just using Greedy with depth. In a plot, the difference in score after the greedy algorithm and after the HillClimber/Simulated Annealing algorithm is shown. If the experiment is now run, Greedy combined with Simulated Annealing will be performed. If you want to combine Greedy with HillClimber, you should use the commented code. This experiment can be run by the following code:

```
python -m code.experiments.greedy_vs_simulated_annealing_or_hillclimber
```

#### Greedy with depth vs Important Parts
This experiment compares the score distribution of running the Greedy with depth algorithm with split size 5 and running the Important Parts algorithm when using splitting on size with a size of 5. This experiment can be run by the following code:

```
python -m code.experiments.greedy_vs_important_parts
```

#### Greedy with depth experiment
There are a 2 important input variables for testing greedy. The splits and the before. The script of the experiment consists of a dictionary with the value for split as the key and the number of iterations as values. Because for the larger splits, the time will increase, using the dictionary you can adjust the amount of iterations so it will not run too long. In the script it also states that for every split, the 'before' will take the value for all numbers lower than the split number. 

To do the experiment we have conducted, the first split has 2000 iterations. Since the split 2 will do the number of iterations twice (for every 'before' possible), give this a number of iterations of 1000. This was done until split 5. Now all splits have the same number of iterations. The experiment plots a boxplot were you can see what the score distribution is of chosing a split and also gives the time it takes in a dataframe. 
This experiment can be run by the following code:

```
python -m code.experiments.greedy_with_depth
```

#### Simulated Annealing parameters
To find the best parameters (temperature scheme and starting temperature) for the Simulated Annealing addition to the Hill-climber, an experiment was created to compare the results for these parameters.  
This script executes a number of simulated annealing run for all different parameter options, then takes the average scores of these runs per parameter settings and plots them. The graph then shows that a linear temperature scheme is far from ideal, which is why we use the exponential temperature scheme. The difference in starting temperature is less clear, but we have chosen to use 10 for this.

This experiment can be run by the following code:

```
python -m code.experiments.simulated_annealing_params
```

#### Hill-Climber vs Simulated Annealing
This experiment shows a distribution of the scores achieved by both the Hill-Climber with, and the Hill-Climber without the Simulated Annealing option. However, it also shows the difference in the amount of bonds changed per iteration. The Hill-Climber without sim. annealing constantly changes a minimum of n specified bonds. The Hill-Climber with sim. annealing reduces the amount of n during the run.
In the graph created by this experiment, the distribution of the scores per (starting) amount of n is shown.

It turns out that the Hill-climber without sim. annealing profits most from n=1, which seems reasonable as a higher n makes it more likely to get stuck while improving the protein. The Hill-Climber with sim. annealing seems to profit most from starting with n=8. That also seems reasonable, since a bigger n at the start could lead to bigger improvements and since the sim. annealing reduces the n throughout the run, it is less likely to get stuck than the Hill-Climber without the sim. annealing. This experiment can be run by the following code:

```
python -m code.experiments.hill_climber_vs_sim_annealing
```

### Structure
This list contains the structure of this repository and what different folders or files contain:
- **/code**: contains all code of the project
    - **code/algorithms**: contains the code for the algorithms
    - **code/classes**: contains the two classes that are necessary for this case
    - **code/visualization**: contains the code to visualize the final fold
    - **code/experiments**: contains the code to run different experiments
- **/results**: contains all visualizations of the results of the experiments
## Authors
- Alaya Storm
- Luna van der Velden
- Sofie Perizonius
