### Overview
Code associated with the manuscript "Mutation accumulation in exponentially growing
populations". Preprint version obtainable at https://arxiv.org/abs/2208.02088 . 

The aim of this repository is to provide code to be able to numerically evaluate the main results and to reproduce all figures in manuscript.

## Numerical evaluate main results

To reproduce manuscript images, the code below should be run in terminal with working directory as accumulate_nmutations.
Figure 1c (cell numbers over time): 
```
python code/multitype_bp_draw.py 
```
Script multitype_bp_draw.py uses precalculated data in /results/simout. However editing script and setting redoSims='T' will result in simulation data being recreated.
