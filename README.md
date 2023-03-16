### Overview
Code associated with the manuscript "Mutation accumulation in exponentially growing
populations". Preprint version obtainable at https://arxiv.org/abs/2208.02088 . 

Usage: 
```
git clone https://github.com/MichaelDNicholson/accumulate_nmutations.git
cd accumulate_nmutations
```

To reproduce manuscript images, the code below should be run in terminal with working directory as accumulate_nmutations.
Figure 1c (cell numbers over time): 
```
python code/multitype_bp_draw.py 

```
Script multitype_bp_draw.py uses precalculated data in /results/simout. However editing script and setting redoSims='T' will result in simulation data being recreated.
