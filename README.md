### Overview
Code associated with the manuscript "Mutation accumulation in exponentially growing
populations". Preprint version obtainable at https://arxiv.org/abs/2208.02088 . 

The aim of this repository is to provide code to be able to numerically evaluate the main results and to reproduce figures in manuscript.

### Numerically evaluate main results
The main results from the paper are as follows (please refer to paper for explanation of notation). The number of cells of type $n$ at time $t$ is 

```math
Z_n(t)\approx V_n t^{r_n-1}e^{\delta_n t},
```
where $V_n$ is Mittag-Leffler distributed with tail parameter  $\lambda_1/\delta_n$, and scale parameter $\omega_n$. 

In R, to numerically evaluating the scale parameter $\omega_n$ from a set of birth/death/mutation parameters, the function get_scale_tail_vec in /code/functionsCoreAccumulateNmuts.R can be used. This numerical value of $\omega_n$ is required, e.g. for evaluating $\mathbb{P}(Z_n(t)>k) \approx \mathbb{P}(V_n > k t^{1-r_n}e^{-\delta_n t})$. 

With $\omega_n$ in hand, the second main result can be numerically evaluated, which is that the arrival time for type $n+1$ cells is approximately 
```math 
\mathbb{P}(\tau_{n+1} >t) \approx \left[1+ \exp\left(\lambda_1 (t-t_{1/2}^{(n+1)})\right)\right]^{-1},
```
with $t_{1/2}^{(n+1)}=\frac{1}{\delta_n}\log\frac{\delta_n}{\omega_n \nu_n [\delta_n^{-1}\log(\nu_n^{-1})]^{r_n-1}}$.


### Reproduce figures
To reproduce manuscript images, the code below should be run in terminal with working directory as accumulate_nmutations.
Figure 1c (cell numbers over time): 
```
python code/multitype_bp_draw.py 
```
Script multitype_bp_draw.py uses precalculated data in /results/simout. However editing script and setting redoSims='T' will result in simulation data being recreated.
