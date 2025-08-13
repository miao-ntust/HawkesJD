# Exact Analytical Option Pricing Formula Under a Jump-Diffusion Model with Hawkes Intensity
This repository contains MATLAB code and figures accompanying our study on option pricing under a jump-diffusion model with Hawkes jump intensity. For further inquiries, please contact the corresponding author.

## Hawkes Process

The Hawkes process intensity is given by

$$
\lambda_t = \gamma + (\lambda_0 - \gamma) e^{-\beta t} + \sum_{k} Y_k e^{-\beta(t - t_k)},
$$

where:
- $\gamma$ — baseline intensity  
- $\lambda_0$ — initial intensity  
- $Y_k \sim \mathrm{Exp}(\alpha)$ — intensity jumps with mean $1/\alpha$  
- $\beta$ — exponential decay rate

At each jump time $t_k$, the intensity $\lambda_t$ increases by $Y_k$, and each increment decays exponentially at rate $\beta$.

---

### Function: `hawkesSim.m` and `hawkesSimC.m`

#### `hawkesSim.m` — Unconditional Simulation
```matlab
Pi = hawkesSim(r, a, b, t, paths)

