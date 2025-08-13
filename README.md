# Exact Analytical Option Pricing Formula Under a Jump-Diffusion Model with Hawkes Intensity

This repository contains the implementation and documentation for the analytical option pricing framework under a jump-diffusion model with Hawkes intensity. For inquiries regarding this research, please contact the corresponding author.

## Hawkes Process

The intensity function of a Hawkes process is defined as:

$$
\lambda_t = \gamma + (\lambda_0-\gamma)e^{-\beta t} + \sum Y_k e^{-\beta(t-t_k)}
$$

where $\gamma$ represents the baseline intensity and $\lambda_0$ denotes the initial intensity. At each jump occurrence time $t_k$ (for $k = 1,2,\ldots$), the intensity process $\lambda_t$ experiences an upward jump of magnitude $Y_k$. These jump magnitudes follow an exponential distribution with cumulative distribution function $G(y) = 1-e^{-\alpha y}$ and expected value $1/\alpha$. Each historical jump contributes to the current intensity through an exponential decay mechanism with decay rate $\beta$.

### Functions: `hawkesSim.m` and `hawkesSimC.m`

These functions perform Monte Carlo simulation of Hawkes processes for both unconditional and conditional scenarios.

**Usage:**

```matlab
Pi = hawkesSim(r, a, b, t, paths)
```

```matlab
Pi = hawkesSimC(r, lambda0, a, b, t, paths)
```

**Parameters:**
- `r`: baseline intensity parameter
- `a`: jump intensity rate
- `b`: decay rate parameter
- `t`: time vector (row vector)
- `lambda0`: initial intensity value (for conditional simulation only)
- `paths`: number of simulation paths

**Returns:**
- `Pi`: simulated process values with dimensions (N, t)

### Function: `hawkesPMF.m`

This function computes the probability mass function (PMF) of the Hawkes process, specifically $\mathrm{P}(N_t = n)$, using a closed-form analytical expression.

**Usage:**

```matlab
[probUC, probC] = hawkesPMF(ga, alpha, beta, t, lambda0, N)
```

**Parameters:**
- `ga`: $\gamma$ (baseline intensity parameter)
- `alpha`: $\alpha$ (exponential distribution parameter)
- `beta`: $\beta$ (intensity decay rate)
- `t`: time horizon
- `lambda0`: $\lambda_0$ (initial intensity value)
- `N`: maximum value of $n$ for computation

**Returns:**
- Column vector of length $N+1$ (corresponding to $n = 0, 1, \ldots, N$)
  - `probUC`: unconditional probabilities
  - `probC`: conditional probabilities
 
**Computational Constraints:**
- For $N = 70$: requires $\beta \leq 5$
- For $N = 50$: requires $\beta \leq 9$
- No constraint on $\alpha \cdot \beta$, however large $\beta$ values may cause numerical instability

## Return Distribution under HJD

The logarithmic return under the jump-diffusion model is expressed as:

$$ R_T = \ln\left(\frac{S_T}{S_0}\right) = \left(r - \frac{\sigma^2}{2} - \nu_T e^{\theta + \frac{\delta^2}{2}}\right)T + \sigma W_T + \sum_{j=1}^{N_T} \ln Y_j $$

Let $R_n$ denote the conditional return $R_T$ given that $N_T = n$. The distribution of $R_n$ follows a normal distribution:

$$ R_n \sim \mathcal{N}(a_n, b_n^2) $$

with moment parameters:

$$ a_n = \left(r - \frac{\sigma^2}{2} - \nu_T e^{\theta + \frac{\delta^2}{2}}\right)T + n \theta, \quad b_n^2 = \sigma^2 T + n \delta^2 $$

The probability density function (PDF) of $R_T$ is characterized as a mixture of normal distributions:

$$ f_{R_T}(y) = \sum_{n=0}^{\infty} P(N_T = n) \frac{1}{\sqrt{2\pi} b_n} \exp\left(-\frac{(y - a_n)^2}{2 b_n^2}\right) $$

### Function: `returnMJD.m`

This function calculates both the probability density function (PDF) and cumulative distribution function (CDF) of logarithmic returns under the Merton jump-diffusion model with Hawkes intensity.

**Usage:**
```matlab
[f, F] = returnMJD(y, r, q, sigma, t, gamma, delta, jumpProb)
```

**Parameters:**
- `y`: logarithmic return values for evaluation
- `r`: risk-free interest rate
- `q`: dividend yield rate
- `sigma`: volatility of the underlying asset
- `t`: time to maturity
- `gamma`: mean jump size parameter
- `delta`: standard deviation of jump sizes
- `jumpProb`: probability distribution of $N_T = n$

**Returns:**
- `f`: PDF values $f(y)$ evaluated at each point in `y`
- `F`: CDF values $F(y)$ evaluated at each point in `y`

## Option Pricing under HJD

The European call option price under the Hawkes Jump-Diffusion (HJD) model is given by:

$$
C_{\mathrm{HJD}}(S_0, K, T, r, q, \sigma;\alpha,\beta,\gamma;\theta,\delta) = \sum_{n=0}^\infty \mathrm{P} (N_T = n) e^{-(r - \tilde{r_n}) T} {C_{\mathrm{BS}}(S_0, K, T, \tilde{r}_n, q, \tilde{\sigma}_n)}
$$

where the adjusted parameters are defined as:

$$
\tilde{r}_n = \left(r-\frac{1}{2}\sigma^2 + \frac{1}{2}\tilde{\sigma}_n^2\right) + \frac{n\theta-\eta}{T}, \qquad \tilde{\sigma}_n^2 = \sigma^2 + \frac{n\delta^2}{T}
$$

$$
\eta = \ln \left(\sum_{n=0}^\infty \left[\mathrm{P}(N_T = n) \exp\left(n \cdot \theta +\frac{n \cdot \delta^2}{2}\right)\right]\right)
$$

### Function: `callJD.m`

This function computes the European call option price under the jump-diffusion framework.

**Usage:**

```matlab
price = callJD(S0, K, r, q, sigma, t, theta, delta, jumpProb)
```

**Parameters:**
- `S0`: current price of the underlying asset
- `K`: strike price of the option
- `r`: risk-free interest rate
- `q`: dividend yield rate
- `sigma`: volatility of the underlying asset
- `t`: time to maturity
- `theta`: mean of logarithmic jump sizes
- `delta`: standard deviation of logarithmic jump sizes
- `jumpProb`: probability distribution of $N_T = n$

**Returns:**
- `price`: European call option price
