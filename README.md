# Exact Analytical Option Pricing Formula Under a Jump-Diffusion Model with Hawkes Intensity
This repository includes code and figures for the study mentioned above. For additional inquiries, please contact the corresponding author.
## Hawkes Process

The intensity of a Hawkes process is specified as

$$
\lambda_t = \gamma + (\lambda_0-\gamma)e^{-\beta t} + \sum Y_k e^{-\beta(t-t_k)}
$$


where $\gamma$ is the baseline intensity and $\lambda_0$ is the initial intensity. When jumps happen at time $t_k, k = 1,2,\cdots$, the intensity process $\lambda_t$ jumps up by $Y_k$, which follows an exponential distribution with distribution function $G(y) = 1-e^{-\alpha y}$ and mean $1/\alpha$. Each up-jump happening in the recent past decays exponentially at a rate $\beta$.

### Function: `hawkesPMF.m`

This function calculates the probability mass function (PMF) of the Hawkes process, $\mathrm{P}(N_t = n)$.

**Usage:**

```matlab
[probUC, probC] = hawkesPMF(ga, alpha, beta, t, lambda0, N)
```

**Input:**
- `ga`: $\gamma$ (baseline intensity)
- `alpha`: $\alpha$ (parameter of the exponential distribution)
- `beta`: $\beta$ (decay rate of intensity)
- `t`: time
- `lambda0`: $\lambda_0$ (initial intensity)
- `N`: maximum value of $n$ to calculate

**Output:**
- A column vector of size $n+1$ (for $n = 0$ to $N$)
  - `probUC`: unconditional probability
  - `probC`: conditional probability
 
**Parameter Constraints**
- For $N = 70$, $\beta \leq 5$
- For $N = 50$, $\beta \leq 9$
- No restriction on $\alpha \cdot \beta$, but large $\beta$ may cause computational issues


## Retrun Distribution of HJD
The logarithmic return of the jump-diffusion model can be expressed as
\begin{equation}
    R_T = \ln\left(\frac{S_T}{S_0}\right) = \left(r-\frac{\sigma^2}{2}-\nu_T(e^{\theta+\frac{\delta^2}{2}})\right)T + \sigma W_T + \sum_{j=1}^{N_T}\ln Y_j.
\end{equation}
Let \( R_n \) denote the return \( R_T \) conditional on \( N_T = n \), then \( R_n \sim \mathrm{N}(a_n, b_n^2) \) with
\begin{eqnarray}
    a_n &=& \left(r-\frac{\sigma^2}{2}-  \nu_T(e^{\theta+\frac{\delta^2}{2}})   \right)T + n\theta, \nonumber\\[5pt]
    b_n^2 &=& \sigma^2 T + n\delta^2. \nonumber
\end{eqnarray}
The PDF of \( R_T \) takes a mixture form of normal distributions as shown below:
\begin{eqnarray}
    f_{R_T}(y) &=& \sum_{n=0}^\infty \mathrm{P}(N_T=n)\frac{1}{\sqrt{2\pi}b_n}e^{\frac{-(y-a_n)^2}{2b_n^2}}.\nonumber
\end{eqnarray}


## Option Price of HJD

The call option price under the Hawkes Jump-Diffusion (HJD) model is given by

$$
C_{\mathrm{HJD}}(S_0, K, T, r, q, \sigma;\alpha,\beta,\gamma;\theta,\delta) = \sum_{n=0}^\infty \mathrm{P} (N_T = n) e^{-(r - \tilde{r_n}) T} {C_{\mathrm{BS}}(S_0, K, T, \tilde{r}_n, q, \tilde{\sigma}_n)}
$$

where

$$
\tilde{r}_n = \left(r-\frac{1}{2}\sigma^2 + \frac{1}{2}\tilde{\sigma}_n^2\right) + \frac{n\theta-\eta}{T}, \qquad \tilde{\sigma}_n^2 = \sigma^2 + \frac{n\delta^2}{T}
$$

$$
\eta = \ln \left(\sum_{n=0}^\infty \left[\mathrm{P}(N_T = n) \exp\left(n \cdot \theta +\frac{n \cdot \delta^2}{2}\right)\right]\right)
$$

### Function: `callJD.m`

This function calculates the call option price under the Jump-Diffusion model.

**Usage:**

```matlab
price = callJD(S0, K, r, q, sigma, t, theta, delta, jumpProb)
```

**Input:**
- `S0`: current price of the underlying asset
- `K`: strike price
- `r`: risk-free interest rate
- `q`: dividend yield
- `sigma`: volatility of the underlying asset
- `t`: time to maturity
- `theta`: mean of jump sizes
- `delta`: standard deviation of jump sizes
- `jumpProb`: probability of $N_T = n$

**Output:**
- `price`: call option price
