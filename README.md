# Stochastic Asset Pricing Models in C++

This repository contains C++ implementations of advanced stochastic processes used in quantitative finance. Compared with a normal distribution, the
logarithm of a price process with jumps is often leptokurtotic, meaning that it
has a high peak and heavy tails, resembling the price movements in actual markets. Models implemented in this repository:
1. Merton Jump-Diffusion Model
2. Variance Gamma Process
3. Normal Inverse Gamma Process

---
## Mathematical Theory

### 1. Merton Jump-Diffusion (MJD)

The MJD model addresses the continuity limitation of Geometric Brownian Motion by adding a Poisson-driven jump component.

* **Stochastic Differential Equation (SDE):**

    $$\frac{dS_t}{S_{t-}} = (r - \lambda \kappa) dt + \sigma dW_t + d\left( \sum_{i=1}^{N_t} (Y_i - 1)   \right)$$


  - $r$: The risk-free rate.
  - $\lambda$: The jump intensity (average number of jumps per year).
  - $\kappa$: The expected relative jump size, $E[Y-1]$.
  - $Y_i$: The $i$-th jump size.
  - $N_t$: A poisson process.

* **The Solution:**
To simulate the price at time $t+\Delta t$, we use the solution to the SDE:


$$S_{t+\Delta t} = S_t \cdot \exp\left( (r - \lambda \kappa - \frac{1}{2}\sigma^2)\Delta t + \sigma\sqrt{\Delta t}Z \right) \cdot \prod_{i=N_t+1}^{N_{t+\Delta t}} Y_i$$

* **Logic:** The price is multiplied by each jump magnitude $Y_i$ that occurs during the interval $\Delta t$.

---

### 2. Variance Gamma (VG) Process

The core idea of the VG process is that the "economic clock" does not run at a constant speed. Market activity accelerates during periods of high information flow and slows down during lulls.

Mathematically, we model this by taking a Brownian Motion with drift and evaluating it at a random time-change

* **The SDE:**


$$\frac{dS_t}{S_{t-}} = (r + \omega) dt + \theta dG_t + \sigma dW_{G_t}$$

* **The Solution (Multiplicative Form):**


    $$S_{t+\Delta t} = S_t \cdot \exp\left( (r + \omega) \Delta t + \theta \Delta G + \sigma \sqrt{\Delta G} Z \right)$$

  - $\Delta G$: An increment of the Gamma process $G \sim \text{Gamma}(\frac{\Delta t}{\beta}, \beta)$.
  - $\omega$: The drift correction $\frac{1}{\beta} \ln(1 - \theta \beta - \sigma^2 \beta / 2)$.

* **Logic:** Here, the "jump" behavior is continuous. Every time step is a "jump" of random size governed by the Gamma subordinator.

### 3. Normal Inverse Gaussian (NIG) Process
The NIG process is a pure-jump Lévy process with infinite activity and infinite variation. It is constructed by subordinating a Brownian motion with drift to an Inverse Gaussian process, allowing it to capture pronounced skewness and kurtosis in asset returns.

* **The representation (subordination form):**
  $$
  X_t = \beta I_t + \sqrt{I_t} \, W_{I_t}
  $$
  where $W_t$ is standard Brownian motion and $I_t$ is an Inverse Gaussian Lévy process with mean $t$ and shape parameter related to $\gamma$.

* **Increment per time step (simulation form):**
  $$
  \Delta X = \beta Y + \sqrt{Y} \, Z
  $$
  - $Y \sim \text{IG}(\delta \Delta t, \gamma)$  (inverse Gaussian with mean $\Delta t$ when parametrized appropriately)
  - $Z \sim \mathcal{N}(0,1)$

* **The exponential Lévy model (multiplicative form for asset price):**
  $$
  S_{t+\Delta t} = S_t \cdot \exp\left( (r + \omega) \Delta t + \Delta X \right)
  $$
  where
  - $r$ is the short rate (risk-free rate),
  - $\omega$ is the **drift correction** (martingale adjustment / compensator) ensuring the discounted asset price is a martingale under the risk-neutral measure.

* **Drift correction $\omega$ (common parametrization):**
  $$
  \omega = -\delta \left( \gamma - \sqrt{\gamma^2 - 2\beta - 1} \right)
  $$

* **Logic:** The Inverse Gaussian subordinator $I_t$ (or its increments $Y$) plays the role of a random "business time" or activity clock. Because the IG process has infinite activity near zero, the NIG process exhibits infinitely many small jumps in any interval, resulting in paths of **infinite variation** (more irregular / wiggly than VG). This makes NIG particularly suitable for modeling assets with heavy tails, skewness, and high-frequency jump behavior.


## Why These Models Matter

Standard Black-Scholes Geometric Brownian Motion (GBM) assumes returns are normally distributed, which significantly underestimates the probability of market crashes.

* **Tail Risk:** Jump models explicitly price the risk of Black Swan events.
* **Skewness/Kurtosis:** Lévy processes like VG and NIG capture the leptokurtotic nature of actual price movements.

---