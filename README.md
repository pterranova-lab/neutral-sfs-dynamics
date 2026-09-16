# neutral-sfs-dynamics
A computational model of time-dependent site-frequency spectrum dynamics under neutral Wright–Fisher drift and recurrent mutation.

# neutral-sfs-dynamics

A computational model of time-dependent site-frequency spectrum dynamics under neutral genetic drift and recurrent mutation.

## Overview

The Site Frequency Spectrum (SFS) is a summary statistic describing the frequencies of mutations observed in a population.

For a sample of size $n$, define

$$
\xi_j = \text{number of mutations observed at allele frequency } j,
$$

where $j = 1,\dots,n-1$.

The SFS can therefore be represented as the vector

$$
\boldsymbol{\xi}
=
(\xi_1,\xi_2,\dots,\xi_{n-1}).
$$

Under a constant population size and neutral evolution, mutation and genetic drift produce a characteristic equilibrium SFS. In the standard neutral model, the expected spectrum has the form

$$
\xi_j = \frac{2n\mu}{j},
$$

where $\mu$ is the per-site mutation rate. Thus,

$$
\boldsymbol{\xi}
=
\left(
2n\mu,\,
n\mu,\,
\frac{2}{3}n\mu,\,
\dots,\,
\frac{2n\mu}{n-1}
\right).
$$

This project asks a forward-modeling question:

> Given a specified evolutionary process, how does the SFS change through time before reaching its equilibrium distribution?

Rather than beginning with an observed SFS and attempting to infer the evolutionary history that produced it, the model begins with an initial SFS and explicitly evolves it under a specified mutation-drift process.

## Mathematical Model

The SFS is modeled as a continuous-time dynamical system,

$$
\frac{d\boldsymbol{\xi}}{dt}
=
G\boldsymbol{\xi}+\boldsymbol{m},
$$

where:

- $\boldsymbol{\xi}$ is the vector of mutation counts across frequency classes
- $G$ is the drift matrix
- $\boldsymbol{m}$ represents the introduction of new mutations

For a discrete time step $\Delta t$, the corresponding Euler update is

$$
\boldsymbol{\xi}_{t+\Delta t}
=
\boldsymbol{\xi}_t
+
\left(
G\boldsymbol{\xi}_t+\boldsymbol{m}
\right)\Delta t.
$$

At equilibrium,

$$
\frac{d\boldsymbol{\xi}}{dt}=0,
$$

so

$$
G\boldsymbol{\xi}_{\mathrm{eq}}
=
-\boldsymbol{m}.
$$

New mutations are introduced into the lowest-frequency class. In this model,

$$
\boldsymbol{m}
=
(2n\mu,0,\dots,0).
$$

## Constructing the Drift Matrix

The equilibrium condition alone does not uniquely determine $G$. There are $n-1$ equilibrium equations but $(n-1)^2$ entries in an unrestricted $(n-1)\times(n-1)$ matrix.

A specific drift matrix can therefore be constructed by imposing additional structure on the evolutionary process.

This implementation uses a nearest-neighbor transition structure motivated by the Moran model. For frequency class $j$,

$$
G_{j,j}
=
-\frac{2j(n-j)}{n},
$$

with neighboring transition terms

$$
G_{j,j+1}
=
\frac{(j+1)(n-j-1)}{n}
$$

and

$$
G_{j,j-1}
=
\frac{(j-1)(n-j+1)}{n}.
$$

All other entries are zero.

This produces a tridiagonal drift matrix in which mutations move between neighboring allele-frequency classes through genetic drift.

## Analytic Solution

Because $G$ and $\boldsymbol{m}$ are constant, the differential equation has the closed-form solution

$$
\boldsymbol{\xi}(t)
=
e^{Gt}
\left(
\boldsymbol{\xi}_0-\boldsymbol{\xi}_{\mathrm{eq}}
\right)
+
\boldsymbol{\xi}_{\mathrm{eq}},
$$

where

$$
\boldsymbol{\xi}_{\mathrm{eq}}
=
-G^{-1}\boldsymbol{m}.
$$

The analytic implementation evaluates this solution using the matrix exponential.

This allows the SFS to be evaluated directly at selected time points rather than integrating the system one time step at a time.

## Numerical Implementation

Two implementations are provided.

### Euler Integration

`simulate_sfs()` evolves the SFS numerically using the discrete Euler update:

```r
SFS <- SFS + (G %*% SFS + m) * dt
