# neutral-sfs-dynamics
A computational model of time-dependent site-frequency spectrum dynamics under neutral Wright–Fisher drift and recurrent mutation.

# Time-Dependent Site-Frequency Spectrum Dynamics

A computational model for exploring how a site-frequency spectrum (SFS) evolves through time under neutral genetic drift and recurrent mutation using a Wright–Fisher framework.

## Overview

The site-frequency spectrum describes the number of mutations observed at different allele frequencies within a sampled population.

This project models the SFS as a dynamical system derived from the neutral Wright–Fisher model. Rather than starting with an observed SFS and attempting to infer the evolutionary process that produced it, the model starts with a specified evolutionary process and examines how the frequency distribution changes over time.

The central dynamical system is

$$
\frac{d\xi}{dt} = G\xi + m
$$

where:

- $\xi$ is the vector of mutation counts across allele-frequency classes
- $G$ is the transition matrix describing neutral Wright–Fisher genetic drift
- $m$ is the mutation input vector

For constant $G$ and $m$, the system has the analytic solution

$$
\xi(t) =
e^{Gt}(\xi_0-\xi_{\mathrm{eq}})
+\xi_{\mathrm{eq}}
$$

where the equilibrium distribution satisfies

$$
G\xi_{\mathrm{eq}} + m = 0
$$

and therefore

$$
\xi_{\mathrm{eq}} = -G^{-1}m.
$$

The model can therefore be used to examine how an arbitrary initial SFS approaches the neutral equilibrium expected under the Wright–Fisher process.

## Wright–Fisher Model

The underlying population-genetic framework is the neutral Wright–Fisher model, with genetic drift represented as transitions between neighboring allele-frequency classes.

For a population represented by $n$ sampled chromosomes, the model contains $n-1$ segregating-frequency classes. The drift matrix $G$ is constructed from the Wright–Fisher transition structure.

Mutation is represented by an input vector in which new mutations enter the lowest-frequency class.

Under the assumptions of the neutral model, the equilibrium has the classical inverse-frequency form

$$
\xi_j \propto \frac{1}{j},
$$

corresponding to the familiar power-law structure of the neutral SFS.

## Implementations

Two implementations are provided.

### Numerical Euler Integration

`simulate_sfs()` evolves the SFS using a discrete Euler update:

```r
SFS <- SFS + (G %*% SFS + m) * dt
