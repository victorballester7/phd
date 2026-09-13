# Frequency in Arnoldi method

## Frequency explanation

Once you use Arnoldi method to find the dominant eigenvalue you get an eigenvector and associate eigenvalue. Remember that the Modified arnoldi actually computes the eigenvalues $\mu$, which are linked to the ones ($\lambda$) of the linear operator $L$ in LNS: $u_t = L u$
$$
\mu = \exp(\lambda * T)
$$
where $T$ is the time between iterations. So when we do $log(\mu)$ we need to choose a branch of the complex logarithm. Generally we have
$$
\log(\mu) = \log|\mu| + i (\arg(\mu) + 2k pi) = \lambda * T = (\sigma + i \omega) * T
$$
where $k$ is an integer. So
$$
\sigma = \frac{\log|\mu|}{T}
$$
ALWAYS, but
$$
\omega = \frac{\arg(\mu) + 2k pi}{T}
$$
depends on the choice of $k$. This means that the frequency that we get from Arnoldi is not unique, but depends on the choice of $k$. To get the correct frequency we need to see the physical problem and estimate the correct $k$ (from the computed $\arg(\mu)$ and $T$). If we want to have the right frequency with k = 0, then the time T that we should use in Arnoldi should be such that $\omega * T$ is in the range $[-pi, pi] = [0, pi]$ (because the frequency is positive). 




