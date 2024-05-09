# Ising Model

## Theory

### Diagonalization method

The code is available in [diagonalization.py](./diagonalization.py).

We assume a lattice of $N_i \times N_j$ (with $N_i$ rows and $N_j$ columns). Each site has a spin, $s_{i j}$. Thus, we define an operator $\mathbf{A}_{i j}$ that when be applied on the state of this lattice, $\left| s_{00} s_{01} \ldots s_{N_i-1, N_j-1}  \right>$, gives us the spin of that specific site, $s_{i j}$. In other words, we want to create an eigenvalue problem.

$$
    \mathbf{A}_{i j} \left| s_{00} s_{01} \ldots s_{N_i-1, N_j-1}  \right> = s_{i j} \left| s_{00} s_{01} \ldots s_{N_i-1, N_j-1}  \right>
$$

In order to create these operators, we used the following formula:

$$
    \mathbf{A}_{i j} = \underbrace{\overbrace{\mathbf{I} \otimes \ldots \otimes\mathbf{I}}^{N_j i + j\ \text{operators}} \otimes \sigma^z \otimes \overbrace{\mathbf{I} \otimes \ldots \otimes\mathbf{I}}^{N_j N_i - N_j i - j - 1\ \text{operators}}}_{N_i N_j\ \text{operators}}, \quad i,j = 0, 1, 2, \ldots
$$

where $\otimes$ means tensor product and $\sigma^z$ is the third Pauli matrix, however, because we want to manipulate the result as a matrix, the resultant tensor should be flattened.

The energy corresponding to $s_{i j}$ is

$$
    H_{i j} = \epsilon s_{i j} + J \sum_{n} s_{i j} s_n,
$$

where $\epsilon$ is the energy resulted from an external magnetic field, $J$ is called the exchange energy, and $n$ refers to the neighboring sites. Similarly, the total energy of the two-dimensional lattice can be found by summing over all of the sites:

$$
    \mathbf{H} = \epsilon \sum_{i, j = 0}^{N_i - 1, N_j - 1} \mathbf{A}_{i j} + \frac{1}{2} J \sum_{i, j = 0}^{N_i - 1, N_j - 1} \mathbf{A}_{i j} (\mathbf{A}_{i + 1, j} + \mathbf{A}_{i - 1, j} + \mathbf{A}_{i, j + 1} + \mathbf{A}_{i, j - 1})
$$

The fraction $1/2$ has emerged because in its absence, the energy between two neighboring sites will be summed twice. It is worth noting that we only have created the hamiltonian, so when $\mathbf{H}$ be applied on the state of the system, it will give the total energy.

$$
    \mathbf{H} \left| s_{00} s_{01} \ldots s_{N_i-1, N_j-1}  \right> = E \left| s_{00} s_{01} \ldots s_{N_i-1, N_j-1}  \right>
$$

The least amount of energy can be found once the matrix $\mathbf{H}$ is diagonalized. The elements on its main diagonal are the eigenvalues or the energies of the system, and their corresponding eigenstates will be the column vector made out of the diagonalized matrix. Once the ground-state energy, $E_0$, is found, we will have its eigenstate $\left| \psi_0 \right>$.

Not all the values on the main diagonal are needed. Only the smallest ones are required. They why of this will be explained later.

The magnetization operator, $\mathbf{M}$, can be defined as a sum over all of the $\mathbf{A}_{i j}$ operators,

$$
    \mathbf{M} = \sum_{i, j = 0}^{N_i-1, N_j-1} \mathbf{A}_{i j},
$$

and as we have $\left| \psi_0 \right>$, the total magnetization at the ground state can be computed by sandwiching the $\mathbf{M}$ operator.

$$
    M_0 = \left< \psi_0 \right| \mathbf{M} \left| \psi_0 \right>
$$

In order to find the magnetization for higher temperatures, we can use the partition function, $Z$.

$$
    Z = \sum_n e^{-E_n / k_B T},
$$

where $E_n$ is one of the eigenvalues of the hamiltonian matrix, $k_B$ is Boltzman's constant, and $T$ is the temperature of the system.

The system might have numerous eigenvalues with large positive amounts. Consequently, their participation in the partition function will be negligible. They can be, therefore, omitted to reduce the runtime.

The helmholtz energy of the system, $F$, can be calculated as it is shown below.

$$
    F = - k_B T \ln Z
$$

The heat capacity of it is

$$
    C_V = k_B T \left( 2 \left. \frac{\partial \ln Z}{\partial T} \right|_V + T \left. \frac{\partial^2 \ln Z}{\partial T^2} \right|_V \right).
$$

However, in order to compute them numerically, we require the following approximations:

$$
    \frac{\partial \ln Z}{\partial T} = \lim_{\delta T \rightarrow 0} \frac{\ln Z(T + \delta T) - \ln Z(T)}{\delta T}
$$

$$
    \frac{\partial^2 \ln Z}{\partial^2 T} = \lim_{\delta T \rightarrow 0} \frac{\frac{\ln Z(T + 2\delta T) - \ln Z(T + \delta T)}{\delta T} - \frac{\ln Z(T + \delta T) - \ln Z(T)}{\delta T}}{\delta T}
$$

Finally, the total magnetization at the arbitrary temperature $T$ is

$$
    M = - \left. \frac{\partial F}{\partial \epsilon} \right|_T = \lim_{\delta \epsilon \rightarrow 0} \frac{F_\epsilon - F_{\epsilon + \delta \epsilon}}{\delta \epsilon}.
$$

### Monte Carlo method

The code is available in [montecarlo.py](./montecarlo.py).

## TODO

- [ ] Complete [Monte Carlo method](#monte-carlo-method) section.
- [ ] Add coding notes/experience.
- [ ] Check the correctness of the [diagonalization.py](./diagonalization.py) results in a small lattice.
- [ ] Add the [References](#references) section.