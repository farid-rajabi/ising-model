import numpy as np
import matplotlib.pyplot as plt


class MonteCarlo:
    """
    MonteCarlo
    ==========

    Monte Carlo class assists with applying the Metropolis algorithm on
    the given lattice.
    """

    def __init__(
        self,
        lattice: np.ndarray,
        J: (int | float | np.ndarray) = 1,
        epsilon: float = 0,
        k_B: float = 1,
        T: float = 0
    ) -> None:
        self.lattice = lattice
        self.dimen = lattice.shape
        self.J = J
        self.epsilon = epsilon
        self.k_B = k_B
        self.T = T

    def summary(self) -> tuple[int, int, int]:
        """Return a summary about the spins in the given lattice.

        Returns
        -------
            tuple[int, int, int]
                The number of positive spins, the number of negative
                spins, and the total spin.
        """
        return (self.lattice > 0).sum(), (self.lattice < 0).sum(), \
            self.lattice.sum()

    def description(self) -> None:
        """Print a description about the spins in the given lattice.

        The output contains:
            - the size of the lattice and its dimensions,
            - the number of positive ones and negative ones,
            - the total spin of the lattice,
            - and an illustration of the spin of the sites.
        """
        print('Num. of spins: {:.0f} ({:.0f}, {:.0f})'.format(
            self.lattice.size, self.dimen[0], self.dimen[1]))
        print('Num. of positive spins: {:.0f}'.format(
            (self.lattice > 0).sum()))
        print('Num. of negative spins: {:.0f}'.format(
            (self.lattice < 0).sum()))
        print('Total spin: {:.0f}'.format(self.lattice.sum()))
        plt.imshow(X=self.lattice, cmap='bwr')
        plt.show()

    def E_flip(
        self,
        coord: tuple[int, int],
        site_spin: (int | float),
        neighbor_spin: np.ndarray
    ) -> (int | float):
        """Return the energy required to flip the spin.

        The energy required of the spin in its current state is:
            :math:`E_1 = J Σ s_{ij} s_{ij neighbor + ε s_{ij}`.
        But if it be flipped, the energy becomes
            :math:`E_2 = - E_1`.
        The flipped energy is
            :math:`E_flip = E_2 - E_1 = -2 E_1`.

        Parameters
        ----------
        site_spin : int or float
            The spin of the site for which E_flip is being calculated.

        neighbor_spin : ndarray
            The list-like object containing the spin of the neighboring
            sites.

        Returns
        -------
        int or float
            The flipping energy of the spin.
        """
        return -2 * (site_spin * (self.J * neighbor_spin).sum() +
                     self.epsilon * site_spin)

    def flip_decider(
        self,
        E_flip: (int | float)
    ) -> bool:
        """Returns the boolean specifying whether the spin should be
        flipped.

        Parameters
        ----------
        E_flip : int or float
            The energy required to flip the spin.
            It can be calculated for any site using the `E_flip` method.

        Returns
        -------
        bool
            Whether the spin should be flipped or not.
        """
        if E_flip < 0:
            return True
        elif np.random.ranf() <= np.exp(-E_flip / (self.k_B * self.T)):
            return True
        else:
            return False

    def sweep(
        self,
        BC: str
    ):
        """Sweep the whole lattice once.
        """
        # ! USING match RESULTS IN ERROR
        if BC == 'OBC':
            for n_i in range(self.dimen[0]):
                for n_j in range(self.dimen[1]):
                    spin_ij = self.lattice[n_i, n_j]
                    spin_neighbor = []
                    # Left neighbor
                    if n_j == 0:
                        # Set the spin of the left neighbor equal to zero,
                        # because we want to impact from beyond the closed
                        # boundaries.
                        spin_neighbor.append(0)
                    else:
                        spin_neighbor.append(self.lattice[n_i, n_j-1])
                    # Right neighbor
                    if n_j == (self.dimen[1] - 1):
                        spin_neighbor.append(0)
                    else:
                        spin_neighbor.append(self.lattice[n_i, n_j+1])
                    # Up neighbor
                    if n_i == 0:
                        spin_neighbor.append(0)
                    else:
                        spin_neighbor.append(self.lattice[n_i-1, n_j])
                    # Down neighbor
                    if n_i == (self.dimen[0] - 1):
                        spin_neighbor.append(0)
                    else:
                        spin_neighbor.append(self.lattice[n_i+1, n_j])
                    # Deciding whether to flip or not
                    flip_energy = self.E_flip((n_i, n_j),
                                              spin_ij,
                                              np.array(spin_neighbor))
                    if self.flip_decider(flip_energy):
                        self.lattice[n_i, n_j] *= -1
        elif BC == 'PBC':
            for n_i in range(self.dimen[0]):
                for n_j in range(self.dimen[1]):
                    spin_ij = self.lattice[n_i, n_j]
                    spin_neighbor = []
                    # Left neighbor
                    spin_neighbor.append(self.lattice[n_i, n_j-1])
                    # Right neighbor
                    try:
                        spin_neighbor.append(self.lattice[n_i, n_j+1])
                    except:
                        spin_neighbor.append(self.lattice[n_i, 0])
                    # Up neighbor
                    spin_neighbor.append(self.lattice[n_i-1, n_j])
                    # Down neighbor
                    try:
                        spin_neighbor.append(self.lattice[n_i+1, n_j])
                    except:
                        spin_neighbor.append(self.lattice[0, n_j])
                    # Deciding whether to flip or not
                    flip_energy = self.E_flip((n_i, n_j),
                                              spin_ij,
                                              np.array(spin_neighbor))
                    if self.flip_decider(flip_energy):
                        self.lattice[n_i, n_j] *= -1
        else:
            raise Exception('`BC` can get values of `OBC` or `PBC`.')

    def magnetization(
        self,
        n: int = 1,
        BC: str = 'PBC'
    ) -> float:
        """Return the average magnetization of the lattice over `n`
        sweeps.

        Parameters
        ----------
        n : int, optional # TODO: check optional parameter doc format.
            The number of the repetition of the sweep.

        Returns
        -------
        float

        Raises
        ------
        Exception
        """
        n = int(n)
        if n < 1:
            raise Exception(
                '\'n\' should be a positive integer larger than 0.')
        M = 0
        for _ in range(n):
            self.sweep(BC)
            M += self.summary()[2]
        return M / n
