import numpy as np
import scipy.sparse as sparse
import sympy


class Diagonalization:

    def __init__(
        self,
        lattice_shape: tuple[int, int],
        k_B: (int | float) = 1,
        J: (int | float) = 1,
        epsilon: (int | float) = 0,
        T: (int | float) = 0,
        delta_epsilon: float = 0.001,
        delta_T: float = 0.001

    ):
        self.lattice_shape = lattice_shape
        self.N_i = lattice_shape[0]
        self.N_j = lattice_shape[1]
        self.k_B = k_B
        self.J = J
        self.epsilon = epsilon
        self.T = T
        self.delta_epsilon = delta_epsilon
        self.delta_T = delta_T

        self.eig_vals, self.E_0, self.psi_0 = self.eigens()

    def n_identities(self, n: int):
        """Return a matrix made of the tensor product of n identity
        matrices.

        Parameters
        ----------
        n : int
            The number of the identity matrices. The least value can be
            `0`.

        Returns
        -------
        Any
            An integers or an identity matrix of shape `(2**n, 2**n)`.
        """
        return sparse.identity(2**n) if (n > 0) else 1

    def A_operators(self) -> np.ndarray:
        """Return the A operators of the lattice.

        The A operator of the site located in `(i, j)` can be created as
        follows:
            :math:`A_{ij} = I ⊗...⊗ I ⊗ σ^z ⊗ I ⊗...⊗ I`
        where :math:`i, j = 0, 1, ...`. The number of `I` matrices
        before :math:`σ^z` is :math:`N_j i + j`, while there are in
        total :math:`N_i N_j` matrices multiplied.

        Parameters
        ----------
        N_i : int
            The number of the rows.

        N_j : int
            The number of columns.

        Returns
        -------
        ndarray
            An array of shape `(N_i, N_j)`.
        """
        pauli_z = sparse.diags_array([1, -1])
        # A list to store the operators A
        A_ops = np.empty(self.lattice_shape,
                         dtype=sparse._coo.coo_matrix)
        for i in range(self.N_i):
            for j in range(self.N_j):
                A_ij = sparse.kron(sparse.kron(self.n_identities((self.N_j*i)+j),
                                               pauli_z),
                                   self.n_identities((self.N_j*self.N_i)
                                                     - (self.N_j*i)-j-1))
                A_ops[i, j] = A_ij
        return A_ops

    def hamiltonian(self):  # TODO: the class of output.
        """Return the hamiltonian of the system.

        The hamiltonian can be calculated as follows:
            :math:`H = ε Σ A_{ij} + (1/2) J Σ A_{ij} (A_{i+1,j} \
                + A_{i-1,j} + A_{i,j+1} + A_{i,j-1})`

        Parameters
        ----------
        N_i : int
            The number of rows.

        N_j : int
            The number of columns.

        A_ops : ndarray
            An array of A operators.

        epsilon : int or float
            The magnetic field exerted on the system.

        J : int or float
            The exchange energy.

        Returns
        -------
        csr_matrix # TODO: the class of output.
            An sparse matrix representing the Hamiltonian of shape
            `(2**(N_i*N_j), 2**(N_i*N_j))`.

        Raises
        ------
        Exception
            The neighboring site on the boundaries of the lattice was
            not found (or incorrect index).
        """
        A_ops = self.A_operators()
        # Creating the Hamiltonian
        H = sparse.dia_matrix((2**(self.N_i*self.N_j),
                               2**(self.N_i*self.N_j)),
                              dtype=float)
        # There used to be a BC_type variable used in a match clause,
        # but it had been causing issue.
        for i in range(self.N_i):
            for j in range(self.N_j):
                H += self.epsilon * A_ops[i, j]
                # Neighbors on the same row
                for index in [j+1, j-1]:
                    # If index is (i, -1), there will be no problem
                    try:
                        H += 0.5 * self.J * (A_ops[i, j] @ A_ops[i, index])
                    # IndexError means that the index is out of bound,
                    # so index 0 should be used.
                    except IndexError:
                        try:
                            H += 0.5 * self.J * (A_ops[i, j] @ A_ops[i, 0])
                        except:
                            raise Exception('failed to find element \
                                ({:.0f},{:.0f}) in A_ops'.format(i, index))
                # Neighbors on the same column
                for index in [i+1, i-1]:
                    # If index is (-1, j), there will be no problem
                    try:
                        H += 0.5 * self.J * (A_ops[i, j] @ A_ops[index, j])
                    # IndexError means that the index is out of bound,
                    # so index 0 should be used.
                    except IndexError:
                        try:
                            H += 0.5 * self.J * (A_ops[i, j] @ A_ops[0, j])
                        except:
                            raise Exception('failed to find element \
                                ({:.0f},{:.0f}) in A_ops'.format(index, j))
        return H

    def eigens(self):  # TODO: find the class of H_op.
        """Return half of the eigenvalues, the least eigenvalue (E_0)
        and its corresponding eigenvector.

        Parameters
        ----------
        H_op : csr_matrix # TODO: find the class.
            The hamiltonian of the system.
            If `None` be used, the hamiltonian of the system itself will
            be computed and used.

        Returns
        -------
        tuple[ndarray, float, ndarray]
            The eigenvalues, the least eigenvalue E_0, and its eigenvector.
        """
        # Getting the smallest eigenvalues and their corresponding
        # eigenvectors.
        # sparse.linalg.eigs(H_op, k=H_op.shape[0]//2, which='SR') could
        # be used but the results were inconsistent.
        # k=H_op.shape[0]//2 returns half of the eigenvalues of the system.
        eig_vals, eig_vecs = np.linalg.eig(self.hamiltonian().todense())
        eig_vals = np.real(eig_vals)
        # Turning column presentation of eig_vecs into row presentation.
        eig_vecs = eig_vecs.T
        # Removing unnecessary computed values.
        eig_vecs[abs(eig_vecs) < 1e-10] = 0
        eig_vals = np.round(eig_vals, decimals=10)  # TODO: Is it necessary?
        # Finding the lowest energy.
        E_0 = np.min(eig_vals)
        # Finding the state with the lowest energy.
        # Turning the row presentation to column presentation.
        psi_0 = eig_vecs[eig_vals == E_0].T
        # Checking the dimensions of psi_0; in case the degeneracy have
        # resulted in more than one eigenstate.
        _, psi_0_col = psi_0.shape
        if psi_0_col > 1:
            # Setting the first state as the ground state.
            psi_0 = psi_0[:, 0].reshape(-1, 1)
        return eig_vals, E_0, psi_0

    def magnetization_0(self) -> float:
        """Return the magnetization of the ground state at :math:`T = 0`.

        The magnetization operator can be created by summing the
        A operators.
            :math:`M = Σ A_{ij}, M_0 = <ψ_0| M |ψ_0>`
        The eigenvalue of the ground state can be calculated using it,
        as it follows:
            :math:`M_0 = <ψ_0| M |ψ_0>`

        Parameters
        ----------
        A_ops : ndarray
            An array of the A operators of the system.
            If `None` be used, the A operators of the system will be
            computed and used.

        psi_0 : ndarray
            The column eigenvector of the ground state.
            If `None` be used, the ground state of the system will be
            computed and used.

        Returns
        -------
        float
            Magnetization at :math:`T = 0`.
        """
        A_ops = self.A_operators()
        # A_ops.sum() is the magnetization operator
        M = self.psi_0.conjugate().T @ A_ops.sum() @ self.psi_0
        # Retrieving tne value from the M array
        if M.shape == (1, 1):
            M = np.real(M[0, 0])
        return M

    def partition(self, T: (float | int | None) = None):
        """Return the value of the partition function of a system.

        :math:`Z = Σ exp(-E_n / k_B T)`

        Parameters
        ----------
        eig_vals : array_like
            An array of eigenvalues.

        T : float
            Temperature.
            If `None`, the initial value will be used.

        Returns
        -------
        Any # TODO: find the type.
            The numerical value of the partition function.
        """
        if T == None:
            T = self.T
        Z = 0
        for energy in self.eig_vals:
            Z += sympy.exp(-energy / (self.k_B * T))  # type: ignore
        return Z

    def helmholtz(self) -> float:
        """Return the value of Helmholtz energy.

        :math:`F = -k_B T ln(Z)`

        Parameters
        ----------
        k_B :float
            Boltzman's constant.
        T : float
            Temperature.
        Z : Any # TODO: find the type.
            The value of the partition function.

        Returns
        -------
        float
            The Helmholtz energy.
        """
        return -self.k_B * self.T * float(sympy.log(self.partition()))

    def d_ln_Z_dT(self, T: (float | int | None) = None) -> float:
        """Return the numerical value of :math:`(d/dT)ln(Z)`.

        :math:`(∂/∂T) ln(Z) = lim_{δT→0} {ln(Z[T+δT]) - ln(Z[T])} / {δT}`

        Parameters
        ----------
        eig_vals (ndarray): eigenvalues (energy levels) of the system
        k_B (float): Boltzman's constant
        T (float): temperature
        delta_T (float): differentiation step

        Returns
        -------
        float
            The derivative of :math:`ln(Z)`, :math:`(d/dT)ln(Z)`.
        """
        if T == None:
            T = self.T
        return float(sympy.log(self.partition(T=T+self.delta_T))
                     - sympy.log(self.partition(T=T))) / self.delta_T  # type: ignore

    def d2_ln_Z_d2T(self, T: (float | int | None) = None) -> float:
        """Return the numerical value of :math:`(d^2/dT^2)ln(Z)`.

        Parameters
        ----------
        eig_vals : ndarray
            The eigenvalues (energy levels) of the system.
        T (float): temperature

        Returns
        -------
        float
            The second derivative of :math:`ln(Z)`, :math:`(d^2/dT^2)ln(Z)`.
        """
        if T == None:
            T = self.T
        return (self.d_ln_Z_dT(T=T+self.delta_T) - self.d_ln_Z_dT(T=T)) / self.delta_T

    def heat_cap_v(
        self,
        eig_vals: (np.ndarray | list | None) = None,
        T: (float | int | None) = None
    ) -> float:
        """Return the heat capacity of the system.

        :math:`C_v = k_B T (2 [∂ ln(Z) / ∂T]_V + T [∂^2 ln(Z) / ∂T^2]_V)`

        Parameters
        ----------
        eig_vals (ndarray): eigenvalues (energy levels) of the system
        k_B (float): Boltzman's constant
        T (float): temperature
        delta_T (float): differentiation step

        Returns
        -------
        float
            The heat capacity of the system at the temperature T.
        """
        if eig_vals == None:
            eig_vals, _, _ = self.eigens()
        if T == None:
            T = self.T
        return self.k_B * T * (2 * self.d_ln_Z_dT(T=T) + T * self.d2_ln_Z_d2T(T=T))

    def magnetization_T(self) -> float:
        """Return the magnetization of the system at temperature T.

        :math:`M = - [∂F/∂ε]_T = lim_{δT→0} {[F(ε) - F(ε+δε)] / δε}`

        Parameters
        ----------
        N_i (int): number of rows
        N_j (int): number of columns
        k_B (float): Boltzman's constant
        epsilon (float): intensity of the magnetic field
        delta_epsilon (float): differentiation step
        J (float): exchange energy
        T (float): temperature

        Returns
        -------
        float
            The magnetization of the system at the temperature T.
        """
        lattice_2 = Diagonalization(lattice_shape=self.lattice_shape,
                                    k_B=self.k_B,
                                    J=self.J,
                                    epsilon=self.epsilon+self.delta_epsilon,
                                    T=self.T,
                                    delta_epsilon=self.delta_epsilon,
                                    delta_T=self.delta_T)
        eig_vals_2 = lattice_2.eig_vals
        # Passing eigenvalues to the partition function and calculating F
        F_1 = self.helmholtz()
        F_2 = lattice_2.helmholtz()
        return (F_1 - F_2) / self.delta_epsilon
