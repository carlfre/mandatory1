import math

import numpy as np
import sympy as sp
import scipy.sparse as sparse
from scipy.sparse.linalg import spsolve
from scipy.interpolate import interpn

x, y = sp.symbols('x,y')

class Poisson2D:
    r"""Solve Poisson's equation in 2D::

        \nabla^2 u(x, y) = f(x, y), in [0, L]^2

    where L is the length of the domain in both x and y directions.
    Dirichlet boundary conditions are used for the entire boundary.
    The Dirichlet values depend on the chosen manufactured solution.

    """

    def __init__(self, L, ue):
        """Initialize Poisson solver for the method of manufactured solutions

        Parameters
        ----------
        L : number
            The length of the domain in both x and y directions
        ue : Sympy function
            The analytical solution used with the method of manufactured solutions.
            ue is used to compute the right hand side function f.
        """
        self.L = L
        self.ue = ue
        self.f = sp.diff(self.ue, x, 2)+sp.diff(self.ue, y, 2)

    def create_mesh(self, N):
        """Create 2D mesh and store in self.xij and self.yij"""
        self.xij, self.yij = np.meshgrid(
            np.linspace(0, self.L, N+1),
            np.linspace(0, self.L, N+1)
        )
        self.N = N
        self.h = self.L / self.N

    def D2(self):
        """Return second order differentiation matrix"""
        N, h = self.N, self.h
        D = sparse.diags([1, -2, 1], [-1, 0, 1], (N+1, N+1), 'csr')
        D[0, :4] = 2, -5, 4, -1
        D[-1, -4:] = -1, 4, -5, 2
        return (D / h**2)

    def laplace(self):
        """Return vectorized Laplace operator"""
        N = self.N
        D = self.D2() # D2x equals D2y in this case.
        term1 = sparse.kron(D, sparse.eye(N + 1))
        term2 = sparse.kron(sparse.eye(N + 1), D)
        return term1 + term2

    def get_boundary_indices(self) -> np.array:
        """Return indices of vectorized matrix that belongs to the boundary"""
        N = self.N
        indices = np.array(
            [(0, i) for i in range(N)]
            + [(i, N) for i in range(N)]
            + [(N, i) for i in range(1, N+1)]
            + [(i, 0) for i in range(1, N+1)]
        )
        return np.ravel_multi_index((indices[:, 0], indices[:, 1]), (N+1, N+1))

    def assemble(self) -> tuple[np.ndarray, np.ndarray]:
        """Return assembled matrix A and right hand side vector b"""
        x_flattened, y_flattened = self.xij.flatten(), self.yij.flatten()
        f_func_vectorized = sp.lambdify((x, y), self.f, 'numpy')

        A = self.laplace()
        b = f_func_vectorized(x_flattened, y_flattened)

        #insert boundary conditions
        boundary_indices = self.get_boundary_indices()
        A[boundary_indices, :] = 0
        A[boundary_indices, boundary_indices] = 1

        ue_func_vectorized = sp.lambdify((x, y), self.ue, 'numpy')
        b[boundary_indices] = ue_func_vectorized(x_flattened[boundary_indices], y_flattened[boundary_indices])
        return A, b
    
    def l2_error(self, u: np.ndarray) -> float:
        """Return l2-error norm"""
        ue_func = sp.lambdify((x, y), self.ue, 'numpy')
        ue_evaluated = ue_func(self.xij, self.yij)
        return np.sqrt(self.h**2 * np.sum(((u - ue_evaluated)**2).flatten()))

    def __call__(self, N: int) -> np.ndarray:
        """Solve Poisson's equation.

        Parameters
        ----------
        N : int
            The number of uniform intervals in each direction

        Returns
        -------
        The solution as a Numpy array

        """
        self.create_mesh(N)
        A, b = self.assemble()
        self.U = spsolve(A, b.flatten()).reshape((N+1, N+1))
        return self.U

    def convergence_rates(self, m: int = 6) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        """Compute convergence rates for a range of discretizations

        Parameters
        ----------
        m : int
            The number of discretization levels to use

        Returns
        -------
        3-tuple of arrays. The arrays represent:
            0: the orders
            1: the l2-errors
            2: the mesh sizes
        """
        E = []
        h = []
        N0 = 8
        for m in range(m):
            u = self(N0)
            E.append(self.l2_error(u))
            h.append(self.h)
            N0 *= 2
        r = [np.log(E[i-1]/E[i])/np.log(h[i-1]/h[i]) for i in range(1, m+1, 1)]
        return r, np.array(E), np.array(h)

    def eval(self, x: float, y: float) -> float:
        """Return u(x, y) using interpn over the local square"""
        x_indices = np.array([math.floor(x / self.h), math.floor(x / self.h) + 1])
        y_indices = np.array([math.floor(y / self.h), math.floor(y / self.h) + 1]) 

        x_vals = x_indices * self.h
        y_vals = y_indices * self.h

        interpolation_points = (x_vals, y_vals)
        U_vals = self.U[y_indices.reshape(-1, 1), x_indices].T
        xy = np.array([[x, y]])
        interpolated_val = interpn(interpolation_points, U_vals, xy)[0]
        return interpolated_val


def test_convergence_poisson2d():
    # This exact solution is NOT zero on the entire boundary
    ue = sp.exp(sp.cos(4*sp.pi*x)*sp.sin(2*sp.pi*y))
    sol = Poisson2D(1, ue)
    r, E, h = sol.convergence_rates()
    assert abs(r[-1]-2) < 1e-2


def test_interpolation():
    ue = sp.exp(sp.cos(4*sp.pi*x)*sp.sin(2*sp.pi*y))
    sol = Poisson2D(1, ue)
    U = sol(100)
    assert abs(sol.eval(0.52, 0.63) - ue.subs({x: 0.52, y: 0.63}).n()) < 1e-3
    assert abs(sol.eval(sol.h/2, 1-sol.h/2) - ue.subs({x: sol.h/2, y: 1-sol.h/2}).n()) < 1e-3

