import math

import numpy as np
import sympy as sp
import scipy.sparse as sparse
import matplotlib.pyplot as plt
from matplotlib import cm
import matplotlib.animation as animation

x, y, t = sp.symbols('x,y,t')

class Wave2D:

    def create_mesh(self, N, sparse=False):
        """Create 2D mesh and store in self.xij and self.yij"""
        self.xij, self.yij = np.meshgrid(
            np.linspace(0, 1, N+1),
            np.linspace(0, 1, N+1),
            indexing='ij',
            sparse=sparse
        )
        self.N = N
        self.h = 1 / self.N

    def D2(self, N):
        """Return second order differentiation matrix"""
        D = sparse.diags([1, -2, 1], [-1, 0, 1], (N+1, N+1), 'lil')
        D[0, :4] = 2, -5, 4, -1
        D[-1, -4:] = -1, 4, -5, 2
        return D

    @property
    def w(self):
        """Return the dispersion coefficient"""
        m = math.sqrt(self.mx**2 + self.my**2)
        k = m * math.pi
        return self.c * k        

    def ue(self, x, y, t):
        """Return the exact standing wave"""
        mx, my, w = self.mx, self.my, self.w
        return np.sin(mx*np.pi*x)*np.sin(my*np.pi*y)*np.cos(w*t)

    def initialize(self, N, mx, my, c, dt, h, sparse_mesh=False):
        r"""Initialize the solution at $U^{n}$ and $U^{n-1}$

        Parameters
        ----------
        N : int
            The number of uniform intervals in each direction
        mx, my : int
            Parameters for the standing wave
        """
        self.create_mesh(N, sparse=sparse_mesh)
        
        U0 = self.ue(self.xij, self.yij, 0)
        D = self.D2(N) / h**2
        U1 = U0 +  (c * dt)**2 / 2 * (D @ U0 + U0 @ D.T)
        U1 = self.apply_bcs(U1)
        return U1, U0

    
    def l2_error(self, u, t0):
        """Return l2-error norm

        Parameters
        ----------
        u : array
            The solution mesh function
        t0 : number
            The time of the comparison
        """
        ue_evaluated = self.ue(self.xij, self.yij, t0)
        return math.sqrt(self.h**2 * np.sum((u - ue_evaluated).flatten()**2))

    def apply_bcs(self, U):
        box = np.zeros((self.N+1, self.N+1))
        box[1:-1, 1:-1] = 1
        return U * box

    
    def evolve(self, Un, Unm1, c, dt, N, h):
        D = self.D2(N) / h**2
        return 2 * Un - Unm1 + (c * dt)**2 * (D @ Un + Un @ D.T)

    def __call__(self, N, Nt, cfl=0.5, c=1.0, mx=3, my=3, store_data=-1):
        """Solve the wave equation

        Parameters
        ----------
        N : int
            The number of uniform intervals in each direction
        Nt : int
            Number of time steps
        cfl : number
            The CFL number
        c : number
            The wave speed
        mx, my : int
            Parameters for the standing wave
        store_data : int
            Store the solution every store_data time step
            Note that if store_data is -1 then you should return the l2-error
            instead of data for plotting. This is used in `convergence_rates`.

        Returns
        -------
        If store_data > 0, then return a dictionary with key, value = timestep, solution
        If store_data == -1, then return the two-tuple (h, l2-error)
        """
        h = 1 / N
        dt = cfl * h / c
        self.c, self.dt = c, dt
        self.mx, self.my = mx, my
        
        U1, U0 = self.initialize(N, mx, my, c, dt, h)
        if store_data == 1:
            data = {0: U0, 1: U1}
        elif store_data > 1:
            data = {0: U0}

        Un, Unm1 = U1, U0
        errors = [self.l2_error(Ui, ti) for Ui, ti in zip([Unm1, Un], [0, dt])]
        for n in range(1, Nt):
            Unp1 = self.evolve(Un, Unm1, c, dt, N, h)
            Unp1 = self.apply_bcs(Unp1)
            errors.append(self.l2_error(Unp1, (n+1)*dt))

            if store_data > 0 and (n+1) % store_data == 0:
                data[n+1] = Unp1
            Un, Unm1 = Unp1, Un


        if store_data > 0:
            return data
        else:
            return self.h, errors         

    def convergence_rates(self, m=4, cfl=0.1, Nt=10, mx=3, my=3):
        """Compute convergence rates for a range of discretizations

        Parameters
        ----------
        m : int
            The number of discretizations to use
        cfl : number
            The CFL number
        Nt : int
            The number of time steps to take
        mx, my : int
            Parameters for the standing wave

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
            dx, err = self(N0, Nt, cfl=cfl, mx=mx, my=my, store_data=-1)
            E.append(err[-1])
            h.append(dx)
            N0 *= 2
            Nt *= 2
        r = [np.log(E[i-1]/E[i])/np.log(h[i-1]/h[i]) for i in range(1, m+1, 1)]
        return r, np.array(E), np.array(h)

class Wave2D_Neumann(Wave2D):

    def D2(self, N):
        """Return second order differentiation matrix with Neumann boundary conditions"""
        D = sparse.diags([1, -2, 1], [-1, 0, 1], shape=(N+1, N+1), format='lil')

        # we use finite difference u'(0) \approx (u_1 - u_(-1)) /(2h), and set u'(0) = 0
        # This gives u(-1) = u(1). Inserting into finite difference for u''(0) gives
        D[0, 0:2] = -2, 2

        # Same for the right boundary
        D[N, N-1:] = 2, -2

        return D

    def apply_bcs(self, U):
        """boundary conditions are encoded in the matrix D2. This is just the identity."""
        return U

    def ue(self, x, y, t):
        mx, my, w = self.mx, self.my, self.w
        return np.cos(np.pi * mx * x ) * np.cos( np.pi * my * y) * np.cos(w*t) 


def test_convergence_wave2d():
    sol = Wave2D()
    r, E, h = sol.convergence_rates(mx=2, my=3)
    assert abs(r[-1]-2) < 1e-2


def test_convergence_wave2d_neumann():
    solN = Wave2D_Neumann()
    r, E, h = solN.convergence_rates(mx=2, my=3)
    assert abs(r[-1]-2) < 0.05


def test_exact_wave2d():
    sol = Wave2D()
    r, E, h = sol.convergence_rates(mx=2, my=2, cfl = 1/np.sqrt(2))
    for Ei in E:
        assert Ei < 1e-12

    solN = Wave2D_Neumann()
    rN, EN, hN = solN.convergence_rates(mx=2, my=2, cfl = 1/np.sqrt(2))
    for Ei in EN:
        assert Ei < 1e-12




def generate_animation(use_neumann_bc=False):

    if use_neumann_bc:
        sol = Wave2D_Neumann()
        name = "Neumann"
    else:
        sol = Wave2D()
        name = "Dirichlet"

    N = 50
    Nt = 100
    cfl = 1/ np.sqrt(2)

    my = 2
    mx = 2
    store_data = 1

    data = sol(N, Nt, cfl=cfl, mx=mx, my=my, store_data=store_data)

    xij, yij = sol.xij, sol.yij

    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.set_zlim(-1, 1)
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.set_zlabel('u')
    

    frames = sorted(data.keys())

    surf = [ax.plot_surface(xij, yij, data[frames[0]], cmap='viridis')]

    def update_plot(frame_number):
        surf[0].remove()
        surf[0] = ax.plot_surface(xij, yij, data[frame_number], cmap='viridis')
        ax.set_title(name)
        return surf[0],

    ani = animation.FuncAnimation(
        fig, update_plot, frames=frames, interval=50, blit=False
    )

    ani.save('wave_eq_animation_' + name + '.gif', writer='ffmpeg', fps=10)
    plt.show()

if __name__ == "__main__":
    generate_animation(False)
    generate_animation(True)
