# Spinor-Evolution
This repository contains code for numerically simulating time evolution of a gaussian spinor in a barrier potential

### Theory :
---

The Hamiltonian of the system containing spin orbit coupling is as follows

$$ \hat{H} = -\dfrac{\hbar^2}{2m}\dfrac{\hat{d}^2}{dz^2} + \hat{V}(z) + \dfrac{\hbar^2}{4m^2c^2}\dfrac{d\hat{V}(z)}{dz}\hat{\vec{s}}\times\hat{\vec{p}}\cdot\nabla V $$

If we expand spin using Pauli Matrices $s=\hbar\mathbf{\sigma}$

$$ \hat{H} = -\dfrac{\hbar^2}{2m}\dfrac{\hat{d}^2}{dz^2} + \hat{V}(z) + \dfrac{\hbar^2}{4m^2c^2}\dfrac{d\hat{V}(z)}{dz}\hat{K} $$

The $\hat{K}$ matrix defined has the following form :

$$ \hat{K} = \begin{pmatrix} 0 & k^+ \\ k^- & 0 \end{pmatrix} $$

Here $k^+=k_x+ik_y$ and $k^-=k_x-ik_y$. Then they move to define the spinor wavefunction as, 

$$ \psi(z, t) = \begin{pmatrix} \phi_{\uparrow}(z,t) \\ \phi_{\downarrow}(z,t) \end{pmatrix} $$

### Implementation :
---

If we follow the central difference method for dealing with the partial derivatives in the Hamiltonian terms, then we end up with a recursion relation, which relates the spinor wavefunction $\mathbf{\psi}^n(z)$ with wavefunction at the next time step $\mathbf{\psi}^{n+1}(z)$.

I have provided the form of central difference method in the context of eigen-value equation $H\psi = E\psi$. Here I have used the $j$ as the spatial index and $n$ as the the temporal index.

$$\left( \dfrac{\partial^2}{\partial^2z} - V_j\right)\psi^n_j = \dfrac{1}{\epsilon^2} (\psi^n_{j+1} - 2\psi^n_{j} + \psi^n_{j-1}) - V_j\psi^n_j$$

Some constants that are used are defined here, $\epsilon$ and $\delta$ are the space and the time discretization parameters which are referred to as <b>del_z</b> and <b>del_t</b> in the program.
 
$$\Gamma = \dfrac{4m\epsilon^2}{\hbar^2} \qquad \lambda = \dfrac{2m\epsilon^2}{\hbar\delta} \qquad \rho = \dfrac{\epsilon^2}{2mc^2} $$

The implementation of the discretized time evolution equation is done by using the split operator method which ensures the preservation of norm of the wavefunction at all time steps, the form is given below.

$$ e^{-i\delta \hat{H}} = \dfrac{\hat{I}-\dfrac{i\delta \hat{H}}{2}}{\hat{I}+\dfrac{i\delta \hat{H}}{2}} $$

