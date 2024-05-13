# Spinor-Evolution
It contains python code for numerically simulation time evolution of a gaussian spinor in a barrier potential
---
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
