import sympy as sp
from sympy import physics as spp
from sympy.physics import vector as sppv


# Set manufactured solution
vmax = sp.symbols('vmax')

t = sp.symbols('t')

R = sppv.ReferenceFrame('R')

x = R[0]

y = R[1]

u = [0, 0]

u[0] = vmax*sp.sin(2.*sp.pi*y)

u[1] = 0.

p = -u[0]**2


# Print in muparser format
print("\nManufactured solution:")

print(("Function expression = "+str(u[0])+"; "+str(u[1])+"; "+str(p)).replace('**', '^').replace('R_', ''))


# Derive manufactured source term
# @todo For now I'm not taking the symmetric part of the stress tensor, because this is complicating the symbolic implementation.
mu = sp.symbols('mu')

gamma = sp.symbols('gamma')

grad_p = sppv.gradient(p, R).to_matrix(R)

f0 = sp.diff(u[0], t) + sppv.dot(u[0]*R.x + u[1]*R.y, sppv.gradient(u[0], R)) - sppv.divergence(mu*sppv.gradient(u[0], R), R) + grad_p[0]

f1 = sp.diff(u[1], t) + sppv.dot(u[0]*R.x + u[1]*R.y, sppv.gradient(u[1], R)) - sppv.divergence(mu*sppv.gradient(u[1], R), R) + grad_p[1]

f2 = sppv.divergence(u[0]*R.x + u[1]*R.y, R) + gamma*p


# Print in muparser format
print("\nDerived manufactured source:")

print(("Function expression = "+str(f0)+"; "+str(f1)+"; "+str(f2)).replace('**', '^').replace('R_', ''))