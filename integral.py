#caso axial
from sympy import *
from sympy.matrices import Matrix

x = S("x")
L = S("L")
δv_0 = S("δv_0")
δv_3 = S("δv_3")

δv_1 = S("δv_1")
δv_2 = S("δv_2")
δv_4 = S("δv_4")
δv_5 = S("δv_5")

qx =  S("qx")
qy = S("qy")

densidad = 2500 #kg/m3
area = 1
#axial shape functions
N0 = (1-x/L)
N3 = x/L

#bending shape functions
N1 = 1 - 3*x**2/L**2 + 2*x**3/L**3
N2 = x - 2*x**2/L + x**3/L**2
N4 = 3*x**2/L**2 - 2*x**3/L**3
N5 = -x**2/L + x**3/L**2

ux_axial_function = (N0*δv_0 + N3*δv_3)*qx
uy_bending_function = (N1*δv_1 + N2*δv_2 + N4+δv_4 + N5*δv_5)*qy

axial_integral = integrate(ux_axial_function,(x,0,L))
bending_integral = integrate(uy_bending_function,(x,0,L))



print(pretty(axial_integral))
print("\n")
print(pretty(bending_integral))














































