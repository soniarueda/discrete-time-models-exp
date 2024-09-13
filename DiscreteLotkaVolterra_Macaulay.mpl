#INSTRUCTIONS: Read DiffResM13.txt
# 
# Symbolic-numeric algorithm for parameter estimation in discrete-time models with exp
# 
# Example 12: (Discrete Lotka-Volterra models)
# 
# Case n=2, N=2
# 
# Dynamic Taylor

# Degree of Taylor expansion N=2.

N:=2; 
epoly := proc(f)
   eval(convert(taylor(exp(x), x, N+1), polynom), x = f);
end proc:
epoly(y+cpa);
# Output of Algorithm 3 for Lotka-Volterra model, Taylor expansion at 0.
# 
F1:=x1[t+1]- x1[t]*epoly(r1-a11*x1[t]-a12*x2[t]);
F2:=x2[t+1]- x2[t]*epoly(r2-a21*x1[t]-a22*x2[t]);
# Shift to obtain a square polynomial system. Case of disjoint parameter sets
# 
# Equation (34) for i=1 and i=2.
# 
with(LREtools):
sF1:=subs(t=0,[F1,seq(shift(F1,t,i),i=1..2)]); 
indets(F1);
sF2:=subs(t=0,[F2,seq(shift(F2,t,i),i=1..2)]);
indets(F2);
# Macaulay resultant, equation (37) and (38).
F10:=expand(sF1[1]);
F11:=expand(sF1[2]);
F12:=expand(sF1[3]);
indets([F10,F11,F12]);
# Compute the Macaulay resultant of F10, F11, F12 with respect to a11, a12. 
# 
#read(`..DiffResM13.txt`):
# 
# The commands in the DR package will compute Macaulay resultants in the case of order zero polynomials, list of orders [0,0,0] in this occasion.
# 
Y:=[alg_pars([a11,a12],[0,0,0])];#algebraic polynomials of orders [0,0,0]
R12:=diffres_matrices([F10,F11,F12],[a11,a12],[0,0,0],[],Y);
# 
# R(r1) is the quotient of determinants det(M0)/det(A).
# 
# 
with(linalg):
M0:=evalm(R12[1]);
A:=evalm(R12[2]);
dA:=factor(det(A));
# We compute next equation (38), the leading coefficient of R(r1), the coefficient of r1^8.
# 
# 
[coeffs(det(M0),r1,'m1')]:m1;
# 
d:=degree(det(M0),r1);
LC:=factor(coeff(det(M0),r1^d));
LC/dA;
# 
# Finally Gamma8 equals
# 
Gamma8:=LC/dA;
# Macaulay resultant, equation (36).
# Polynomials in  (35) for N=2 and i=1, that is the homogeneization of F10 and F11, at h=0.
# 
F10h:=numer(subs({a11=a11/h,a12=a12/h},F10));
bF10:=subs(h=0,F10h);
F11h:=numer(subs({a11=a11/h,a12=a12/h},F11));
bF11:=subs(h=0,F11h);
# 
# The Macauly resultant is a univariate resultant is this case
# 
factor(resultant(subs(a11=1,bF10),subs(a11=1,bF11),a12));
# Jacobian and equation (39)
# 
# The Jacobian of F10, F11, F12 with respect to all the parameters r1, a11, a12
# 
with(VectorCalculus):
cF1:=[F10,F11,F12];
J:=Jacobian(cF1,[r1,a11,a12]);
# Independent term of the determinant of the Jacobian, equation (39)
# 
factor(subs({r1=0,a11=0,a12=0},det(J))); 
