
# 
# DISCRETE LOTKA-VOLTERRA MODEL
# 
# n=2, N=3=m-1
# 
restart:
# Dynamic Taylor
m:=3;# N=m-1 degree of Taylor approximation 
epoly := proc(f)
   eval(convert(taylor(exp(x), x, m), polynom), x = f);
end proc:
epoly(y+cpa);
vr1:=1;vr2:=1;
va11:=1/2;va12:=3/2;va21:=3/2; va22:=1/2
;


X100 := 2.00:
X200 := 3.00:

data_gen := proc(X10, X20)
  local X11, X21:
  X11 := X10*exp(vr1-va11*X10-va12*X20);
  X21 := X20*exp(vr2-va21*X10-va22*X20);
  return X11, X21
end proc:
DD := (data_gen@@5)(X100, X200);

x1:=X100:x2:=X200: N:=19:
for i from 1 to N do DD:=(data_gen@@i)(X100, X200); x1:=x1,DD[1];x2:=x2,DD[2]; od:
x1:=[x1];x2:=[x2];

f1:=X1[t+1]- X1[t]*epoly(r1-a11*X1[t]-a12*X2[t]);
f2:=X2[t+1]- X2[t]*epoly(r2-a21*X1[t]-a22*X2[t]);

NULL;
# Disjoint parameter sets
with(LREtools):
#######################################################
;
#x1,x2 known.
;
F1:=subs(t=0,[f1,seq(shift(f1,t,i),i=1..2)]); F2:=subs(t=0,[f2,seq(shift(f2,t,i),i=1..2)]);
indets(F1);
indets(F2);
NULL;
# Macaulay resultant of F10, F11, F12
read(`C:/Users/srued/OneDrive - Universidad Politécnica de Madrid/Sonia/Software/PROFESIONAL/INVESTIGACION/diff alg/DR/DiffResM13.txt`):
with(linalg):
F10:=expand(F1[1]);
F11:=expand(F1[2]);
F12:=expand(F1[3]);
indets([F10,F11,F12]);
Y:=[alg_pars([a11,a12],[0,0,0])];#algebraic polynomials of orders [0,0,0]
R12:=diffres_matrices([F10,F11,F12],[a11,a12],[0,0,0],[],Y);
M11:=evalm(R12[1]);
M12:=evalm(R12[2]);
#Macaulay Resultant:=det(M11)/det(M12);

dM12:=factor(det(M12));
cc1:=[coeffs(det(M11),r1,'m1')]:m1;
LC:=factor(cc1[nops(cc1)]);
#Leading coefficient in r1 of the Macaulay resultant:
LC/dM12;
NULL;
# Macaulay resultant of barF10 and barF11
# 
F1;
F10h:=numer(subs({a11=a11/h,a12=a12/h},F1[1]));
bF10:=subs(h=0,F10h);
F11h:=numer(subs({a11=a11/h,a12=a12/h},F1[2]));
bF11:=subs(h=0,F11h);
factor(resultant(subs(a11=1,bF10),subs(a11=1,bF11),a12));

NULL;
# Jacobian
# 
with(VectorCalculus):
cF1:=[F1[1],F1[2],F1[3]];
J:=Jacobian(cF1,[r1,a11,a12]);
factor(det(J));
ccJ:=[coeffs(det(J),[r1,a11,a12],'mJ')]:mJ;
factor(ccJ[nops(ccJ)]);
factor(ccJ[4]);
NULL;
# 
# 
# 
NULL;
NULL;
