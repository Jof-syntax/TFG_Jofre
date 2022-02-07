%% check topology derivative
clc;

mu = sym('mu','real');
mu2 = sym('mu2','real'); %[mu]
v = sym('v','real');
lam = sym('lam','real');
lam2 = sym('lam2','real'); %[lam]

%Elas3D
m1n = 15*mu*mu2*(v-1);
m1d = 15*mu*(1-v)+2*mu2*(5*v-4);
m1 = m1n/m1d;

m2n = lam2*(15*mu*lam*(1-v) + 2*lam*mu2*(5*v-4)) - 2*mu2*(lam*mu2-5*mu*v*lam2);
m2d = 5*mu2*(3*mu*lam*(1-v)-3*mu*lam*lam2-lam*mu2*(1-2*v));
m2 = m2n/m2d;

eq1FT = 4/3*m1*2;
eq1ST = 4/3*m1*m2;
disp('ELAS3D =')
eq1FT = simplify(eq1FT);
pretty(eq1FT)
eq1ST = simplify(eq1ST);
pretty(eq1ST)


% 5021

an = -5*mu*v*lam2+lam*mu2;
ad = 15*lam*mu*(1-v);
a = an/ad;

bn = 15*mu*(1-v)-2*mu2*(4-5*v);
bd = 15*mu*(1-v);
b = bn/bd;


eq2FT = -2/(3*b)*2*mu2;
eq2ST = -2/(3*b)*(lam2*b-2*mu2*a)/(3*a+b);
disp('5021 =')
eq2FT = simplify(eq2FT);
pretty(eq2FT)
eq2ST = simplify(eq2ST);
pretty(eq2ST)

%%
disp('Proportion = ')

eq1 = eq2FT/eq1FT;
eq1 = simplify(eq1);
pretty(eq1)
eq2 = eq2ST/eq1ST;
eq2 = simplify(eq2);
pretty(eq2)
