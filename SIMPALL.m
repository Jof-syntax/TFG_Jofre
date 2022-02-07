function SIMPALL(E1,E0,nu1,nu0);

%% SIMPALL(1,0.01,1/3,1/3)
mu=@(E,nu) E/(2*(1+nu));
kappa = @(E,nu) E/(2*(1-nu)); %valid en 3D?

mu0 = mu(E0,nu0); %mu
mu1 = mu(E1,nu1); %muPreDef

kappa0 = kappa( E0 , nu0 ); %kappa
kappa1 = kappa( E1 , nu1 ); %kPreDef

%etamu = @(mu, kappa)  mu*(9*kappa+8*mu)/(6*kappa+12*mu); %3D
etamu = @(mu, kappa) kappa*mu/(2*mu+kappa); % 2D
 
%etakappa = @(mu, kappa) 4*mu/3;

A = @(f0, f1, df0, df1) -(-f0^2+2*f0*f1-f1^2+df0*df1) /(df1-f1+f0);
B = @(f0, f1, df0, df1) -(df0*f1-df0*df1+df1*f0-2*f0*f1+2*f0^2)/(df1-f1+f0);
C = @(f0)               f0;
D = @(f0, f1, df0, df1) -(df0+df1+2*f0-2*f1)/(df1-f1+f0);
 

 f = @(rho, A, B, C, D) (A.*rho.^2+B.*rho+C)./(D.*rho+1); %SIMPALL
 df = @(rho, A, B, C, D) (A*D.*rho.^2+2*A.*rho+B-C*D)./(D.*rho+1).^2;
 
 
df0 = @(eta, f0, f1) (f0+eta)*(f1-f0)/(f1+eta);
df1 = @(eta, f0, f1) (f1+eta)*(f1-f0)/(f0+eta);



F0 = mu0;
F1 = mu1;
Eta0 = etamu(mu0, kappa0);
Eta1 = etamu(mu1, kappa1);
DF0 = df0(Eta0, F0, F1);
DF1 = df1(Eta1, F0, F1);

AA = A(F0, F1, DF0, DF1);
BB = B(F0, F1, DF0, DF1);
CC = C(F0);
DD = D(F0, F1, DF0, DF1);


plot(linspace(0,1,100),f(linspace(0,1,100), AA, BB, CC, DD), 'r');
hold on

SIMP = @(rho) (1-rho.^3).*F0+rho.^3.*F1;

plot(linspace(0,1,100),SIMP(linspace(0,1,100)), 'b');




