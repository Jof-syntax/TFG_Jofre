function SIMPALL(E1,E0,nu1,nu0);

%% SIMPALL(1,0.01,1/3,1/3)
mu=@(E,nu) E/(2*(1+nu));
kappa = @(E,nu) E/(2*(1-nu)); %valid en 3D?

mu0 = mu(E0,nu0); %mu
mu1 = mu(E1,nu1); %muPreDef

kappa0 = kappa( E0 , nu0 ); %kappa
kappa1 = kappa( E1 , nu1 ); %kPreDef

etamu = @(mu, kappa, muPreDef )  (mu-muPreDef)*(9*kappa+8*mu)/(5*muPreDef*(3*kappa+4*mu));
etakappa = @(mu,kappa, kPreDef) 4*mu*(kappa-kPreDef)/(kPreDef*(3*kappa+4*mu));


A = @(f0, f1, df0, df1) (-df1*df0+(f1-f0)^2)/(df1+f1-f0);
B = @(f0, f1, df0, df1) (2*f0*(f1-f0)+df1*df0-df1*f0-f1*df0)/(df1+f1-f0);
C = @(f0)               f0;
D = @(f0, f1, df0, df1) (2*(f1-f0)-(df1-df0))/(df1+f1-f0);
 

 f = @(rho, A, B, C, D) (A.*rho.^2+B.*rho+C)./(D.*rho+1);
 df = @(rho, A, B, C, D) (A*D.*rho.^2+2*A.*rho+B-C*D)./(D.*rho+1).^2;
 
 
df0 = @(eta0, f0, f1) (f0+eta0)*(f1-f0)/(f1+eta0);
df1 = @(eta1, f0, f1) (f1+eta1)*(f1-f0)/(f0+eta1);

%% mu case

F0 = mu0;
F1 = mu1;
Eta0 = etamu(mu0, kappa0, mu1);
Eta1 = etamu(mu1, kappa1, mu0);
DF0 = df0(Eta0, F0, F1);
DF1 = df1(Eta1, F0, F1);

AA = A(F0, F1, DF0, DF1);
BB = B(F0, F1, DF0, DF1);
CC = C(F0);
DD = D(F0, F1, DF0, DF1);


plot(linspace(0,1,100),f(linspace(0,1,100), AA, BB, CC, DD));



