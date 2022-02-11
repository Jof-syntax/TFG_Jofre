function SIMPALL(E1,E0,nu1,nu0, dimension)


%% mirar si expressio de kappa al inici es correcta


%% HS valides en 3D?


%% SIMPALL(1,0.01,1/3,1/3, N)         A
%% SIMPALL(1,0.001, 0, 1/3, N)        B
%% SIMPALL(1,0.01, -0.5, 1/3, N)      C
%% SIMPALL(1,0.001, -0.75, 1/3, N)    D
%% SIMPALL(1,0.001, -0.9, 1/3, N)     E
%% Code
mu=@(E,nu) E/(2*(1+nu));
kappa = @(E,nu) E/(3*(1-2*nu)); %expressio correcta?

mu0 = mu(E0,nu0); %mu
mu1 = mu(E1,nu1); %muPreDef

kappa0 = kappa( E0 , nu0 ); %kappa
kappa1 = kappa( E1 , nu1 ); %kPreDef

if dimension == 2
    etamu = @(mu, kappa) kappa*mu/(2*mu+kappa); % 2D
    etakappa = @(mu, kappa) mu; %2D
elseif dimension  == 3
    etamu = @(mu, kappa)  mu*(9*kappa+8*mu)/(6*kappa+12*mu); %3D
    etakappa = @(mu, kappa) 4*mu/3; %3D
end

A = @(f0, f1, df0, df1) -(-f0^2+2*f0*f1-f1^2+df0*df1) /(df1-f1+f0);
B = @(f0, f1, df0, df1) -(df0*f1-df0*df1+df1*f0-2*f0*f1+2*f0^2)/(df1-f1+f0);
C = @(f0)               f0;
D = @(f0, f1, df0, df1) -(df0+df1+2*f0-2*f1)/(df1-f1+f0);

df0 = @(eta, f0, f1) (f0+eta)*(f1-f0)/(f1+eta);
df1 = @(eta, f0, f1) (f1+eta)*(f1-f0)/(f0+eta);

%% Case mu

F0 = mu0;
F1 = mu1;
Eta0 = etamu(mu0, kappa0);
Eta1 = etamu(mu1, kappa1);
DF0 = df0(Eta0, F0, F1);
DF1 = df1(Eta1, F0, F1);

data = computeData(F0, F1, DF0, DF1, Eta0, Eta1, A, B, C, D);
computePlot('Shear modulus', dimension, data);

%% Case kappa
F0 = kappa0;
F1 = kappa1;
Eta0 = etakappa(mu0, kappa);
Eta1 = etakappa(mu1, kappa);
DF0 = df0(Eta0, F0, F1);
DF1 = df1(Eta1, F0, F1);
data = computeData(F0, F1, DF0, DF1, Eta0, Eta1, A, B, C, D);
computePlot('Bulk modulus', dimension, data);

end

function data = computeData(F0, F1, DF0, DF1, Eta0, Eta1, A, B, C, D)
data.AA = A(F0, F1, DF0, DF1);
data.BB = B(F0, F1, DF0, DF1);
data.CC = C(F0);
data.DD = D(F0, F1, DF0, DF1);
data.F0 = F0;
data.F1 = F1;
data.Eta0 = Eta0;
data.Eta1 = Eta1;
end

function computePlot(Title, dimension, data)
AA = data.AA;
BB = data.BB;
CC = data.CC;
DD = data.DD;
F0 = data.F0;
F1 = data.F1;
Eta0 = data.Eta0;
Eta1 = data.Eta1;
figure;
hold on
f = @(rho, A, B, C, D) (A.*rho.^2+B.*rho+C)./(D.*rho+1); %SIMPALL
%df = @(rho, A, B, C, D) (A*D.*rho.^2+2*A.*rho+B-C*D)./(D.*rho+1).^2;
plot(linspace(0,1,100),f(linspace(0,1,100), AA, BB, CC, DD), 'r', 'LineWidth',2); %SIMP-ALL plot
SIMP = @(rho) (1-rho.^3).*F0+rho.^3.*F1; %SIMP
plot(linspace(0,1,100),SIMP(linspace(0,1,100)), 'b', 'LineWidth',2); %SIMP plot
if dimension == 2
    fHsLb = @(rho) F0.*(1-rho)+F1.*rho-((1-rho).*rho*(F1-F0)^2)./(F0.*rho+F1.*(1-rho)+Eta0); %LB
    fHsUb = @(rho) F0.*(1-rho)+F1.*rho-((1-rho).*rho*(F1-F0)^2)./(F0.*rho+F1.*(1-rho)+Eta1); %UB
    plot(linspace(0,1,100),fHsLb(linspace(0,1,100)), 'k', 'LineWidth',2); %LB plot
    plot(linspace(0,1,100),fHsUb(linspace(0,1,100)), 'k', 'LineWidth',2); %UB plot
elseif dimension == 3
    fHsLb = @(rho) F0+(rho-1)./(1/(F0-F1)-rho./(F0+Eta0));
    fHsUb = @(rho) F1+rho./(1/(F0-F1)+(rho-1)./(F1+Eta1));
    plot(linspace(0,1,100),fHsLb(linspace(0,1,100)), 'k', 'LineWidth',2); %LB plot
    plot(linspace(0,1,100),fHsUb(linspace(0,1,100)), 'k', 'LineWidth',2); %UB plot
end
set(gca,'fontsize',16)
set(gca,'XTick',0:0.5:1)
title(Title,'fontsize',20)
xlabel('Density','fontsize',20)
legend({'SIMP-ALL','SIMP', 'HS-LB','HS-UB'},'location','best','fontsize',20)
end
