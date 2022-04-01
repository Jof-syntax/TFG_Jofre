%% Obtain critic cases SIMPALL 3D

function DataStoraged = Critic3dSIMPALL()
clear all;
clc;
    DataStoraged = zeros(1.5*100,4);
    NumberOfCriticCases = 1;
    E1 = 1;
    E0 = 0.01;
    % nu1 determined by means of the 'for'
    nu0 =  1/3;
    N = 3;
    
    for i = -1:0.01:0.5
        a = Comparator(E1,E0,i,nu0, N);
        if (a.checker == true)
            DataStoraged(NumberOfCriticCases, 1) = E1;
            DataStoraged(NumberOfCriticCases, 2) = E0;
            DataStoraged(NumberOfCriticCases, 3) = i; %nu1
            DataStoraged(NumberOfCriticCases, 4) = nu0;
            NumberOfCriticCases = NumberOfCriticCases+1;
        end
    end
    disp('The number of critic cases when the SIMP is outside the HS bounds are: ');
    disp(NumberOfCriticCases);
end

function Data = Comparator(E1,E0,nu1,nu0, N)
Data.checker = false; %true if constraint is found

mu=@(E,nu) E/(2*(1+nu));
kappa = @(E,nu, N) (N*nu*E+E-2*E*nu)/(N-nu*N);

mu0 = mu(E0,nu0);
mu1 = mu(E1,nu1);
kappa0 = kappa(E0, nu0, N);
kappa1 = kappa(E1, nu1, N);

etaMu = @(mu, kappa, N)  -(mu*(4*mu - kappa*N^2 - 2*mu*N^2 + 2*mu*N))/(2*N*(kappa + 2*mu));
etaKappa = @(mu, kappa, N) (2*mu*(N - 1))/N;

%% Case mu
muData = zeros(3,3);
muDataCounter = 1;
F0 = mu0;
F1 = mu1;
Eta0 = etaMu(mu0, kappa0, N);
Eta1 = etaMu(mu1, kappa1, N);
Simp = @(rho) (1-rho.^3).*F0+rho.^3.*F1; %SIMP
fHs = @(rho, Eta) F0.*(1-rho)+F1.*rho-((1-rho).*rho*(F1-F0)^2)./(F0.*rho+F1.*(1-rho)+Eta);
for i = 0.25:0.25:0.75  %only checks 3 points
    muData(1, muDataCounter) = fHs(i, Eta1);    % HS UB
    muData(2, muDataCounter) = fHs(i, Eta0);	% HS LB
    muData(3, muDataCounter) = Simp(i);         % SIMP
    muDataCounter = muDataCounter + 1;
end

%% Case kappa
kappaData = zeros(3,3);
kappaDataCounter = 1;
F0 = kappa0;
F1 = kappa1;
Eta0 = etaKappa(mu0, kappa0, N);
Eta1 = etaKappa(mu1, kappa1, N);
Simp = @(rho) (1-rho.^3).*F0+rho.^3.*F1; %SIMP
fHs = @(rho, Eta) F0.*(1-rho)+F1.*rho-((1-rho).*rho*(F1-F0)^2)./(F0.*rho+F1.*(1-rho)+Eta);
    for i = 0.25:0.25:0.75  %only checks 3 points
        kappaData(1, kappaDataCounter) = fHs(i, Eta1);    % HS UB
        kappaData(2, kappaDataCounter) = fHs(i, Eta0);	% HS LB
        kappaData(3, kappaDataCounter) = Simp(i);         % SIMP
        kappaDataCounter = kappaDataCounter + 1;
    end
% Check constraint Simp above or lower the HS bounds
    if (muData(1,1) <= muData(3,1) || muData(2,1) >= muData(3,1) || kappaData(1,1) <= kappaData(3,1) || kappaData(2,1)>=kappaData(3,1))
        Data.checker    = true;
    elseif (muData(1,2) <= muData(3,2) || muData(2,2) >= muData(3,2) || kappaData(1,2) <= kappaData(3,2) || kappaData(2,2)>=kappaData(3,2))
        Data.checker    = true;
    elseif (muData(1,3) <= muData(3,3) || muData(2,3) >= muData(3,3) || kappaData(1,3) <= kappaData(3,3) || kappaData(2,3)>=kappaData(3,3))
        Data.checker    = true;
    end
end



