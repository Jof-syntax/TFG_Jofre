function HashinShtrikmanBounds();
clc;
% 3 Dimension case
N = 3;
porosity = sym('porosity','real');
rho = sym('rho','real');
muUB = sym('muUB','real');
muLB = sym('muLB','real');
kappaUB = sym('kUB','real');
kappaLB = sym('kLB','real');
mu1 = sym('mu1','real');
mu0 = sym('mu0','real');
k1 = sym('k1','real');
k0 = sym('k0','real');
lambda1 = k1-mu1*2/N;
lambda0 = k0-mu0*2/N;

muUBf = solve(porosity/(2*(mu1-muUB)) == 1/(2*(mu1-mu0)) + ((1-porosity)*(N-1)*(k1+2*mu1))/((N^2+N-2)*mu1*(2*mu1+lambda1))   , muUB);
muUBf = subs(muUBf, porosity, 1-rho);
muUBf = simplify(muUBf);
disp('muUB = ');
pretty(muUBf);

muLBf = solve((1-porosity)/(2*(muLB-mu0))== 1/(2*(mu1-mu0)) + porosity*(N-1)*(k0+2*mu0)/((N^2+N-2)*mu0*(2*mu0+lambda0))   , muLB);
muLBf = subs(muLBf, porosity, 1-rho);
muLBf = simplify(muLBf);
disp('muLB = ');
pretty(muLBf);

kappaUBf = solve( porosity/(k1-kappaUB) == 1/(k1-k0)-(1-porosity)/(2*mu1+lambda1) , kappaUB);
kappaUBf = subs(kappaUBf, porosity, 1-rho);
kappaUBf = simplify(kappaUBf);
disp('kUB = ');
pretty(kappaUBf);

kappaLBf = solve( (1-porosity)/(kappaLB-k0) == 1/(k1-k0)+porosity/(2*mu0+lambda0)  , kappaLB);
kappaLBf = subs(kappaLBf, porosity, 1-rho);
kappaLBf = simplify(kappaLBf);
disp('kLB = ');
pretty(kappaLBf);

F1 = sym('F1','real');
F0 = sym('F0','real');
Eta1 = sym('Eta1','real');
Eta0 = sym('Eta0','real');

% term to add to F1*rho + F0*(1-rho) + ...
fsh1 = (-(rho).*(F0-F1)*(F0+Eta0) + (- F1*rho + F0*rho )* (F1*(1-rho)+F0*rho+Eta0))./(F1*(1-rho)+F0*rho+Eta0);
fsh1 = simplify(fsh1);
pretty(fsh1);
% term to add to F1*rho + F0*(1-rho) + ...
fsh2 = ((1-rho).*(F0-F1)*(F1+Eta1)+(F1-F1*rho-F0*(1-rho))*(F1*(1-rho)+F0*rho+Eta1))./(F1*(1-rho)+F0*rho+Eta1);
fsh2 = simplify(fsh2);
pretty(fsh2);

%% N dimension

N = sym('N','real');
porosity = sym('porosity','real');
rho = sym('rho','real');
muUB = sym('muUB','real');
muLB = sym('muLB','real');
kappaUB = sym('kUB','real');
kappaLB = sym('kLB','real');
F1k = sym('F1k','real');
F0k = sym('F0k','real');
F1mu = sym('F1mu','real');
F0mu = sym('F0mu','real');
lambda1 = F1k-F1mu*2/N;
lambda0 = F0k-F0mu*2/N;
Eta = sym('Eta','real');

muUBf = solve(porosity/(2*(F1mu-muUB)) == 1/(2*(F1mu-F0mu)) + ((1-porosity)*(N-1)*(F1k+2*F1mu))/((N^2+N-2)*F1mu*(2*F1mu+lambda1))   , muUB);
muUBf = subs(muUBf, porosity, 1-rho);
muUBf = simplify(muUBf);
disp('muUB = ');
pretty(muUBf);

muLBf = solve((1-porosity)/(2*(muLB-F0mu))== 1/(2*(F1mu-F0mu)) + porosity*(N-1)*(F0k+2*F0mu)/((N^2+N-2)*F0mu*(2*F0mu+lambda0))   , muLB);
muLBf = subs(muLBf, porosity, 1-rho);
muLBf = simplify(muLBf);
disp('muLB = ');
pretty(muLBf);

kappaUBf = solve( porosity/(F1k-kappaUB) == 1/(F1k-F0k)-(1-porosity)/(2*F1mu+lambda1) , kappaUB);
kappaUBf = subs(kappaUBf, porosity, 1-rho);
kappaUBf = simplify(kappaUBf);
disp('kUB = ');
pretty(kappaUBf);

kappaLBf = solve( (1-porosity)/(kappaLB-F0k) == 1/(F1k-F0k)+porosity/(2*F0mu+lambda0)  , kappaLB);
kappaLBf = subs(kappaLBf, porosity, 1-rho);
kappaLBf = simplify(kappaLBf);
disp('kLB = ');
pretty(kappaLBf);

ETAk = solve(kappaLBf==F0k - rho/(1/(F0k - F1k) + (rho - 1)/(F0k +Eta)), Eta);
ETAk = simplify(ETAk);
disp('ETA k =');
pretty(ETAk);

ETAmu = solve((2*N*(F0k + 2*F0mu)*(N - 1))/(F0mu*(N^2 + N - 2)*(F0k*N - 2*F0mu + 2*F0mu*N))==1/(F0mu+Eta), Eta);
ETAmu = simplify(ETAmu);
disp('ETA mu =');
pretty(ETAmu);

