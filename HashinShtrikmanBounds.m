function HashinShtrikmanBounds();
clc;

N = 3;
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


muUB = solve(rho/(2*(mu1-muUB))== 1/(2*(mu1-mu0)) + (1-rho)*(N-1)*(k1+2*mu1)/((N^2+N-2)*mu1*(2*mu1+lambda1))   , muUB);
muUB = simplify(muUB);
disp('muUB = ');
pretty(muUB);

muLB = solve((1-rho)/(2*(muLB-mu0))== 1/(2*(mu1-mu0)) + rho*(N-1)*(k0+2*mu0)/((N^2+N-2)*mu0*(2*mu0+lambda0))   , muLB);
muLB = simplify(muLB);
disp('muLB = ');
pretty(muLB);


kappaUB = solve(  rho/(k1-kappaUB) == 1/(k1-k0)+(1-rho)/(2*mu1+lambda1) , kappaUB);
kappaUB = simplify(kappaUB);
disp('kUB = ');
pretty(kappaUB);

kappaLB = solve( (1-rho)/(kappaLB-k0) == 1/(k1-k0)+rho/(2*mu0+lambda0)  , kappaLB);
kappaLB = simplify(kappaLB);
disp('kLB = ');
pretty(kappaLB);

