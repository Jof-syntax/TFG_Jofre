function EtaComputer2D();
clc;

mu = sym('mu','real');
muPreDef = sym('muPreDef','real'); 

kappa = sym('k','real');
kPreDef = sym('kPreDef','real');

v = @(ki, mui) (ki-mui)/(ki+mui);

E = @(ki, mui) (4*ki*mui)/(ki+mui);

alpha = (1+v(kappa, mu))/(1-v(kappa, mu));
beta = (3-v(kappa, mu))/(1+v(kappa, mu));
gamma = E(kPreDef, muPreDef)/E(kappa, mu);
tau1 = (1+v(kPreDef, muPreDef))/(1+v(kappa, mu));
tau2 = (1-v(kPreDef, muPreDef))/(1-v(kappa, mu));
tau3 = (v(kPreDef, muPreDef)*(3*v(kappa, mu)-4)+1)/(v(kappa, mu)*(3*v(kappa, mu)-4)+1);

P1n = (1+beta)*(tau1-gamma);
P1d = beta*gamma+tau1;
P1 = P1n/P1d;
P1 = simplify(P1);

P2n = (alpha-beta)*(gamma*(gamma-2*tau3)+tau1*tau2);
P2d = 2*(beta*gamma+tau1)*(alpha*gamma+tau2);
P2 = P2n/P2d;
P2 = simplify(P2);

Pmu = P1;
Pmu = simplify(Pmu);
Pk = 2*P2+P1;
Pk = simplify(Pk);

qmu = Pmu/(mu-muPreDef);
qmu = simplify(qmu);

qk = Pk/(kappa-kPreDef);
qk = simplify(qk);

% Obtain Etas
etaM = (qmu*mu*muPreDef-mu)/(1-qmu*mu);
etaM = simplify(etaM);
etaK = (qk*kappa*kPreDef-kappa)/(1-qk*kappa);
etaK = simplify(etaK);

disp('etaM = ');
pretty(etaM);
disp('etaK = ');
pretty(etaK);

% Obtain dk(kappa, mu) & dmu(kappa, mu)

dmu = qmu*mu*(muPreDef-mu);
dmu = simplify(dmu);
dk = qk*kappa*(kPreDef-kappa);
dk = simplify(dk);
disp('dmu(kappa, mu) = ');
pretty(dmu);

disp('dk(kappa, mu) = ');
pretty(dk);

% Obtain dk(etaK, etaMu) & dmu(etaK, etaMu)

etaKappa = sym('etaKappa','real');
etaMu = sym('etaMu','real');

muAsEtaK = solve(etaK == etaKappa, mu);
muAsEta = simplify(muAsEtaK);

kappaAsEtaM = solve(etaM == etaMu, kappa, 'ReturnConditions',true); %Conditions
kappaAsEtaM = subs(kappaAsEtaM, mu, muAsEtaK);
kappaAsEta = simplify(kappaAsEtaM.k);



dmu = subs(dmu, kappa, kappaAsEta);
dmu = subs(dmu, mu, muAsEta);
dmu = simplify(dmu);

dk = subs(dk, kappa, kappaAsEta);
dk = subs(dk, mu, muAsEta);
dk = simplify(dk);
dk = simplify(dk);



