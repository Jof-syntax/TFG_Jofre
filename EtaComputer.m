function EtaComputer();
clc;

mu = sym('mu','real');
muPreDef = sym('muPreDef','real'); 
mu2 = mu - muPreDef;

lambda = sym('lam','real');
lamPreDef = sym('lamPreDef','real'); 
lam2 = lambda - lamPreDef;

v = sym('v','real');

kappa = sym('k','real');
kPreDef = sym('kPreDef','real');

%Elas3D definitions
m1n = 15*mu*mu2*(v-1);
m1d = 15*mu*(1-v)+2*mu2*(5*v-4);
m1 = m1n/m1d;

m2n = lam2*(15*mu*lambda*(1-v) + 2*lambda*mu2*(5*v-4)) - 2*mu2*(lambda*mu2-5*mu*v*lam2);
m2d = 5*mu2*(3*mu*lambda*(1-v)-3*mu*v*lam2-lambda*mu2*(1-2*v));
m2 = m2n/m2d;

% Change of variables

m1 = subs(m1, v, lambda/(2*(lambda+mu)));
m1 = subs(m1, lambda, kappa-2/3*mu);
m1 = subs(m1, lamPreDef, kPreDef-2/3*muPreDef);
m1 = simplify(m1);

m2 = subs(m2, v, lambda/(2*(lambda+mu)));
m2 = subs(m2, lambda, kappa-2/3*mu);
m2 = subs(m2, lamPreDef, kPreDef-2/3*muPreDef);
m2 = simplify(m2);

% New definitions

dmu = m1;
dmu = simplify(dmu);
dk = m1*m2+m1*2/3;
dk = simplify(dk);

qmu = dmu/(mu*(muPreDef-mu));
qmu = simplify(qmu);

qk = dk/(kappa*(kPreDef-kappa));
qk = simplify(qk);


etaM = (qmu*mu*muPreDef-mu)/(1-qmu*mu);
etaM = simplify(etaM);
pretty(etaM);
etaK = (qk*kappa*kPreDef-kappa)/(1-qk*kappa);
etaK = simplify(etaK);
pretty(etaK);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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




