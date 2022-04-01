
%% Derivation of HS generalized equation
function DerivationHS()


f0 = sym('f0','real');
f1 = sym('f1','real');
eta = sym('eta','real');
rho = sym('rho','real');

fHS = f0*(1-rho)+f1*rho-(rho*(1-rho)*(f0-f1)^2)/(eta+f0*rho+f1*(1-rho));
dfHS = diff(fHS, rho);
dfHS = simplify(dfHS);
pretty(dfHS);
latex(dfHS)

dfHS0 = subs(dfHS, rho, 0);
dfHS0 = simplify(dfHS0);
pretty(dfHS0);
latex(dfHS0)

dfHS1 = subs(dfHS, rho, 1);
dfHS1 = simplify(dfHS1);
pretty(dfHS1);
latex(dfHS1)
end