function SIMPALLcoefficients();
clc;

%Contraints
f0 = sym('f0','real');
f1 = sym('f1','real');
df0 = sym('df0','real');
df1 = sym('df1','real');
%Coefficients
A = sym('A','real');
B = sym('B','real');
C = sym('C','real');
D = sym('D','real');

 f = @(rho) (A.*rho.^2+B.*rho+C)./(D.*rho+1); %SIMPALL function
 df = @(rho) (A*D.*rho.^2+2*A.*rho+B-C*D)./(D.*rho+1).^2; % derivation checked
 
 [A, B, C, D] = solve(f(0)==f0, df(0)==df0, f(1)==f1, df(1)==df1, A, B, C, D);
 A = simplify(A);
 B = simplify(B);
 C = simplify(C);
 D = simplify(D);
 
  disp('A = ');
pretty(A);
 disp('B = ');
pretty(B);
 disp('C = ');
pretty(C);
 disp('D = ');
pretty(D);
 
 