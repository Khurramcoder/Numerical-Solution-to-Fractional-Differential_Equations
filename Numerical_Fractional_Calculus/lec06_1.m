%Baleanu, D., Jajarmi, A., & Hajipour, M. (2018).
%On the nonlinear dynamical systems within the generalized fractional derivatives with Mittag-Leffler kernel.
%Nonlinear dynamics, 94(1), 397-414.
% The numerical coded in this m file is given in eqs.(33)-(34) in the above
% research paper

clc; clear; close all;
% Initial Conditions
t(1)=0; y(1)=0;
% Inputs
alpha=3/4;ABC=1;h= 0.1; tfinal=1; t=t(1):h:tfinal;
% Number of Iterations
N=ceil((tfinal-t(1))/h);
% The given fractional-orderODE under the ABC operator
f = @(t,y) t.^2;
Exact = (1-alpha)/ABC*t.^2+(2/(gamma(alpha)*ABC*(alpha^2+3*alpha+2)))*t.^(alpha+2);
% ABC Fractional Euler Method starts
tic;
for n=1: N
    j=1:n;
y(n+1)=y(1)+((1-alpha)*f(t(n),y(n)))/(ABC)+...
((alpha.*h^alpha)./(ABC.*gamma (alpha+1))).*...
sum(((n-j+1).^alpha-(n-j).^alpha).*f(t(j),y(j)));
t(n+1)=t(n)+h;
end
toc;
Errors=abs(Exact-y);
Max_Error = max(Errors);Last_Error = Errors(end);
Ave_Error=mean(Errors);Norm_Error=norm(Errors);
ERRORS=[Max_Error Last_Error;Ave_Error Norm_Error]