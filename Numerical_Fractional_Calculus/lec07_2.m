%Jajarmi, A., & Baleanu, D.(2018). A new fractional analysis on the interaction of HIV
%with CD4+ T-cells, Chaos, Solitons & Fractals, 113,221-229.
%The numerical method coded in this m file is given in Eq(23)in the above
%research paper

clc;clear;close all;
%Initial Conditions
t(1)=0; y(1)=0;
%Inputs
alpha=0.874;ABC=1;h=0.1; tfinal=1; t=t(1):h:tfinal;
%Number of Iterations
N=ceil((tfinal-t(1))/h);
%The given fractional-order ODE under the CF operator
f=@(t,y)t.^2;
Exact=(1-alpha)*t.^2+(alpha/3)*t.^3;
%%Caputo-Fabrizio Euler Method Starts
tic;
for n=1:N
    j=1:n;
    y(n+1)=y(1)+((1-alpha)*f(t(n),y(n)))+(alpha.*h).*sum(f(t(j),y(j))); 
    t(n+1)=t(n)+h;
end
toc;
Errors=abs(Exact-y);
Max_Error=max(Errors);Last_Error=Errors(end);
Ave_Error=mean(Errors);Norm_Error=norm(Errors);
ERRORS=[Max_Error Last_Error;Ave_Error Norm_Error]