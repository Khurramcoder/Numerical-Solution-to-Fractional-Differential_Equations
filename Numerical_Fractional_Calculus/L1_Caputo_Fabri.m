% Akman,T.,Yildiz,B.,& Baleanu,D.(2018).
% Numerical discretization of Caputo-Fabrizio Derivative.
% Computational & Applied Mathematics,37(3),3307-3333.


% Numerical L1 Method to approximate the Caputo Fabrizio Derivative

%Example f(x)=cos(4x) over [0,2].

clc;clear;close all;


%Inputs

h=0.1;
x(1)=0;
xlast=2;
x=x(1):h:xlast;
n=ceil((xlast-x(1))/h);


alpha =0.1; 

M=1;w=4;

f=@(x)cos(w*x);


% The Algorithm
A=M/(h*alpha);

lam=-alpha/(1-alpha);

k=1:n;
L1CF=A.*sum((f(k.*h)-f(h.*(k-1))).*(exp(lam.*h.*(n-k))-exp(lam.*h.*(n-k+1))));



% Exact Solution
syms x


Exact=(M/(1-alpha))*int(diff(f,x)*exp(lam*(xlast-x)),x,0,2);

Exact=double(Exact);
Error=abs(Exact-L1CF);
disp('   Steps    Stepsize    Exact    Approximate    Error');


disp('-----------------------------------------------------');


Result = [n   h   Exact    L1CF      Error]