% Akman,T.,Yildiz,B.,& Baleanu,D.(2018).
% Numerical discretization of Caputo-Fabrizio Derivative.
% Computational & Applied Mathematics,37(3),3307-3333.


% Numerical L1-2 Method to approximate the Caputo Fabrizio-Central Derivative

%Example f(x)=cos(4x) over [0,2].

clc;clear;close all;


%Inputs

h=0.1;
x(1)=0;
xlast=2;
x=x(1):h:xlast;
k=ceil((xlast-x(1))/h);


alpha =0.1; 

M=1;w=4;
lam=-alpha/(1-alpha);
f=@(x)cos(w*x);
% The Algorithm
j=1:k;


B=(2/(h))*(f((j-1)*h).*(1/2-2.*j))+(4/h).*(f((j-1/2)*h).*(2*j-1))+(2./(h)).*(f(j*h).*(3/2-2*j));
R1=exp(lam*(k*h-j*h));
R2=exp(lam*(k*h-(j-1)*h));
G=(1/h^2)*(4*f((j-1)*h)-8*f((j-1/2)*h)+4*f(j*h));
P=j.*h.*R1-(j-1).*h.*R2-((1-alpha)/alpha)*R1+((1-alpha)/alpha)*R2;

L12CF=(1/alpha).*sum(B.*(R1-R2)+G.*P);




%Preparation for the Exact Solution

syms x

Exact=(M/(1-alpha))*int(diff(f,x)*exp(lam*(xlast-x)),x,0,2);

Exact=double(Exact);

Error = abs(Exact-L12CF);

disp('   Steps    Stepsize    Exact    Approximate    Error');


disp('-----------------------------------------------------');


Result = [k   h   Exact    L12CF      Error]