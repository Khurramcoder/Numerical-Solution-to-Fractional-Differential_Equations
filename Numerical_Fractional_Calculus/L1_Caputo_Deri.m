% Numerical L1 Method for Caputo Derivative


% Example f(x)=sin(x) over [0,1]

clc;clear;close all;


%Inputs

h=0.1;
x(1)=0;
xlast=1;
x=x(1):h:xlast;
n=ceil((xlast-x(1))/h);


alpha =0.5; 

f=@(x)sin(x);


% The Algorithm

k=1:n-1;

A=(h.^(-alpha))./(gamma(2-alpha));


Bnk=((n-k).^(1-alpha)-(n-k-1).^(1-alpha));

L1_Caputo_Der=A*sum(Bnk.*(f(x(k+1))-f(x(k))));

% Exact Solution

x=xlast;

Exact=(cos(x)*fresnelc(sqrt(x)*sqrt(2)/sqrt(pi))+sin(x)*fresnels(sqrt(x)*sqrt(2)/sqrt(pi)))*sqrt(2);

Error=abs(Exact-L1_Caputo_Der);

disp('   Steps    Stepsize    Exact    Approximate    Error');


disp('-----------------------------------------------------');


Result=[n   h   Exact    L1_Caputo_Der      Error]
