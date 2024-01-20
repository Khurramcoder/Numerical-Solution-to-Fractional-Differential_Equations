% Numerical L2C Method to approximate the left Caputo Derivative


% Example f(x)=sin(x) over [0,1] & alpha=1.5

clc;clear;close all;


%Inputs

h=0.1;
x(1)=0;
xlast=1;
x=x(1):h:xlast;
n=ceil((xlast-x(1))/h);


alpha =1.5; 

f=@(x)sin(x);


% The Algorithm
A=(h.^(-alpha))./(2*gamma(3-alpha));


k=3:n;




%Bnk=((n-k).^(1-alpha)-(n-k-1).^(1-alpha));

L2C_left_Caputo_Der=A.*sum(((n-k+1).^(2-alpha)-(n-k).^(2-alpha)).*(f(x(k+1))-f(x(k))-f(x(k-1))+f(x(k-2))));

% Exact Solution

x=xlast;

Exact=(cos(x)*fresnelc(sqrt(x)*sqrt(2)/sqrt(pi))-sin(x)*fresnels(sqrt(x)*sqrt(2)/sqrt(pi)))*sqrt(2);

Error=abs(Exact-L2C_left_Caputo_Der);

disp('   Steps    Stepsize    Exact    Approximate    Error');


disp('-----------------------------------------------------');


Result=[n   h   Exact    L2C_left_Caputo_Der      Error]