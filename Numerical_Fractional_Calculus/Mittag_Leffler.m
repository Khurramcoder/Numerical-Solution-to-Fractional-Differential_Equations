% Mittag Leffler Function
clc; clear; close all;
syms x k
alpha=0.95;

beta=0.87;

gma=0.68;

mittag_series=(1./gamma(gma)).*symsum(gamma(gma+k).*...
    (x.^k)./(factorial(k).*gamma(alpha*k+1)),k,0,Inf);


P=11;a=0;b=1;

X=linspace(a,b,P);

for i=1:P
    mittag(i,1)=double(subs(mittag_series,x,X(i)));
end

xvalues=X;
mittagvalues=mittag';
A=[xvalues;mittagvalues]'
plot(X,mittag,'linewidth',2)

