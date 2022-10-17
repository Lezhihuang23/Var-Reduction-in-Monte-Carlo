% Vanilla Monte Carlo Simulation 
clear all 
format long 
r=0.02; 
S0=100; 
sig=0.2; 
K1=60; 
K2=50; 
T=1; 
N=1000000; 
E=randn(1,N); 
ST=S0*exp((r-0.5*sig^2)*T + sig*sqrt(T).*E); 
optionprice= exp(-r*T)*(max(K1-ST,0)-max(K2-ST,0)); 
MC=cumsum(optionprice)./[1:N]; 
se=sqrt( sum ((optionprice-MC(N)).^2) / ((N-1)*N)); 
fprintf('1c) No. of Path= %.0f\n',N) 
fprintf('MC estimate= %f\n',MC(N)) 
fprintf('Standard Error %f',se) 
plot(MC(1:50000)) 
title('1c) Convergence Diagram') 
xlabel('No. of Path') 
ylabel('Option price') 

% Monte Carlo Simulation 
% antithetic variable technique 
clear all 
format long 
r=0.02; 
S0=100; 
sig=0.2; 
K1=60; 
K2=50; 
T=1; 
N=1000000; 
E1=randn(1,N/2); 
E=[E1,-E1]; 
ST=S0*exp((r-0.5*sig^2)*T + sig*sqrt(T).*E); 
optionprice= exp(-r*T)*(max(K1-ST,0)-max(K2-ST,0)); 
MC=cumsum(optionprice)./[1:N]; 
se=sqrt( sum ((optionprice-MC(N)).^2) / ((N-1)*N)); 
fprintf('1d) No. of Path= %.0f\n',N) 
fprintf('MC estimate= %f\n',MC(N)) 
fprintf('Standard Error %f',se) 
plot(MC(1:50000)) 
title('1d) Convergence Diagram') 
xlabel('No. of Path') 
ylabel('Option price')

% Monte Carlo Simulation 
%Moment Matching 
clear all 
format long 
r=0.02; 
S0=100; 
sig=0.2; 
K1=60; 
K2=50; 
T=1; 
N=1000000; 
E1=randn(1,N); 
mu=sum(E1)/N; 
std=sqrt(sum((E1-mu).^2)/(N-1)); 
E=(E1-mu)/std; 
ST=S0*exp((r-0.5*sig^2)*T + sig*sqrt(T).*E); 
optionprice= exp(-r*T)*(max(K1-ST,0)-max(K2-ST,0)); 
MC=cumsum(optionprice)./[1:N]; 
se=sqrt( sum ((optionprice-MC(N)).^2) / ((N-1)*N)); 
fprintf('1e) No. of Path= %.0f\n',N) 
fprintf('MC estimate= %f\n',MC(N)) 
fprintf('Standard Error %.10f',se) 
plot(MC(1:50000)) 
title('1e) Convergence Diagram') 
xlabel('No. of Path') 
ylabel('Option price') 

% Vanilla Monte Carlo Simulation 
%important sampling technique 
clear all 
format long 
r=0.02; 
S0=100; 
sig=0.2; 
K1=60; 
K2=50; 
T=1; 
alpha=0.1-5*log(55/100); 
N=1000000; 
n=randn(1,N); 
ST=S0*exp((r-0.5*sig^2)*T + sig*sqrt(T).*(n-alpha)); 
optionprice= exp(-r*T)*(max(K1-ST,0)-max(K2-ST,0)).*exp(alpha*n-0.5*alpha^2); 
MC=cumsum(optionprice)./[1:N]; 
se=sqrt( sum ((optionprice-MC(N)).^2) / ((N-1)*N)); 
fprintf('1f) No. of Path= %.0f\n',N) 
fprintf('MC estimate= %f\n',MC(N))
fprintf('Standard Error %.10f',se) 
plot(MC(1:50000)) 
title('1f) Convergence Diagram') 
xlabel('No. of Path') 
ylabel('Option price') 

% Monte Carlo Simulation 
% Important Sampling with antithetic variable technique 
clear all 
format long 
r=0.02; 
S0=100; 
sig=0.2; 
K1=60; 
K2=50; 
T=1; 
alpha=0.1-5*log(55/100); 
N=1000000; 
E1=randn(1,N/2); 
n=[E1,-E1]; 
ST=S0*exp((r-0.5*sig^2)*T + sig*sqrt(T).*(n-alpha)); 
optionprice= exp(-r*T)*(max(K1-ST,0)-max(K2-ST,0)).*exp(alpha*n-0.5*alpha^2); 
MC=cumsum(optionprice)./[1:N]; 
se=sqrt( sum ((optionprice-MC(N)).^2) / ((N-1)*N)); 
fprintf('1g) No. of Path= %.0f\n',N) 
fprintf('MC estimate= %f\n',MC(N)) 
fprintf('Standard Error %.10f',se) 
plot(MC(1:50000)) 
title('1g) Convergence Diagram') 
xlabel('No. of Path') 
ylabel('Option price') 

% Monte Carlo Simulation 
% Important Sampling with moment matching technique 
clear all 
format long 
r=0.02; 
S0=100; 
sig=0.2; 
K1=60; 
K2=50; 
T=1; 
alpha=0.1-5*log(55/100); 
N=1000000; 
E1=randn(1,N); 
mu=sum(E1)/N; 
std=sqrt(sum((E1-mu).^2)/(N-1)); 
n=(E1-mu)/std; 
ST=S0*exp((r-0.5*sig^2)*T + sig*sqrt(T).*(n-alpha)); 
optionprice= exp(-r*T)*(max(K1-ST,0)-max(K2-ST,0)).*exp(alpha*n-0.5*alpha^2); 
MC=cumsum(optionprice)./[1:N]; 
se=sqrt( sum ((optionprice-MC(N)).^2) / ((N-1)*N)); 
fprintf('1h) No. of Path= %.0f\n',N) 
fprintf('MC estimate= %f\n',MC(N)) 
fprintf('Standard Error %.10f',se) 
plot(MC(1:50000)) 
title('1h) Convergence Diagram') 
xlabel('No. of Path') 
ylabel('Option price')