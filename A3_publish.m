clear
close all
clc

%% Kalman filter creation
load('dataset0.mat');

[xUpd, xPred] = KalmanFilterFunc(0.1,1,1,1,data);

figure(1)
plot(dtime,data);
hold on
plot(dtime,xUpd);
plot(dtime,xPred);
xlabel("Time");
ylabel("Data");
legend("Market data", "Filtered estimates", "Prediction");
title("Kalman filter dataset0, not optimised");
grid on

%% Parameter estimation dataset0
clear

load("dataset0.mat");

options = optimset('MaxFunEvals',2000);

theta0 = [0.1 5 5 0.1];
x0 = 0;
p0 = 0;

fun = @(theta)max_like1(theta, x0, p0);

optParameters0 = fminsearch(fun,theta0,options)

optArguments0 = num2cell(optParameters0);

[xUpdOpt0, xPredOpt0] = KalmanFilterFunc(optArguments0{:},data);

figure(2)
plot(dtime,data);
hold on
plot(dtime,xUpdOpt0);
plot(dtime,xPredOpt0);
xlabel("Time");
ylabel("Data");
legend("Market data", "Filtered estimates", "Prediction");
title("Kalman filter dataset0, optimised");
grid on

%% Residual dataset0
theta = optParameters0;
load('dataset0.mat');

x0 = 0;
p0 = 0;

x1(1) = theta(1)*x0;
p1(1) = theta(1)*p0*theta(1)'+theta(3);

sk(1) = theta(2)*p1(1)*theta(2)'+theta(4);
L(1) =-0.5*log(2*pi)-0.5*log(abs(sk(1)))-0.5*(data(1)-theta(2)*x1(1))'*sk(1)^-1*(data(1)-theta(2)*x1(1));
res(1) = (data(1)-theta(2)*x1(1))/sk(1);

for i=1:1:length(data)-1
    K(i) = p1(i)*theta(2)'*((theta(2)*p1(i)*theta(2)'+theta(4))^(-1));
    
    xm(i) = x1(i)+K(i)*(data(i)-theta(2)*x1(i));
    pm(i) = (1-K(i)*theta(2))*p1(i);
    
    x1(i+1) = theta(1)*xm(i);
    p1(i+1) = theta(1)*pm(i)*theta(1)'+theta(3);
    sk(i+1) = theta(2)*p1(i+1)*theta(2)'+theta(4);
    res(i+1) = (data(i+1)-theta(2)*x1(i+1))/sk(i+1);
end

figure(3)
x = -2.5:0.1:2.5;
histu = hist(res,x);
bar(x,histu);
hold on
norm = 20*normpdf(x,0,1);
plot(x,norm,'LineWidth',2);
xlabel('Residual');
ylabel('Frequency');
legend("Residual", "20*N(0,1)");
title("Residual for dataset0");

%% Parameter estimation dataset1
clear

load("dataset1.mat");
t = datenum(dtime);
T0 = t(1);
t = t-T0;

options = optimset('MaxFunEvals',2000);

theta0 = [2 1 0.2 0.2];
x0 = 0;
p0 = 0;

fun = @(theta)max_like2(theta, x0, p0);

optParameters1 = fminsearch(fun,theta0,options)

optArguments1 = num2cell(optParameters1);

[xUpdOpt1, xPredOpt1] = KalmanFilterFunc(optArguments1{:},data);

figure(4)
plot(t,data);
hold on
plot(t,xUpdOpt1);
plot(t,xPredOpt1);
xlabel("Time");
ylabel("Data");
legend("Market data", "Filtered estimates", "Prediction");
title("Kalman filter dataset1, optimised");
grid on

%% Residual dataset1
theta = optParameters1;
load('dataset1.mat');

t = datenum(dtime);
T0 = t(1);
t = t-T0;

x0 = 0;
p0 = 0;

x1(1) = theta(1)*x0;
p1(1) = theta(1)*p0*theta(1)'+theta(3);

sk(1) = theta(2)*p1(1)*theta(2)'+theta(4);
L(1) =-0.5*log(2*pi)-0.5*log(abs(sk(1)))-0.5*(data(1)-theta(2)*x1(1))'*sk(1)^-1*(data(1)-theta(2)*x1(1));
res(1) = (data(1)-theta(2)*x1(1))/sk(1);

for i=1:1:length(data)-1
    K(i) = p1(i)*theta(2)'*((theta(2)*p1(i)*theta(2)'+theta(4))^(-1));
    
    xm(i) = x1(i)+K(i)*(data(i)-theta(2)*x1(i));
    pm(i) = (1-K(i)*theta(2))*p1(i);
    
    x1(i+1) = theta(1)*xm(i);
    p1(i+1) = theta(1)*pm(i)*theta(1)'+theta(3);
    sk(i+1) = theta(2)*p1(i+1)*theta(2)'+theta(4);
    res(i+1) = (data(i+1)-theta(2)*x1(i+1))/sk(i+1);
end

figure(5)
x = -2.5:0.1:2.5;
histu = hist(res,x);
bar(x,histu);
hold on
norm = 20*normpdf(x,0,1);
plot(x,norm,'LineWidth',2);
xlabel('Residual');
ylabel('Frequency');
legend("Residual", "20*N(0,1)");
title("Residual dataset1");
