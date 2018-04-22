%% Assignment 3 - Kalman filter for parameter estimation
%
% Written by Peili Guo (peili.guo.7645@student.uu.se) and Francisco José
% Peralta Alguacil (franciscojose.peraltaalguacil.0481@student.uu.se)
% This report is for Computational Finance: Calibration and Estimation
% Assignment 3. 
%
% In this project we apply kalman filter to pricing models and optimise it
% to find the corresponding set of parameters. 
%
% First, we build the Kalman filter and observe how it works with
% prediction and data update. and try optimization of parameters with 10
% interations of "guessing optimizer".
%
% In the second part, the maximum likelihood was computed and by calling fminsearch
% we find the set of parameter that maximize the likelihood of the
% filtering processs. And this methods was also applied to ABB stock data
% on 2015-02-05.
%
% In the last part, we plot the residual and visualize it with histogram to
% check if it is normal distributed. 

clear
close all
clc

%% Kalman filter creation
load('dataset0.mat');

[xUpd, xPred] = KalmanFilterFunc(0.1,1,1,1,data);

figure(1)
plot(dtime,data,'bo--');
hold on
plot(dtime,xUpd,'r*--');
plot(dtime,xPred,'kx-');
xlabel("Time");
ylabel("Data");
legend({"Market data", "Filtered estimates", "Prediction"},'FontSize',12);
title("Kalman filter dataset0, not optimised",'FontSize', 16);
grid on

%% 10 interations for guessing optimizer 

load('dataset0.mat');

%iteration 1 
ite = [1 1 1 1];
[xUpd, xPred] = KalmanFilterFunc(ite(1),ite(2),ite(3),ite(4),data);

f1 = figure('position', [0, 0, 600, 500]);
plot(dtime,data,'bo--',dtime,xUpd,'r*--',dtime,xPred,'kx-');
xlabel("Time");
ylabel("Data");
legend({"Market data", "Filtered estimates", "Prediction"},'FontSize',12);
str = sprintf('Fk = %s, Hk = %s, Qk = %s, Rk = %s',string(ite(1)),string(ite(2)),string(ite(3)),string(ite(4)));
title({"Kalman filter dataset0 with parameter 1";str},'FontSize', 16);
grid on

%iteration 2
ite = [0.3 1 0.5 1];
[xUpd, xPred] = KalmanFilterFunc(ite(1),ite(2),ite(3),ite(4),data);

f2 = figure('position', [0, 0, 600, 500]);
plot(dtime,data,'bo--',dtime,xUpd,'r*--',dtime,xPred,'kx-');
xlabel("Time");
ylabel("Data");
legend({"Market data", "Filtered estimates", "Prediction"},'FontSize',12);
str = sprintf('Fk = %s, Hk = %s, Qk = %s, Rk = %s',string(ite(1)),string(ite(2)),string(ite(3)),string(ite(4)));
title({"Kalman filter dataset0 with parameter 2";str},'FontSize', 16);
grid on

%iteration 3
ite = [0.3 1 0.5 2];
[xUpd, xPred] = KalmanFilterFunc(ite(1),ite(2),ite(3),ite(4),data);

f3 = figure('position', [0, 0, 600, 500]);
plot(dtime,data,'bo--',dtime,xUpd,'r*--',dtime,xPred,'kx-');
xlabel("Time");
ylabel("Data");
legend({"Market data", "Filtered estimates", "Prediction"},'FontSize',12);
str = sprintf('Fk = %s, Hk = %s, Qk = %s, Rk = %s',string(ite(1)),string(ite(2)),string(ite(3)),string(ite(4)));
title({"Kalman filter dataset0 with parameter 3";str},'FontSize', 16);
grid on

%iteration 4
ite = [0.3 5 0.5 2];
[xUpd, xPred] = KalmanFilterFunc(ite(1),ite(2),ite(3),ite(4),data);

f4 = figure('position', [0, 0, 600, 500]);
plot(dtime,data,'bo--',dtime,xUpd,'r*--',dtime,xPred,'kx-');
xlabel("Time");
ylabel("Data");
legend({"Market data", "Filtered estimates", "Prediction"},'FontSize',12);
str = sprintf('Fk = %s, Hk = %s, Qk = %s, Rk = %s',string(ite(1)),string(ite(2)),string(ite(3)),string(ite(4)));
title({"Kalman filter dataset0 with parameter 4";str},'FontSize', 16);
grid on

%iteration 5
ite = [2 5 2.5 2];
[xUpd, xPred] = KalmanFilterFunc(ite(1),ite(2),ite(3),ite(4),data);

f5 = figure('position', [0, 0, 600, 500]);
plot(dtime,data,'bo--',dtime,xUpd,'r*--',dtime,xPred,'kx-');
xlabel("Time");
ylabel("Data");
legend({"Market data", "Filtered estimates", "Prediction"},'FontSize',12);
str = sprintf('Fk = %s, Hk = %s, Qk = %s, Rk = %s',string(ite(1)),string(ite(2)),string(ite(3)),string(ite(4)));
title({"Kalman filter dataset0 with parameter 5";str},'FontSize', 16);
grid on


%iteration 6
ite = [2 5 2.5 0.25];
[xUpd, xPred] = KalmanFilterFunc(ite(1),ite(2),ite(3),ite(4),data);

f6 = figure('position', [0, 0, 600, 500]);
plot(dtime,data,'bo--',dtime,xUpd,'r*--',dtime,xPred,'kx-');
xlabel("Time");
ylabel("Data");
legend({"Market data", "Filtered estimates", "Prediction"},'FontSize',12);
str = sprintf('Fk = %s, Hk = %s, Qk = %s, Rk = %s',string(ite(1)),string(ite(2)),string(ite(3)),string(ite(4)));
title({"Kalman filter dataset0 with parameter 6";str},'FontSize', 16);
grid on


%iteration 7
ite = [2 1 0.5 0.25];
[xUpd, xPred] = KalmanFilterFunc(ite(1),ite(2),ite(3),ite(4),data);

f7 = figure('position', [0, 0, 600, 500]);
plot(dtime,data,'bo--',dtime,xUpd,'r*--',dtime,xPred,'kx-');
xlabel("Time");
ylabel("Data");
legend({"Market data", "Filtered estimates", "Prediction"},'FontSize',12);
str = sprintf('Fk = %s, Hk = %s, Qk = %s, Rk = %s',string(ite(1)),string(ite(2)),string(ite(3)),string(ite(4)));
title({"Kalman filter dataset0 with parameter 7";str},'FontSize', 16);
grid on

%iteration 8
ite = [1.2 1.3 0.5 0.25];
[xUpd, xPred] = KalmanFilterFunc(ite(1),ite(2),ite(3),ite(4),data);

f8 = figure('position', [0, 0, 600, 500]);
plot(dtime,data,'bo--',dtime,xUpd,'r*--',dtime,xPred,'kx-');
xlabel("Time");
ylabel("Data");
legend({"Market data", "Filtered estimates", "Prediction"},'FontSize',12);
str = sprintf('Fk = %s, Hk = %s, Qk = %s, Rk = %s',string(ite(1)),string(ite(2)),string(ite(3)),string(ite(4)));
title({"Kalman filter dataset0 with parameter 8";str},'FontSize', 16);
grid on


%iteration 9
ite = [0.65 1.15 0.5 0.25];
[xUpd, xPred] = KalmanFilterFunc(ite(1),ite(2),ite(3),ite(4),data);

f9 = figure('position', [0, 0, 600, 500]);
plot(dtime,data,'bo--',dtime,xUpd,'r*--',dtime,xPred,'kx-');
xlabel("Time");
ylabel("Data");
legend({"Market data", "Filtered estimates", "Prediction"},'FontSize',12);
str = sprintf('Fk = %s, Hk = %s, Qk = %s, Rk = %s',string(ite(1)),string(ite(2)),string(ite(3)),string(ite(4)));
title({"Kalman filter dataset0 with parameter 9";str},'FontSize', 16);
grid on



%iteration 10
ite = [0.65 1.15 0.65 0.05];
[xUpd, xPred] = KalmanFilterFunc(ite(1),ite(2),ite(3),ite(4),data);

f10 = figure('position', [0, 0, 600, 500]);
plot(dtime,data,'bo--',dtime,xUpd,'r*--',dtime,xPred,'kx-');
xlabel("Time");
ylabel("Data");
legend({"Market data", "Filtered estimates", "Prediction"},'FontSize',12);
str = sprintf('Fk = %s, Hk = %s, Qk = %s, Rk = %s',string(ite(1)),string(ite(2)),string(ite(3)),string(ite(4)));
title({"Kalman filter dataset0 with parameter 10";str},'FontSize', 16);
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
plot(dtime,data,'bo--');
hold on
plot(dtime,xUpdOpt0,'r*--');
plot(dtime,xPredOpt0,'kx-');
xlabel("Time");
ylabel("Data");
legend({"Market data", "Filtered estimates", "Prediction"},'FontSize',12);
title("Kalman filter dataset0, optimised",'FontSize',16);
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
legend({"Residual", "20*N(0,1)"},'FontSize',12);
title("Residual for dataset0",'FontSize',16);

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
plot(t,data,'bo--');
hold on
plot(t,xUpdOpt1,'r*--');
plot(t,xPredOpt1,'kx-');
xlabel("Time");
ylabel("Data");
legend({"Market data", "Filtered estimates", "Prediction"},'FontSize',12);
title("Kalman filter dataset1, optimised",'FontSize',16);
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
legend({"Residual", "20*N(0,1)"},'FontSize',12);
title("Residual dataset1",'FontSize',16);

%% Refletions

% In the multi parameter optimisation process, fminsearch are likely to
% find the local minimum. It improved slightly of the predictions. 
% The optimization worked better on dataset1 with the ABB stock prices. We
% also plot the residual to see if it is normally distributed. For
% dataset0, we observe large residuals, which might be related to the Q4
% report release. 


%% Functions called in the scripts

dbtype('KalmanFilterFunc.m');
dbtype('max_like1.m');
dbtype('max_like2.m');

