clear
close all
clc

%% Kalman filter
load('dataset0.mat');

Fk = 0.1;
Hk = 1;
Qk = 1;
Rk = 1;

x0 = 0;
P0 = 0;

xPred = zeros(length(data),1);
PPred = zeros(length(data),1);

xUpd = zeros(length(data),1);
PUpd = zeros(length(data),1);

xPred(1) = x0;
PPred(1) = P0;

K = PPred(1)*Hk/(Hk*PPred(1)*Hk+Rk);

xUpd(1) = xUpd(1)+K*(data(1)-Hk*xUpd(1));
PUpd(1) = (1-K*Hk)*PUpd(1);

for i=1:length(data)-1
    xPred(i+1) = Fk*xUpd(i);
    PPred(i+1) = Fk*PUpd(i)*Fk+Qk;
    
    K = PPred(i+1)*Hk/(Hk*PPred(i+1)*Hk+Rk);
    
    xUpd(i+1) = xPred(i+1)+K*(data(i+1)-Hk*xPred(i+1));
    PUpd(i+1) = (1-K*Hk)*PPred(i+1);
end

figure(1)
plot(dtime,data);
hold on
plot(dtime,xUpd);
xlabel("Time");
ylabel("Data");
legend("Market data", "Filtered estimates");
grid on

%% Optimization
x0 = [0.1 5 5 0.1];
optParameters = fminsearch(@objectiveFunc,x0)