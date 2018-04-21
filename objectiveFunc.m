function logLikelihood = objectiveFunc(x)

    load('dataset0.mat');

    x0 = 0;
    P0 = 0;

    xPred = zeros(length(data),1);
    PPred = zeros(length(data),1);

    xUpd = zeros(length(data),1);
    PUpd = zeros(length(data),1);

    xPred(1) = x0;
    PPred(1) = P0;

    K = PPred(1)*x(2)/(x(2)*PPred(1)*x(2)+x(4));

    xUpd(1) = xUpd(1)+K*(data(1)-x(2)*xUpd(1));
    PUpd(1) = (1-K*x(2))*PUpd(1);

    for i=1:length(data)-1
        xPred(i+1) = x(1)*xUpd(i);
        PPred(i+1) = x(1)*PUpd(i)*x(1)+x(3);

        K = PPred(i+1)*x(2)/(x(2)*PPred(i+1)*x(2)+x(4));

        xUpd(i+1) = xPred(i+1)+K*(data(i+1)-x(2)*xPred(i+1));
        PUpd(i+1) = (1-K*x(2))*PPred(i+1);
    end

    logp = zeros(length(data),1);
    
    for k=1:length(dtime)
        Sk = x(2)*PPred(k)*x(2)+x(4);
        logp(k) = -0.5*log(2*pi)-0.5*log(abs(Sk))-0.5*(data(k)-x(2)*xPred(k))*(1/Sk)*(data(k)-x(2)*xPred(k));
    end
    logLikelihood = 1/sum(logp);
end