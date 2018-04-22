function y = max_like2(theta,x0,p0)
    load('dataset1.mat');

    t = datenum(dtime);
    T0 = t(1);
    t = t-T0;


    x1(1) = theta(1)*x0;
    p1(1) = theta(1)*p0*theta(1)'+theta(3);


    sk(1) = theta(2)*p1(1)*theta(2)'+theta(4);
    L(1) = -0.5*log(2*pi)-0.5*log(abs(sk(1)))-0.5*(data(1)-theta(2)*x1(1))'*sk(1)^-1*(data(1)-theta(2)*x1(1));

    for i=1:1:length(data)-1

        K(i) = p1(i)*theta(2)'*((theta(2)*p1(i)*theta(2)'+theta(4))^(-1));

        xm(i) = x1(i)+K(i)*(data(i)-theta(2)*x1(i));
        pm(i) = (1-K(i)*theta(2))*p1(i);

        x1(i+1) = theta(1)*xm(i);
        p1(i+1) = theta(1)*pm(i)*theta(1)'+theta(3);
        sk(i+1) = theta(2)*p1(i+1)*theta(2)'+theta(4);
        L(i+1) = -0.5*log(2*pi)-0.5*log(abs(sk(i+1)))-0.5*(data(i+1)-theta(2)*x1(i+1))'*(sk(i+1)^(-1))*(data(i+1)-theta(2)*x1(i+1));
    end

    y = -sum(L);
end
 