clear all;

%theta = [theta(1) theta(2), theta(3),theta(4)]
%theta = [Fk Hk Qk Rk]
%theta = [0.9530    0.0139  118.9087    1.0583];
theta =  [0.2752    1.4364    0.4473    0.2721];
load('dataset0.mat');
x0 = 0;
p0 = 0;

x1(1) = theta(1) * x0;
p1(1) = theta(1) * p0 * theta(1)' +theta(3);

sk(1) = theta(2)*p1(1)*theta(2)' + theta(4);
L(1) = -0.5 * log(2*pi) - 0.5 * log(abs(sk(1))) - 0.5 * (data(1) - theta(2) * x1(1))' * sk(1)^-1 * (data(1) - theta(2) * x1(1));
res(1) = (data(1) - theta(2) * x1(1))/sk(1);
for i = 1:1:length(data)-1
  
    K(i) = p1(i) * theta(2)' * ((theta(2) * p1(i) * theta(2)' + theta(4))^(-1));
    
    xm(i) = x1(i) + K(i) * (data(i) - theta(2)*x1(i));
    pm(i) = (1 - K(i)*theta(2))*p1(i);
    
        
    x1(i+1) = theta(1) * xm(i);
    p1(i+1) = theta(1) * pm(i) * theta(1)' + theta(3);
    sk(i+1) = theta(2) * p1(i+1) *theta(2)' + theta(4);
    res(i+1) = (data(i+1) - theta(2) * x1(i+1))/sk(i+1);
    
        

end
histu = hist(res,[-2.5:0.1:2.5]);
x = -2.5:0.1:2.5;
bar(x,histu);
xlabel('residule');
ylabel('frequency');


 