clear all;
options = optimset('MaxFunEvals',250000);



%theta0 = [0.1 5 5 0.1];
x0 = 0;
p0 = 0;

fun = @(theta)max_like1(theta, x0, p0);


for i = 1:1:100
    theta0 = 10 * rand(1,4);
    result(i,:) = fminsearch(fun,theta0,options);
    y(i) = max_like1(result(i,:),x0,p0);
    i
end
y = y';
new_result = [result y];