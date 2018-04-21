clear all;
options = optimset('MaxFunEvals',20000);

theta0 = [0.1 5 5 0.1];
x0 = 0;
p0 = 0;

fun = @(theta)max_like1(theta, x0, p0);

result = fminsearch(fun, theta0);
