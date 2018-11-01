% HW3 
% PB16000953 ∆Ó¡’∑Â
% 2018.10.22
x = 0.25:0.25:2.5;  % xdata
y = [1.284, 1.648, 2.117, 2.718, 3.427, 2.798, 3.534, 4.456, 5.465, 5.894]; % ydata
% solve Ax=b
A = [sum(sin(x).^2),sum(sin(x).*cos(x));sum(cos(x).*sin(x)),sum(cos(x).^2)];
b = [sum(y.*sin(x));sum(y.*cos(x))];
X = inv(A)*b;
% draw img
scatter(x,y);
hold on;
res = X(1)*sin(x)+X(2)*cos(x);
plot(x,res);
% res
fprintf('a=%.15f, b=%.15f, error = %.15f\n', X(1), X(2), sum((y-res).^2));
