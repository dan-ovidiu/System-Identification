clear
 load("y.mat");
t = 0 : 0.01 : 2;
u = @(t)(0.1*(t > 0) + 0.9*(t >= 0.40 & t < 0.42) + 0.9*(t >= 1.20 & t < 1.22) + 0.9*(t >= 1.70 & t < 1.72));
plot(u(t))
% plot(t,u(t))
% hold on
% y = DCMRun.run(u(t));
% plot(y);
%%
t = 0 : 0.01 : 2;

% gasim T,K pe datele de identificare (primul impuls)
 % plot(t(1:65),y(1:65)) 
% valT = 0.368*(max(y(1:60)) - 270.69)% valoarea unde il gasim pe T pe
% datele de identificare => x = 0.51;
y0 = 0.488;
T = (0.51 - y0);

% plot(t,y)
% figure;
yss = 270.69;
uss = 0.1;
K = yss / uss;

A = -1/T;
B = K/T;
C = 1;
D = 0;

H = ss(A,B,C,D);
output = lsim(H,u(t),t,yss);

% plot(t,u(t));
% figure;
plot(t,y);
hold on;
plot(t,output);
legend('y','aprox');

calculateMSE(output,y)


function MSE = calculateMSE(yAprox, y)
N = size(y,2);
e = yAprox' - y;
MSE = sum(e.^2) / N;
end
