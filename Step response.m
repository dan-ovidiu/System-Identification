clear
load("y.mat");
t = 0 : 0.01 : 3;

u = @(t)0.7*(t >= 0.3 & t <= 1) + 0.4*(t >= 1.3 & t <= 2) + (-0.5)*(t >= 2.3 & t <= 3);

% y = DCMRun.run(u(t));
plot(y)
%%
tId = 0 : 0.01 : 3;
xId = u(tId);
yId = y;

% plot(xId);
% % figure;
plot(yId);
hold on;

uss = 0.7;
u0 = 0;

yss = 2169.085; % media valorilor citite de pe datele de id.
y0 = 0;

K = (yss - y0) / (uss - u0);
T = 0.037;

H = tf(K,[T,1]);
output = lsim(H,xId,tId);
plot(output);