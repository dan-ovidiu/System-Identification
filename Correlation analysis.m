clear
t = 0 : 0.01 : 8;
u = zeros(1,800);
randomVector = unifrnd(-0.7,0.7,[499,1]);
contor = 1;
for i = 1 : 600
    if(i < 30)
        u(i) = 0;
    else
        if(i < 529)
            u(i) = randomVector(contor);
            contor = contor + 1;
        else
            u(i) = 0;
        end
    end
end
u(600:670) = 0.2;
% y = DCMRun.run(u);
%%
close, clear
load("y.mat");
load("u.mat");

uId = u(31:529);
yId = y(31:529);

uVal = u(600:670);
yVal = y(600:670);

uDet = detrend(uId);
yDet = detrend(yId);

M = 80;
N = size(uDet,2);

Y = calcMatrixY(uDet,yDet,N);

MSE = zeros(1,M);
for m = 1 : M
    phi = calcPhi(m,N,uDet);
    H = phi \ Y;
    
    yIdHat = conv(H,uDet);
    yIdHat = yIdHat(1:length(yId));
    MSE(m) = calculateMSE(yIdHat,yId');
end
plot(MSE)
title('MSE');

minValue = min(MSE);
bestM = find(MSE == minValue);
fprintf("Minimum MSE value is %d and the best M is %d\n",minValue, bestM);

% calculam phi si H pentru cel mai bun M si generam graficele
phi = calcPhi(bestM,N,uDet);
H = phi \ Y;

yIdHat = conv(H,uDet);
figure;
plot(yDet), grid, hold on
plot(yIdHat);
title('Date identificare pentru M = 50');
legend('yId','yIdHat');

yValHat = conv(H,uVal);
figure;
plot(yVal), grid, hold on
plot(yValHat);
title('Date validare pentru M = 50');
legend('yVal','yValHat');

figure;
plot(H), grid
title('Raspunsul la impuls pentru M = 50');

% functions

% ru correlation
function ru = calcCorrelationU(N,tau,u)
ru = 0;
for k = 1 : N - tau
    ru = ru + u(k + tau) * u(k);
end
ru = ru/N;
end

% ryu correlation
function ryu = calcCorrelationUY(N,tau,y,u)
ryu = 0;
for k = 1 : N - tau
    ryu = ryu + y(k + tau) * u(k);
end
ryu = ryu/N;
end

% Y
function Y = calcMatrixY(u,y,N)
Y = zeros(N,1);
for tau = 0 : N - 1
    Y(tau + 1) = calcCorrelationUY(N,tau,y,u);
end
end

% phi
function phi = calcPhi(M,N,u)
phi = zeros(N,M);

% generam vectorul ru de N elemente
ru = zeros(1,N);
for tau = 0 : N - 1
    ru(tau + 1) = calcCorrelationU(N,tau,u);
end

for n = 1 : N
    for m = 1 : M
        phi(n,m) = ru(abs(n-m) + 1);
    end
end
end

% MSE
function MSE = calculateMSE(yAprox, y)
N = size(y,2);
e = yAprox' - y;
MSE = sum(e.^2) / N;
end