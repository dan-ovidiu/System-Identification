load("lab6_8.mat")

uId = id.InputData{1};
yId = id.OutputData{1};

uVal = val.InputData{1};
yVal = val.OutputData{1};

N = size(yId,1);
na = 8;
nb = 10;
phiId = calcPhi(uId,yId,na,nb);
thetta = phiId \ yId;

phiVal = calcPhi(uVal,yVal,na,nb);
yAproxVal = phiVal * thetta;

% predictie
plot(yVal)
hold on
plot(yAproxVal)
title("Predictie");
legend('real','approximated')

% simulare
ySim = calcYSim(uVal,na,nb,thetta);
figure;
plot(yVal)
hold on
plot(ySim)
title("Simulare")
legend('real','simulated')

function phi = calcPhi(u,y,na,nb)
N = size(y,1);
phi1 = zeros(N,na);
phi2 = zeros(N,nb);

% phi1 (na x N)
for i = 1 : N
    for j = 1 : na
        if(i - j <= 0)
            phi1(i,j) = 0;
        else
            phi1(i,j) = -y(i - j);
        end
    end
end

% phi2 (nb x N)
for i = 1 : N
    for j = 1 : nb
        if(i - j <= 0)
            phi2(i,j) = 0;
        else
            phi2(i,j) = u(i - j);
        end
    end
end
phi = [phi1 phi2];
end

function ySim = calcYSim(u,na,nb,thetta)
N = size(u,1);
ySim = zeros(N,1);
for i = 1 : N
    for k1 = 1 : na
        if(i - k1 > 0)
            ySim(i) = ySim(i) - ySim(i - k1) * thetta(k1);
        end
    end

    for k2 = 1 : nb
        if(i - k2 > 0)
            ySim(i) = ySim(i) + u(i - k2) * thetta(na + k2);
        end
    end
end
end
