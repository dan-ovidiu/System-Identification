clear
load('lab2_10.mat')

%plot(id.X,id.Y);
%hold on;

r = 29; % grad polinom
MSE = zeros(1,r);
for n = 1 : r
    xId = id.X;
    yId = id.Y;

    N1 = size(xId);
    N1 = N1(2);
    phi = zeros(N1,n);

    for i = 1 : N1
        for j = 1 : n
            phi(i, j) = xId(i)^(j-1);
        end
    end

    yTr = yId';
    thetta = phi \ yTr;
    yAproxId = phi * thetta;  % aproximarea pe datele de identificare

    %plot(xId, yAproxId)
    %xlabel('x cu datele identificare')
    %ylabel('y aprox. cu datele identificare')
    %legend('beforeAprox','afterAproxId');

    %pas2
    xVal = val.X;
    yVal = val.Y;
    N2 = size(xVal);
    N2 = N2(2);

    phiVal = zeros(N2,n);
    for i = 1 : N2
        for j = 1 : n
            phiVal(i, j) = xVal(i)^(j-1);
        end
    end

    yAproxVal = phiVal * thetta; % aproximarea pe datele de validare

    %plot(xVal,yAproxVal);
    %xlabel('x pe datele de validare');
    %ylabel('y aprox. pe datele de validare');
    %legend('beforeAprox','afterAproxVal');

    %pas3: Calcul MSE
    MSE(n) = calculateMSE(yAproxVal,yVal);
end

plot(MSE);
xlabel('Polynom grade');
ylabel('MSE value for each polynom');
title('MSE graphic(1 < n < 30)');

m = min(MSE);
bestGrade = find(MSE == m);
fprintf("Minimum MSE value is %d and the best polynom grade is %d.\n",m,bestGrade);

function phi = calcPhi(N,n,x)
   phi = zeros(N,n);
   for i = 1 : N
        for j = 1 : n
            phi(i, j) = x(i)^(j-1);
        end
    end
end


function MSE = calculateMSE(yAprox, y)
N = size(y,2);
e = yAprox' - y;
MSE = sum(e.^2) / N;
end








