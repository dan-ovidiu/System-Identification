% generate SPAB
m = 10;
N = 200;
a = -0.7;
b = 0.7;

SPAB = generateSPAB(N,m,a,b);

u = zeros(1,100);
u(20:90) = 0.4;

input = [zeros(1,30),SPAB,zeros(1,30),u,zeros(1,30)];

u_id = input(1:250);
u_val = input(251:390);

% y = DCMRun.run(input);
%%
close all
load("y.mat");
y_id = y(1:250);
y_val = y(251:390);

alpha = 0.1;
thetta0 = [0.01;100];
delta = 1e-5;
lmax = 10000;
nk = 3;

thetta = implementAlg(input,y,nk,alpha,thetta0,delta,lmax);
f = thetta(1);
b = thetta(2);
Ts = 0.01;

A = 1;
C = 1;
D = 1;
F = [1 f];
B = [zeros(1,nk) b];

mOE = idpoly(A,B,C,D,F,0,Ts);

data_id = iddata(y_id(:),u_id(:),Ts);
data_val = iddata(y_val(:),u_val(:),Ts);

figure, compare(data_id,mOE);
figure, compare(data_val,mOE);

% plot(data_id)
% figure
% plot(data_val)

function thetta_res = implementAlg(u,y,nk,alpha,thetta0,delta,lmax)
N = length(u);

l = 0;
thetta = thetta0;

while(1)
    l = l + 1;
    if l == lmax
        break;
    end

    e = zeros(1,N);
    dE_dThetta = zeros(2,N);

    for k = 1 : nk
        e(k) = y(k);
        dE_dThetta(:,k) = [0 0];
    end

    for k = nk + 1 : N
        f = thetta(1,l);
        b = thetta(2,l);
        e(k) = y(k) + f * y(k - 1) - b * u(k - nk) - f * e(k - 1);
        dE_dThetta(1,k) = y(k - 1) - e(k - 1) - f * dE_dThetta(1,k - 1);
        dE_dThetta(2,k) = -u(k - nk) - f * dE_dThetta(2,k - 1);
    end

    % gradient
    dV_dThetta = 0;
    for k = 1 : N
        dV_dThetta = dV_dThetta + e(k) * dE_dThetta(:,k);
    end
    dV_dThetta = 2 / (N - nk) * dV_dThetta;
    % disp(dV_dThetta)

    % hessian
    H = 0;
    for k = 1 : N
        H = H + dE_dThetta(:,k) * dE_dThetta(:,k)';
    end
    H = 2 / (N - nk) * H;

    thetta(:,l+1) = thetta(:,l) - alpha * (H \ dV_dThetta);

    if norm(thetta(l+1) - thetta(l)) <= delta
        break;
    end
end
    thetta_res = thetta(:,l);
end


function SPAB = generateSPAB(N,m,a,b)
SPAB = zeros(1,N);
A = generateMatrixA(m);
x = zeros(m,N);

x(1,1) = 1;
x(1,3) = 1;
for k = 1 : N
    x(:,k+1) = mod(A*x(:,k),2);
end

for k = 1 : N
    SPAB(k) = a + (b-a)*x(m,k);
end

end

function A = generateMatrixA(m)
A = zeros(m,m);
a = zeros(1,m);
if(m == 3)
    a(1) = 1;
    a(3) = 1;
else
    if(m == 4)
        a(1) = 1;
        a(4) = 1;
    else
        if(m == 5)
            a(2) = 1;
            a(5) = 1;
        else
            if(m == 6)
                a(1) = 1;
                a(6) = 1;
            else
                if(m == 7)
                    a(1) = 1;
                    a(7) = 1;
                else
                    if(m == 8)
                        a(1) = 1;
                        a(2) = 1;
                        a(7) = 1;
                        a(8) = 1;
                    else
                        if(m == 9)
                            a(4) = 1;
                            a(9) = 1;
                        else
                            if(m == 10)
                                a(3) = 1;
                                a(10) = 1;
                            end
                        end
                    end
                end
            end
        end
    end
end
I = eye(m - 1);
emptyColumn = zeros(m-1, 1);
A = [a;I emptyColumn];
end