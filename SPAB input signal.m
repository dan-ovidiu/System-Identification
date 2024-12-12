N = 200;
a = -0.7;
b = 0.7;

m1 = 3;
SPAB1 = generateSPAB(N,m1,a,b);

m2 = 10;
SPAB2 = generateSPAB(N,m2,a,b);

u = zeros(1,100);
u(20:90) = 0.4;

input = [zeros(1,30),SPAB1,zeros(1,30),SPAB2,zeros(1,30),u];

u_id1 = input(1:250);
u_id2 = input(251:480);
u_val = input(481:590);

% y = DCMRun.run(input);
%%
close all
load("y.mat")
y_id1 = y(1:250);
y_id2 = y(251:480);
y_val = y(481:590);

Ts = 0.01;
id1 = iddata(y_id1(:),u_id1(:),Ts);
id2 = iddata(y_id2(:),u_id2(:),Ts);

na_id1 = 3;
nb_id1 = 3;
nk_id1 = 1;

model1 = arx(id1,[na_id1,nb_id1,nk_id1]);
% compare(model1,id1)

na_id2 = 3;
nb_id2 = 3;
nk_id2 = 1;

model2 = arx(id2,[na_id2,nb_id2,nk_id2]);
figure
% compare(model2,id2)

val = iddata(y_val(:),u_val(:),Ts);
figure
compare(val,model1,model2);

% figure
% compare(val,model2);

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