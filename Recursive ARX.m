N = 200;

% id
interval = [-0.8, 0.8];
u_id = idinput(N,'PRBS',[],interval);

% val
u_val = zeros(1,110);
u_val(20:90) = 0.3;

% y_val = zeros(1,length(u_val));
% motor = DCMRun.start('Ts',10e-3);
% for k = 1 : length(u_val)
%   y_val(k) = motor.step(u_val(k));
%   motor.wait();
% end
% plot(y_val)
%%
% plot(y_val)
load("y_val.mat");
load("thetta_RARX.mat");

N = 200;
na = 6;
nb = 6;
Ts = 10e-3;
thetta0 = zeros(na + nb, 1);
Pinv0 = 1000 * eye(na + nb,na + nb);

y_id = zeros(1,length(u_id));
% thetta_RARX = ARX_rec(u_id,y_id,N,na,nb,thetta0,Pinv0);
thetta = thetta_RARX';
A = [1 thetta(200, 1:na)];
B = [0 thetta(200, na + 1 : end)];
model = idpoly(A,B,[],[],[],0,Ts);
data_val = iddata(y_val(:),u_val(:),Ts);
compare(model,data_val);

function thetta_ARX_rec = ARX_rec(u,y,N,na,nb,thetta0,Pinv0)
thetta_ARX_rec(:,1) = thetta0;
current_P = zeros(na + nb, na + nb);
last_P = Pinv0;
motor = DCMRun.start('Ts',10e-3);
for k = 1 : N
  y(k) = motor.step(u(k));
  phi = calc_phi(u,y,k,na,nb);
  % phi
  % size(phi)
  % size(thetta_ARX_rec(:,k))
  e_pred(k) = y(k) - phi * thetta_ARX_rec(:,k);
  current_P = last_P - (last_P * phi' * phi * last_P)/(1 + phi * last_P * phi');
  % current_P
  last_P = current_P;
  % size(phi)
  % size(current_P)
  W = current_P * phi';
  thetta_ARX_rec(:,k + 1) = thetta_ARX_rec(:,k) + W * e_pred(k);
  % thetta_ARX_rec(:,k + 1) = thetta_ARX_rec(:,k);
  motor.wait();
end
end

function phi = calc_phi(input,output,k,na,nb)
phi = zeros(1,na + nb);
for k1 = 1 : na
    if(k - k1 > 0)
        phi(k1) = -output(k - k1);
    end
end
for k2 = 1 : nb
    if(k - k2 > 0)
        phi(k2 + na) = input(k - k2);
    end
end
end