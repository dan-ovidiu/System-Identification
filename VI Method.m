load("lab9_3.mat")
close all
input_id = id.InputData;
output_id = id.OutputData;
Ts_id = id.Ts;
    
input_val = val.InputData;
output_val = val.OutputData;
Ts_val = val.Ts;

na = n + 1;
nb = n + 1;
phiId = calcPhi(input_id,output_id,na,nb);
thetta = phiId \ output_id;

data_id = iddata(output_id(:),input_id(:),Ts_id);
data_val = iddata(output_val(:),input_val(:),Ts_val);

nk = 1;
model = arx(data_id,[na nb nk]);

ySim = compare(data_val,model).OutputData;
phi_sim = calcPhi(input_val, ySim, na, nb);

y_VI = generateVI(input_val,output_val,ySim,phi_sim,na,nb);

model_VI = iddata(y_VI(:),input_val(:),Ts_val);

y_VI_simple = generateVI_simple(input_val,output_val,phi_sim,na,nb);
model_VI_simple = iddata(y_VI_simple(:),output_val(:),Ts_val);

compare(data_val,model_VI,model_VI_simple);
title("VI methods");

model_ARX = iddata(ySim(:),input_val(:),Ts_val);
figure, compare(data_val, model);
title("ARX method");

function phi = calcPhi(u,y,na,nb)
N = size(y,1);
phi1 = zeros(N,na);
phi2 = zeros(N,nb);

for i = 1 : N
    for j = 1 : na
        if(i - j <= 0)
            phi1(i,j) = 0;
        else
            phi1(i,j) = -y(i - j);
        end
    end
end

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


function y_VI = generateVI(u,y,ySim,phi_sim,na,nb)
N = size(ySim,1);
Z = zeros(N,na+nb);
phi_tilda=zeros(na+nb,na+nb);
Y_tilda=zeros(na+nb,1);

for k = 1 : N
    for i = 1 : na
        if(k - i > 0)
           Z(k,i) = -ySim(k - i);
        end
    end
    for j = 1 : nb
        if(k - j > 0)
            Z(k,j + na) = u(k - j);
        end
    end
    phi_tilda = phi_tilda + (Z(k,:)' * phi_sim(k,:));
    Y_tilda = Y_tilda + (Z(k,:)' * y(k));
end
phi_tilda = phi_tilda / N;
Y_tilda = Y_tilda / N;

new_theta = phi_tilda \ Y_tilda;
y_VI = phi_sim * new_theta;
end



function y_VI_simple = generateVI_simple(u,y,phi_sim,na,nb)
N = size(y,1);
Z = zeros(N,na+nb);
phi_tilda=zeros(na+nb,na+nb);
Y_tilda=zeros(na+nb,1);

for k = 1 : N
    for i = 1 : na
        if(k - nb - i > 0)
           Z(k,i) = u(k - nb - i);
        end
    end
    for j = 1 : nb
        if(k - j > 0)
            Z(k,j + na) = u(k - j);
        end
    end
    phi_tilda = phi_tilda + (Z(k,:)' * phi_sim(k,:));
    Y_tilda = Y_tilda + (Z(k,:)' * y(k));
end
phi_tilda = phi_tilda / N;
Y_tilda = Y_tilda / N;

new_theta = phi_tilda \ Y_tilda;
y_VI_simple = phi_sim * new_theta;
end


