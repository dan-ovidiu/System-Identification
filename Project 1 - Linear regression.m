close all; clear; clc
load("proj_fit_14.mat")

% id
m = 16;
X1_id = id.X{1};
X2_id = id.X{2};
Y_id = id.Y;
N1_id = id.dims(1);
N2_id = id.dims(2);

phi_id = calc_phi(m, X1_id, X2_id, N1_id, N2_id);
Y_id_reshape = reshape(Y_id, [1, N1_id*N2_id]);
theta = phi_id \ Y_id_reshape';
Y_id_aprox = phi_id * theta;

Y_id_aprox_reshape = reshape(Y_id_aprox, [N1_id, N2_id]);
figure; mesh(X1_id, X2_id, Y_id); hold on; mesh(X1_id, X2_id, Y_id_aprox_reshape, 'EdgeColor', 'r')
xlabel("X1 id"); ylabel("X2 id"); zlabel("Y id")
MSE_id = calc_best_MSE(X1_id, X2_id, N1_id, N2_id, Y_id, 50, "identification");

% val
X1_val = val.X{1};
X2_val = val.X{2};
Y_val = val.Y;
N1_val = val.dims(1);
N2_val = val.dims(2);

phi_val = calc_phi(m, X1_val, X2_val, N1_val, N2_val);
% Y_val_reshape = reshape(Y_val, [1, N1_val*N2_val]);
Y_val_aprox = phi_val * theta;

Y_val_aprox_reshape = reshape(Y_val_aprox, [N1_val, N2_val]);
figure; mesh(X1_val, X2_val, Y_val); hold on; mesh(X1_val, X2_val, Y_val_aprox_reshape, 'EdgeColor', 'r')
xlabel("X1 val"); ylabel("X2 val"); zlabel("Y val");
MSE_val = calc_best_MSE(X1_val, X2_val, N1_val, N2_val, Y_val, 50, "validation");

MSEs = abs(MSE_val - MSE_id);
figure; plot(MSEs); title("abs(MSE_v_a_l - MSE_i_d)"); xlabel("grad pol"); ylabel("error"); grid on
fprintf('Best pol grad overall: %d, error: %d\n', find(MSEs==min(MSEs)), min(MSEs));

function phi = calc_phi(pol_grad, X1, X2, N1, N2)
    for j=1:N2
        for i=1:N1
           idx = N1*(j-1) + i;
           term_idx = 1;

           for p1 = 0:pol_grad
               for p2 = 0:pol_grad
                    if (p1+p2 > pol_grad)
                        continue
                    end

                    phi(idx, term_idx) = X1(i)^p1 * X2(j)^p2;
                    term_idx = term_idx + 1;
               end
           end
        end
    end
end


function MSE = calc_MSE(y_aprox, Y)
    N = size(Y, 2);
    epsilon = y_aprox' - Y;
    MSE = sum(epsilon.^2) / N;
end


function MSEs = calc_best_MSE(X1, X2, N1, N2, Y, pol_grad, type)
    MSEs = [];
    for m=1:pol_grad
        phi = calc_phi(m, X1, X2, N1, N2);
        Y_reshape = reshape(Y, [1, N1*N2]);
        theta = phi \ Y_reshape';
        Y_aprox = phi * theta;

        MSEs = [MSEs calc_MSE(Y_aprox, Y_reshape)];
    end

    figure; plot(MSEs); title(sprintf("MSE %s", type)); xlabel("grad pol"); ylabel("error"); grid on
    fprintf('Best pol grad %s: %d, error: %d\n', type, find(MSEs==min(MSEs)), min(MSEs));
end
