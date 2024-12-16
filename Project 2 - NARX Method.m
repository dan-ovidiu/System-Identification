clc, clear, close all

load("iddata-08.mat")

t_id = id_array(:,1);
u_id = id_array(:,2);
y_id = id_array(:,3);

t_val = val_array(:,1);
u_val = val_array(:,2);
y_val = val_array(:,3);

calc_NARX(4, 5, y_id, u_id, y_val, u_val);


function calc_NARX(m_max, n_max, y_id, u_id, y_val, u_val)
    best_mse = inf; 
    best_theta = 0;
    best_y_sim = 0;
    best_m = 0;
    best_n = 0;

    for m=1:m_max
        for n=1:n_max
            phi_id = calc_phi(u_id, y_id, m, n, n);
            
            lambda = 1;
            theta = inv(phi_id' * phi_id + lambda * eye(size(phi_id, 2))) * phi_id' * y_id;

            y_sim = calc_Ys(u_val, theta, m, n, n);

            MSE(m, n) = calc_MSE(y_sim, y_val);

            if (best_mse > MSE(m, n))
                best_theta = theta;
                best_y_sim = y_sim;
                best_m = m;
                best_n = n;
                best_mse = MSE(m, n);
            end
        end
    end

    phi_val = calc_phi(u_val,  y_val, best_m, best_n, best_n);
    y_val_pred = phi_val * best_theta;
    
    figure, plot(y_val, 'r'), hold on, plot(y_val_pred, "b--"), title("y_v_a_l & y_p_r_e_d"), legend("y_v_a_l", "y_p_r_e_d")
    figure, plot(y_val, 'r'), hold on, plot(best_y_sim, "b--"), title("y_v_a_l & y_s_i_m"), legend("y_v_a_l", "y_s_i_m")
    figure, imagesc(MSE), colorbar, title('MSE'), xlabel('na = nb'), ylabel('m')
    figure, surf(MSE), shading interp, title('MSE surface'), xlabel('na = nb'), ylabel('m'), zlabel('mse value')

    fprintf("best m: %d, best na=nb: %d, mse value: %d\n", best_m, best_n, best_mse)
end


function phi = calc_phi(u, y, m, na, nb)
    y_size = size(y, 1);
    d_k = zeros(y_size, na+nb);

    pows_len = size(d_k, 2);
    pow_combinations = find_combinations(m, pows_len);
    pow_combinations_len = size(pow_combinations, 1);

    phi = zeros(y_size, 1);

    for k=1:y_size
        for i=1:na
            if (k-i > 0)
                d_k(k, i) = y(k-i);
            end
        end
    
        for j=1:nb
            if (k-j > 0)
                d_k(k, na+j) = u(k-j);
            end
        end

        for i=1:pow_combinations_len
            phi(k, i) = prod(d_k(k, :) .^ pow_combinations(i, :));
        end
    end
end


function Ys = calc_Ys(u, theta, m, na, nb)
    u_size = size(u, 1);
    d_k = zeros(u_size, na+nb);

    pows_len = size(d_k, 2);
    pow_combinations = find_combinations(m, pows_len);
    pow_combinations_len = size(pow_combinations, 1);

    Ys = zeros(u_size, 1);

    for k=1:u_size
        for i=1:na
            if (k-i > 0)
                d_k(k, i) = Ys(k-i);
            end
        end

        for j=1:nb
            if (k-j > 0)
                d_k(k, na+j) = u(k-j);
            end
        end

        for i=1:pow_combinations_len
            Ys(k) = Ys(k) + (prod(d_k(k, :) .^ pow_combinations(i, :))) * theta(i);
        end
    end
end


function valid_combinations = find_combinations(max_sum, num_numbers)
    numbers = 0:max_sum;
    valid_combinations = generate_combinations([], numbers, num_numbers, max_sum);
end


function combinations = generate_combinations(current, numbers, num_numbers, max_sum)
    if num_numbers == 0
        if sum(current) <= max_sum
            combinations = current;
        else
            combinations = [];
        end
        return;
    end

    combinations = [];

    for i = 1:length(numbers)
        next = [current, numbers(i)];
        new_combinations = generate_combinations(next, numbers, num_numbers - 1, max_sum);
        if ~isempty(new_combinations)
            combinations = [combinations; new_combinations];
        end
    end
end


function MSE = calc_MSE(y_aprox, y)
    N = size(y, 2);
    epsilon = y_aprox - y;
    MSE = sum(epsilon.^2) / N;
end


function MSEs = calc_best_MSE_val(X1_id, X2_id, N1_id, N2_id, Y_id, X1_val, X2_val, N1_val, N2_val, Y_val, pol_grad, type)
    MSEs = [];
    for m=1:pol_grad
        phi_id = calc_phi(m, X1_id, X2_id, N1_id, N2_id);
        Y_id_reshape = reshape(Y_id, [1, N1_id*N2_id]);
        theta = phi_id \ Y_id_reshape';

        phi_val = calc_phi(m, X1_val, X2_val, N1_val, N2_val);
        Y_val_reshape = reshape(Y_val, [1, N1_val*N2_val]);
        Y_aprox = phi_val * theta;

        MSEs = [MSEs calc_MSE(Y_aprox, Y_val_reshape)];
    end

    figure; plot(MSEs); title(sprintf("MSE %s", type)); xlabel("grad pol"); ylabel("error"); grid on
    fprintf('Best pol grad %s: %d, error: %d\n', type, find(MSEs==min(MSEs)), min(MSEs));
end
