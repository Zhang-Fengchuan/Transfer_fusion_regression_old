function [simulation_size, X_target, Y_target, sample_size, p, beta_target_real] = Data_generation_setting_2(target_if, h)
if target_if == 1
    % target study data generation
    rng(1,'philox'); % setting the seed
    p = 20; % the dimension of coefficient 25
    simulation_size = 100; % the number of repeat %!!!!!!!!!!!!!!!!!!!!!!!!
    sample_size = 150; % the sample size 400
    rho = 0.6; % autocorrelation coefficient of covariance
    sigma = zeros(p,p)*nan; % create an empty covariance matrix
    % create an covariance matrix with autocorrelation
    for i = 1 : p
        for j = 1 : p
            sigma(i, j) = rho^(abs(i-j));
        end
    end
    [sigma_vec,sigma_value] = eig(sigma); % Perform eigen decomposition on the covariance matrix
    sigma_sef = sigma_vec*sqrt(sigma_value); % Decompose the covariance matrix into the product of an orthogonal matrix and its transpose
    X_target = zeros(sample_size, p, simulation_size); % Create matrix x (100 simulations, each with 60 vector samples).
    for i = 1:simulation_size
        for j = 1:sample_size
            X_target(j,:,i) = (sigma_sef*randn(p, 1))';
        end
    end
    mu = 1;
    beta_target_real1 = repmat([mu*ones(p,1)], floor((1/3)*sample_size), 1);
    beta_target_real2 = repmat([-0.5*mu*ones(0.5*p,1);0.5*mu*ones(0.5*p,1)], (sample_size-2*floor((1/3)*sample_size)), 1);
    beta_target_real3 = repmat([-mu*ones(p,1)], floor((1/3)*sample_size), 1);
    beta_target_real = [beta_target_real1; beta_target_real2; beta_target_real3];
    epsilon = 0.2 * randn(sample_size, simulation_size); % create the error epsilon
    Y_target = zeros(sample_size, simulation_size)*nan; % create the initial response variable
    for i = 1: simulation_size
        for j = 1: sample_size
            Y_target(j, i) = X_target(j,:,i)*beta_target_real(((j-1)*p+1):j*p, 1) + epsilon(j, i);
        end
    end
else
    % auxiliary study data generation
    rng(2,'philox'); % setting the seed
    p = 20; % the dimension of coefficient
    simulation_size = 100; % the number of repeat %!!!!!!!!!!!!!!!!!!!!!!!!
    sample_size = 500; % the sample size
    h = h/p;
    sigma_B = h*eye(p);
    [sigma_vec,sigma_value] = eig(sigma_B); % Perform eigen decomposition on the covariance matrix
    sigma_sef = sigma_vec*sqrt(sigma_value); % Decompose the covariance matrix into the product of an orthogonal matrix and its transpose
    error_B = (sigma_sef*randn(p, 1));
    mu = 1;
    beta_auxiliary_real1 = repmat(([mu*ones(p,1)]+error_B), floor((1/3)*sample_size), 1);
    beta_auxiliary_real2 = repmat(([-0.5*mu*ones(0.5*p,1);0.5*mu*ones(0.5*p,1)]+error_B),...
                                    (sample_size-2*floor((1/3)*sample_size)), 1);
    beta_auxiliary_real3 = repmat(([-mu*ones(p,1)]+error_B), floor((1/3)*sample_size), 1);
    %生成辅助域四的样本X_Auxiliary4
    rho = 0.6; % autocorrelation coefficient of covariance
    sigma = zeros(p,p)*nan; % create an empty covariance matrix
    % create an covariance matrix with autocorrelation
    for i = 1 : p
        for j = 1 : p
            sigma(i, j) = rho^(abs(i-j));
        end
    end
    [sigma_vec,sigma_value] = eig(sigma); % Perform eigen decomposition on the covariance matrix
    sigma_sef = sigma_vec*sqrt(sigma_value); % Decompose the covariance matrix into the product of an orthogonal matrix and its transpose
    X_Auxiliary = zeros(sample_size, p, simulation_size); % Create matrix x (100 simulations, each with 60 vector samples).
    for i = 1:simulation_size
        for j = 1:sample_size
            X_Auxiliary(j,:,i) = (sigma_sef*randn(p, 1))';
        end
    end
    %生成辅助域四的样本Y_Auxiliary4
    beta_Auxiliary_real = [beta_auxiliary_real1; beta_auxiliary_real2; beta_auxiliary_real3];
    epsilon = 0.2 * randn(sample_size, simulation_size); % create the error epsilon
    Y_Auxiliary = zeros(sample_size, simulation_size)*nan; % create the initial response variable
    for i = 1: simulation_size
        for j = 1: sample_size
            Y_Auxiliary(j, i) = X_Auxiliary(j,:,i)*beta_Auxiliary_real(((j-1)*p+1):j*p, 1) + epsilon(j, i);
        end
    end
    X_target = X_Auxiliary;
    Y_target = Y_Auxiliary;
    beta_target_real = beta_Auxiliary_real;
end
end