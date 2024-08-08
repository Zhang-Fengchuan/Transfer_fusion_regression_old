function [simulation_size, X_target, Y_target, sample_size, p, beta_target_real] = Data_generation_setting_1(which_index, h)
if which_index == 0
    rng(1,'philox'); % setting the seed
    p = 20; % the dimension of coefficient 25
    simulation_size = 100; % the number of repeat
    sample_size = 150; % the sample size 400
    rho = 0.6; % autocorrelation coefficient of covariance
    sigma = zeros(p,p)*nan; % create an empty covariance matrix
    % create an covariance matrix with autocorrelation
    for i = 1 : p
        for j = 1 : p
            sigma(i, j) = rho^(abs(i-j));
            %sigma(i, j) = 1;
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
    % 参数高维非稀疏设置 p=25
    % beta_target_real1 = repmat([mu*ones(p,1)], floor((1/3)*sample_size), 1);
    % beta_target_real2 = repmat([-0.5*mu*ones(0.4*p,1);0.5*mu*ones(0.6*p,1)], (sample_size-2*floor((1/3)*sample_size)), 1);
    % beta_target_real3 = repmat([-mu*ones(p,1)], floor((1/3)*sample_size), 1);
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
    % Extract the data needed for the i-th simulation
    yangbencishu = 1;
    y = Y_target(:,yangbencishu); %(n*1)
    x = X_target(:,:,yangbencishu); % (n*p)
elseif which_index == 1
    % 辅助域零参数设置
    rng(2,'philox'); % setting the seed
    p = 20; % the dimension of coefficient
    simulation_size = 100; % the number of repeat 
    sample_size = 500; % the sample size
    % h0 = 20;
    % h0 = h0/p;
    h = 5;
    sigma_A = h*eye(p);
    [sigma_vec,sigma_value] = eig(sigma_A); % Perform eigen decomposition on the covariance matrix
    sigma_sef = sigma_vec*sqrt(sigma_value); % Decompose the covariance matrix into the product of an orthogonal matrix and its transpose
    error_A = (sigma_sef*randn(p, 1));

    h_out1 = 5;
    sigma_out = h_out1*eye(p);
    [sigma_vec,sigma_value] = eig(sigma_out); % Perform eigen decomposition on the covariance matrix
    sigma_sef = sigma_vec*sqrt(sigma_value); % Decompose the covariance matrix into the product of an orthogonal matrix and its transpose
    error_out = (sigma_sef*randn(p, 1));
    mu = 1;
    beta_auxiliary_real1 = repmat((2.5*[ones(p,1)]), floor((1/2)*sample_size), 1);
    beta_auxiliary_real2 = repmat((-2.5*[ones(p,1)]), floor((1/2)*sample_size), 1);
    %生成辅助域零的样本X_Auxiliary0
    rho = 0.6; % autocorrelation coefficient of covariance
    sigma = zeros(p,p)*nan; % create an empty covariance matrix
    % create an covariance matrix with autocorrelation
    for i = 1 : p
        for j = 1 : p
            sigma(i, j) = rho^(abs(i-j));
            %sigma(i, j) = 1;
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
    %生成辅助域一的样本Y_Auxiliary0
    beta_Auxiliary_real = [beta_auxiliary_real1; beta_auxiliary_real2];
    epsilon = 0.2 * randn(sample_size, simulation_size); % create the error epsilon
    Y_Auxiliary = zeros(sample_size, simulation_size)*nan; % create the initial response variable
    for i = 1: simulation_size
        for j = 1: sample_size
            Y_Auxiliary(j, i) = X_Auxiliary(j,:,i)*beta_Auxiliary_real(((j-1)*p+1):j*p, 1) + epsilon(j, i);
        end
    end
    % Extract the data needed for the i-th simulation
    yangbencishu = 1;
    y = Y_Auxiliary(:,yangbencishu); %(n*1)
    x = X_Auxiliary(:,:,yangbencishu); % (n*p)
    X_target = X_Auxiliary;
    Y_target = Y_Auxiliary;
    beta_target_real = beta_Auxiliary_real;
elseif which_index == 2
    % 辅助域一参数设置
    rng(3,'philox'); % setting the seed
    p = 20; % the dimension of coefficient
    simulation_size = 100; % the number of repeat 
    sample_size = 500; % the sample size
    h = h/p;
    sigma_A = h*eye(p);
    [sigma_vec,sigma_value] = eig(sigma_A); % Perform eigen decomposition on the covariance matrix
    sigma_sef = sigma_vec*sqrt(sigma_value); % Decompose the covariance matrix into the product of an orthogonal matrix and its transpose
    error_A = (sigma_sef*randn(p, 1));

    h_out1 = 5;
    sigma_out = h_out1*eye(p);
    [sigma_vec,sigma_value] = eig(sigma_out); % Perform eigen decomposition on the covariance matrix
    sigma_sef = sigma_vec*sqrt(sigma_value); % Decompose the covariance matrix into the product of an orthogonal matrix and its transpose
    error_out = (sigma_sef*randn(p, 1));
    mu = 1;
    beta_auxiliary_real1 = repmat(([mu*ones(p,1)]+error_A), floor((1/3)*sample_size), 1);
    beta_auxiliary_real2 = repmat(([-mu*ones(p,1)]+error_A), (sample_size-2*floor((1/3)*sample_size)), 1);
    beta_auxiliary_real3 = repmat((2.5*[ones(p,1)]), floor((1/3)*sample_size), 1);

    % beta_auxiliary_real1 = repmat(([mu*ones(p,1)]+error_A), floor((1/3)*sample_size), 1);
    % beta_auxiliary_real2 = repmat(([-0.5*mu*ones(0.4*p,1);0.5*mu*ones(0.6*p,1)]+error_A), (sample_size-2*floor((1/3)*sample_size)), 1);
    % beta_auxiliary_real3 = repmat((5*[ones(p,1)]), floor((1/3)*sample_size), 1);
    %生成辅助域一的样本X_Auxiliary1
    rho = 0.6; % autocorrelation coefficient of covariance
    sigma = zeros(p,p)*nan; % create an empty covariance matrix
    % create an covariance matrix with autocorrelation
    for i = 1 : p
        for j = 1 : p
            sigma(i, j) = rho^(abs(i-j));
            %sigma(i, j) = 1;
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

    %生成辅助域一的样本Y_Auxiliary1
    beta_Auxiliary_real = [beta_auxiliary_real1; beta_auxiliary_real2; beta_auxiliary_real3];
    epsilon = 0.2 * randn(sample_size, simulation_size); % create the error epsilon
    Y_Auxiliary = zeros(sample_size, simulation_size)*nan; % create the initial response variable
    for i = 1: simulation_size
        for j = 1: sample_size
            Y_Auxiliary(j, i) = X_Auxiliary(j,:,i)*beta_Auxiliary_real(((j-1)*p+1):j*p, 1) + epsilon(j, i);
        end
    end
    % Extract the data needed for the i-th simulation
    yangbencishu = 1;
    y = Y_Auxiliary(:,yangbencishu); %(n*1)
    x = X_Auxiliary(:,:,yangbencishu); % (n*p)
    X_target = X_Auxiliary;
    Y_target = Y_Auxiliary;
    beta_target_real = beta_Auxiliary_real;
elseif which_index == 3
    rng(4,'philox'); % setting the seed
    p = 20; % the dimension of coefficient
    simulation_size = 100; % the number of repeat
    sample_size = 500; % the sample size
    h = h/p;
    sigma_A = h*eye(p);
    [sigma_vec,sigma_value] = eig(sigma_A); % Perform eigen decomposition on the covariance matrix
    sigma_sef = sigma_vec*sqrt(sigma_value); % Decompose the covariance matrix into the product of an orthogonal matrix and its transpose
    error_A = (sigma_sef*randn(p, 1));

    h_out2 = 6;
    sigma_out = h_out2*eye(p);
    [sigma_vec,sigma_value] = eig(sigma_out); % Perform eigen decomposition on the covariance matrix
    sigma_sef = sigma_vec*sqrt(sigma_value); % Decompose the covariance matrix into the product of an orthogonal matrix and its transpose
    error_out = (sigma_sef*randn(p, 1));
    mu = 1;

    % beta_target_real1 = repmat([mu*ones(8,1);zeros(12,1)], floor((1/3)*sample_size), 1);
    % beta_target_real2 = repmat([-mu*ones(4,1);mu*ones(4,1);zeros(12,1)], (sample_size-2*floor((1/3)*sample_size)), 1);
    % beta_target_real3 = repmat([-mu*ones(8,1);zeros(12,1)], floor((1/3)*sample_size), 1);

    beta_auxiliary_real1 = repmat(([-0.5*mu*ones(0.5*p,1);0.5*mu*ones(0.5*p,1)]+error_A), floor((1/3)*sample_size), 1);
    beta_auxiliary_real2 = repmat(([-mu*ones(p,1)]+error_A), (sample_size-2*floor((1/3)*sample_size)), 1);
    beta_auxiliary_real3 = repmat((2.5*[ones(p,1)]), floor((1/3)*sample_size), 1);

    %生成辅助域一的样本X_Auxiliary2
    rho = 0.6; % autocorrelation coefficient of covariance
    sigma = zeros(p,p)*nan; % create an empty covariance matrix
    % create an covariance matrix with autocorrelation
    for i = 1 : p
        for j = 1 : p
            sigma(i, j) = rho^(abs(i-j));
            %sigma(i, j) = 1;
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

    %生成辅助域二的样本Y_Auxiliary2
    beta_Auxiliary_real = [beta_auxiliary_real1; beta_auxiliary_real2; beta_auxiliary_real3];
    epsilon = 0.2 * randn(sample_size, simulation_size); % create the error epsilon
    Y_Auxiliary = zeros(sample_size, simulation_size)*nan; % create the initial response variable
    for i = 1: simulation_size
        for j = 1: sample_size
            Y_Auxiliary(j, i) = X_Auxiliary(j,:,i)*beta_Auxiliary_real(((j-1)*p+1):j*p, 1) + epsilon(j, i);
        end
    end
    % Extract the data needed for the i-th simulation
    yangbencishu = 1;
    y = Y_Auxiliary(:,yangbencishu); %(n*1)
    x = X_Auxiliary(:,:,yangbencishu); % (n*p)
    X_target = X_Auxiliary;
    Y_target = Y_Auxiliary;
    beta_target_real = beta_Auxiliary_real;
elseif which_index == 4
    % 辅助域三参数设置
    rng(8,'philox'); % setting the seed
    p = 20; % the dimension of coefficient
    simulation_size = 100; % the number of repeat
    sample_size = 500; % the sample size
    h = h/p;
    sigma_B = h*eye(p);
    [sigma_vec,sigma_value] = eig(sigma_B); % Perform eigen decomposition on the covariance matrix
    sigma_sef = sigma_vec*sqrt(sigma_value); % Decompose the covariance matrix into the product of an orthogonal matrix and its transpose
    error_B = (sigma_sef*randn(p, 1));

    h_out3 = 7;
    sigma_out = h_out3*eye(p);
    [sigma_vec,sigma_value] = eig(sigma_out); % Perform eigen decomposition on the covariance matrix
    sigma_sef = sigma_vec*sqrt(sigma_value); % Decompose the covariance matrix into the product of an orthogonal matrix and its transpose
    error_out3_1 = (sigma_sef*randn(p, 1));
    error_out3_2 = (sigma_sef*randn(p, 1));
    mu = 1;
    % beta_target_real1 = repmat([mu*ones(8,1);zeros(12,1)], floor((1/3)*sample_size), 1);
    % beta_target_real2 = repmat([-mu*ones(4,1);mu*ones(4,1);zeros(12,1)], (sample_size-2*floor((1/3)*sample_size)), 1);
    % beta_target_real3 = repmat([-mu*ones(8,1);zeros(12,1)], floor((1/3)*sample_size), 1);

    beta_auxiliary_real1 = repmat(([mu*ones(p,1)]+error_B), floor((1/3)*sample_size), 1);
    beta_auxiliary_real2 = repmat(([-0.5*mu*ones(0.5*p,1);0.5*mu*ones(0.5*p,1)]+error_B), (sample_size-2*floor((1/3)*sample_size)), 1);
    beta_auxiliary_real3 = repmat((-2.5*[ones(p,1)]), floor((1/3)*sample_size), 1);
    %生成辅助域三的样本X_Auxiliary3
    rho = 0.6; % autocorrelation coefficient of covariance
    sigma = zeros(p,p)*nan; % create an empty covariance matrix
    % create an covariance matrix with autocorrelation
    for i = 1 : p
        for j = 1 : p
            sigma(i, j) = rho^(abs(i-j));
            %sigma(i, j) = 1;
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

    %生成辅助域三的样本Y_Auxiliary3
    beta_Auxiliary_real = [beta_auxiliary_real1; beta_auxiliary_real2; beta_auxiliary_real3];
    epsilon = 0.2 * randn(sample_size, simulation_size); % create the error epsilon
    Y_Auxiliary = zeros(sample_size, simulation_size)*nan; % create the initial response variable
    for i = 1: simulation_size
        for j = 1: sample_size
            Y_Auxiliary(j, i) = X_Auxiliary(j,:,i)*beta_Auxiliary_real(((j-1)*p+1):j*p, 1) + epsilon(j, i);
        end
    end
    % Extract the data needed for the i-th simulation
    yangbencishu = 1;
    y = Y_Auxiliary(:,yangbencishu); %(n*1)
    x = X_Auxiliary(:,:,yangbencishu); % (n*p)
    X_target = X_Auxiliary;
    Y_target = Y_Auxiliary;
    beta_target_real = beta_Auxiliary_real;
end
% elseif which_index == 5
%     % 辅助域四参数设置
%     rng(6,'philox'); % setting the seed
%     p = 20; % the dimension of coefficient
%     simulation_size = 100; % the number of repeat
%     sample_size = 500; % the sample size
%     h = h/p;
%     sigma_B = h*eye(p);
%     [sigma_vec,sigma_value] = eig(sigma_B); % Perform eigen decomposition on the covariance matrix
%     sigma_sef = sigma_vec*sqrt(sigma_value); % Decompose the covariance matrix into the product of an orthogonal matrix and its transpose
%     error_B = (sigma_sef*randn(p, 1));
% 
%     h_out4 = 8;
%     sigma_out = h_out4*eye(p);
%     [sigma_vec,sigma_value] = eig(sigma_out); % Perform eigen decomposition on the covariance matrix
%     sigma_sef = sigma_vec*sqrt(sigma_value); % Decompose the covariance matrix into the product of an orthogonal matrix and its transpose
%     error_out4_1 = (sigma_sef*randn(p, 1));
%     error_out4_2 = (sigma_sef*randn(p, 1));
%     mu = 1;
% 
%     % beta_target_real1 = repmat([mu*ones(8,1);zeros(12,1)], floor((1/3)*sample_size), 1);
%     % beta_target_real2 = repmat([-mu*ones(4,1);mu*ones(4,1);zeros(12,1)], (sample_size-2*floor((1/3)*sample_size)), 1);
%     % beta_target_real3 = repmat([-mu*ones(8,1);zeros(12,1)], floor((1/3)*sample_size), 1);
% 
%     beta_auxiliary_real1 = repmat(([mu*ones(p,1)]+error_B), floor((1/3)*sample_size), 1);
%     beta_auxiliary_real2 = repmat(([-0.5*mu*ones(0.5*p,1);0.5*mu*ones(0.5*p,1)]+error_B), (sample_size-2*floor((1/3)*sample_size)), 1);
%     beta_auxiliary_real3 = repmat(([-mu*ones(p,1)]+error_B), floor((1/3)*sample_size), 1);
%     %生成辅助域四的样本X_Auxiliary4
%     rho = 0.6; % autocorrelation coefficient of covariance
%     sigma = zeros(p,p)*nan; % create an empty covariance matrix
%     % create an covariance matrix with autocorrelation
%     for i = 1 : p
%         for j = 1 : p
%             sigma(i, j) = rho^(abs(i-j));
%             %sigma(i, j) = 1;
%         end
%     end
%     [sigma_vec,sigma_value] = eig(sigma); % Perform eigen decomposition on the covariance matrix
%     sigma_sef = sigma_vec*sqrt(sigma_value); % Decompose the covariance matrix into the product of an orthogonal matrix and its transpose
%     X_Auxiliary = zeros(sample_size, p, simulation_size); % Create matrix x (100 simulations, each with 60 vector samples).
%     for i = 1:simulation_size
%         for j = 1:sample_size
%             X_Auxiliary(j,:,i) = (sigma_sef*randn(p, 1))';
%         end
%     end
% 
%     %生成辅助域四的样本Y_Auxiliary4
%     beta_Auxiliary_real = [beta_auxiliary_real1; beta_auxiliary_real2; beta_auxiliary_real3];
%     epsilon = 0.2 * randn(sample_size, simulation_size); % create the error epsilon
%     Y_Auxiliary = zeros(sample_size, simulation_size)*nan; % create the initial response variable
%     for i = 1: simulation_size
%         for j = 1: sample_size
%             Y_Auxiliary(j, i) = X_Auxiliary(j,:,i)*beta_Auxiliary_real(((j-1)*p+1):j*p, 1) + epsilon(j, i);
%         end
%     end
%     % Extract the data needed for the i-th simulation
%     yangbencishu = 1;
%     y = Y_Auxiliary(:,yangbencishu); %(n*1)
%     x = X_Auxiliary(:,:,yangbencishu); % (n*p)
%     X_target = X_Auxiliary;
%     Y_target = Y_Auxiliary;
%     beta_target_real = beta_Auxiliary_real;
% elseif which_index == 6
%    % 辅助域四参数设置
%     rng(7,'philox'); % setting the seed
%     p = 20; % the dimension of coefficient
%     simulation_size = 100; % the number of repeat
%     sample_size = 500; % the sample size
%     h = h/p;
%     sigma_B = h*eye(p);
%     [sigma_vec,sigma_value] = eig(sigma_B); % Perform eigen decomposition on the covariance matrix
%     sigma_sef = sigma_vec*sqrt(sigma_value); % Decompose the covariance matrix into the product of an orthogonal matrix and its transpose
%     error_B = (sigma_sef*randn(p, 1));
% 
%     h_out4 = 8;
%     sigma_out = h_out4*eye(p);
%     [sigma_vec,sigma_value] = eig(sigma_out); % Perform eigen decomposition on the covariance matrix
%     sigma_sef = sigma_vec*sqrt(sigma_value); % Decompose the covariance matrix into the product of an orthogonal matrix and its transpose
%     error_out4_1 = (sigma_sef*randn(p, 1));
%     error_out4_2 = (sigma_sef*randn(p, 1));
%     mu = 1;
% 
%     % beta_target_real1 = repmat([mu*ones(8,1);zeros(12,1)], floor((1/3)*sample_size), 1);
%     % beta_target_real2 = repmat([-mu*ones(4,1);mu*ones(4,1);zeros(12,1)], (sample_size-2*floor((1/3)*sample_size)), 1);
%     % beta_target_real3 = repmat([-mu*ones(8,1);zeros(12,1)], floor((1/3)*sample_size), 1);
% 
%     beta_auxiliary_real1 = repmat(([mu*ones(p,1)]+error_B), floor((1/3)*sample_size), 1);
%     beta_auxiliary_real2 = repmat(([-0.5*mu*ones(0.5*p,1);0.5*mu*ones(0.5*p,1)]+error_B), (sample_size-2*floor((1/3)*sample_size)), 1);
%     beta_auxiliary_real3 = repmat(([-mu*ones(p,1)]+error_B), floor((1/3)*sample_size), 1);
%     %生成辅助域四的样本X_Auxiliary4
%     rho = 0.6; % autocorrelation coefficient of covariance
%     sigma = zeros(p,p)*nan; % create an empty covariance matrix
%     % create an covariance matrix with autocorrelation
%     for i = 1 : p
%         for j = 1 : p
%             sigma(i, j) = rho^(abs(i-j));
%             %sigma(i, j) = 1;
%         end
%     end
%     [sigma_vec,sigma_value] = eig(sigma); % Perform eigen decomposition on the covariance matrix
%     sigma_sef = sigma_vec*sqrt(sigma_value); % Decompose the covariance matrix into the product of an orthogonal matrix and its transpose
%     X_Auxiliary = zeros(sample_size, p, simulation_size); % Create matrix x (100 simulations, each with 60 vector samples).
%     for i = 1:simulation_size
%         for j = 1:sample_size
%             X_Auxiliary(j,:,i) = (sigma_sef*randn(p, 1))';
%         end
%     end
% 
%     %生成辅助域四的样本Y_Auxiliary4
%     beta_Auxiliary_real = [beta_auxiliary_real1; beta_auxiliary_real2; beta_auxiliary_real3];
%     epsilon = 0.2 * randn(sample_size, simulation_size); % create the error epsilon
%     Y_Auxiliary = zeros(sample_size, simulation_size)*nan; % create the initial response variable
%     for i = 1: simulation_size
%         for j = 1: sample_size
%             Y_Auxiliary(j, i) = X_Auxiliary(j,:,i)*beta_Auxiliary_real(((j-1)*p+1):j*p, 1) + epsilon(j, i);
%         end
%     end
%     % Extract the data needed for the i-th simulation
%     yangbencishu = 1;
%     y = Y_Auxiliary(:,yangbencishu); %(n*1)
%     x = X_Auxiliary(:,:,yangbencishu); % (n*p)
%     X_target = X_Auxiliary;
%     Y_target = Y_Auxiliary;
%     beta_target_real = beta_Auxiliary_real;
% end
end