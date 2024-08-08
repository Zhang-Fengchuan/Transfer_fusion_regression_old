function [result] = residual_single_regression(x, y, sample_size, p, beta_real, k_min, k_max, num_time)
result = struct;
result.per_if = [];
result.sc = [];
result.k = [];
result.mse_beta = [];
result.beta_back = [];
% 首先利用所有样本进行一次矩阵回归得到残差
index_store = cell(1, 1);
index_store{1,1} = 1:sample_size;
[beta_back] = matrix_homogeneity_regression(x, y, index_store, sample_size, p);
%one_vec_row = ones(p,1);
% one_beta = kron(eye(sample_size),one_vec_row');
% D_beta = one_beta*diag(beta_back);
x_cell = mat2cell(x,ones(1,sample_size),p);
X = blkdiag(x_cell{:});
% residual = y-D_beta*X*gamma_back;
residual = y-X*beta_back;
% 获得单次simulation的k_means的分组数和分类情况
[class, k_best] = residual_k_means_sc(residual, k_min, k_max, num_time, sample_size);
index_store = cell(k_best, 1);
for i = 1: k_best
    index_store{i,1} = find(class==i);
end
% 真实参数下的分组情况
real_class = zeros(sample_size, sample_size);
for i = 1: (sample_size - 1)
    for j = (i+1): sample_size
        ifsame = norm([beta_real(((i-1)*p+1):(i*p), 1)]-[beta_real(((j-1)*p+1):(j*p), 1)]);
        if ifsame ~= 0
            real_class(i, j) = 1;
            real_class(j, i) = 1;
        end
    end
end
% 判断基于响应变量k_means聚类的结果与真实系数类数是否相同
[~, Ireal, ~] = unique(real_class,'rows');
k_real_num = length(Ireal);
per_if = 0;
if k_real_num == k_best
    per_if = 1;
end
% 计算一致性指数SC
est_class = ones(sample_size, sample_size);
est_class(sample_size, sample_size) = 0;
for i = 1: (sample_size-1)
    for j = (i+1): sample_size
        for k = 1:k_best
            if any(index_store{k}(:) == i)&&any(index_store{k}(:) == j)
                est_class(i, j) = 0;
                est_class(j, i) = 0;
            end
        end
    end
    est_class(i, i) = 0;
end
count = 0;
for i = 1: (sample_size-1)
    for j = (i+1): sample_size
        if real_class(i,j)==est_class(i,j)
            count = count + 1;
        end
    end
end
sc = count/nchoosek(sample_size,2);
% 对每个分组分别用矩阵回归线性模型进行拟合
[beta_back] = matrix_homogeneity_regression(x, y, index_store, sample_size, p);
% 计算估计参数的MSE
mse_beta = sum((beta_back-beta_real).^2)/(sample_size*p);
result.per_if = per_if;
result.sc = sc;
result.k = k_best;
result.mse_beta = mse_beta;
result.beta_back = beta_back;
end