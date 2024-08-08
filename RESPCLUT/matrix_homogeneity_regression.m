function [beta_back] = matrix_homogeneity_regression(x, y, partitions, sample_size, p)
num_partitions = length(partitions);
coefficients = zeros(p,num_partitions);
beta_back = zeros(sample_size*p,1);
beta_back0 = zeros(sample_size*p,1);
residual_coefficients = 1;
%初始值随机化机制
beta0 = rand(p,1);
for i = 1:num_partitions
   coefficients(:,i) = beta0; 
end
coefficients0 = coefficients;
%对每个分组的数据进行矩阵回归
for i = 1: num_partitions
    iter_initial_in = 0;
    n_part = length(partitions{i});
    x_part = x(partitions{i},:)-mean(x(partitions{i},:));
    y_part = y(partitions{i});
    beta0_part = coefficients0((1:p),i);
    beta_part = (x_part'*x_part)\(x_part'*y_part);
    coefficients((1:p),i) = beta_part;
end
for j = 1:num_partitions
   for k = 1:length(partitions{j})
       beta_back((((partitions{j}(k)-1)*p)+1):(partitions{j}(k)*p),:) = coefficients((1:p),j);
   end
end
end