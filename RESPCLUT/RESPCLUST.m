function [Result_table, Result_table_summary] = RESPCLUST(simulation_size, Data_x, y_real, sample_size, p, beta_real, k_min, k_max, num_time)
%-------------------------------------函数功能----------------------------------------%
% 根据响应变量对样本进行k_means聚类，聚类个数由轮廓系数来进行选择，对聚类后的模型进行
% 矩阵变量线性回归，比较sc值，均方误差等参数
%-------------------------------------使用方法----------------------------------------%
% [Result_table, Result_table_summary] = RESPCLUST(simulation_size, X_target, Y_target,...
% sample_size, p, beta_target_real, [], [], [])
%----------------------------------需要的前置函数-------------------------------------%
% Matlab Function:         respclust_single_regression     
%                          respclust_k_means_sc
%-----------------------------------输出变量说明--------------------------------------%
% Result_table:            每次simulation结果的表格
% Result_table_summary:    所有simulation的评价结果取均值或标准差的表格汇总
%-----------------------------------输入变量说明--------------------------------------%
% simulation_size:         模拟重复的次数
% Data_x:                  矩阵数据(row_size*col_size*sample_size*simulation_size)
% y_real:                  连续响应变量
% sample_size:             样本量
% col_size:                矩阵变量的列数
% row_size:                矩阵变量的行数
% beta_real:               真实的行回归系数
% gamma_real:              真实的列回归系数
% eps_initial:             矩阵回归模型的对数似然的收敛容差
% iter_max_initial_in      矩阵回归模型迭代的最大步数
% k_min                    kmeans聚类搜索个数的最小类数
% k_max                    kmeans聚类搜索个数的最大类数
% num_time                 kmeans聚类搜索个数计算sc随机重复的最大个数
% 解决多次模拟和单次模拟的维度不一致
if simulation_size == 1
    Data_x = reshape(Data_x, [sample_size,p,1]);
    y_real = reshape(y_real, [sample_size, 1]);
end
% 设置默认参数
if isempty(beta_real)
    beta_real = zeros(sample_size*row_size,1);
end
if isempty(k_min)
    k_min = 1;
end
if isempty(k_max)
    k_max = 10;
end
if isempty(num_time)
    num_time = 10;
end
% 设置存储的结构体Result
Result = struct;
Result(simulation_size).per_if = [];
Result(simulation_size).sc = [];
Result(simulation_size).k = [];
Result(simulation_size).mse_beta = [];
Result(simulation_size).beta_back = [];
% 抽出第i次模拟的样本
for i = 1:simulation_size
    y = y_real(:, i);
    x = Data_x(:,:,i);
    [result] = respclust_single_regression(x, y, sample_size, p, beta_real, k_min, k_max, num_time);
    Result(i).per_if = result.per_if;
    Result(i).sc = result.sc;
    Result(i).k = result.k;
    Result(i).mse_beta = result.mse_beta;
    Result(i).beta_back = result.beta_back;
    fprintf("已完成基于响应变量聚类的异质矩阵回归模型第%.1f %%\n",(i/simulation_size)*100)
end
Result_table = struct2table(Result,"AsArray",true);
per = sum(Result_table.per_if)/simulation_size;
K_mean = mean(Result_table.k);
K_sd = std(Result_table.k);
SC_mean = mean(Result_table.sc);
SC_sd  = std(Result_table.sc);
MSE_beta_mean = mean(Result_table.mse_beta);
MSE_beta_sd = std(Result_table.mse_beta);
Result_table_summary = [K_mean,K_sd,per,SC_mean,SC_sd,MSE_beta_mean,MSE_beta_sd];
Result_table_summary = array2table(Result_table_summary, 'VariableNames', {'K_mean','K_sd','per','SC_mean','SC_sd','MSE_beta_mean','MSE_beta_sd'});
disp(Result_table);
disp(Result_table_summary);
%save('D:\MATLAB_Document2\ADMM_BCD_new\多方法结果\RESPCLUST\非平衡设计\mu=1\B2\Result_table',"Result_table")
%save('D:\MATLAB_Document2\ADMM_BCD_new\多方法结果\RESPCLUST\非平衡设计\mu=1\B2\Result_table_summary',"Result_table_summary")
end















