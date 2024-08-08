%%

for which_index = 0:2 %%%%%%%%%%%%%%%%%%%%%%%可修改 0:辅助域的个数，0是目标域本身
%-----------------------------辅助域参数估计-------------------------------------%
[simulation_size, X_target, Y_target, sample_size, p, beta_target_real] = Data_generation_real_data(which_index);
Aux = Auxiliary_heterogeneity_regression();
if which_index == 0
   num_partitions_list = 3:4; %2:5
else
   num_partitions_list = 5:10; %51:100
end
high_if = 1;
random_number_list = 51:100;
[Results,results,Results_opt,results_opt,...
   initial,Results_list_opt,Results_single_opt,Class_summary] =  Aux.Auxiliary_heterogeneous_regression(simulation_size,[],X_target,Y_target,sample_size,p,beta_target_real,...
   num_partitions_list,[],[],[],[],random_number_list,[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],high_if,[],[],[],[],[],[],[],[],[],[]);
save(fullfile('C:\Users\张丰川\Documents\MATLAB\Transfer_fusion_regression\real_data\Result\auxiliary', ['results_' num2str(which_index) '.mat']), 'Results', 'results',...
    'Results_opt', 'results_opt', 'initial', 'Results_list_opt', 'Results_single_opt', 'Class_summary');
end







%%
%-----------------------------利用辅助域参数进行迁移估计-------------------------------------%
multi_if = 0; %%%%%%%%%%%%%%%%%%%%%%%可修改
which_index = 2;%%%%%%%%%%%%%%%%%%%%%%%辅助域的个数
num_partitions_list = 4; %5
random_number_list = 51:100; %51:100
high_if = 1;
[simulation_size, X_target, Y_target, sample_size, p, beta_target_real] = Data_generation_real_data(0);
 Trans = Transfer_heterogeneity_regression_pro();
 Auxiliary_dateset_number = which_index;
 Auxiliary = struct;
 Auxiliary(Auxiliary_dateset_number).position_provide = [];
for i = 1:which_index
% 本地地址
 Auxiliary(i).position_provide = 'C:\Users\张丰川\Documents\MATLAB\Transfer_fusion_regression\real_data\Result\auxiliary';
end
 [Results,results,Results_opt,results_opt,...
     initial,Results_list_opt,Results_single_opt,Class_summary] = Trans.Transfer_heterogeneous_regression_pro(Auxiliary_dateset_number,...
     Auxiliary,multi_if,simulation_size,[],X_target,Y_target,sample_size,p,beta_target_real,num_partitions_list,[],...
     [],[],[],random_number_list,[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],high_if,[],[],[],[],[],[],[],[],[],[],[]);
%----------------------------------------------------------------------------------------------------------------%
save(fullfile('C:\Users\张丰川\Documents\MATLAB\Transfer_fusion_regression\real_data\Result\transfer',...
    ['results_' num2str(multi_if) '.mat']), 'Results', 'results', 'Results_opt', 'results_opt', 'initial', 'Results_list_opt', 'Results_single_opt', 'Class_summary');







%%
%-----------------------------利用辅助域参数进行迁移估计-------------------------------------%
multi_if = 1; %%%%%%%%%%%%%%%%%%%%%%%可修改
which_index = 2; %%%%%%%%%%%%%%%%%%%%%%%辅助域的个数
[simulation_size, X_target, Y_target, sample_size, p, beta_target_real] = Data_generation_real_data(0);
 Trans = Transfer_heterogeneity_regression_pro();
 Auxiliary_dateset_number = which_index;
 Auxiliary = struct;
 Auxiliary(Auxiliary_dateset_number).position_provide = [];
for i = 1:which_index
% 本地地址
 Auxiliary(i).position_provide = 'C:\Users\张丰川\Documents\MATLAB\Transfer_fusion_regression\real_data\Result\auxiliary';
end
 [Results,results,Results_opt,results_opt,...
     initial,Results_list_opt,Results_single_opt,Class_summary] = Trans.Transfer_heterogeneous_regression_pro(Auxiliary_dateset_number,...
     Auxiliary,multi_if,simulation_size,[],X_target,Y_target,sample_size,p,beta_target_real,[],[],...
     [],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[]);
%----------------------------------------------------------------------------------------------------------------%
save(fullfile('C:\Users\张丰川\Documents\MATLAB\Transfer_fusion_regression\real_data\Result\transfer',...
    ['results_' num2str(multi_if) '.mat']), 'Results', 'results', 'Results_opt', 'results_opt', 'initial', 'Results_list_opt', 'Results_single_opt', 'Class_summary');




%%
% 提取X和y数据
X = Class_summary.X{1};
y = Class_summary.y{1};

% 拟合线性回归模型
model = fitlm(X, y);

% 显示回归结果
disp(model);

% 提取回归系数和p值
coefficients = model.Coefficients.Estimate;
pValues = model.Coefficients.pValue;

% 显示回归系数和p值
disp('Coefficients:');
disp(coefficients);
disp('p-values:');
disp(pValues);

%%
%R^2 同质时
Y = readtable("C:/Users/张丰川/Documents/MATLAB/Transfer_fusion_regression/real_data/data/Data_target.csv").Y;
X = readtable("C:/Users/张丰川/Documents/MATLAB/Transfer_fusion_regression/real_data/data/Data_target.csv");
X_target = table2array(X(:,{'Data_age', 'Data_sex', 'Data_site_detail', 'Data_summary_stage', 'Data_surg',...
                        'Data_surg_rad_seq', 'Data_radiation', 'Data_chemotherapy', 'Data_systemic_surg_seq',...
                        'Data_first_if', 'Data_total_number', 'Data_marital_status', 'Data_median_income'}));
%X2_target = X(:,{'Data_ER','Data_PR','Data_Her2','Data_grade_clinical','Data_grade_summary','Data_grade_pathological',...
%        'Data_type','Data_subtype','Data_hist_behavior'});
model = fitlm(X_target, Y);
beta_f = model.Coefficients.Estimate;
r_2_f = 1-(mean((Y-(X_target*beta_f((2:end),1))).^2)/mean((Y-mean(Y)).^2))
%%
%R^2 迁移前
load('C:\Users\张丰川\Documents\MATLAB\Transfer_fusion_regression\real_data\Result\auxiliary\results_0.mat')
beta_s1 = Class_summary.coef{1, 1};
beta_s2 = Class_summary.coef{1, 2};
X1 = Class_summary.X{1, 1};
X2 = Class_summary.X{1, 2};
Y1 = Class_summary.y{1, 1};
Y2 = Class_summary.y{1, 2};
%model = fitlm(X1, Y1);
%model = fitlm(X2, Y2);
beta_s1 = inv(X1'*X1)*(X1'*Y1);
beta_s2 = inv(X2'*X2)*(X2'*Y2);
r_2_s = 1-(sum((Y1-(X1*beta_s1)).^2) + sum((Y2-(X2*beta_s2)).^2))/(sum((Y1-mean(Y1)).^2) + sum((Y2-mean(Y2)).^2))
r1 = 1-sum((Y1-(X1*beta_s1)).^2)/sum((Y1-mean(Y1)).^2); 
r2 = 1-sum((Y2-(X2*beta_s2)).^2)/sum((Y2-mean(Y2)).^2);

%%
%R^2 迁移后 single
load('C:\Users\张丰川\Documents\MATLAB\Transfer_fusion_regression\real_data\Result\transfer\results_0.mat')
beta_s1 = Class_summary.coef{1, 1};
beta_s2 = Class_summary.coef{1, 2};
beta_s3 = Class_summary.coef{1, 3};
X1 = Class_summary.X{1, 1};
X2 = Class_summary.X{1, 2};
X3 = Class_summary.X{1, 3};
Y1 = Class_summary.y{1, 1};
Y2 = Class_summary.y{1, 2};
Y3 = Class_summary.y{1, 3};
%model = fitlm(X1, Y1);
%model = fitlm(X2, Y2);
beta_s1 = inv(X1'*X1)*(X1'*Y1);
beta_s2 = inv(X2'*X2)*(X2'*Y2);
beta_s3 = inv(X3'*X3)*(X3'*Y3);
r_2_s = 1-(sum((Y1-(X1*beta_s1)).^2) + sum((Y2-(X2*beta_s2)).^2) + sum((Y3-(X3*beta_s3)).^2))/(sum((Y1-mean(Y1)).^2) + sum((Y2-mean(Y2)).^2) + sum((Y3-mean(Y3)).^2))
r1 = 1-sum((Y1-(X1*beta_s1)).^2)/sum((Y1-mean(Y1)).^2); 
r2 = 1-sum((Y2-(X2*beta_s2)).^2)/sum((Y2-mean(Y2)).^2);
r3 = 1-sum((Y3-(X3*beta_s3)).^2)/sum((Y3-mean(Y3)).^2);
%%
%R^2 迁移后 multiple
load('C:\Users\张丰川\Documents\MATLAB\Transfer_fusion_regression\real_data\Result\transfer\results_0.mat')
beta_s1 = Class_summary.coef{1, 1};
beta_s2 = Class_summary.coef{1, 2};
X1 = Class_summary.X{1, 1};
X2 = Class_summary.X{1, 2};
Y1 = Class_summary.y{1, 1};
Y2 = Class_summary.y{1, 2};
r_2_s = 1-(sum((Y1-(X1*beta_s1)).^2) + sum((Y2-(X2*beta_s2)).^2))/(sum((Y1-mean(Y1)).^2) + sum((Y2-mean(Y2)).^2))
%%










