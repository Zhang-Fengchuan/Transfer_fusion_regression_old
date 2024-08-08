function [simulation_size, X_target, Y_target, sample_size, p, beta_target_real] = Data_generation_real_data(which_index)
if which_index == 0
    Y = readtable("C:/Users/张丰川/Documents/MATLAB/Transfer_fusion_regression/real_data/data/Data_target.csv").Y;
    X = readtable("C:/Users/张丰川/Documents/MATLAB/Transfer_fusion_regression/real_data/data/Data_target.csv");
    % X_target = table2array(X(:,{'Data_age', 'Data_sex', 'Data_behavior', 'Data_site_detail', 'Data_summary_stage',...
        % 'Data_surg', 'Data_surg_rad_seq', 'Data_radiation', 'Data_chemotherapy', 'Data_systemic_surg_seq',...
        % 'Data_time_diag_treat', 'Data_first_if', 'Data_total_number', 'Data_marital_status', 'Data_median_income'}));
    X_target = table2array(X(:,{'Data_surg', 'Data_radiation', 'Data_chemotherapy', 'Data_surg_rad_seq',...
                        'Data_age', 'Data_sex', 'Data_site_detail', 'Data_summary_stage', 'Data_time_diag_treat',...
                        'Data_first_if', 'Data_total_number', 'Data_marital_status', 'Data_median_income'}));
    X2_target = X(:,{'Data_ER','Data_PR','Data_Her2','Data_grade_clinical','Data_grade_summary','Data_grade_pathological',...
        'Data_type','Data_subtype','Data_hist_behavior'});
    Y_target = Y;
    [sample_size, p] = size(X_target);
    simulation_size = 1; % the number of repeat
    mu = 1;
    beta_target_real1 = repmat([mu*ones(p,1)], floor((1/3)*sample_size), 1);
    beta_target_real2 = repmat([-0.5*mu*ones(floor(0.5*p),1);0.5*mu*ones((p-floor(0.5*p)),1)], (sample_size-2*floor((1/3)*sample_size)), 1);
    beta_target_real3 = repmat([-mu*ones(p,1)], floor((1/3)*sample_size), 1);
    beta_target_real = [beta_target_real1; beta_target_real2; beta_target_real3];
elseif which_index == 1
    Y = readtable("C:/Users/张丰川/Documents/MATLAB/Transfer_fusion_regression/real_data/data/Data_1.csv").Y;
    X = readtable("C:/Users/张丰川/Documents/MATLAB/Transfer_fusion_regression/real_data/data/Data_1.csv");
    % X_target = table2array(X(:,{'Data_age', 'Data_sex', 'Data_behavior', 'Data_site_detail', 'Data_summary_stage',...
    %     'Data_surg', 'Data_surg_rad_seq', 'Data_radiation', 'Data_chemotherapy', 'Data_systemic_surg_seq',...
    %     'Data_time_diag_treat', 'Data_first_if', 'Data_total_number', 'Data_marital_status', 'Data_median_income'}));
   X_target = table2array(X(:,{'Data_surg',  'Data_radiation', 'Data_chemotherapy', 'Data_surg_rad_seq',...
                        'Data_age', 'Data_sex', 'Data_site_detail', 'Data_summary_stage', 'Data_time_diag_treat',...
                        'Data_first_if', 'Data_total_number', 'Data_marital_status', 'Data_median_income'}));
    X2_target = X(:,{'Data_ER','Data_PR','Data_Her2','Data_grade_clinical','Data_grade_summary','Data_grade_pathological',...
        'Data_type','Data_subtype','Data_hist_behavior'});
    Y_target = Y;
    [sample_size, p] = size(X_target);
    simulation_size = 1; % the number of repeat
    mu = 1;
    beta_target_real1 = repmat([mu*ones(p,1)], floor((1/3)*sample_size), 1);
    beta_target_real2 = repmat([-0.5*mu*ones(floor(0.5*p),1);0.5*mu*ones((p-floor(0.5*p)),1)], (sample_size-2*floor((1/3)*sample_size)), 1);
    beta_target_real3 = repmat([-mu*ones(p,1)], floor((1/3)*sample_size), 1);
    beta_target_real = [beta_target_real1; beta_target_real2; beta_target_real3];
elseif which_index == 2
    Y = readtable("C:/Users/张丰川/Documents/MATLAB/Transfer_fusion_regression/real_data/data/Data_2.csv").Y;
    X = readtable("C:/Users/张丰川/Documents/MATLAB/Transfer_fusion_regression/real_data/data/Data_2.csv"); 
    % X_target = table2array(X(:,{'Data_age', 'Data_sex', 'Data_behavior', 'Data_site_detail', 'Data_summary_stage',...
    %     'Data_surg', 'Data_surg_rad_seq', 'Data_radiation', 'Data_chemotherapy', 'Data_systemic_surg_seq',...
    %     'Data_time_diag_treat', 'Data_first_if', 'Data_total_number', 'Data_marital_status', 'Data_median_income'}));
    X_target = table2array(X(:,{'Data_surg',  'Data_radiation', 'Data_chemotherapy', 'Data_surg_rad_seq',...
                        'Data_age', 'Data_sex', 'Data_site_detail', 'Data_summary_stage', 'Data_time_diag_treat',...
                        'Data_first_if', 'Data_total_number', 'Data_marital_status', 'Data_median_income'}));
    X2_target = X(:,{'Data_ER','Data_PR','Data_Her2','Data_grade_clinical','Data_grade_summary','Data_grade_pathological',...
        'Data_type','Data_subtype','Data_hist_behavior'});
    Y_target = Y;
    [sample_size, p] = size(X_target);
    simulation_size = 1; % the number of repeat
    mu = 1;
    beta_target_real1 = repmat([mu*ones(p,1)], floor((1/3)*sample_size), 1);
    beta_target_real2 = repmat([-0.5*mu*ones(floor(0.5*p),1);0.5*mu*ones((p-floor(0.5*p)),1)], (sample_size-2*floor((1/3)*sample_size)), 1);
    beta_target_real3 = repmat([-mu*ones(p,1)], floor((1/3)*sample_size), 1);
    beta_target_real = [beta_target_real1; beta_target_real2; beta_target_real3];
elseif which_index == 3
   Y = readtable("C:/Users/张丰川/Documents/MATLAB/Transfer_fusion_regression/real_data/data/Data_3.csv").Y;
    X = readtable("C:/Users/张丰川/Documents/MATLAB/Transfer_fusion_regression/real_data/data/Data_3.csv");
    % X_target = table2array(X(:,{'Data_age', 'Data_sex', 'Data_behavior', 'Data_site_detail', 'Data_summary_stage',...
    %     'Data_surg', 'Data_surg_rad_seq', 'Data_radiation', 'Data_chemotherapy', 'Data_systemic_surg_seq',...
    %     'Data_time_diag_treat', 'Data_first_if', 'Data_total_number', 'Data_marital_status', 'Data_median_income'}));
    X_target = table2array(X(:,{'Data_surg',  'Data_radiation', 'Data_chemotherapy','Data_surg_rad_seq',...
                        'Data_age', 'Data_sex', 'Data_site_detail', 'Data_summary_stage', 'Data_time_diag_treat',...
                        'Data_first_if', 'Data_total_number', 'Data_marital_status', 'Data_median_income'}));
    X2_target = X(:,{'Data_ER','Data_PR','Data_Her2','Data_grade_clinical','Data_grade_summary','Data_grade_pathological',...
        'Data_type','Data_subtype','Data_hist_behavior'});
    Y_target = Y;
    [sample_size, p] = size(X_target);
    simulation_size = 1; % the number of repeat
    mu = 1;
    beta_target_real1 = repmat([mu*ones(p,1)], floor((1/3)*sample_size), 1);
    beta_target_real2 = repmat([-0.5*mu*ones(floor(0.5*p),1);0.5*mu*ones((p-floor(0.5*p)),1)], (sample_size-2*floor((1/3)*sample_size)), 1);
    beta_target_real3 = repmat([-mu*ones(p,1)], floor((1/3)*sample_size), 1);
    beta_target_real = [beta_target_real1; beta_target_real2; beta_target_real3];
end
end