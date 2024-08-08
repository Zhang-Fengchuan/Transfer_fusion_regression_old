classdef Initial_regression
    methods
        % obj = Initial_regression();
        % [Result_table, Result_table_summary] = obj.INITIAL(simulation_size, X_target, Y_target, sample_size, p, beta_target_real, [], [],...
        %                           [],[],[],[],[]);
        function [Result_table, Result_table_summary] = INITIAL(obj,simulation_size, X_target, Y_target, sample_size, p, beta_real,...
                num_partitions_list, eps_out, iter_max_initial_out,...
                combine_principle_list, split_principle_list, random_number_list, disp_if, c,  high_if, CV_number, K_require, K_require_if)
            %-------------------------------------Functionality of the function----------------------------------------------------%
            % The heterogeneous matrix regression model using the initial value search algorithm returns the number of clusters,
            % SC value, mean squared error, and other parameters, as well as the estimated initial values (can be simulated multiple times).
            %------------------------------------------Usage instructions---------------------------------------------------------%
            % [Result_table, Result_table_summary] = INITIAL(simulation_size, X_target, Y_target, sample_size, ...
            %                                p, beta_target_real, [], [], [], [], [], [], [], [], [], [], [], [])
            %-------------------------------------Required prerequisite functions-------------------------------------------------%
            % Matlab Function:         k_means_regression_initial
            %                          k_means_regression_tree_single
            %                          k_means_regression_cv
            %                          k_means_regression
            %-----------------------------------Description of output variables--------------------------------------------------%
            % Result_table:            The table of results for each simulation.
            % Result_table_summary:    The summary table of the mean or standard deviation of evaluation results for all simulations.
            %-----------------------------------Description of input variables.--------------------------------------------------%
            % num_partitions_list       Divide the data into several parts, parameter matrix.
            % iter_max_initial_in       The maximum number of iterations for alternating optimization during each homogeneous model regression.
            % eps_out                   The coefficient difference precision and convergence precision between
            %                           the n-th and (n-1)-th models after sample re-partitioning.
            % iter_max_initial_out      The maximum number of iterations for sample re-partitioning.
            % sample_size               sample size
            % p                         coefficient dimension
            % x                         The explanatory variables in a single simulation. [sample_size*p]
            % y                         The response variables in a single simulation. [sample_size*1]
            % random_number_list        The random seed parameter matrix.
            % combine_principle_list    The parameter matrix for the minimum element count merging criterion.
            % split_principle_list      The parameter matrix for the maximum element count splitting criterion.
            % disp_if                   Whether to output progress information.
            % beta_real:                real coefficient
            % Resolve dimension inconsistencies between multiple simulations and a single simulation.
            % % % if simulation_size == 1
            % % %     Data_x = reshape(Data_x, [row_size,col_size,sample_size,1]);
            % % %     y_real = reshape(y_real, [sample_size, 1]);
            % % % end
            % 设置默认参数
            if isempty(eps_out)
                eps_out = 1e-3;
                %初值算法收敛条件1e-8
            end
            if isempty(iter_max_initial_out)
                iter_max_initial_out = 50;
            end
            if isempty(num_partitions_list)
                %num_partitions_list = min(10, floor(sqrt(sample_size))):10;
                %num_partitions_list = floor(sqrt(sample_size));
                num_partitions_list = 5:10;
            end
            if isempty(combine_principle_list)
                combine_principle_list = linspace(0.1, 0.1, 1)*sample_size;
            end
            if isempty(split_principle_list)
                split_principle_list = linspace(1, 1, 1)*sample_size;
            end
            if isempty(random_number_list)
                random_number_list = 1:50;
            end
            if isempty(disp_if)
                disp_if = 1;
            end
            if isempty(beta_real)
                beta_real = zeros(sample_size*p,1);
            end
            if isempty(high_if)
                high_if = 0;
            end
            if isempty(CV_number)
                CV_number = 5;
            end
            if isempty(c)
                c = 1;
            end
            if isempty(K_require)
                K_require = floor(sqrt(sample_size));
            end
            if isempty(K_require_if)
                K_require_if = 1;
            end
            % if isempty(high_constant)
            %     high_constant = 1;
            % end
            % 设置存储的结构体Result
            Result = struct;
            Result(simulation_size).per_if = [];
            Result(simulation_size).sc = [];
            Result(simulation_size).k = [];
            Result(simulation_size).mse_beta = [];
            Result(simulation_size).beta_back = [];
            Result(simulation_size).convergence = [];
            Result(simulation_size).partitions = [];
            Result(simulation_size).coefficients = [];
            Result(simulation_size).num_partitions_best = [];
            Result(simulation_size).combine_principle_best = [];
            Result(simulation_size).split_principle_best = [];
            Result(simulation_size).random_number_best = [];
            % 抽出第i次模拟的样本
            parfor z = 1:simulation_size
                y = Y_target(:, z);
                x = X_target(:,:,z);
                [~, partitions, result] = obj.k_means_regression_initial(num_partitions_list, eps_out,...
                    iter_max_initial_out, combine_principle_list, split_principle_list, random_number_list,...
                    sample_size, p, x, y, disp_if, beta_real, c, high_if, CV_number, K_require, K_require_if);
                Result(z).per_if = result.per_if;
                Result(z).sc = result.sc;
                Result(z).k = result.k;
                Result(z).mse_beta = result.mse_beta;
                Result(z).beta_back = result.beta_back;
                Result(z).convergence = result.convergence;
                Result(z).partitions = partitions;
                Result(z).coefficients = result.coefficients;
                Result(z).num_partitions_best = result.num_partitions_best;
                Result(z).combine_principle_best = result.combine_principle_best;
                Result(z).split_principle_best = result.split_principle_best;
                Result(z).random_number_best = result.random_number_best;
                fprintf("已完成初值搜索算法的异质回归模型第%.1f %%\n",(z/simulation_size)*100)
            end
            Result_table = struct2table(Result);
            per = sum(Result_table.per_if)/simulation_size;
            K_mean = mean(Result_table.k);
            K_sd = std(Result_table.k);
            SC_mean = mean(Result_table.sc);
            SC_sd  = std(Result_table.sc);
            MSE_beta_mean = mean(Result_table.mse_beta);
            MSE_beta_sd = std(Result_table.mse_beta);
            convergence_mean = mean(Result_table.convergence);
            convergence_sd = std(Result_table.convergence);
            Result_table_summary = [K_mean,K_sd,per,SC_mean,SC_sd,MSE_beta_mean,MSE_beta_sd,convergence_mean,convergence_sd];
            Result_table_summary = array2table(Result_table_summary, 'VariableNames', {'K_mean','K_sd','per','SC_mean','SC_sd','MSE_beta_mean','MSE_beta_sd',...
                'convergence_mean','convergence_sd'});
            disp(Result_table);
            disp(Result_table_summary);
            save('D:\Transfer-learning-heterogeneity-analysis\Result\Initial_target\Result_table',"Result_table")
            save('D:\Transfer-learning-heterogeneity-analysis\Result\Initial_target\Result_table_summary',"Result_table_summary")
        end


























        function [beta0, partitions, result] = k_means_regression_initial(obj, num_partitions_list, eps_out,...
                iter_max_initial_out, combine_principle_list,...
                split_principle_list, random_number_list,...
                sample_size, p, x, y, disp_if, beta_real, c, high_if, CV_number, K_require, K_require_if)
            %-----------------------------------------Functionality of the function-------------------------------------------%
            % The heterogeneous matrix regression model using the initial value search algorithm returns the number of clusters,
            % SC value, mean squared error, and other parameters, as well as the estimated initial values.
            %----------------------------------------------Usage instructions-------------------------------------------------%
            % [beta0, gamma0, partitions, result] = k_means_regression_initial([], [], [],...
            %                                                    [], [], [], sample_size, p,...
            %                                                    x, y, [], [], [], [], [], [], [])
            %-----------------------------------------Required prerequisite functions-----------------------------------------%
            % Matlab Function:         k_means_regression_tree_single
            %                          k_means_regression_cv
            %                          k_means_regression
            %----------------------------------------Description of output variables------------------------------------------%
            % beta0:                   The beta estimates returned by the heterogeneous matrix regression model
            %                          using the initial value search algorithm in a single simulation.
            % partitions:              The grouping results returned by the heterogeneous matrix regression model
            %                          using the initial value search algorithm in a single simulation.
            % result:                  The results structure returned in a single simulation,
            %                          including k, SC, optimal CV parameters, parameter estimates, etc.
            %----------------------------------------Description of input variables-------------------------------------------%
            % num_partitions_list       Divide the data into several parts, parameter matrix.
            % eps_initial               The convergence precision of each model during each homogeneous model regression.
            % eps_out                   The coefficient difference precision and convergence precision between
            %                           the n-th and (n-1)-th models after sample re-partitioning.
            % iter_max_initial_out      The maximum number of iterations for sample re-partitioning.
            % sample_size               sample size
            % p                         coefficient dimension
            % x                         The explanatory variables in a single simulation [row_size*col_size*sample_size]
            % y                         The response variables in a single simulation [sample_size*1]
            % random_number_list        The random seed parameter matrix.
            % combine_principle_list    The parameter matrix for the minimum element count merging criterion.
            % split_principle_list      The parameter matrix for the maximum element count splitting criterion.
            % disp_if                   Whether to output progress information.
            % beta_real:                True regression coefficients.
            if isempty(eps_out)
                eps_out = 1e-3;
                %初值算法收敛条件1e-8
            end
            if isempty(iter_max_initial_out)
                iter_max_initial_out = floor(sqrt(sample_size));
            end
            if isempty(num_partitions_list)
                %num_partitions_list = min(10, floor(sqrt(sample_size))):10;
                %num_partitions_list = floor(sqrt(sample_size));
                num_partitions_list = 5:10;
            end
            if isempty(combine_principle_list)
                combine_principle_list = linspace(0.1,0.1,1)*sample_size;
            end
            if isempty(split_principle_list)
                split_principle_list = 1*sample_size;
            end
            if isempty(random_number_list)
                random_number_list = 1:50;
            end
            if isempty(disp_if)
                disp_if = 1;
            end
            if isempty(beta_real)
                beta_real = zeros(sample_size*p,1);
            end
            if isempty(high_if)
                high_if = 0;
            end
            if isempty(CV_number)
                CV_number = 5;
            end
            if isempty(c)
                c = 1;
            end
            if isempty(K_require)
                K_require = floor(sqrt(sample_size));
            end
            if isempty(K_require_if)
                K_require_if = 1;
            end
            result = struct;
            result.per_if = [];
            result.sc = [];
            result.k = [];
            result.mse_beta = [];
            result.beta_back = [];
            result.convergence = [];
            result.coefficients = [];
            result.num_partitions_best = [];
            result.combine_principle_best = [];
            result.split_principle_best = [];
            result.random_number_best = [];
            [num_partitions_best,combine_principle_best,split_principle_best,random_number_best, cv_matrix] = obj.k_means_regression_cv(num_partitions_list,eps_out,...
                iter_max_initial_out,combine_principle_list,split_principle_list,random_number_list,sample_size,p,x,y,disp_if,c,high_if,CV_number,K_require_if);
            [beta0,coefficients,partitions,combine,residual_coefficients] = obj.k_means_regression(num_partitions_best,eps_out,...
                iter_max_initial_out,combine_principle_best,split_principle_best,sample_size,p,x,y,random_number_best, high_if, CV_number,K_require_if);

            % % %
            % % % if K_require_if == 1
            % % %     partitions = cell(1, 1);
            % % %     partitions{1} = 1:sample_size;
            % % %     coefficients = zeros(p, 1);
            % % %     residuals = zeros(sample_size,1);
            % % %     sample_size_sum = sample_size;
            % % %     iter = 0;
            % % %     while((length(partitions) < K_require) && (iter < iter_max_initial_out))
            % % %         num = length(partitions);
            % % %         partitions = [partitions,cell(1,num)];
            % % %         coefficients = [coefficients, zeros(p,num)];
            % % %         residuals = [residuals, zeros(sample_size,num)];
            % % %         random_number_best = random_number_best + 1;
            % % %         num_partitions_best = 10;
            % % %         for i = 1:num
            % % %             partitions_old = {partitions{i}};
            % % %             [beta0,coefficients_new,...
            % % %                 partitions_new,combine,residual_coefficients] = obj.k_means_regression_tree_single(2,eps_out,...
            % % %                 iter_max_initial_out,(combine_principle_best/sample_size_sum),p,x,y,random_number_best,coefficients(:,i),partitions_old,sample_size_sum,high_if,CV_number,0);
            % % %             % [~,coefficients_new,...
            % % %             %     partitions_new,~,~] = k_means_regression(obj,2,eps_out,...
            % % %             %     iter_max_initial_out,combine_principle_best,split_principle_best,length(partitions{i}),p,x,y,random_number_best, high_if, CV_number, K_require_if);
            % % %             partitions{i} = partitions_new{1};
            % % %             partitions{num+i} = partitions_new{2};
            % % %             coefficients(:,i) = coefficients_new(:,1);
            % % %             coefficients(:,(i+num)) = coefficients_new(:,2);
            % % %             for t = 1:sample_size
            % % %                 residuals(t, i) = abs(y(t)-(x(t,:)-mean(x(partitions{i},:)))*coefficients(:,i)-mean(y(partitions{i})));
            % % %                 residuals(t, i+num) = abs(y(t)-(x(t,:)-mean(x(partitions{(i+num)},:)))*coefficients(:,(i+num))-mean(y(partitions{(i+num)})));
            % % %             end
            % % %
            % % %         end
            % % %         combine = zeros(length(partitions),1);
            % % %         for j = 1:length(partitions)
            % % %             combine(j,1) = length(partitions{j});
            % % %         end
            % % %         combine_index = min(combine);
            % % %         combine_min_index = find(combine==min(combine));
            % % %         while (combine_index < combine_principle_best)
            % % %             num_partitions = length(partitions) - length(combine_min_index);
            % % %             residuals(:,combine_min_index) = [];
            % % %             coefficients(:,combine_min_index) = [];
            % % %             partitions_mid = partitions;
            % % %             partitions(combine_min_index) = [];
            % % %             [~, min_col_indices] = min(residuals, [], 2);
            % % %             for k = 1:length(combine_min_index)
            % % %                 for t = 1:length(partitions_mid{combine_min_index(k)})
            % % %                     for j = 1:length(coefficients)
            % % %                         if min_col_indices(partitions_mid{combine_min_index(k)}(t))==j
            % % %                             partitions{j} = [partitions{j},partitions_mid{combine_min_index(k)}(t)];
            % % %                         end
            % % %                     end
            % % %                 end
            % % %             end
            % % %             combine = zeros(length(partitions),1);
            % % %             for j = 1:num_partitions
            % % %                 combine(j,1) = length(partitions{j});
            % % %             end
            % % %             combine_index = min(combine);
            % % %             combine_min_index = find(combine==min(combine));
            % % %         end
            % % %         iter = iter + 1;
            % % %         % partitions
            % % %         % length(partitions)
            % % %     end
            % % % end

            if K_require_if == 1
                while(length(partitions)<K_require)
                    split_max_index = find(combine==max(combine));
                    split_max_index = split_max_index(1);
                    % [beta_back,coefficients,...
                    % partitions,combine,residual_coefficients] = obj.k_means_regression_tree_single(2,eps_out,...
                    % iter_max_initial_out,combine_principle,p,x,y,random_number,beta0,partitions_old,sample_size_sum,high_if,CV_number)
                    coefficients = [coefficients,zeros(p,1)];
                    partitions = [partitions;cell(1,1)];
                    partitions_old = cell(1,1);
                    partitions_old{1} = partitions{split_max_index};
                    sample_size_sum = sample_size;
                    [~,coefficients_split,...
                        partitions_split,~,residual_coefficients] = obj.k_means_regression_tree_single(2,eps_out,...
                        iter_max_initial_out, 0, p,x,y,random_number_best,coefficients(:,split_max_index),partitions_old,sample_size_sum,high_if,CV_number,K_require_if);
                    partitions{split_max_index} = partitions_split{1};
                    partitions{end} = partitions_split{2};
                    coefficients(:,split_max_index) = coefficients_split(:,1);
                    coefficients(:,end) = coefficients_split(:,2);
                    combine = zeros(length(partitions),1);
                    for j = 1:length(partitions)
                        combine(j,1) = length(partitions{j});
                    end
                end
                for j = 1:length(partitions)
                    for k = 1:length(partitions{j})
                        beta0((((partitions{j}(k)-1)*p)+1):(partitions{j}(k)*p),:) = coefficients((1:p),j);
                    end
                end
            end

            % 获得单次simulation的k_means的分组数和分类情况
            k_best = length(partitions);
            index_store = cell(k_best, 1);
            for i = 1: k_best
                index_store{i,1} = partitions{i}';
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

            % 计算估计参数的MSE
            mse_beta = sum((beta0-beta_real).^2)/(sample_size*p);
            result.per_if = per_if;
            result.sc = sc;
            result.k = k_best;
            result.mse_beta = mse_beta;
            result.beta_back = beta0;
            result.convergence = residual_coefficients;
            result.coefficients = coefficients;
            result.num_partitions_best = num_partitions_best;
            result.combine_principle_best = combine_principle_best/sample_size;
            result.split_principle_best = split_principle_best/sample_size;
            result.random_number_best = random_number_best;
        end

































        function [beta_back,coefficients,...
                partitions,combine,residual_coefficients] = k_means_regression(obj,num_partitions,eps_out,...
                iter_max_initial_out,combine_principle,split_principle,sample_size,p,x,y,random_number, high_if, CV_number, K_require_if)
            %------------------------------------------Functionality of the function---------------------------------------------%
            % (Residual Iteration Model) Solve the initial value results of subgroup analysis
            %  using the residual iteration algorithm given the input parameters.

            %-----------------------------------------Description of output variables---------------------------------------------%
            % beta_back             Row vector values returned by the residual iteration model. [p*sample_size] matrix
            % coefficients          The regression coefficients of the two models returned by the residual iteration model.
            %                       [(row_size+col_size)*The number of groups automatically determined by the model] matrix
            % partitions            The group indices of the two models returned by the residual iteration model (relative to the original x and y group indices).
            %                       [The number of groups automatically determined by the model*1] cell
            % combine               The sample sizes corresponding to the two models returned by the residual iteration model.
            %                       [The number of groups automatically determined by the model*1] matrix
            % residual_coefficients The coefficient difference between the n-th and (n-1)-th models after sample re-partitioning.
            %                       [1*1] double
            %-----------------------------------------Description of input variables---------------------------------------------%
            % num_partitions        Divide the data into several parts.
            % eps_out               The coefficient difference precision between the n-th and (n-1)-th model.
            %                       after sample re-partitioning, and the convergence precision.
            % iter_max_initial_out  The maximum number of iterations for sample re-partitioning.
            % sample_size           sample size.
            % p                     coefficient dimension.
            % x                     The explanatory variables in a single simulation. [sample_size*p]
            % y                     The response variables in a single simulation. [sample_size*1]
            % random_number         Random seed number to ensure reproducibility of results.
            % combine_principle     c*sample_size; The minimum element count merging criterion.
            % split_principle       c*sample_size; Criterion for the maximum number of elements per group.
            % if isempty(high_constant)
            %     high_constant = 1.5;
            % end
            % if isempty(high_if)
            %     high_if = double((sample_size/(num_partitions*p)) >= high_constant);
            % end
            %panduan = (sample_size/(num_partitions*p));
            rng(random_number);
            index = randperm(sample_size);
            partition_lengths = [floor(sample_size / num_partitions) * ones(1, num_partitions - 1), sample_size - floor(sample_size / num_partitions) * (num_partitions - 1)];
            partitions = cell(num_partitions,1);
            start_idx = 1;
            for i = 1:num_partitions
                end_idx = start_idx + partition_lengths(i) - 1;
                partitions{i} = index(start_idx:end_idx);
                start_idx = end_idx + 1;
            end
            for i = 1:num_partitions
                partitions{i} = sort(partitions{i});
            end
            partitions0 = partitions;
            coefficients = zeros(p, num_partitions);
            residuals = zeros(sample_size,num_partitions);
            beta_back = zeros(sample_size*p,1);
            beta_back0 = zeros(sample_size*p,1);
            residual_coefficients = 1;
            beta0 = rand(p,1);
            for i = 1:num_partitions
                coefficients(:,i) = beta0;
            end
            coefficients0 = coefficients;
            iter_initial_out = 0;
            if high_if == 0
                while(residual_coefficients > eps_out&&iter_initial_out<=iter_max_initial_out)
                    warning('off', 'all');
                    for i = 1:num_partitions
                        n_part = length(partitions{i});
                        % if n_part < (0.05*sample_size)
                        %     continue;
                        % end
                        x_part = x(partitions{i},:) - mean(x(partitions{i},:));
                        y_part = y(partitions{i}) - mean(y(partitions{i}));
                        beta_part = (x_part'*x_part)\(x_part'*y_part);
                        coefficients((1:p),i) = beta_part;
                        for k = 1:sample_size
                            residuals(k,i) = abs(y(k)-(x(k,:)-mean(x(partitions{i},:)))*beta_part-mean(y_part));
                        end
                    end
                    partitions = cell(num_partitions,1);
                    [~, min_col_indices] = min(residuals, [], 2);
                    for k = 1:sample_size
                        for j = 1:num_partitions
                            if min_col_indices(k)==j
                                partitions{j} = [partitions{j},k];
                            end
                        end
                    end
                    combine = zeros(num_partitions,1);
                    for j = 1:num_partitions
                        combine(j,1) = length(partitions{j});
                    end
                    combine_index = min(combine);
                    combine_min_index = find(combine==min(combine));
                    while (combine_index < combine_principle)
                        num_partitions = num_partitions - length(combine_min_index);
                        residuals(:,combine_min_index) = [];
                        coefficients(:,combine_min_index) = [];
                        partitions = cell(num_partitions,1);
                        [~, min_col_indices] = min(residuals, [], 2);
                        for k = 1:sample_size
                            for j = 1:num_partitions
                                if min_col_indices(k)==j
                                    partitions{j} = [partitions{j},k];
                                end
                            end
                        end
                        combine = zeros(num_partitions,1);
                        for j = 1:num_partitions
                            combine(j,1) = length(partitions{j});
                        end
                        combine_index = min(combine);
                        combine_min_index = find(combine==min(combine));
                    end
                    split_index = max(combine);
                    split_max_index = find(combine==max(combine));
                    while (split_index > split_principle)
                        num_partitions = num_partitions + length(split_max_index);
                        residuals = [residuals,zeros(sample_size,length(split_max_index))];
                        coefficients = [coefficients,zeros(p,length(split_max_index))];
                        partitions = [partitions;cell(length(split_max_index),1)];
                        for k = 1:length(split_max_index)
                            partitions_old = cell(1,1);
                            partitions_old{1} = partitions{split_max_index(k)};
                            sample_size_sum = sample_size;
                            [beta_back,coefficients_split,...
                                partitions_split,~,~] = obj.k_means_regression_tree_single_mex(2,eps_out,...
                                iter_max_initial_out,combine_principle,p,x,y,random_number,beta0,partitions_old,sample_size_sum);
                            partitions{split_max_index(k)} = partitions_split{1};
                            partitions{num_partitions-length(split_max_index)+k} = partitions_split{2};
                            coefficients(:,split_max_index(k)) = coefficients_split(:,1);
                            coefficients(:,(num_partitions-length(split_max_index)+k)) = coefficients_split(:,2);
                        end
                        combine = zeros(num_partitions,1);
                        for j = 1:num_partitions
                            combine(j,1) = length(partitions{j});
                        end
                        split_index = max(combine);
                        split_max_index = find(combine==max(combine));
                    end
                    iter_initial_out = iter_initial_out + 1;
                    for j = 1:num_partitions
                        for k = 1:length(partitions{j})
                            beta_back((((partitions{j}(k)-1)*p)+1):(partitions{j}(k)*p),:) = coefficients((1:p),j);
                        end
                    end
                    residual_coefficients = norm(beta_back-beta_back0);
                    beta_back0 = beta_back;
                end
            elseif high_if == 1
                while(residual_coefficients > eps_out&&iter_initial_out<=iter_max_initial_out)
                    warning('off', 'all');
                    for i = 1:num_partitions
                        n_part = length(partitions{i});
                        if n_part < (0.05*sample_size)
                            continue;
                        end
                        x_part = x(partitions{i},:) - mean(x(partitions{i},:));
                        y_part = y(partitions{i}) - mean(y(partitions{i}));
                        % [B, FitInfo] = lasso(x_part, y_part,"CV", "resubstitution");
                        % beta_part = B(:, FitInfo.MSE == min(FitInfo.MSE));

                        % [B, FitInfo] = lasso(x_part, y_part,"CV", CV_number);
                        % beta_part = B(:, FitInfo.IndexMinMSE);

                        options = struct('alpha', 1, 'nlambda', 100, 'nfold', CV_number); % Lasso回归的设置
                        options.lambda = cvglmnet(x_part, y_part, 'gaussian', options).lambda_min;
                        beta_part = glmnet(x_part, y_part, [], options).beta;

                        coefficients((1:p),i) = beta_part;
                        for k = 1:sample_size
                            residuals(k,i) = abs(y(k)-(x(k,:)-mean(x(partitions{i},:)))*beta_part-mean(y_part));
                        end
                    end
                    partitions = cell(num_partitions,1);
                    [~, min_col_indices] = min(residuals, [], 2);
                    for k = 1:sample_size
                        for j = 1:num_partitions
                            if min_col_indices(k)==j
                                partitions{j} = [partitions{j},k];
                            end
                        end
                    end
                    combine = zeros(num_partitions,1);
                    for j = 1:num_partitions
                        combine(j,1) = length(partitions{j});
                    end
                    combine_index = min(combine);
                    combine_min_index = find(combine==min(combine));
                    while (combine_index < combine_principle)
                        num_partitions = num_partitions - length(combine_min_index);
                        residuals(:,combine_min_index) = [];
                        coefficients(:,combine_min_index) = [];
                        partitions = cell(num_partitions,1);
                        [~, min_col_indices] = min(residuals, [], 2);
                        for k = 1:sample_size
                            for j = 1:num_partitions
                                if min_col_indices(k)==j
                                    partitions{j} = [partitions{j},k];
                                end
                            end
                        end
                        combine = zeros(num_partitions,1);
                        for j = 1:num_partitions
                            combine(j,1) = length(partitions{j});
                        end
                        combine_index = min(combine);
                        combine_min_index = find(combine==min(combine));
                    end
                    split_index = max(combine);
                    split_max_index = find(combine==max(combine));
                    while (split_index > split_principle)
                        num_partitions = num_partitions + length(split_max_index);
                        residuals = [residuals,zeros(sample_size,length(split_max_index))];
                        coefficients = [coefficients,zeros(p,length(split_max_index))];
                        partitions = [partitions;cell(length(split_max_index),1)];
                        for k = 1:length(split_max_index)
                            partitions_old = cell(1,1);
                            partitions_old{1} = partitions{split_max_index(k)};
                            sample_size_sum = sample_size;
                            [beta_back,coefficients_split,...
                                partitions_split,~,~] = obj.k_means_regression_tree_single(2,eps_out,...
                                iter_max_initial_out,combine_principle,p,x,y,random_number,beta0,partitions_old,sample_size_sum, 1, CV_number, K_require_if);
                            partitions{split_max_index(k)} = partitions_split{1};
                            partitions{num_partitions-length(split_max_index)+k} = partitions_split{2};
                            coefficients(:,split_max_index(k)) = coefficients_split(:,1);
                            coefficients(:,(num_partitions-length(split_max_index)+k)) = coefficients_split(:,2);
                        end
                        combine = zeros(num_partitions,1);
                        for j = 1:num_partitions
                            combine(j,1) = length(partitions{j});
                        end
                        split_index = max(combine);
                        split_max_index = find(combine==max(combine));
                    end
                    iter_initial_out = iter_initial_out + 1;
                    for j = 1:num_partitions
                        for k = 1:length(partitions{j})
                            beta_back((((partitions{j}(k)-1)*p)+1):(partitions{j}(k)*p),:) = coefficients((1:p),j);
                        end
                    end
                    residual_coefficients = norm(beta_back-beta_back0);
                    beta_back0 = beta_back;
                end
            end
        end




















        function [beta_back,coefficients,...
                partitions,combine,residual_coefficients] = k_means_regression_tree_single(obj,num_partitions,eps_out,...
                iter_max_initial_out,combine_principle,p,x,y,random_number,beta0,partitions_old,sample_size_sum,high_if,CV_number,K_require_if)
            %------------------------------------------Functionality of the function---------------------------------------------%
            % For any given partitions_old (a set of index values), binary splitting can be performed, and
            % the residual iteration algorithm can be used to solve it.
            % Obtain the split groups' coefficients and the corresponding group indices partitions.
            % It should be noted that the required format for partitions_old is:
            % partitions_old = cell(1,1);partitions_old{1} = partitions{k};
            %-----------------------------------------Description of output variables---------------------------------------------%
            % beta_back             a row vector returned by a binary tree model.  [p*sample_size] matrix
            % coefficients          the regression coefficients of the two models returned by a binary tree model. [p*2] matrix
            % combine               the sample sizes corresponding to the two models returned by a binary tree model. [2*1] matrix
            % partitions            the group indices of the two models returned by a binary tree model.
            %                       (relative to the original x and y group indices). [2*1] cell
            %-----------------------------------------Description of input variables---------------------------------------------%
            % num_partitions        Divide the data into several parts, with 2 as the default.
            % eps_out               The coefficient difference precision between the n-th and (n-1)-th model.
            %                       after sample re-partitioning, and the convergence precision.
            % iter_max_initial_out  The maximum number of iterations for sample re-partitioning.
            % sample_size           sample size.
            % p                     coefficient dimension.
            % x                     The explanatory variables in a single simulation. [sample_size*p]
            % y                     The response variables in a single simulation. [sample_size*1]
            % random_number         Random seed number to ensure reproducibility of results.
            % beta0                 The initial value of beta for each homogeneous model calculation.
            % combine_principle     c*sample_size; The minimum element count merging criterion.
            % sample_size_sum       Total number of samples.
            % partitions_old        Indices to be split.
            rng(random_number);
            sample_size = length(partitions_old{1});
            combine_principle = combine_principle*sample_size_sum;
            index = randperm(sample_size);
            partition_lengths = [floor(sample_size / num_partitions) * ones(1, num_partitions - 1), sample_size - floor(sample_size / num_partitions) * (num_partitions - 1)];
            partitions = cell(num_partitions,1);
            start_idx = 1;
            for i = 1:num_partitions
                end_idx = start_idx + partition_lengths(i) - 1;
                partitions{i} = partitions_old{1}(index(start_idx:end_idx));
                start_idx = end_idx + 1;
            end
            for i = 1:num_partitions
                partitions{i} = sort(partitions{i});
            end
            partitions0 = partitions;
            coefficients = zeros(p,num_partitions);
            residuals = zeros(sample_size,num_partitions);
            beta_back = zeros(sample_size_sum*p,1);
            beta_back0 = zeros(sample_size_sum*p,1);
            residual_coefficients = 1;
            for i = 1:num_partitions
                coefficients(:,i) = beta0;
            end
            coefficients0 = coefficients;
            iter_initial_out = 0;
            if high_if == 0
                while(residual_coefficients > eps_out&&iter_initial_out<=iter_max_initial_out)
                    warning('off', 'all');
                    for i = 1:num_partitions
                        n_part = length(partitions{i});
                        if n_part < (0.05*sample_size_sum)
                            continue;
                        end
                        x_part = x(partitions{i},:)-mean(x(partitions{i},:));
                        y_part = y(partitions{i})- mean(y(partitions{i}));
                        beta_part = (x_part'*x_part)\(x_part'*y_part);
                        coefficients((1:p),i) = beta_part;
                        for k = 1:sample_size
                            residuals(k,i) = abs(y(partitions_old{1}(k))-(x(partitions_old{1}(k),:)-mean(x(partitions{i},:)))*beta_part-mean(y_part));%残差存储矩阵
                        end
                    end
                    partitions = cell(num_partitions,1);
                    for i = 1:num_partitions
                        partitions{i,1} = 0;
                    end
                    [~, min_col_indices] = min(residuals, [], 2);
                    for k = 1:sample_size
                        for j = 1:num_partitions
                            if min_col_indices(k)==j
                                partitions{j} = [partitions{j},partitions_old{1}(k)];
                            end
                        end
                    end

                    for j = 1:num_partitions
                        partitions{j} = partitions{j}(2:end);
                    end
                    combine = zeros(num_partitions,1);
                    for j = 1:num_partitions
                        combine(j,1) = length(partitions{j});
                    end
                    combine_index = min(combine);
                    combine_min_index = find(combine==min(combine));
                    iter_initial_out = iter_initial_out + 1;
                    for j = 1:num_partitions
                        for k = 1:length(partitions{j})
                            beta_back((((partitions{j}(k)-1)*p)+1):(partitions{j}(k)*p),:) = coefficients((1:p),j);
                        end
                    end
                    residual_coefficients = norm(beta_back-beta_back0);
                    beta_back0 = beta_back;

                end
            elseif high_if == 1
                while(residual_coefficients > eps_out&&iter_initial_out<=iter_max_initial_out)
                    warning('off', 'all');
                    for i = 1:num_partitions
                        n_part = length(partitions{i});
                        if n_part < (0.05*sample_size_sum)
                            continue;
                        end
                        x_part = x(partitions{i},:)-mean(x(partitions{i},:));
                        y_part = y(partitions{i})- mean(y(partitions{i}));
                        % [B, FitInfo] = lasso(x_part, y_part,"CV", "resubstitution");
                        % beta_part = B(:, FitInfo.MSE == min(FitInfo.MSE));

                        % [B, FitInfo] = lasso(x_part, y_part,"CV", CV_number);
                        % beta_part = B(:, FitInfo.IndexMinMSE);

                        options = struct('alpha', 1, 'nlambda', 100, 'nfold', CV_number); % Lasso回归的设置
                        options.lambda = cvglmnet(x_part, y_part, 'gaussian', options).lambda_min;
                        beta_part = glmnet(x_part, y_part, [], options).beta;
                        coefficients((1:p),i) = beta_part;
                        for k = 1:sample_size
                            residuals(k,i) = abs(y(partitions_old{1}(k))-(x(partitions_old{1}(k),:)-mean(x(partitions{i},:)))*beta_part-mean(y_part));%残差存储矩阵
                        end
                    end


                    if K_require_if ~= 1
                        partitions = cell(num_partitions,1);
                        for i = 1:num_partitions
                            partitions{i,1} = 0;
                        end
                        [~, min_col_indices] = min(residuals, [], 2);
                        for k = 1:sample_size
                            for j = 1:num_partitions
                                if min_col_indices(k)==j
                                    partitions{j} = [partitions{j},partitions_old{1}(k)];
                                end
                            end
                        end
                        for j = 1:num_partitions
                            partitions{j} = partitions{j}(2:end);
                        end
                        combine = zeros(num_partitions,1);
                        for j = 1:num_partitions
                            combine(j,1) = length(partitions{j});
                        end
                        combine_index = min(combine);
                        combine_min_index = find(combine==min(combine));
                    end


                    iter_initial_out = iter_initial_out + 1;
                    for j = 1:num_partitions
                        for k = 1:length(partitions{j})
                            beta_back((((partitions{j}(k)-1)*p)+1):(partitions{j}(k)*p),:) = coefficients((1:p),j);
                        end
                    end
                    residual_coefficients = norm(beta_back-beta_back0);
                    beta_back0 = beta_back;

                end
            end

            combine = zeros(num_partitions,1);
            for j = 1:num_partitions
                combine(j,1) = length(partitions{j});
            end
        end



















        function [beta_back,coefficients,...
                partitions,combine,residual_coefficients] = k_means_regression_tree_single_mex(obj,num_partitions,eps_out,...
                iter_max_initial_out,combine_principle,p,x,y,random_number,beta0,partitions_old,sample_size_sum)

            [beta_back,coefficients,...
                partitions,combine,residual_coefficients] = k_means_regression_tree_single_mex(num_partitions,eps_out,...
                iter_max_initial_out,combine_principle,p,x,y,random_number,beta0,partitions_old,sample_size_sum);
        end















        function [num_partitions_best,combine_principle_best,split_principle_best,random_number_best,cv_matrix] = k_means_regression_cv(obj,num_partitions_list,eps_out,...
                iter_max_initial_out,combine_principle_list,split_principle_list,random_number_list,sample_size,p,x,y,disp_if,c,high_if,CV_number,K_require_if)
            %------------------------------------------Functionality of the function-----------------------------------------------%
            % Select the optimal parameters for the cross-validation model of the residual iteration model.
            %-----------------------------------------Description of output variables---------------------------------------------%
            % num_partitions_best      The optimal number of folds obtained through cross-validation.
            % combine_principle_best   The optimal merging criterion number obtained through cross-validation.
            % split_principle_best     The optimal splitting criterion number obtained through cross-validation.
            % random_number_best       The optimal random seed number obtained through cross-validation.
            % cv_matrix                The error tensor corresponding to the parameters obtained through cross-validation.
            %-----------------------------------------Description of input variables---------------------------------------------%
            % num_partitions_list       Divide the data into several parts, parameter matrix.
            % eps_initial               The convergence precision of each model during each homogeneous model regression.
            % iter_max_initial_in       The maximum number of iterations for alternating optimization during each homogeneous model regression.
            % eps_out                   The coefficient difference precision and convergence precision between the n-th and (n-1)-th models after sample re-partitioning.
            % iter_max_initial_out      The maximum number of iterations for sample re-partitioning.
            % sample_size               Sample size.
            % row_size                  Row vector dimensions.
            % col_size                  Column vector dimensions.
            % x                         The explanatory variables in a single simulation. [row_size*col_size*sample_size]
            % y                         The response variables in a single simulation. [sample_size*1]
            % random_number_list        The random seed parameter matrix.
            % combine_principle_list    The parameter matrix for the minimum element count merging criterion.
            % split_principle_list      The parameter matrix for the maximum element count splitting criterion.
            cv_matrix = inf*ones(length(num_partitions_list),length(combine_principle_list),length(split_principle_list),length(random_number_list));
            for kk = 1:length(random_number_list)
                for ii = 1:length(num_partitions_list)
                    for jj = 1:length(combine_principle_list)
                        for qq = 1:length(split_principle_list)
                            % if isempty(high_constant)
                            %     high_constant = 1.5;
                            % end
                            % if isempty(high_if)
                            %     high_if = double((sample_size/(num_partitions_list(ii)*p)) >= high_constant);
                            % end
                            % panduan = sample_size/(p*num_partitions_list(ii));
                            warning('off', 'all');
                            rng(random_number_list(kk));
                            num_partitions = num_partitions_list(ii);
                            combine_principle = combine_principle_list(jj);
                            split_principle = split_principle_list(qq);
                            index = randperm(sample_size);
                            partition_lengths = [floor(sample_size / num_partitions) * ones(1, num_partitions - 1), sample_size - floor(sample_size / num_partitions) * (num_partitions - 1)];
                            partitions = cell(num_partitions,1);
                            start_idx = 1;
                            for i = 1:num_partitions
                                end_idx = start_idx + partition_lengths(i) - 1;
                                partitions{i} = index(start_idx:end_idx);
                                start_idx = end_idx + 1;
                            end
                            for i = 1:num_partitions
                                partitions{i} = sort(partitions{i});
                            end
                            partitions0 = partitions;
                            coefficients = zeros(p,num_partitions);
                            residuals = zeros(sample_size,num_partitions);
                            beta_back = zeros(sample_size*p,1);
                            beta_back0 = zeros(sample_size*p,1);
                            residual_coefficients = 1;
                            beta0 = rand(p,1);
                            for i = 1:num_partitions
                                coefficients(:,i) = beta0;
                            end
                            coefficients0 = coefficients;
                            iter_initial_out = 0;

                            if high_if == 0
                                while(residual_coefficients > eps_out&&iter_initial_out<=iter_max_initial_out)
                                    for i = 1:num_partitions
                                        iter_initial_in = 0;
                                        n_part = length(partitions{i});
                                        if n_part < (0.05*sample_size)
                                            continue;
                                        end
                                        x_part = x(partitions{i},:)-mean(x(partitions{i},:));
                                        y_part = y(partitions{i})- mean(y(partitions{i}));
                                        beta_part = (x_part'*x_part)\(x_part'*y_part);
                                        coefficients((1:p),i) = beta_part;
                                        for k = 1:sample_size
                                            residuals(k,i) = abs(y(k)-(x(k,:)-mean(x(partitions{i},:)))*beta_part-mean(y_part));
                                        end

                                    end
                                    partitions = cell(num_partitions,1);
                                    [~, min_col_indices] = min(residuals, [], 2);
                                    for k = 1:sample_size
                                        for j = 1:num_partitions
                                            if min_col_indices(k)==j
                                                partitions{j} = [partitions{j},k];
                                            end
                                        end
                                    end
                                    combine = zeros(num_partitions,1);
                                    for j = 1:num_partitions
                                        combine(j,1) = length(partitions{j});
                                    end
                                    combine_index = min(combine);
                                    combine_min_index = find(combine==min(combine));
                                    while (combine_index < combine_principle)
                                        num_partitions = num_partitions - length(combine_min_index);
                                        residuals(:,combine_min_index) = [];
                                        coefficients(:,combine_min_index) = [];
                                        partitions = cell(num_partitions,1);
                                        [~, min_col_indices] = min(residuals, [], 2);
                                        for k = 1:sample_size
                                            for j = 1:num_partitions
                                                if min_col_indices(k)==j
                                                    partitions{j} = [partitions{j},k];
                                                end
                                            end
                                        end
                                        combine = zeros(num_partitions,1);
                                        for j = 1:num_partitions
                                            combine(j,1) = length(partitions{j});
                                        end

                                        combine_index = min(combine);
                                        combine_min_index = find(combine==min(combine));
                                    end
                                    split_index = max(combine);
                                    split_max_index = find(combine==max(combine));
                                    while (split_index > split_principle)
                                        num_partitions = num_partitions + length(split_max_index);
                                        residuals = [residuals,zeros(sample_size,length(split_max_index))];
                                        coefficients = [coefficients,zeros((p),length(split_max_index))];
                                        partitions = [partitions;cell(length(split_max_index),1)];
                                        for k = 1:length(split_max_index)
                                            partitions_old = cell(1,1);
                                            partitions_old{1} = partitions{split_max_index(k)};
                                            sample_size_sum = sample_size;
                                            random_number = random_number_list(kk);
                                            [beta_back,coefficients_split,...
                                                partitions_split,~,~] = obj.k_means_regression_tree_single_mex(2,eps_out,...
                                                iter_max_initial_out,combine_principle,p,x,y,random_number,beta0,partitions_old,sample_size_sum);
                                            partitions{split_max_index(k)} = partitions_split{1};
                                            partitions{num_partitions-length(split_max_index)+k} = partitions_split{2};
                                            coefficients(:,split_max_index(k)) = coefficients_split(:,1);
                                            coefficients(:,(num_partitions-length(split_max_index)+k)) = coefficients_split(:,2);
                                        end
                                        combine = zeros(num_partitions,1);
                                        for j = 1:num_partitions
                                            combine(j,1) = length(partitions{j});
                                        end
                                        split_index = max(combine);
                                        split_max_index = find(combine==max(combine));
                                    end
                                    iter_initial_out = iter_initial_out + 1;
                                    for j = 1:num_partitions
                                        for k = 1:length(partitions{j})
                                            beta_back((((partitions{j}(k)-1)*p)+1):(partitions{j}(k)*p),:) = coefficients((1:p),j);
                                        end
                                    end
                                    residual_coefficients = norm(beta_back-beta_back0);
                                    beta_back0 = beta_back;
                                end
                                try
                                    cv_matrix(ii,jj,qq,kk) = 0;
                                    for k = 1:sample_size
                                        cv_matrix(ii,jj,qq,kk) = cv_matrix(ii,jj,qq,kk) + (1/sample_size)*(y(k)-(x(k,:)-mean(x(partitions{find(cellfun(@(x) any(x == k), partitions), 1)},:)))*beta_back(((k-1)*p+1):(k*p),1)...
                                            -mean(y(partitions{find(cellfun(@(x) any(x == k), partitions), 1)}))).^2;

                                    end
                                    % if K_require_if ~= 1
                                    %     cv_matrix(ii,jj,qq,kk) = log(cv_matrix(ii,jj,qq,kk)) + c*log(sample_size*p)*((log(sample_size)/sample_size)*length(partitions))*p; %c=2
                                    % else
                                    %     cv_matrix(ii,jj,qq,kk) = log(cv_matrix(ii,jj,qq,kk)) + ((log(sample_size)/sample_size)*length(partitions))*p;
                                    % end
                                    cv_matrix(ii,jj,qq,kk) = log(cv_matrix(ii,jj,qq,kk)) + c*log(sample_size*p)*((log(sample_size)/sample_size)*length(partitions))*p; %c=2
                                    cv_matrix(ii,jj,qq,kk) = log(cv_matrix(ii,jj,qq,kk)) + ((log(sample_size)/sample_size)*length(partitions))*p;

                                    %cv_matrix(ii,jj,qq,kk) = log(cv_matrix(ii,jj,qq,kk)) + c*log(sample_size*p)*((log(sample_size)/sample_size)*length(partitions));
                                    %c=5
                                catch
                                    cv_matrix(ii,jj,qq,kk) = inf;
                                    continue;
                                end
                            elseif high_if == 1
                                while(residual_coefficients > eps_out&&iter_initial_out<=iter_max_initial_out)
                                    for i = 1:num_partitions
                                        iter_initial_in = 0;
                                        n_part = length(partitions{i});
                                        if n_part < (0.05*sample_size)
                                            continue;
                                        end
                                        x_part = x(partitions{i},:)-mean(x(partitions{i},:));
                                        y_part = y(partitions{i})- mean(y(partitions{i}));
                                        % [B, FitInfo] = lasso(x_part, y_part,"CV", "resubstitution");
                                        % beta_part = B(:, FitInfo.MSE == min(FitInfo.MSE));

                                        % [B, FitInfo] = lasso(x_part, y_part,"CV", CV_number);
                                        % beta_part = B(:, FitInfo.IndexMinMSE);

                                        options = struct('alpha', 1, 'nlambda', 100, 'nfold', CV_number); % Lasso回归的设置
                                        options.lambda = cvglmnet(x_part, y_part, 'gaussian', options).lambda_min;
                                        beta_part = glmnet(x_part, y_part, [], options).beta;
                                        coefficients((1:p),i) = beta_part;
                                        for k = 1:sample_size
                                            residuals(k,i) = abs(y(k)-(x(k,:)-mean(x(partitions{i},:)))*beta_part-mean(y_part));
                                        end

                                    end
                                    partitions = cell(num_partitions,1);
                                    [~, min_col_indices] = min(residuals, [], 2);
                                    for k = 1:sample_size
                                        for j = 1:num_partitions
                                            if min_col_indices(k)==j
                                                partitions{j} = [partitions{j},k];
                                            end
                                        end
                                    end
                                    combine = zeros(num_partitions,1);
                                    for j = 1:num_partitions
                                        combine(j,1) = length(partitions{j});
                                    end
                                    combine_index = min(combine);
                                    combine_min_index = find(combine==min(combine));
                                    while (combine_index < combine_principle)
                                        num_partitions = num_partitions - length(combine_min_index);
                                        residuals(:,combine_min_index) = [];
                                        coefficients(:,combine_min_index) = [];
                                        partitions = cell(num_partitions,1);
                                        [~, min_col_indices] = min(residuals, [], 2);
                                        for k = 1:sample_size
                                            for j = 1:num_partitions
                                                if min_col_indices(k)==j
                                                    partitions{j} = [partitions{j},k];
                                                end
                                            end
                                        end
                                        combine = zeros(num_partitions,1);
                                        for j = 1:num_partitions
                                            combine(j,1) = length(partitions{j});
                                        end

                                        combine_index = min(combine);
                                        combine_min_index = find(combine==min(combine));
                                    end
                                    split_index = max(combine);
                                    split_max_index = find(combine==max(combine));
                                    while (split_index > split_principle)
                                        num_partitions = num_partitions + length(split_max_index);
                                        residuals = [residuals,zeros(sample_size,length(split_max_index))];
                                        coefficients = [coefficients,zeros((p),length(split_max_index))];
                                        partitions = [partitions;cell(length(split_max_index),1)];
                                        for k = 1:length(split_max_index)
                                            partitions_old = cell(1,1);
                                            partitions_old{1} = partitions{split_max_index(k)};
                                            sample_size_sum = sample_size;
                                            random_number = random_number_list(kk);
                                            [beta_back,coefficients_split,...
                                                partitions_split,~,~] = obj.k_means_regression_tree_single(2,eps_out,...
                                                iter_max_initial_out,combine_principle,p,x,y,random_number,beta0,partitions_old,sample_size_sum, 1, CV_number, K_require_if);
                                            partitions{split_max_index(k)} = partitions_split{1};
                                            partitions{num_partitions-length(split_max_index)+k} = partitions_split{2};
                                            coefficients(:,split_max_index(k)) = coefficients_split(:,1);
                                            coefficients(:,(num_partitions-length(split_max_index)+k)) = coefficients_split(:,2);
                                        end
                                        combine = zeros(num_partitions,1);
                                        for j = 1:num_partitions
                                            combine(j,1) = length(partitions{j});
                                        end
                                        split_index = max(combine);
                                        split_max_index = find(combine==max(combine));
                                    end
                                    iter_initial_out = iter_initial_out + 1;
                                    for j = 1:num_partitions
                                        for k = 1:length(partitions{j})
                                            beta_back((((partitions{j}(k)-1)*p)+1):(partitions{j}(k)*p),:) = coefficients((1:p),j);
                                        end
                                    end
                                    residual_coefficients = norm(beta_back-beta_back0);
                                    beta_back0 = beta_back;
                                end
                                try
                                    cv_matrix(ii,jj,qq,kk) = 0;
                                    for k = 1:sample_size
                                        cv_matrix(ii,jj,qq,kk) = cv_matrix(ii,jj,qq,kk) + (1/sample_size)*(y(k)-(x(k,:)-mean(x(partitions{find(cellfun(@(x) any(x == k), partitions), 1)},:)))*beta_back(((k-1)*p+1):(k*p),1)...
                                            -mean(y(partitions{find(cellfun(@(x) any(x == k), partitions), 1)}))).^2;

                                    end
                                    % if K_require_if ~= 1
                                    %     cv_matrix(ii,jj,qq,kk) = log(cv_matrix(ii,jj,qq,kk)) + c*log(sample_size*p)*((log(sample_size)/sample_size)*length(partitions))*p; %c=2
                                    % else
                                    %     cv_matrix(ii,jj,qq,kk) = log(cv_matrix(ii,jj,qq,kk)) + ((log(sample_size)/sample_size)*length(partitions))*p;
                                    % end
                                    cv_matrix(ii,jj,qq,kk) = log(cv_matrix(ii,jj,qq,kk)) + c*log(sample_size*p)*((log(sample_size)/sample_size)*length(partitions))*p; %c=2
                                    %cv_matrix(ii,jj,qq,kk) = log(cv_matrix(ii,jj,qq,kk)) + ((log(sample_size)/sample_size)*length(partitions))*p;

                                    %cv_matrix(ii,jj,qq,kk) = log(cv_matrix(ii,jj,qq,kk)) + c*log(sample_size*p)*((log(sample_size)/sample_size)*length(partitions));
                                    %c=5
                                catch
                                    cv_matrix(ii,jj,qq,kk) = inf;
                                    continue;
                                end
                            end
                        end
                    end
                    if disp_if >0
                        fprintf('Initial value solving process completed %.1f%%\n',(100*kk/length(random_number_list)))
                    end
                end
            end

            [minValue, linearIndex] = min(cv_matrix(:));
            [row, col, layer, layer2] = ind2sub(size(cv_matrix), linearIndex);
            num_partitions_best = num_partitions_list(row);
            combine_principle_best = combine_principle_list(col);
            split_principle_best = split_principle_list(layer);
            random_number_best = random_number_list(layer2);
        end
    end
end

