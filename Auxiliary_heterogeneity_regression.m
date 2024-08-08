classdef Auxiliary_heterogeneity_regression
    methods
        function [Data] = admm(obj, sample_size, p, iter_max1, eps1,...
                beta0, H_p, X, y, k, a, lambda, c, v, kesai,...
                v_vec, v_vec_bar, kesai_vec, beta_real, HH_p, E, iter_max2,eps2,iter_max_in,eps_in,...
                    eps_in_in,iter_max_in_in,full_if)
            %-------------------------------------------函数功能-----------------------------------------------%
            % 使用ADMM算法求解单个lambda对应的目标函数
            %-----------------------------------------输出变量说明---------------------------------------------%
            % Data.lambda           返回当前lambda
            % Data.beta             返回估计出的参数beta  [row_size*sample_size] matrix
            % Data.class            返回亚组划分的结果，每一行代表一个组，这一行为这一组的原始索引值
            % Data.index            返回亚组划分的结果，每一行的第一列代表样本索引，第二列代表组号
            % Data.subgroup_number  返回识别的亚族数
            % Data.per_if           返回亚组数准确识别标记，若为1则为亚组数识别准确，为0则不准确
            % Data.mse_beta         估计出的beta的均方误差
            % Data.sc               估计出的一致性指数
            % Data.BIC              返回BIC值
            % Data.BIC_part1        BIC值的第一部分
            % Data.BIC_part2        BIC值的第二部分
            % Data.convergence      第 n 步与第 n-1 步ADMM交替求解的原始残差的二范数（考虑维度平均意义下的差值二范数）
            % Data.iter             ADMM迭代步数
            %-----------------------------------------输入变量说明---------------------------------------------%
            % sample_size           sample_size
            % p                     coefficient dimension
            % iter_max1             ADMM最大交替迭代步数
            % eps1                  ADMM原始残差二范数的容差
            % beta0                 初值beta0
            % H_p                   H_p维度为（组合数*p）*(样本数*p) matrix
            % X                     解释变量数据 （样本数）*(样本数*p) matrix
            % y                     连续相应变量 （样本数*1）matrix
            % a                     MCP惩罚函数里的控制凹凸性的参数
            % k                     ADMM算法的惩罚参数
            % v                     beta组合差信息储存矩阵（对偶变量）
            % lambda                惩罚调节参数
            % c                     BIC准则第二部分的权重参数
            % kesai　　　　　　　　　关于beta的增广lagrange乘子
            % v_vec                 beta组合差信息储存矩阵（对偶变量）按组合顺序拉成向量
            % kesai_vec             beta的增广lagrange乘子按组合顺序拉成向量
            % beta_real             真实的beta
            % HH_p                  H_p'*H_p
            beta = beta0;
            iter = 0;
            residual_primal = 1;
            per_if = 0;
            sc_hat_store = zeros(sample_size,sample_size);
            sc_real_store = zeros(sample_size,sample_size);
            sc_store = zeros(sample_size,sample_size);
            max_time = 60; % 设置最大执行时间为 60 秒
            start_time = tic; % 开始计时

            
            for i = 1:(sample_size-1)%#######
                v(i,i,:) = 0;%#######
                for j = i+1:sample_size%#######
                    v(i,j,:) = beta(((i-1)*p+1):(i*p)) - beta(((j-1)*p+1):(j*p));%#######
                    v(j,i,:) = v(i,j,:);%#######
                end%#######
            end%#######
            [Class_iter,~,~,~,~] = obj.Cluster_iter_pro(sample_size, p, v, X, y, beta, beta0);%#######
            while ((residual_primal > eps1)&&(iter<=iter_max1))
                %---------------------------------------------------------------------------------------------%
                %----------------------------------------ADMM第一步-------------------------------------------%
                %---------------------------------------------------------------------------------------------%
                %beta = (X'*X + k*HH_p)\(X'*y+k*reshape(((reshape(v_vec, p, sample_size*(sample_size-1)*0.5)-...
                %    (1/k)*reshape(kesai_vec, p, sample_size*(sample_size-1)*0.5))*E), [], 1));

                [beta,~] = obj.admm_part(sample_size,p,iter_max2,eps2,iter_max_in,eps_in,...
                    eps_in_in,iter_max_in_in,beta,X,y,a,kesai,Class_iter);
                % admm_part(obj, sample_size,p,iter_max2,eps2,iter_max_in,eps_in,...
                % eps_in_in,iter_max_in_in,beta,X,y,a,kesai,Class_iter)
                % [beta_back,convergence] = admm_part(obj, sample_size,p,iter_max2,eps2,iter_max_in,eps_in,...
                % eps_in_in,iter_max_in_in,beta,X,y,a,kesai,Class_iter)
                %---------------------------------------------------------------------------------------------%
                %----------------------------------------ADMM第二步-------------------------------------------%
                %---------------------------------------------------------------------------------------------%
                %%v_vec = zeros(sample_size*(sample_size-1)*0.5*row_size,1);
                %%kesai_vec = zeros((sample_size*(sample_size-1)*0.5*row_size),1);
                %v_vector中关于beta第ij系数差的范围(((((2*sample_size-i)*(i-1)/2+(j-i))-1)*row_size)+1):((((2*sample_size-i)*(i-1)/2)+(j-i))*row_size)
                %%v_vec_bar = zeros(sample_size*(sample_size-1)*0.5*row_size,1);
                %v(i,j,:) = beta(((i-1)*row_size+1):(i*row_size)) - beta(((j-1)*row_size+1):(j*row_size));
                %beta(((i-1)*row_size+1):(i*row_size),:)为第i个beta
                for i = 1:(sample_size-1)
                    for j = (i+1):sample_size
                        v_vec_bar((((((2*sample_size-i)*(i-1)/2+(j-i))-1)*p)+1):((((2*sample_size-i)*(i-1)/2)+(j-i))*p))=...
                            beta(((i-1)*p+1):(i*p))-beta(((j-1)*p+1):(j*p))+...
                            +kesai_vec((((((2*sample_size-i)*(i-1)/2+(j-i))-1)*p)+1):((((2*sample_size-i)*(i-1)/2)+(j-i))*p))/k;

                        v_bar_norm = norm([v_vec_bar((((((2*sample_size-i)*(i-1)/2+(j-i))-1)*p)+1):((((2*sample_size-i)*(i-1)/2)+(j-i))*p))]);
                        if v_bar_norm>(a*lambda)
                            v_vec((((((2*sample_size-i)*(i-1)/2+(j-i))-1)*p)+1):((((2*sample_size-i)*(i-1)/2)+(j-i))*p))=...
                                v_vec_bar((((((2*sample_size-i)*(i-1)/2+(j-i))-1)*p)+1):((((2*sample_size-i)*(i-1)/2)+(j-i))*p));
                        else
                            v_vec((((((2*sample_size-i)*(i-1)/2+(j-i))-1)*p)+1):((((2*sample_size-i)*(i-1)/2)+(j-i))*p))=...
                                (max((1-lambda/(k*v_bar_norm)),0)/(1-1/(a*k)))*v_vec_bar((((((2*sample_size-i)*(i-1)/2+(j-i))-1)*p)+1):((((2*sample_size-i)*(i-1)/2)+(j-i))*p));
                        end
                        kesai_vec((((((2*sample_size-i)*(i-1)/2+(j-i))-1)*p)+1):((((2*sample_size-i)*(i-1)/2)+(j-i))*p))=...
                            kesai_vec((((((2*sample_size-i)*(i-1)/2+(j-i))-1)*p)+1):((((2*sample_size-i)*(i-1)/2)+(j-i))*p))+...
                            k*(beta(((i-1)*p+1):(i*p))-beta(((j-1)*p+1):(j*p))-v_vec((((((2*sample_size-i)*(i-1)/2+(j-i))-1)*p)+1):((((2*sample_size-i)*(i-1)/2)+(j-i))*p)));
                    end
                end
                %--------------------------------------------更新v,w,kesai,yita---------------------------------------------%
                for i = 1:(sample_size-1)
                    v(i,i,:) = 0;
                    for j = (i+1):(sample_size)
                        v(i,j,:) = v_vec((((((2*sample_size-i)*(i-1)/2+(j-i))-1)*p)+1):((((2*sample_size-i)*(i-1)/2)+(j-i))*p));
                        kesai(i,j,:) = kesai_vec((((((2*sample_size-i)*(i-1)/2+(j-i))-1)*p)+1):((((2*sample_size-i)*(i-1)/2)+(j-i))*p));
                        v(j,i,:) = v(i,j,:);
                    end
                end
                %--------------------------------------------系数差和原始残差计算-----------------------------------------------%
                    % [combine_number,~] = size(E);
                    % residual_vec = zeros(combine_number*p,1);
                    % for i = 1:combine_number
                    %     residual_vec(((i-1)*p+1):(i*p),1) = kron(E(i,:),eye(p))*beta;
                    % end
                    % residual_primal = norm((residual_vec-v_vec)/sqrt(0.5*sample_size*(sample_size-1)*p));
                    residual_primal = obj.Residual_primal(full_if,E,beta,v_vec,sample_size,p,H_p);
                    % if full_if == 1
                    %     [combine_number,~] = size(E);
                    %     residual_vec = zeros(combine_number*(p),1);
                    %     for i = 1:combine_number
                    %         residual_vec(((i-1)*p+1):(i*p),1) = kron(E(i,:),eye(p))*beta;
                    %     end
                    %     residual_primal = norm((residual_vec-v_vec)/sqrt(0.5*sample_size*(sample_size-1)*p));
                    % else
                    %     residual_primal = norm((H_p*beta-v_vec)/sqrt(0.5*sample_size*(sample_size-1)*p));
                    % end
                    iter = iter + 1;

                    % [Class_iter,~,~,~,~] = obj_cluster.Cluster_iter_pro(sample_size, p, v, X, y, beta, beta0);%#######
                    % for i = 1 : Class_iter.subgroup_number%#######
                    %     Class_iter.coef{i} = beta(((Class_iter.index{i}(1)-1)*p+1):((Class_iter.index{i}(1))*p));%#######
                    % end%#######
                    % for t = 1:sample_size%#######
                    %     panduan = inf*ones(1, Class_iter.subgroup_number);%#######
                    %     %f_store = obj.store2(Class_iter, v, kesai, t);
                    %     which_old = find(cellfun(@(c) ismember(t, c), Class_iter.index));%#######
                    %     for i = 1 : Class_iter.subgroup_number%#######
                    %         panduan(1, i) = 0.5*(y(t) - X(t,((t-1)*p+1):(t*p))*Class_iter.coef{i}).^2;%#######
                    %         %+ mcp_function(norm(Class_iter.coef{i}-W.coef(:,which_aux)),a,lambda_2);
                    %         %+ f_store(i).store;
                    %     end%#######
                    %     [~, which] = min(panduan(1, :));%#######
                    %     beta((((t-1)*p+1):(t*p))) = Class_iter.coef{which};%#######
                    %     if which_old ~= which%#######
                    %         Class_iter.index{which} = sort([Class_iter.index{which},t]);%#######
                    %         Class_iter.index{which_old}(Class_iter.index{which_old} == t) = [];%#######
                    %     end%#######
                    % 
                    %     Class_iter.subgroup_number = sum(double(~cellfun(@isempty, Class_iter.index)));%#######
                    %     Class_iter.coef(cellfun(@isempty, Class_iter.index)) = [];%#######
                    %     Class_iter.index(cellfun(@isempty, Class_iter.index)) = [];%#######
                    %     Class_iter.subgroup_number_size = cell(1, Class_iter.subgroup_number);%#######
                    %     Class_iter.X = cell(1,Class_iter.subgroup_number);%#######
                    %     Class_iter.y = cell(1,Class_iter.subgroup_number);%#######
                    %     for i = 1 : Class_iter.subgroup_number%#######
                    %         Class_iter.subgroup_number_size{i} = length(Class_iter.index{i}); %#######
                    %         Class_iter.X{i} = zeros(Class_iter.subgroup_number_size{i},p); %#######
                    %         Class_iter.y{i} = y(Class_iter.index{i},1); %#######
                    %         for j = 1 : Class_iter.subgroup_number_size{i} %#######
                    %             Class_iter.X{i}(j,:) = X(Class_iter.index{i}(j),(((Class_iter.index{i}(j)-1)*p)+1):(Class_iter.index{i}(j)*p));%#######
                    %         end %#######
                    %     end %#######
                    % end %#######

                [Class_iter,~,~,~,~] = obj.Cluster_iter_pro(sample_size, p, v, X, y, beta, beta0);%#######




                
                elapsed_time = toc(start_time);
                if elapsed_time > max_time
                    %disp('超过给定时间，跳出循环。');
                    break;
                end
            end
            %----------------------------------------------亚组结构识别-------------------------------------------------%
            class_matrix = squeeze(v(:,:,1));
            class_matrix(sample_size,sample_size) = 0;
            for i = 1:(sample_size-1)
                class_matrix(i,i) = 0;
                for j = (i+1):sample_size
                    if class_matrix(i,j)~= 0
                        class_matrix(i,j) = 1;
                        class_matrix(j,i) = 1;
                    else
                        class_matrix(j,i) = 0;
                    end
                end
            end
            [v_class,Iv,~] = unique(class_matrix,'rows');
            class_store = zeros(length(Iv),sample_size);
            for i = 1:length(Iv)
                store_tem = [];
                for j = 1:sample_size
                    if class_matrix(j,:)==v_class(i,:)
                        store_tem = [store_tem,j];
                    end
                end
                class_store(i,(1:length(store_tem))) = store_tem;
            end
            class_store(class_store==0)= nan;
            index_store = [(1:sample_size)',zeros(sample_size,1)];
            for i = 1:sample_size
                for j = 1:length(Iv)
                    if class_matrix(i,:) == v_class(j,:)
                        index_store(i,2)= j;
                        break
                    end
                end
            end
            %class_store %存储每组下的原始序号（每一行代表一组，元素为原始序号）
            %index_store %存储原始序号对应的组标签（第一列为原始序号，第二列为对应的组标签）
            %----------------------------------------\hat{k}_number是否等于k_real_number--------------------------------------%
            theta_real_re = [reshape(beta_real,p,sample_size)];
            [~,Itheta,~] = unique(theta_real_re','rows');
            k_real_number = length(Itheta);
            k_hat_number = length(Iv);%为估计出的亚组数
            if k_hat_number == k_real_number
                per_if = 1;%判定估计出的亚组数与真实系数亚组数是否相等
            end
            %------------------------------------------真实亚组结构识别-------------------------------------------------------%
            v_real = zeros(sample_size,sample_size,p);
            for i = 1:(sample_size-1)
                for j = (i+1):sample_size
                    v_real(i,j,:) = abs(theta_real_re(:,i)-theta_real_re(:,j));
                end
            end
            class_real_matrix = squeeze(sum(abs(v_real(:,:,:)),3));
            class_real_matrix(sample_size,sample_size) = 0;
            for i = 1:(sample_size-1)
                for j = (i+1):sample_size
                    if class_real_matrix(i,j)~= 0
                        class_real_matrix(i,j) = 1;
                        class_real_matrix(j,i) = 1;
                    end
                end
            end
            [v_real_class,Ireal,~] = unique(class_real_matrix,'rows');
            index_real_store = [(1:sample_size)',zeros(sample_size,1)];
            for i = 1:sample_size
                for j = 1:length(Ireal)
                    if class_real_matrix(i,:) == v_real_class(j,:)
                        index_real_store(i,2)= j;
                        break
                    end
                end
            end
            %------------------------------------------------计算估计参数的mse------------------------------------------------%
            mse_beta = sqrt(sum((beta-beta_real).^2)/(sample_size*p));
            %-----------------------------------------------计算亚组划分一致性------------------------------------------------%
            for i = 1:(sample_size-1)
                for j = (i+1):sample_size
                    if index_store(i,2)==index_store(j,2)
                        sc_hat_store(i,j)=1;
                    end
                end
            end
            for i = 1:(sample_size-1)
                for j = (i+1):sample_size
                    if index_real_store(i,2)==index_real_store(j,2)
                        sc_real_store(i,j)=1;
                    end
                end
            end
            for i = 1:(sample_size-1)
                for j = (i+1):sample_size
                    if sc_hat_store(i,j)==sc_real_store(i,j)
                        sc_store(i,j)=1;
                    end
                end
            end
            sc = sum(sum(sc_store))/nchoosek(sample_size,2);
            %-----------------------------------------------计算预测误差值---------------------------------------------------%
            wucha = log(sum((y-X*beta).^2)/sample_size);
            %------------------------------------------------计算BIC值---------------------------------------------------%
            %bic = log(sum((y-beta'*X*gamma).^2)/sample_size)+...
            %      log(log(sample_size*(row_size+col_size)))*(log(sample_size)/sample_size)*(k_hat_number*(row_size+col_size));

            bic = wucha+c*log(log(sample_size*p))*(log(sample_size)/sample_size)*(k_hat_number);
            %%%%%%%%bic = wucha + c*log(sample_size*p)*(log(sample_size)/sample_size)*(k_hat_number)*p;
            %bic = wucha+log(sample_size*(row_size+col_size))*(log(sample_size)/sample_size)*(k_hat_number*(row_size+col_size));
            bic = wucha + (log(sample_size)/sample_size)*(k_hat_number)*p;
            %-----------------------------------------------计算预测BIC第二项---------------------------------------------------%
            BIC_part2 = c*log(log(sample_size*p))*(log(sample_size)/sample_size)*(k_hat_number);
            %%%%%%%BIC_part2 = c*log(sample_size*p)*(log(sample_size)/sample_size)*(k_hat_number)*p;
            %BIC_part2 = log(sample_size*(row_size+col_size))*(log(sample_size)/sample_size)*(k_hat_number*(row_size+col_size));
            BIC_part2 = (log(sample_size)/sample_size)*(k_hat_number)*p;
            %------------------------------------------------函数输出---------------------------------------------------%
            BIC_part2_back = BIC_part2;
            beta_back = beta;
            class_back = class_store;
            index_back = index_store;
            subgroup_number = k_hat_number;
            convergence = residual_primal;
            BIC_back = bic;
            wucha_back = wucha;
            iter_back = iter-1;
            Data.lambda = lambda;
            Data.beta = beta_back;
            Data.per_if = per_if;
            Data.subgroup_number = subgroup_number;
            Data.sc = sc;
            Data.mse_beta = mse_beta;
            Data.class = class_back;
            Data.index = index_back;
            Data.BIC = BIC_back;
            Data.BIC_part1 = wucha_back;
            Data.BIC_part2 = BIC_part2;
            Data.convergence = convergence;
            Data.iter = iter_back;
        end









        function residual_primal = Residual_primal(obj,full_if,E,beta,v_vec,sample_size,p,H_p)
            if full_if == 1
                [combine_number,~] = size(E);
                residual_vec = zeros(combine_number*(p),1);
                for i = 1:combine_number
                    residual_vec(((i-1)*p+1):(i*p),1) = kron(E(i,:),eye(p))*beta;
                end
                residual_primal = norm((residual_vec-v_vec)/sqrt(0.5*sample_size*(sample_size-1)*p));
            else
                residual_primal = norm((H_p*beta-v_vec)/sqrt(0.5*sample_size*(sample_size-1)*p));
            end
        end





        function [Class_iter,class_matrix,class_store,index_store, Iv] = Cluster_iter_pro(obj, sample_size, p, v, X, y, beta, beta0)
            class_matrix = squeeze(v(:,:,1));
            class_matrix(sample_size,sample_size) = 0;
            for i = 1:(sample_size-1)
                class_matrix(i,i) = 0;
                for j = (i+1):sample_size
                    if class_matrix(i,j)~= 0
                        class_matrix(i,j) = 1;
                        class_matrix(j,i) = 1;
                    else
                        class_matrix(j,i) = 0;
                    end
                end
            end
            [v_class,Iv,~] = unique(class_matrix,'rows');
            class_store = zeros(length(Iv),sample_size);
            for i = 1:length(Iv)
                store_tem = [];
                for j = 1:sample_size
                    if class_matrix(j,:)==v_class(i,:)
                        store_tem = [store_tem, j];
                    end
                end
                class_store(i,(1:length(store_tem))) = store_tem;
            end
            class_store(class_store==0)= nan;
            index_store = [(1:sample_size)',zeros(sample_size,1)];
            for i = 1:sample_size
                for j = 1:length(Iv)
                    if class_matrix(i,:) == v_class(j,:)
                        index_store(i,2)= j;
                        break
                    end
                end
            end
            %class_store %存储每组下的原始序号（每一行代表一组，元素为原始序号）
            %index_store %存储原始序号对应的组标签（第一列为原始序号，第二列为对应的组标签）
            [row, ~] = size(class_store);
            Class_iter = struct();
            Class_iter.subgroup_number = row;
            Class_iter.index = cell(1,row);
            Class_iter.coef = cell(1,row);
            Class_iter.subgroup_sample_size = cell(1,row);
            Class_iter.X = cell(1,row);
            Class_iter.y = cell(1,row);
            for j = 1 : row
                Class_iter.index{1,j} = rmmissing(class_store(j,:));
                Class_iter.subgroup_sample_size{1,j} = length(rmmissing(class_store(j,:)));
                Class_iter.y{1,j} = y(Class_iter.index{1,j},1);
            end
            for j = 1 : row
                Class_iter.coef{1,j} = zeros(p, 1);
                for r = 1:Class_iter.subgroup_sample_size{1,j}
                    Class_iter.coef{1,j} = Class_iter.coef{1,j} + beta(((Class_iter.index{1,j}(r)-1)*p+1):((Class_iter.index{1,j}(r))*p));
                end
                Class_iter.coef{1,j} = (1/Class_iter.subgroup_sample_size{1,j})*Class_iter.coef{1,j};
            end
            % % % for j = 1 : row
            % % %     Class_iter.coef{1,j} = zeros(p, 1);
            % % %     panduan = zeros(p, Class_iter.subgroup_sample_size{1,j});
            % % %     for r = 1:Class_iter.subgroup_sample_size{1,j}
            % % %         %Class_iter.coef{1,j} = Class_iter.coef{1,j} + beta(((Class_iter.index{1,j}(r)-1)*p+1):((Class_iter.index{1,j}(r))*p));
            % % %         panduan(:, r) = beta0(((Class_iter.index{1,j}(r)-1)*p+1):((Class_iter.index{1,j}(r))*p));
            % % %     end
            % % %     [~,Iv2,count] = unique(panduan','rows');
            % % %     Class_iter.coef{1,j} = beta0(((Iv2(mode(count))-1)*p+1):(Iv2(mode(count))*p));
            % % %     %Class_iter.coef{1,j} = (1/Class_iter.subgroup_sample_size{1,j})*Class_iter.coef{1,j};
            % % % end
            for j = 1 : row
                Class_iter.X{1,j} = zeros(Class_iter.subgroup_sample_size{1,j}, p);
                for r = 1 : Class_iter.subgroup_sample_size{1,j}
                    Class_iter.X{1,j}(r,:) = X(Class_iter.index{1,j}(1,r),((Class_iter.index{1,j}(1,r)-1)*p+1):(Class_iter.index{1,j}(1,r)*p));
                end
            end
            % for j = 1 : row
            %     Class_iter.coef{1,j} = zeros(p, 1);
            %     xx = Class_iter.X{1,j} - mean(Class_iter.X{1,j});
            %     yy = Class_iter.y{1,j} - mean(Class_iter.y{1,j});
            %     %options = struct('alpha', 1, 'nlambda', 100, 'nfold', 1); % Lasso回归的设置
            %     %options.lambda = cvglmnet(xx, yy, 'gaussian', options).lambda_min;
            %     %Class_iter.coef{1,j} = glmnet(xx, yy, [], options).beta;
            %     Class_iter.coef{1,j} = (xx'*xx)\(xx'*yy);
            % end

        end



















        function kesai_store = store(obj, Class_iter, kesai)
            %使用时:kesai_store = store(Class_iter, kesai)
            %kesai_store(idx).store
            [~, ~, p] = size(kesai);
            kesai_store = struct();
            kesai_store(length(Class_iter.index)).store = inf*ones(p,1);
            for j = 1:length(Class_iter.index)
                idx = j;
                % idx = find(cellfun(@(c) ismember(element, c), Class_iter.index));
                otherCells = Class_iter.index;
                otherCells(idx) = [];
                otherElements = [otherCells{:}];
                kesai_idx = zeros(p,1);
                for i = 1:length(Class_iter.index{idx})
                    element = Class_iter.index{idx}(i);
                    smallerElements = otherElements(otherElements < element);
                    largerElements = otherElements(otherElements > element);
                    kesai_large = squeeze(kesai(element,largerElements,:))';
                    kesai_small = squeeze(kesai(smallerElements,element,:))';
                    kesai_element = sum(kesai_large,2)-sum(kesai_small,2);
                    kesai_idx = kesai_idx + kesai_element;
                end
                kesai_store(idx).store = kesai_idx;
            end
        end






















        function [beta_back,convergence] = admm_part(obj, sample_size,p,iter_max2,eps2,iter_max_in,eps_in,...
                eps_in_in,iter_max_in_in,beta,X,y,a,kesai,Class_iter)

                beta_back = beta;
                for i = 1 : Class_iter.subgroup_number
                    for j = 1: length(Class_iter.index{i})
                        beta_back(((Class_iter.index{i}(j)-1)*p+1):((Class_iter.index{i}(j))*p)) = Class_iter.coef{i};
                        %beta_back(((Class_iter.index{i}(j)-1)*p+1):((Class_iter.index{i}(j))*p)) = ...
                        %   beta0(((Class_iter.index{i}(j)-1)*p+1):((Class_iter.index{i}(j))*p));
                    end
                end
                %v2 = beta_back;
                iter = 0;
                residual_primal = 1;
                beta_back_old = beta_back;
                %kesai_store = obj.store(Class_iter, kesai); %更新完分组后重新计算储存的kesai
                while ((residual_primal > eps2)&&(iter<=iter_max2))
                    kesai_store = obj.store(Class_iter, kesai); %更新完分组后重新计算储存的kesai
                    for i = 1 : Class_iter.subgroup_number
                        %iter_in = 0;
                        %residual_primal_in = 1;
                        %residual_dual_in = 1;
                        %v2_old = v2;
                        %%%k = 1;
                        kesai_group = kesai_store(i).store;
                        %beta_group = beta_back(((Class_iter.index{i}(1)-1)*p+1):((Class_iter.index{i}(1))*p));
                        %v2_group = v2(((Class_iter.index{i}(1)-1)*p+1):((Class_iter.index{i}(1))*p));
                        %kesai2_group = kesai2(((Class_iter.index{i}(1)-1)*p+1):((Class_iter.index{i}(1))*p));
                        %v2_old_group = v2_old(((Class_iter.index{i}(1)-1)*p+1):((Class_iter.index{i}(1))*p));
                        % while ((residual_primal_in > eps_in)&&(iter_in <= iter_max_in))
                            % %-----------------------S-AMA第一步------------------------------%
                            %beta_group = ((Class_iter.X{i}'*Class_iter.X{i}) + k * eye(p))\...
                            %    (Class_iter.X{i}'*Class_iter.y{i} - kesai_group + k * v2_group - kesai2_group);
                            %beta_group = ((Class_iter.X{i}'*Class_iter.X{i}))\...
                            %   (Class_iter.X{i}'*Class_iter.y{i} - kesai_group - kesai2_group);
                            beta_group = ((Class_iter.X{i}'*Class_iter.X{i}))\...
                               (Class_iter.X{i}'*Class_iter.y{i} - kesai_group);

                            % % beta_back(((i-1)*p+1):(i*p)) = (X(i,((i-1)*p+1):(i*p))'*X(i,((i-1)*p+1):(i*p)) + k*eye(p))\...
                            % %     (X(i,((i-1)*p+1):(i*p))'*y(i,1) + k*v2(((i-1)*p+1):(i*p),1) - kesai2(((i-1)*p+1):(i*p),1));
                            % % %-----------------------S-AMA第二步------------------------------%
                            % % 第i个分量的索引((i-1)*p+1):(i*p)
                            % v2_bar_group = beta_group + (1/k)*kesai2_group;
                            % %v2_bar(((i-1)*p+1):(i*p)) = beta_back(((i-1)*p+1):(i*p))+(1/k)*kesai2(((i-1)*p+1):(i*p));%%
                            % %[~, col] = size(W.coef);
                            % [~, which_aux] = feval(@(x) min(x(1, 2:end)), squareform(pdist([v2_bar_group'; W.coef'])));
                            % v2_bar_norm = norm(v2_bar_group-W.coef(:,which_aux));
                            % if v2_bar_norm >= (a*lambda_2)
                            %     v2_group = v2_bar_group;
                            % else
                            %     v2_group = W.coef(:,which_aux) +...
                            %         max(0,(1-(lambda_2*Class_iter.subgroup_sample_size{i})/(k*v2_bar_norm)))*(1-Class_iter.subgroup_sample_size{i}*(1/(a*k)))*(v2_bar_group-W.coef(:,which_aux));
                            % end
                            % kesai2_group = kesai2_group +...
                            %     k*(beta_group-v2_group);
                        %     residual_primal_in = norm((beta_group-v2_group)/sqrt(p));
                        %     residual_dual_in = norm(k * (v2_group - v2_old_group)/sqrt(p));
                        %     if residual_primal_in > 10 * residual_dual_in
                        %         k = k * 2;
                        %     elseif residual_dual_in > 10 * residual_primal_in
                        %         k = k / 2;
                        %     end
                        %     v2_old_group = v2_group;
                        %     iter_in = iter_in + 1;
                        % end
                        for j = 1: length(Class_iter.index{i})
                            beta_back(((Class_iter.index{i}(j)-1)*p+1):((Class_iter.index{i}(j))*p)) = beta_group;
                            %v2(((Class_iter.index{i}(j)-1)*p+1):((Class_iter.index{i}(j))*p)) = v2_group;
                            %kesai2(((Class_iter.index{i}(j)-1)*p+1):((Class_iter.index{i}(j))*p)) = kesai2_group;
                            %v2_old(((Class_iter.index{i}(j)-1)*p+1):((Class_iter.index{i}(j))*p)) =...
                                %v2(((Class_iter.index{i}(j)-1)*p+1):((Class_iter.index{i}(j))*p));
                        end
                        Class_iter.coef{i} = beta_group;
                    end


                    beta_back_old_in = beta_back;
                    residual_in_in = 1;
                    iter_in_in = 0;
                    while ((residual_in_in>eps_in_in) && (iter_in_in < iter_max_in_in))
                        for t = 1:sample_size
                            panduan = inf*ones(1, Class_iter.subgroup_number);
                            %f_store = obj.store2(Class_iter, v, kesai, t);
                            which_old = find(cellfun(@(c) ismember(t, c), Class_iter.index));
                            for i = 1 : Class_iter.subgroup_number
                                %[~, which_aux] = feval(@(x) min(x(1, 2:end)), squareform(pdist([Class_iter.coef{i}'; W.coef'])));
                                %v2_bar_norm = norm(beta_back(((t-1)*p+1):(t*p))-W.coef(:,which_aux));
                                panduan(1, i) = 0.5*(y(t) - X(t,((t-1)*p+1):(t*p))*Class_iter.coef{i}).^2;
                                %+ mcp_function(norm(Class_iter.coef{i}-W.coef(:,which_aux)),a,lambda_2);
                                %+ f_store(i).store;
                            end
                            [~, which] = min(panduan(1, :));
                            beta_back((((t-1)*p+1):(t*p))) = Class_iter.coef{which};
                            if which_old ~= which
                                Class_iter.index{which} = sort([Class_iter.index{which},t]);
                                Class_iter.index{which_old}(Class_iter.index{which_old} == t) = [];
                            end

                            Class_iter.subgroup_number = sum(double(~cellfun(@isempty, Class_iter.index)));
                            Class_iter.coef(cellfun(@isempty, Class_iter.index)) = [];
                            Class_iter.index(cellfun(@isempty, Class_iter.index)) = [];
                            Class_iter.subgroup_number_size = cell(1, Class_iter.subgroup_number);
                            Class_iter.X = cell(1,Class_iter.subgroup_number);
                            Class_iter.y = cell(1,Class_iter.subgroup_number);
                            for i = 1 : Class_iter.subgroup_number
                                Class_iter.subgroup_number_size{i} = length(Class_iter.index{i});
                                Class_iter.X{i} = zeros(Class_iter.subgroup_number_size{i},p);
                                Class_iter.y{i} = y(Class_iter.index{i},1);
                                for j = 1 : Class_iter.subgroup_number_size{i}
                                    Class_iter.X{i}(j,:) = X(Class_iter.index{i}(j),(((Class_iter.index{i}(j)-1)*p)+1):(Class_iter.index{i}(j)*p));
                                end
                            end
                        end
                        residual_in_in = norm((beta_back-beta_back_old_in)/sqrt(sample_size*p));
                        beta_back_old_in = beta_back;
                        iter_in_in = iter_in_in + 1;
                    end

                    
                    %v2 = beta_back;
                    %kesai2 = zeros(sample_size*p, 1);
                    residual_primal = norm((beta_back-beta_back_old)/sqrt(sample_size*p));
                    beta_back_old = beta_back;
                    iter = iter + 1;
                end
            convergence = residual_primal;
        end














        function [Data] = admm_mex(obj, sample_size, p, iter_max1, eps1,...
                beta0, H_p, X, y, k, a, lambda, c, v, kesai,...
                v_vec, v_vec_bar, kesai_vec, beta_real, HH_p, E)
            [Data] = admm_mex(sample_size, p, iter_max1, eps1,...
                beta0, H_p, X, y, k, a, lambda, c, v, kesai,...
                v_vec, v_vec_bar, kesai_vec, beta_real, HH_p, E);
        end

        function residual_primal = Residual_primal_mex(obj,full_if,E,beta,v_vec,sample_size,p,H_p)
            residual_primal = Residual_primal_mex(full_if,E,beta,v_vec,sample_size,p,H_p);
        end















        function [lambda_left_back,lambda_right_back] = interval_lambda(obj,lambda_start_point,change_constant,sample_size,...
                p,iter_max1_interval,eps1_interval,beta0,H_p,X,y,k,a,c,v,kesai,v_vec,v_vec_bar,kesai_vec,...
                beta_real,min_class_num,HH_p,E,iter_max2,eps2,iter_max_in,eps_in,eps_in_in,iter_max_in_in,full_if)
            %-------------------------------------------函数功能-----------------------------------------------%
            % 使用admm 求解正则化路径的lambda左右端点
            %----------------------------------------需要的前置函数---------------------------------------------%
            % MATLAB Function:      admm
            %-----------------------------------------输出变量说明---------------------------------------------%
            % lambda_left_back      惩罚因子区间左端点 double
            % lambda_right_back     惩罚因子区间右端点 double
            %-----------------------------------------输入变量说明---------------------------------------------%
            % lambda_start_point    lambda搜索起始点
            % change_constant       lambda因子变动系数
            % sample_size           样本量
            % p                     coefficient dimension
            % iter_max1_interval    区间搜索时ADMM最大交替迭代步数
            % eps1_interval         区间搜索时ADMM原始残差二范数的容差
            % beta0                 初值beta0
            % H_p                   H_p维度为（组合数*p）*(样本数*p) matrix
            % X                     解释变量数据 （样本数*p）*(样本数*q) matrix
            % y                     连续相应变量 （样本数*1）matrix
            % a                     MCP惩罚函数里的控制凹凸性的参数
            % k                     ADMM算法的惩罚参数
            % v                     beta组合差信息储存矩阵（对偶变量）
            % c                     BIC准则第二部分的权重参数
            % kesai　　　　　　　　 关于beta的增广lagrange乘子
            % v_vec                 beta组合差信息储存矩阵（对偶变量）按组合顺序拉成向量
            % kesai_vec             beta的增广lagrange乘子按组合顺序拉成向量
            % store_beta            内部最小二乘迭代时更新beta时用到的信息
            % beta_real             真实的beta
            % HH_p                  H_p'*H_p
            % min_class_num         lambda右端点对应的可接受的最小亚组数
            for i = 1:(sample_size-1)%%%%%%%
                v(i,i,:) = 0;%%%%%%%
                for j = i+1:sample_size%%%%%%%
                    v(i,j,:) = beta0(((i-1)*p+1):(i*p)) - beta0(((j-1)*p+1):(j*p));%%%%%%%
                    v(j,i,:) = v(i,j,:);%%%%%%%
                end%%%%%%%
            end%%%%%%%
            [Class_iter,~,~,~,~] = obj.Cluster_iter_pro(sample_size, p, v, X, y, beta0, beta0);%%%%%%%
            left_index = Class_iter.subgroup_number;%%%%%%%
            lambda_left = lambda_start_point;
            lambda_right = lambda_start_point;
            iter_max1 = iter_max1_interval;
            eps1 = eps1_interval;
            subgroup_number_start = obj.admm(sample_size, p, iter_max1, eps1,...
                beta0, H_p, X, y, k, a, lambda_start_point, c, v, kesai,...
                v_vec, v_vec_bar, kesai_vec, beta_real, HH_p, E,iter_max2,eps2,iter_max_in,eps_in,...
                    eps_in_in,iter_max_in_in,full_if).subgroup_number;
            subgroup_number_left = subgroup_number_start;
            subgroup_number_right = subgroup_number_start;
            while subgroup_number_left==sample_size
                lambda_left = lambda_left/change_constant;
                subgroup_number_left = obj.admm(sample_size, p, iter_max1, eps1,...
                    beta0, H_p, X, y, k, a, lambda_left, c, v, kesai,...
                    v_vec, v_vec_bar, kesai_vec, beta_real, HH_p, E, iter_max2,eps2,iter_max_in,eps_in,...
                    eps_in_in,iter_max_in_in,full_if).subgroup_number;
            end
            lambda_left = lambda_left*change_constant;
            while subgroup_number_right == 1
                lambda_right = lambda_right*change_constant;
                subgroup_number_right = obj.admm(sample_size, p, iter_max1, eps1,...
                    beta0, H_p, X, y, k, a, lambda_right, c, v, kesai,...
                    v_vec, v_vec_bar, kesai_vec, beta_real, HH_p, E, iter_max2,eps2,iter_max_in,eps_in,...
                    eps_in_in,iter_max_in_in,full_if).subgroup_number;
            end
            lambda_right = lambda_right/change_constant;
            while subgroup_number_left < left_index
                lambda_left = lambda_left*change_constant;
                subgroup_number_left = obj.admm(sample_size, p, iter_max1, eps1,...
                    beta0, H_p, X, y, k, a, lambda_left, c, v, kesai,...
                    v_vec, v_vec_bar, kesai_vec, beta_real, HH_p, E, iter_max2,eps2,iter_max_in,eps_in,...
                    eps_in_in,iter_max_in_in,full_if).subgroup_number;
                if lambda_left < 0.001
                    lambda_left = 0.001;
                    subgroup_number_left = left_index;
                end
            end
            change_constant_mid = change_constant;
            while subgroup_number_right > min_class_num
                lambda_right = lambda_right/change_constant_mid;
                subgroup_number_right = obj.admm(sample_size, p, iter_max1, eps1,...
                    beta0, H_p, X, y, k, a, lambda_right, c, v, kesai,...
                    v_vec, v_vec_bar, kesai_vec, beta_real, HH_p, E, iter_max2,eps2,iter_max_in,eps_in,...
                    eps_in_in,iter_max_in_in,full_if).subgroup_number;
                % if subgroup_number_right < min_class_num
                %     lambda_right = lambda_right*change_constant_mid;
                %     %change_constant_mid = sqrt(change_constant_mid);
                %     %lambda_right = lambda_right/change_constant_mid;
                %     subgroup_number_right = obj.admm_mex(sample_size, p, iter_max1, eps1,...
                %         beta0, H_p, X, y, k, a, lambda_right, c, v, kesai,...
                %         v_vec, v_vec_bar, kesai_vec, beta_real, HH_p, E).subgroup_number;
                %     if subgroup_number_right > min_class_num
                %         change_constant_mid = sqrt(change_constant_mid);
                %         lambda_right = lambda_right/change_constant_mid;
                %         subgroup_number_right = obj.admm_mex(sample_size, p, iter_max1, eps1,...
                %         beta0, H_p, X, y, k, a, lambda_right, c, v, kesai,...
                %         v_vec, v_vec_bar, kesai_vec, beta_real, HH_p, E).subgroup_number;
                %     end
                % elseif subgroup_number_right > min_class_num
                %     lambda_right = lambda_right/change_constant_mid;
                %     subgroup_number_right = obj.admm_mex(sample_size, p, iter_max1, eps1,...
                %         beta0, H_p, X, y, k, a, lambda_right, c, v, kesai,...
                %         v_vec, v_vec_bar, kesai_vec, beta_real, HH_p, E).subgroup_number;
                %     if subgroup_number_right < min_class_num
                %         change_constant_mid = sqrt(change_constant_mid);
                %         lambda_right = lambda_right*change_constant_mid;
                %         subgroup_number_right = obj.admm_mex(sample_size, p, iter_max1, eps1,...
                %         beta0, H_p, X, y, k, a, lambda_right, c, v, kesai,...
                %         v_vec, v_vec_bar, kesai_vec, beta_real, HH_p, E).subgroup_number;
                %     end
            end

            
            %lambda_right = lambda_right*change_constant;%%%%%%%%%%%%%%%%%%%%%%%%
            lambda_left_back = lambda_left;
            lambda_right_back = lambda_right;
        end




















        function [DATA_struct,data,DATA_opt,data_opt] = admm_lambda(obj,ifauto,lambda_size,lambda_start_point,change_constant,...
                beta0,sample_size,p,iter_max1,iter_max1_interval,eps1,eps1_interval,H_p,X,y,k,a,...
                lambda_left,lambda_right,c,v,kesai,v_vec,v_vec_bar,kesai_vec,beta_real,min_class_num,HH_p,E,iter_max2,eps2,iter_max_in,eps_in,...
                    eps_in_in,iter_max_in_in,full_if)
            %-------------------------------------------函数功能-----------------------------------------------%
            % 使用interval_lambda、admm求解一次模拟中最优lambda对应的估计信息
            %----------------------------------------需要的前置函数---------------------------------------------%
            % MATLAB Function:      admm、interval_lambda
            %-----------------------------------------输出变量说明---------------------------------------------%
            % DATA_struct      　　　一次模拟中每个lambda对应的详细估计信息
            % data                  一次模拟中每个lambda对应的简略估计信息
            % DATA_opt              一次模拟中最优lambda对应的详细估计信息
            % data_opt              一次模拟中最优lambda对应的简略估计信息
            %-----------------------------------------输入变量说明---------------------------------------------%
            % ifauto                是否使用interval_lambda自动搜索lambda,1为使用，0为不使用
            % lambda_left           手动输入lambda区间左端点
            % lambda_right          手动输入lambda区间右端点
            % lambda_start_point    lambda搜索起始点
            % change_constant       lambda因子变动系数
            % log_change_constant   lambda区间对数变换系数
            % sample_size           样本量
            % p                     coefficien dimension
            % iter_max1_interval    区间搜索时ADMM最大交替迭代步数
            % eps1_interval         区间搜索时ADMM原始残差二范数的容差
            % beta0                 初值beta0
            % H_p                   H_p维度为（组合数）*(样本数*p) matrix
            % X                     解释变量数据 （样本数）*(样本数*p) matrix
            % y                     连续相应变量 （样本数*1）matrix
            % a                     MCP惩罚函数里的控制凹凸性的参数
            % k                     ADMM算法的惩罚参数
            % v                     beta组合差信息储存矩阵（对偶变量）
            % c                     BIC准则第二部分的权重参数
            % kesai　　　　　　　　 关于beta的增广lagrange乘子
            % v_vec                 beta组合差信息储存矩阵（对偶变量）按组合顺序拉成向量
            % kesai_vec             beta的增广lagrange乘子按组合顺序拉成向量
            % beta_real             真实的beta
            % HH_p                  H_p'*H_p
            % min_class_num         lambda右端点对应的可接受的最小亚组数
            %----------------------------------------是否选择使用自动搜索lambda区间-------------------------------------------%
            if ifauto == 1
                [lambda_left_back,lambda_right_back] = obj.interval_lambda(lambda_start_point,change_constant,sample_size,...
                    p,iter_max1_interval,eps1_interval,beta0,H_p,X,y,k,a,c,v,kesai,v_vec,v_vec_bar,kesai_vec,...
                    beta_real,min_class_num,HH_p,E, iter_max2,eps2,iter_max_in,eps_in,...
                    eps_in_in,iter_max_in_in,full_if);
                Lambda = linspace(lambda_left_back,lambda_right_back,lambda_size);
            else
                lambda_left_back = lambda_left;
                lambda_right_back = lambda_right;
                Lambda = linspace(lambda_left_back,lambda_right_back,lambda_size);
            end
            disp("lambda区间搜索完成")
            %----------------------------------------设置结构体储存每个lambda对应的结果----------------------------------------%
            DATA = struct;
            DATA(lambda_size).lambda = [];DATA(lambda_size).beta= [];DATA(lambda_size).class= [];DATA(lambda_size).index= [];
            DATA(lambda_size).subgroup_number = [];DATA(lambda_size).convergence = [];DATA(lambda_size).BIC= [];
            DATA(lambda_size).BIC_part1= [];DATA(lambda_size).BIC_part2= [];DATA(lambda_size).iter = [];
            DATA(lambda_size).per_if = [];DATA(lambda_size).mse_beta = [];DATA(lambda_size).sc = [];DATA(lambda_size).Data = [];
            for i = 1:length(Lambda)
                Data = obj.admm(sample_size,p,iter_max1,eps1,...
                    beta0,H_p,X,y,k,a,Lambda(i),c,v,kesai,v_vec, v_vec_bar,kesai_vec,...
                    beta_real,HH_p,E, iter_max2,eps2,iter_max_in,eps_in,...
                    eps_in_in,iter_max_in_in,full_if);

                DATA(i).lambda = Lambda(i);
                DATA(i).beta= Data.beta;
                DATA(i).class= Data.class;
                DATA(i).index= Data.index;
                DATA(i).subgroup_number = Data.subgroup_number;
                DATA(i).per_if = Data.per_if;
                DATA(i).mse_beta = Data.mse_beta;
                DATA(i).sc = Data.sc;
                DATA(i).convergence = Data.convergence;
                DATA(i).iter = Data.iter;
                DATA(i).BIC= Data.BIC;
                DATA(i).BIC_part1 = Data.BIC_part1;
                DATA(i).BIC_part2 = Data.BIC_part2;
                DATA(i).Data = Data;
                fprintf('已完成一次模拟中 %.1f%% 的融合惩罚\n',(100*i/lambda_size))
            end
            %------------------------------------每次模拟的详细信息在DATA_struct中--------------------------------%
            DATA_struct = DATA;
            %------------------------------------每次模拟的简略信息在data中---------------------------------------%
            lambda = [DATA_struct.lambda];
            BIC = [DATA_struct.BIC];
            BIC_part1 = [DATA_struct.BIC_part1];
            BIC_part2 = [DATA_struct.BIC_part2];
            mse_beta = [DATA_struct.mse_beta];
            subgroup_number = [DATA_struct.subgroup_number];
            per_if =[DATA_struct.per_if];
            sc = [DATA_struct.sc];
            convergence =[DATA_struct.convergence];
            iter = [DATA_struct.iter];
            data = [lambda',BIC',BIC_part1',BIC_part2',mse_beta',subgroup_number',per_if',sc',convergence',iter'];
            data = array2table(data, 'VariableNames', {'lambda','BIC','BIC_part1','BIC_part2','mse_beta','subgroup_number','per_if','sc','convergence','iter'});
            DATA_opt = DATA_struct(min(find(BIC == min(BIC))));
            lambda_opt = [DATA_struct(min(find(BIC == min(BIC)))).lambda];
            BIC_opt = [DATA_opt.BIC];
            BIC_part1_opt = [DATA_opt.BIC_part1];
            BIC_part2_opt = [DATA_opt.BIC_part2];
            mse_beta_opt = [DATA_opt.mse_beta];
            subgroup_number_opt = [DATA_opt.subgroup_number];
            per_if_opt =[DATA_opt.per_if];
            sc_opt = [DATA_opt.sc];
            convergence_opt =[DATA_opt.convergence];
            iter_opt = [DATA_opt.iter];
            data_opt = [lambda_opt,BIC_opt,BIC_part1_opt,BIC_part2_opt,mse_beta_opt,...
                subgroup_number_opt,per_if_opt,sc_opt,convergence_opt,iter_opt];
            data_opt = array2table(data_opt, 'VariableNames', {'lambda','BIC','BIC_part1','BIC_part2','mse_beta','subgroup_number','per_if','sc','convergence','iter'});
        end




















        function [Results,results,Results_opt,results_opt,...
                initial,Results_list_opt,Results_single_opt,Class_summary] = Auxiliary_heterogeneous_regression(obj,simulation_size,simulation_index,Data_x,y_real,sample_size,p,beta_real,...
                num_partitions_list,eps_out,iter_max_initial_out,combine_principle_list,split_principle_list,...
                random_number_list,...
                ifauto,lambda_size,lambda_start_point,change_constant,log_change_constant,iter_max1,iter_max1_interval,...
                eps1,eps1_interval,k,a,lambda_left,lambda_right,c,min_class_num,disp_if,initial_if,initial_position, high_if, CV_number, K_require, K_require_if,...
                iter_max2,eps2,iter_max_in,eps_in, eps_in_in,iter_max_in_in,full_if)

            %
            %-----------------------------------输出参数------------------------------------------%
            % Results:                  每次模拟的DATA(每个lambda的详细信息)
            % results:                  每次模拟的data(每个lambda的简略信息)
            % Results_opt:              每次模拟的DATA_opt(最优lambda对应的详细信息)
            % results_opt:              每次模拟的data_opt(最优lambda对应的简略信息)
            % Results_list_opt:         每次模拟的最优lambda对应的简略信息拼成的table
            % Results_single_opt:       综合每次模拟最优lambda对应的评价信息获得的table
            % initial:                  每次模拟初值搜索算法的信息
            %----------------------------------使用方法-------------------------------------------%
            %模拟时必须输入的参数有simulation_size、Data_x、y_real、sample_size、p、beta_real
            % % [Results,results,Results_opt,results_opt,...
            % %              initial,Results_list_opt,Results_single_opt,Class_summary] =  Aux.Auxiliary_heterogeneous_regression(simulation_size,[],X_target,Y_target,sample_size,p,beta_target_real,...
            % %                                                                           [],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[])
            %实际数据时必须输入的参数有simulation_size、Data_x、y_real、sample_size、p
            % 解决多次模拟和单次模拟的维度不一致
            if isempty(simulation_index)
                simulation_index = 1:simulation_size;
            end
            if length(simulation_index) ~= simulation_size
                Data_x = Data_x(1:sample_size,1:p,simulation_index);
                y_real = y_real(1:sample_size,simulation_index);
            end
            if isempty(num_partitions_list)
                %num_partitions_list = min(10, floor(sqrt(sample_size))); %默认初值搜索分割组数区间10
                num_partitions_list =  5:10; %默认初值搜索分割组数区间10
            end
            if isempty(eps_out)
                eps_out = 1e-3; %默认初值算法收敛条件1e-3
            end
            if isempty(iter_max_initial_out)
                iter_max_initial_out = 50; %默认迭代最大步数10
            end
            if isempty(combine_principle_list)
                combine_principle_list = linspace(0.1, 0.1, 1)*sample_size;
            end
            if isempty(split_principle_list)
                split_principle_list = linspace(1, 1, 1)*sample_size;
            end
            if isempty(random_number_list)
                random_number_list = 51:100; %默认随机种子数列表1:100
            end
            if isempty(ifauto)
                ifauto = 1; %默认采取自动搜索区间 1
            end
            if isempty(lambda_size)
                lambda_size = 10; %默认lambda因子个数 25
            end
            if isempty(lambda_start_point)
                lambda_start_point = sqrt(p*log(p)/sample_size); %%lambda因子搜索初始点
            end
            if isempty(change_constant)
                change_constant = 0.5; %lambda因子区间搜索变动系数为0.9
            end
            if isempty(log_change_constant)
                log_change_constant = 1; %lambda因子区间log-exp变换系数
            end
            if isempty(iter_max1)
                iter_max1 = 50; %默认ADMM最大交替迭代步数为50
            end
            if isempty(iter_max1_interval)
                iter_max1_interval = 50; %默认区间搜索时ADMM最大交替迭代步数为50
            end
            if isempty(eps1)
                eps1 = 5e-2; %默认ADMM原始残差二范数的容差为1e-3
            end
            if isempty(eps1_interval)
                eps1_interval = 5e-2; %默认区间搜索时ADMM原始残差二范数的容差为1e-3
            end
            if isempty(beta_real)
                beta_real = zeros(sample_size*p,1); %预测时不输入beta_real时默认的beta_real
            end
            if isempty(k)
                k = 1; %默认ADMM算法的惩罚参数为10
            end
            if isempty(a)
                a = 3; %默认ADMM算法的惩罚参数为3.7
            end
            if isempty(lambda_left)
                lambda_left = 0.005; %手动设置lambda左初始点0.005
            end
            if isempty(lambda_right)
                lambda_right = 3; %手动设置lambda右初始点0.6
            end
            if isempty(c)
                c = 1; %默认BIC准则第二部分的权重参数为2
            end
            if isempty(min_class_num)
                min_class_num = 2; %默认lambda右端点对应的可接受的最小亚组数为2
            end
            if isempty(disp_if)
                disp_if = 1;%默认不显示过程信息
            end
            if isempty(initial_if)
                initial_if = 0; %1:不重新计算初值
            end
            if isempty(initial_position)
                initial_position = 'C:\Users\张丰川\Documents\MATLAB\Transfer-learning-heterogeneity-analysis\Result\Initial_auxiliary\'; %初值结果存放地址
            end
            position = strcat(initial_position, '', 'Result_table.mat');
            if isempty(high_if)
                high_if = 0;
            end
            if isempty(CV_number)
                CV_number = 5;
            end
            if isempty(K_require)
                K_require = floor(sqrt(sample_size));
            end
            if isempty(K_require_if)
                K_require_if = 0;
            end
            if isempty(iter_max2)
                iter_max2 = 50;%%%%%25
            end
            if isempty(eps2)
                eps2 = 1e-3;%%%%%%%%1e-3
            end
            if isempty(iter_max_in)
                iter_max_in = 50;%%%%%%%%%%%
            end
            if isempty(eps_in)
                eps_in = 1e-3;
            end
            if isempty(eps_in_in)
                eps_in_in = 1e-3;
            end
            if isempty(iter_max_in_in)
                iter_max_in_in = 50; %50
            end
            if isempty(full_if)
                full_if = 0;
            end
            Results = struct;
            
            Results(length(simulation_index)).DATA = [];
            results = struct;
            results(length(simulation_index)).data = [];
            Results_opt = struct;%
            Results_opt(length(simulation_index)).DATA = [];
            results_opt = struct;%
            results_opt(length(simulation_index)).data = [];
            Results_list_opt = [];
            Class_summary = struct;
            Class_summary(length(simulation_index)).subgroup_number = [];
            Class_summary(length(simulation_index)).coef = [];
            Class_summary(length(simulation_index)).subgroup_sample_size = [];
            Class_summary(length(simulation_index)).index = [];
            Class_summary(length(simulation_index)).X = [];
            Class_summary(length(simulation_index)).y = [];
            initial = struct;
            initial(length(simulation_index)).beta0 = [];
            initial(length(simulation_index)).partitions = [];
            initial(length(simulation_index)).result = [];
            simulation_size = length(simulation_index);%%%%%%%%%%%%%%%%
            if length(simulation_index) > 1
                parfor local_i = 1:simulation_size  %!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                    %local_i = yangbencishu;
                    obj_initial = Initial_regression();
                    if initial_if == 0
                        y = y_real(:,local_i);%-mean(y_real(:,local_i));%将第yangbencishu次模拟中心化后的y抽出来记为y
                        x = Data_x(:,:,local_i);%-mean(Data_x(:,:,local_i));%将第yangbencishu次模拟的矩阵变量抽出来记为x（维度为n*p）
                        %[beta0,gamma0,coefficients,partitions,combine,residual_coefficients,...
                        %num_partitions_best,combine_principle_best,split_principle_best,...
                        %random_number_best] = k_means_regression_initial(num_partitions_list,eps_initial,eps_out,iter_max_initial_in,...
                        %                                            iter_max_initial_out,combine_principle_list,split_principle_list,random_number_list,sample_size,row_size,col_size,x,y,disp_if);
                        [beta0, partitions, result] = obj_initial.k_means_regression_initial(num_partitions_list, eps_out,...
                            iter_max_initial_out, combine_principle_list,...
                            split_principle_list, random_number_list,...
                            sample_size, p, x, y, disp_if, beta_real, c, high_if, CV_number, K_require, K_require_if);
                        %Partitions(local_i).partitions = partitions;
                        initial(local_i).beta0 = beta0;
                        initial(local_i).partitions = partitions;
                        initial(local_i).result = result;
                    else
                        initial_read = load(position);
                        beta0 = cell2mat(initial_read.Result_table.beta_back(local_i));
                        %paritions = initial_read.Result_table.partitions(local_i);
                        %result = initial_read.Result_table.per_if(local_i);
                        initial(local_i).beta0 = beta0;
                        %initial(local_i).partitions = partitions;
                        %initial(local_i).result = result;
                    end
                    %-----------------------------------------------------------------------------------------------------%
                    %-----------------------------------------------------------------------------------------------------%
                    y = y_real(:,local_i)-mean(y_real(:,local_i));%将第1次模拟中心化后的y抽出来记为y
                    x = Data_x(:,:,local_i)-mean(Data_x(:,:,local_i));%将第1次模拟的矩阵变量抽出来记为x（维度为n*p）
                    %x = x - mean(x,3);
                    beta = beta0;
                    %for i = 1:sample_size %按照可识别性条件对参数进行处理
                    %    identify = sign(gamma(((i-1)*col_size+1)))/sqrt(sum(gamma(((i-1)*col_size+1):(i*col_size)).^2));
                    %    beta(((i-1)*row_size+1):(i*row_size)) = beta(((i-1)*row_size+1):(i*row_size))/identify;
                    %    gamma(((i-1)*col_size+1):(i*col_size)) = gamma(((i-1)*col_size+1):(i*col_size))*identify;
                    %end
                    v = zeros(sample_size, sample_size, p)*nan;
                    for i = 1:(sample_size-1)
                        for j = i+1:sample_size
                            v(i,j,:) = beta(((i-1)*p+1):(i*p)) - beta(((j-1)*p+1):(j*p));
                        end
                    end
                    v_vec = zeros(sample_size*(sample_size-1)*0.5*p,1);
                    v_vec_bar = zeros(sample_size*(sample_size-1)*0.5*p,1);
                    for i = 1:(sample_size-1)
                        for j = (i+1):(sample_size)
                            v_vec((((((2*sample_size-i)*(i-1)/2+(j-i))-1)*p)+1):((((2*sample_size-i)*(i-1)/2)+(j-i))*p)) = v(i,j,:);
                        end
                    end
                    kesai = zeros(sample_size, sample_size, p);
                    kesai_vec = zeros((sample_size*(sample_size-1)*0.5*p), 1);
                    for i = 1:(sample_size-1)
                        for j = (i+1):sample_size
                            kesai_vec((((((2*sample_size-i)*(i-1)/2+(j-i))-1)*p)+1):((((2*sample_size-i)*(i-1)/2)+(j-i))*p)) = kesai(i,j,:);
                        end
                    end
                    E = zeros((sample_size*(sample_size-1)*0.5),sample_size);
                    for i = 1:(sample_size-1)
                        for j = (i+1):sample_size
                            e_i = zeros(sample_size,1);
                            e_i(i) = 1;
                            e_j = zeros(sample_size,1);
                            e_j(j) = 1;
                            E((((2*sample_size-i)*(i-1)/2)+(j-i)),:) = (e_i - e_j)';
                        end
                    end
                    E = sparse(E);
                    if full_if == 1
                        H_p = sparse(zeros(2));
                        HH_p = sparse(zeros(2));
                    else
                        H_p = kron(sparse(E),eye(p));
                        HH_p = kron(sparse(E),eye(p))'*kron(sparse(E),eye(p));
                    end
                    x_cell = mat2cell(x,ones(1,sample_size),p);
                    X = blkdiag(x_cell{:});
                    %v_vector中关于beta第ij系数差的范围(((((2*sample_size-i)*(i-1)/2+(j-i))-1)*row_size)+1):((((2*sample_size-i)*(i-1)/2)+(j-i))*row_size)
                    %w_vector中关于gamma第ij系数差的范围(((((2*sample_size-i)*(i-1)/2+(j-i))-1)*col_size)+1):((((2*sample_size-i)*(i-1)/2)+(j-i))*col_size)
                    %E中第((2*sample_size-i)*(i-1)/2)+(j-i))行是1*sample_size关于ij的指示向量
                    %H_p维度为（p）*(样本数*p)
                    %将sample_size个x[:,:,sample_size]拼成一个大的对角块矩阵X(np*nq)
                    %先按照第三个维度转换成sample_size个元胞
                    %f1 = 0.5*norm(y-D_beta*X*gamma).^2+0.5*k*norm(H_p*beta-v_vec+kesai_vec/k).^2+0.5*k*norm(H_q*gamma-w_vec+yita_vec/k).^2;
                    %-----------------------------------------------------------------------------------------------------%
                    %-----------------------------------------------------------------------------------------------------%
                    [DATA_struct,data,DATA_opt,data_opt] = obj.admm_lambda(ifauto,lambda_size,lambda_start_point,change_constant,...
                        beta0,sample_size,p,iter_max1,iter_max1_interval,eps1,eps1_interval,H_p,X,y,k,a,...
                        lambda_left,lambda_right,c,v,kesai,v_vec,v_vec_bar,kesai_vec,beta_real,min_class_num,HH_p,E,iter_max2,eps2,iter_max_in,eps_in,...
                    eps_in_in,iter_max_in_in,full_if);

                    Results(local_i).DATA = DATA_struct;
                    results(local_i).data = data;
                    Results_opt(local_i).DATA = DATA_opt;
                    results_opt(local_i).data = data_opt;
                    fprintf('已完成第%.1f次的融合惩罚\n',local_i)
                end
            else
                for local_i = 1:simulation_size  %!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                    %local_i = yangbencishu;
                    obj_initial = Initial_regression();
                    if initial_if == 0
                        y = y_real(:,local_i);%-mean(y_real(:,local_i));%将第yangbencishu次模拟中心化后的y抽出来记为y
                        x = Data_x(:,:,local_i);%-mean(Data_x(:,:,local_i));%将第yangbencishu次模拟的矩阵变量抽出来记为x（维度为n*p）
                        %[beta0,gamma0,coefficients,partitions,combine,residual_coefficients,...
                        %num_partitions_best,combine_principle_best,split_principle_best,...
                        %random_number_best] = k_means_regression_initial(num_partitions_list,eps_initial,eps_out,iter_max_initial_in,...
                        %                                            iter_max_initial_out,combine_principle_list,split_principle_list,random_number_list,sample_size,row_size,col_size,x,y,disp_if);
                        [beta0, partitions, result] = obj_initial.k_means_regression_initial(num_partitions_list, eps_out,...
                            iter_max_initial_out, combine_principle_list,...
                            split_principle_list, random_number_list,...
                            sample_size, p, x, y, disp_if, beta_real, c, high_if, CV_number, K_require, K_require_if);
                        %Partitions(local_i).partitions = partitions;
                        initial(local_i).beta0 = beta0;
                        initial(local_i).partitions = partitions;
                        initial(local_i).result = result;
                    else
                        initial_read = load(position);
                        beta0 = cell2mat(initial_read.Result_table.beta_back(local_i));
                        %paritions = initial_read.Result_table.partitions(local_i);
                        %result = initial_read.Result_table.per_if(local_i);
                        initial(local_i).beta0 = beta0;
                        %initial(local_i).partitions = partitions;
                        %initial(local_i).result = result;
                    end
                    %-----------------------------------------------------------------------------------------------------%
                    %-----------------------------------------------------------------------------------------------------%
                    y = y_real(:,local_i)-mean(y_real(:,local_i));%将第1次模拟中心化后的y抽出来记为y
                    x = Data_x(:,:,local_i)-mean(Data_x(:,:,local_i));%将第1次模拟的矩阵变量抽出来记为x（维度为n*p）
                    %x = x - mean(x,3);
                    beta = beta0;
                    %for i = 1:sample_size %按照可识别性条件对参数进行处理
                    %    identify = sign(gamma(((i-1)*col_size+1)))/sqrt(sum(gamma(((i-1)*col_size+1):(i*col_size)).^2));
                    %    beta(((i-1)*row_size+1):(i*row_size)) = beta(((i-1)*row_size+1):(i*row_size))/identify;
                    %    gamma(((i-1)*col_size+1):(i*col_size)) = gamma(((i-1)*col_size+1):(i*col_size))*identify;
                    %end
                    v = zeros(sample_size, sample_size, p)*nan;
                    for i = 1:(sample_size-1)
                        for j = i+1:sample_size
                            v(i,j,:) = beta(((i-1)*p+1):(i*p)) - beta(((j-1)*p+1):(j*p));
                        end
                    end
                    v_vec = zeros(sample_size*(sample_size-1)*0.5*p,1);
                    v_vec_bar = zeros(sample_size*(sample_size-1)*0.5*p,1);
                    for i = 1:(sample_size-1)
                        for j = (i+1):(sample_size)
                            v_vec((((((2*sample_size-i)*(i-1)/2+(j-i))-1)*p)+1):((((2*sample_size-i)*(i-1)/2)+(j-i))*p)) = v(i,j,:);
                        end
                    end
                    kesai = zeros(sample_size, sample_size, p);
                    kesai_vec = zeros((sample_size*(sample_size-1)*0.5*p), 1);
                    for i = 1:(sample_size-1)
                        for j = (i+1):sample_size
                            kesai_vec((((((2*sample_size-i)*(i-1)/2+(j-i))-1)*p)+1):((((2*sample_size-i)*(i-1)/2)+(j-i))*p)) = kesai(i,j,:);
                        end
                    end
                    E = zeros((sample_size*(sample_size-1)*0.5),sample_size);
                    for i = 1:(sample_size-1)
                        for j = (i+1):sample_size
                            e_i = zeros(sample_size,1);
                            e_i(i) = 1;
                            e_j = zeros(sample_size,1);
                            e_j(j) = 1;
                            E((((2*sample_size-i)*(i-1)/2)+(j-i)),:) = (e_i - e_j)';
                        end
                    end
                    E = sparse(E);
                    if full_if == 1
                        H_p = sparse(zeros(2));
                        HH_p = sparse(zeros(2));
                    else
                        H_p = kron(sparse(E),eye(p));
                        HH_p = kron(sparse(E),eye(p))'*kron(sparse(E),eye(p));
                    end
                    x_cell = mat2cell(x,ones(1,sample_size),p);
                    X = blkdiag(x_cell{:});
                    %v_vector中关于beta第ij系数差的范围(((((2*sample_size-i)*(i-1)/2+(j-i))-1)*row_size)+1):((((2*sample_size-i)*(i-1)/2)+(j-i))*row_size)
                    %w_vector中关于gamma第ij系数差的范围(((((2*sample_size-i)*(i-1)/2+(j-i))-1)*col_size)+1):((((2*sample_size-i)*(i-1)/2)+(j-i))*col_size)
                    %E中第((2*sample_size-i)*(i-1)/2)+(j-i))行是1*sample_size关于ij的指示向量
                    %H_p维度为（p）*(样本数*p)
                    %将sample_size个x[:,:,sample_size]拼成一个大的对角块矩阵X(np*nq)
                    %先按照第三个维度转换成sample_size个元胞
                    %f1 = 0.5*norm(y-D_beta*X*gamma).^2+0.5*k*norm(H_p*beta-v_vec+kesai_vec/k).^2+0.5*k*norm(H_q*gamma-w_vec+yita_vec/k).^2;
                    %-----------------------------------------------------------------------------------------------------%
                    %-----------------------------------------------------------------------------------------------------%
                    [DATA_struct,data,DATA_opt,data_opt] = obj.admm_lambda(ifauto,lambda_size,lambda_start_point,change_constant,...
                        beta0,sample_size,p,iter_max1,iter_max1_interval,eps1,eps1_interval,H_p,X,y,k,a,...
                        lambda_left,lambda_right,c,v,kesai,v_vec,v_vec_bar,kesai_vec,beta_real,min_class_num,HH_p,E,iter_max2,eps2,iter_max_in,eps_in,...
                    eps_in_in,iter_max_in_in,full_if);

                    Results(local_i).DATA = DATA_struct;
                    results(local_i).data = data;
                    Results_opt(local_i).DATA = DATA_opt;
                    results_opt(local_i).data = data_opt;
                    fprintf('已完成第%.1f次的融合惩罚\n',local_i)
                end

            end
            % % % test_number = 2000;
            % % % PARAMETER_sum =zeros((row_size+col_size),group_number,simulation_size);
            % % % P_VALUE_sum = zeros((row_size+col_size),group_number,simulation_size);
            % % % for i = 1:simulation_size
            % % %     group_number = Results_opt(i).DATA.subgroup_number;%每一次模拟的亚组数
            % % %     subgroup_practical_number = zeros(1,group_number);
            % % %     PARAMETER_SD = zeros((row_size+col_size),test_number,group_number);%自助法抽样获得的第i次模拟对应每个亚组的估计系数
            % % %     PARAMETER = zeros((row_size+col_size),group_number);%第i次模拟估计出的每个亚组的结果
            % % %     P_VALUE = zeros((row_size+col_size),group_number);
            % % %     for j = 1:group_number
            % % %         index = Results_opt(i).DATA.class(j,1);
            % % %         PARAMETER(:,j) = [Results_opt(i).DATA.beta(((index-1)*row_size+1):index*row_size);Results_opt(i).DATA.gamma(((index-1)*col_size+1):index*col_size)];
            % % %         subgroup_practical_number(1,j) = sum(~isnan(Results_opt(i).DATA.class(j,:)));
            % % %     end
            % % %     for j = 1:group_number
            % % %         %x_test = zeros(row_size,col_size,subgroup_practical_number(1,j),test_number,group_number);
            % % %         %y_test = zeros(subgroup_practical_number(1,j),test_number,group_number);
            % % %         A = 1:sample_size;
            % % %         for t = 1:test_number
            % % %             sample_eachgroup_once = datasample(A, subgroup_practical_number(1,j), 'Replace', true);
            % % %             %x_test(:,:,:,t,j) = x(:,:,sample_eachgroup_once);
            % % %             %y_test(:,t,j) = y(sample_eachgroup_once,1);
            % % %             [result] = oracle_single_regression(x(:,:,sample_eachgroup_once), y(sample_eachgroup_once,1),...
            % % %                 subgroup_practical_number(1,j), row_size, col_size,...
            % % %                 zeros(subgroup_practical_number(1,j)*row_size,1),...
            % % %                 zeros(subgroup_practical_number(1,j)*col_size,1),...
            % % %                 eps_initial, iter_max_initial_in);
            % % %             PARAMETER_SD(:,t,j) = [result.beta_back(1:row_size);result.gamma_back(1:col_size)];
            % % %             fprintf('已完成第%.1d个实验数据的第%.1d组亚组系数的%.1f%%的假设检验,共%.1d个亚组\n',i,j,100*(t/test_number),group_number);
            % % %         end
            % % %         %t_val = PARAMETER(:,j)./sqrt(sum((PARAMETER_SD(:,:,j)-mean(PARAMETER_SD(:,:,j),2)*ones(1,test_number)).^2,2)/(subgroup_practical_number(1,j)));
            % % %         %P_VALUE(:,j) = 2*(1 - tcdf(abs(t_val), (test_number - 1)));
            % % %         %abs(PARAMETER(:,j))*ones(1,test_number);
            % % %         P_VALUE(:,j) = (sum(abs(PARAMETER_SD(:,:,j))>=abs(PARAMETER(:,j))*ones(1,test_number),2)/test_number);
            % % %     end
            % % %     PARAMETER_sum((row_size+col_size),group_number,simulation_size);
            % % %     P_VALUE_sum
            % % % end
            % % % PARAMETER;
            % % % P_VALUE;
            % % % fprintf('已完成假设检验\n');
            % % % fprintf('模型分析完成\n');
            
            for i = 1:length(simulation_index)
                Results_list_opt = [Results_list_opt;results_opt(i).data];
                [row, ~] = size(Results_opt(i).DATA.class);
                Class_summary(i).subgroup_number = row;
                Class_summary(i).coef = cell(1,row);
                Class_summary(i).subgroup_sample_size = cell(1,row);
                Class_summary(i).index = cell(1,row);
                Class_summary(i).X = cell(1,row);
                Class_summary(i).y = cell(1,row);
                for j = 1:row
                    Class_summary(i).index{1,j} = rmmissing(Results_opt(i).DATA.class(j,:));
                    Class_summary(i).subgroup_sample_size{1,j} = length(rmmissing(Results_opt(i).DATA.class(j,:)));
                    Class_summary(i).coef{1,j} = Results_opt(i).DATA.beta((Class_summary(i).index{1,j}(1)-1)*p+1:(Class_summary(i).index{1,j}(1))*p);
                    Class_summary(i).X{1,j} = Data_x(Class_summary(i).index{1,j},:,i);
                    Class_summary(i).y{1,j} = y_real(Class_summary(i).index{1,j},i);
                end
            end
            per = sum(Results_list_opt.per_if)/(length(simulation_index));
            K_mean = mean(Results_list_opt.subgroup_number);
            K_sd = std(Results_list_opt.subgroup_number);
            SC_mean = mean(Results_list_opt.sc);
            SC_sd = std(Results_list_opt.sc);
            MSE_beta_mean = mean(Results_list_opt.mse_beta);
            MSE_beta_sd = std(Results_list_opt.mse_beta);
            convergence_mean = mean(Results_list_opt.convergence);
            convergence_sd = std(Results_list_opt.convergence);
            iter_mean = mean(Results_list_opt.iter);
            iter_sd = std(Results_list_opt.iter);
            Results_single_opt = [K_mean,K_sd,per,SC_mean,SC_sd,MSE_beta_mean,MSE_beta_sd,convergence_mean,convergence_sd,iter_mean,iter_sd];
            Results_single_opt = array2table(Results_single_opt, 'VariableNames', {'K_mean','K_sd','per','SC_mean','SC_sd','MSE_beta_mean','MSE_beta_sd',...
                'convergence_mean','convergence_sd','iter_mean','iter_sd'});
            %  save('C:\Users\张丰川\Documents\MATLAB\Transfer-learning-heterogeneity-analysis\Result\h = 0.1\2th-auxiliary\RResults.mat','Results','results','Results_opt',...
            %  'results_opt','Results_list_opt','Results_single_opt','Class_summary')
            %save('C:\Users\张丰川\Documents\MATLAB\Transfer-learning-heterogeneity-analysis\Setting_1\Result\h = 0.1\single_auxiliary_0\RResults.mat','Results','results','Results_opt',...
            %   'results_opt','Results_list_opt','Results_single_opt','Class_summary')
        end
    end
end
