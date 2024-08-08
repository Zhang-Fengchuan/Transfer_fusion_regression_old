% This script is used for the production of Figures 2-3 and Tables 1,S1,S2,S3 in the manuscript
% Please modify the value of OBJ on Line 18 to produce a certain figure or table
% Re-running the simulation program may take a considerable amount of time
% If you want to directly load the previously generated data, please set load_if = 1

% INPUT ARGUMENT
% Please designate the figure or table to be reproduced
% e.g., 'F2' (Figure 2), 'F3' (Figure 3), 'T1' (Table 1),
% 'TS1' (Table S1), 'TS2' (Table S2) 'TS3' (Table S3) for the manuscript

% Please designate Re-running simulation program (load_if = 0)
% or directly load the previously generated data (load_if = 1)

% Please modify the address of the Transfer_fusion_regression folder
% on the local machine in the NewAddress at lines 20-21.


OBJ='F2';
load_if = 1;
NewAddress = 'C:\Users\张丰川\Documents\MATLAB\Transfer_fusion_regression\Setting_1\Result';
NewAddress2 = 'C:\Users\张丰川\Documents\MATLAB\Transfer_fusion_regression\Setting_2\Result';







% ------------------------------------Produce  Figure 2----------------------------------------
if strcmp(OBJ,'F2') == 1
    if load_if ~= 1
        for which_index = 0:4
            h_index = [0.1, 1, 2.5, 5, 7.5, 10];
            for i = 1:length(h_index)
                h = h_index(1, i);
                % Data generation for 1 target domain and 4 auxilary domains in a loop
                [simulation_size, X_target, Y_target, sample_size, p, beta_target_real]...
                    = Data_generation_setting_1(which_index, h);
                % Get the estimation results for 1 target domain and 4 auxilary domains
                % by using FUSION method
                Aux = Auxiliary_heterogeneity_regression();
                [Results,results,Results_opt,results_opt,...
                    initial,Results_list_opt,Results_single_opt,Class_summary]...
                    =  Aux.Auxiliary_heterogeneous_regression(simulation_size,[],...
                    X_target,Y_target,sample_size,p,beta_target_real,[],[],[],...
                    [],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],...
                    [],[],[],[],[],[],[],[],[],[],[]);
                % Save the estimation results for transfer method
                results_dir = fullfile(NewAddress, ['h=' num2str(h)]);
                save(fullfile(results_dir, ['results_' num2str(which_index) '.mat']), 'Results',...
                    'results', 'Results_opt', 'results_opt', 'initial', 'Results_list_opt',...
                    'Results_single_opt', 'Class_summary');
            end
        end
        multi_if = 0;
        for which_index = 1:4
            h_index = [0.1, 1, 2.5, 5, 7.5, 10];
            for t = 1:length(h_index)
                h = h_index(1, t);
                % Data generation for the target domain
                [simulation_size, X_target, Y_target, sample_size, p, beta_target_real]...
                    = Data_generation_setting_1(0, h);
                % Get the estimation results for target domain using S-TRANS
                % when the signal equal to h and the number of auxilary domains is which_index
                Trans = Transfer_heterogeneity_regression_pro();
                Auxiliary_dateset_number = which_index;
                Auxiliary = struct;
                Auxiliary(Auxiliary_dateset_number).position_provide = [];
                results_dir = fullfile(NewAddress,...
                    ['h=' num2str(h)]);
                for i = 1 : which_index
                    Auxiliary(i).position_provide = results_dir;
                end
                [Results,results,Results_opt,results_opt,...
                    initial,Results_list_opt,Results_single_opt,Class_summary]...
                    = Trans.Transfer_heterogeneous_regression_pro(Auxiliary_dateset_number,...
                    Auxiliary,multi_if,simulation_size,[],X_target,Y_target,sample_size,p,...
                    beta_target_real,[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],...
                    [],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[]);
                save(fullfile(results_dir, ['results_' strjoin(arrayfun(@num2str, 1:which_index, 'UniformOutput', false), '') '_multi_if=0.mat']),...
                    'Results', 'results', 'Results_opt', 'results_opt', 'initial',...
                    'Results_list_opt', 'Results_single_opt', 'Class_summary');
            end
        end
        multi_if = 1;
        for which_index = 1:4
            h_index = [0.1, 1, 2.5, 5, 7.5, 10];
            for t = 1:length(h_index)
                h = h_index(1, t);
                % Data generation for the target domain
                [simulation_size, X_target, Y_target, sample_size, p, beta_target_real]...
                    = Data_generation_setting_1(0, h);
                % Get the estimation results for target domain using M-TRANS
                % when the signal equal to h and the number of auxilary domains is which_index
                Trans = Transfer_heterogeneity_regression_pro();
                Auxiliary_dateset_number = which_index;
                Auxiliary = struct;
                Auxiliary(Auxiliary_dateset_number).position_provide = [];
                results_dir = fullfile(NewAddress,...
                    ['h=' num2str(h)]);
                for i = 1 : which_index
                    Auxiliary(i).position_provide = results_dir;
                end
                [Results,results,Results_opt,results_opt,...
                    initial,Results_list_opt,Results_single_opt,Class_summary]...
                    = Trans.Transfer_heterogeneous_regression_pro(Auxiliary_dateset_number,...
                    Auxiliary,multi_if,simulation_size,[],X_target,Y_target,sample_size,p,...
                    beta_target_real,[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],...
                    [],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[]);
                save(fullfile(results_dir, ['results_' strjoin(arrayfun(@num2str, 1:which_index, 'UniformOutput', false), '') '_multi_if=1.mat']),...
                    'Results', 'results', 'Results_opt', 'results_opt', 'initial',...
                    'Results_list_opt', 'Results_single_opt', 'Class_summary');
            end
        end
    else
        % Load the estimation result of FUSION
        simulation_size = 100;
        max_num_aux = 4;
        h_index = [0.1, 1, 2.5, 5, 7.5, 10];
        for t = 1:length(h_index)
            h = h_index(1, t);
            eval([['fusion' strrep(num2str(h), '.', '') '_msemean'] ' =  zeros(max_num_aux,1);']);
            eval([['fusion' strrep(num2str(h), '.', '') '_scmean'] ' =  zeros(max_num_aux,1);']);
            for which_index = 1:max_num_aux
                results_dir = fullfile(NewAddress,...
                    ['h=' num2str(h)]);
                load(fullfile(results_dir, 'results_0.mat'));
                mse_beta = [];
                for l = 1:simulation_size
                    mse_beta = [mse_beta;sqrt(results_opt(l).data.mse_beta)];
                end
                eval([ ['fusion' strrep(num2str(h), '.', '') '_msemean(which_index,1)'] ' =  mean(mse_beta);']);
                eval([ ['fusion' strrep(num2str(h), '.', '') '_scmean(which_index,1)'] ' =  Results_single_opt.SC_mean;']);
            end
        end
        fusion_msemean = [fusion10_msemean';fusion75_msemean';fusion5_msemean';fusion25_msemean';fusion1_msemean';fusion01_msemean'];
        fusion_scmean = [fusion10_scmean';fusion75_scmean';fusion5_scmean';fusion25_scmean';fusion1_scmean';fusion01_scmean'];
        % Load the estimation result of S-TRANS
        simulation_size = 100;
        max_num_aux = 4;
        h_index = [0.1, 1, 2.5, 5, 7.5, 10];
        for t = 1:length(h_index)
            h = h_index(1, t);
            eval([['s_subtrans' strrep(num2str(h), '.', '') '_msemean'] ' =  zeros(max_num_aux,1);']);
            eval([['s_subtrans' strrep(num2str(h), '.', '') '_scmean'] ' =  zeros(max_num_aux,1);']);
            for which_index = 1:max_num_aux
                results_dir = fullfile(NewAddress,...
                    ['h=' num2str(h)]);
                load(fullfile(results_dir, ['results_' strjoin(arrayfun(@num2str,...
                    1:which_index, 'UniformOutput', false), '') '_multi_if=0.mat']));
                mse_beta = [];
                for l = 1:simulation_size
                    mse_beta = [mse_beta;sqrt(results_opt(l).data.mse_beta)];
                end
                eval([ ['s_subtrans' strrep(num2str(h), '.', '') '_msemean(which_index,1)'] ' =  mean(mse_beta);']);
                eval([ ['s_subtrans' strrep(num2str(h), '.', '') '_scmean(which_index,1)'] ' =  Results_single_opt.SC_mean;']);
            end
        end
        s_subtrans_msemean = [s_subtrans10_msemean';s_subtrans75_msemean';s_subtrans5_msemean';s_subtrans25_msemean';s_subtrans1_msemean';s_subtrans01_msemean'];
        s_subtrans_scmean = [s_subtrans10_scmean';s_subtrans75_scmean';s_subtrans5_scmean';s_subtrans25_scmean';s_subtrans1_scmean';s_subtrans01_scmean'];
        % Load the estimation result of M-TRANS
        simulation_size = 100;
        max_num_aux = 4;
        h_index = [0.1, 1, 2.5, 5, 7.5, 10];
        for t = 1:length(h_index)
            h = h_index(1, t);
            eval([['m_subtrans' strrep(num2str(h), '.', '') '_msemean'] ' =  zeros(max_num_aux,1);']);
            eval([['m_subtrans' strrep(num2str(h), '.', '') '_scmean'] ' =  zeros(max_num_aux,1);']);
            for which_index = 1:max_num_aux
                results_dir = fullfile(NewAddress,...
                    ['h=' num2str(h)]);
                load(fullfile(results_dir, ['results_' strjoin(arrayfun(@num2str,...
                    1:which_index, 'UniformOutput', false), '') '_multi_if=1.mat']));
                mse_beta = [];
                for l = 1:simulation_size
                    mse_beta = [mse_beta;sqrt(results_opt(l).data.mse_beta)];
                end
                eval([ ['m_subtrans' strrep(num2str(h), '.', '') '_msemean(which_index,1)'] ' =  mean(mse_beta);']);
                eval([ ['m_subtrans' strrep(num2str(h), '.', '') '_scmean(which_index,1)'] ' =  Results_single_opt.SC_mean;']);
            end
        end
        m_subtrans_msemean = [m_subtrans10_msemean';m_subtrans75_msemean';m_subtrans5_msemean';m_subtrans25_msemean';m_subtrans1_msemean';m_subtrans01_msemean'];
        m_subtrans_scmean = [m_subtrans10_scmean';m_subtrans75_scmean';m_subtrans5_scmean';m_subtrans25_scmean';m_subtrans1_scmean';m_subtrans01_scmean'];
        figure
        % Plot the left subfigure
        subplot(1,2,1)
        which_index = 1 : 4;
        h_index = [0.1, 1, 2.5, 5, 7.5, 10];
        x = which_index;
        y = flip(10 - h_index);
        [X, Y] = meshgrid(x, y);
        % Set different color mappings and transparency
        C1 = [250,127,121]/256;
        C2 = [255,190,122]/256;
        C3 = [130,176,210]/256;
        data = {m_subtrans_msemean, s_subtrans_msemean, fusion_msemean};
        transparencies = [1, 0.5, 0.3];
        handles = gobjects(1, 3);
        handles(1) = surf(X,Y,data{3}, 'EdgeColor', [C3],'LineWidth',2);
        handles(1).CData = repmat(reshape(repmat(C3, [size(data{1}, 1),1]), [], 1, 3), [1 size(data{1}, 2) 1]);
        alpha(handles(1), transparencies(3));
        hold on;
        handles(2) = surf(X,Y,data{1},'EdgeColor', [C1],'LineWidth',2);
        handles(2).CData = repmat(reshape(repmat(C1, [size(data{2}, 1),1]), [], 1, 3), [1 size(data{2}, 2) 1]);
        alpha(handles(2), transparencies(1));
        hold on;
        handles(3) = surf(X,Y,data{2},'EdgeColor', [C2],'LineWidth',2);
        handles(3).CData = repmat(reshape(repmat(C2, [size(data{3}, 1),1]), [], 1, 3), [1 size(data{3}, 2) 1]);
        alpha(handles(3), transparencies(2));
        hold on;
        lighting phong;
        view(25, 23);
        axis tight;
        grid on;
        pos1 = get(gca, 'Position');
        % Create dummy legend handles for custom legend color blocks
        fake_handles = gobjects(1, 3);
        fake_handles(1) = plot(nan, nan, 'Color', C1, 'LineWidth', 5);
        fake_handles(2) = plot(nan, nan, 'Color', C2, 'LineWidth', 5);
        fake_handles(3) = plot(nan, nan, 'Color', C3, 'LineWidth', 5);
        legend_entries = {'M-SUBTRANS', 'S-SUBTRANS', 'FUSION'};
        % Add axis labels and title
        xlabel('Number of auxiliary studies(K)', 'FontSize', 14, 'FontName', 'Times New Roman');
        ylabel('Signal(h)', 'FontSize', 14, 'FontName', 'Times New Roman');
        zlabel('SC', 'FontSize', 14, 'FontName', 'Times New Roman');
        title('Trend of MSE with the signal h and the number of auxiliary studies K', 'FontSize', 16, 'FontName', 'Times New Roman');
        % Set the axis ticks, font, and size
        ax = gca;
        ax.XAxis.FontSize = 14;
        ax.YAxis.FontSize = 14;
        ax.ZAxis.FontSize = 14;
        ax.XAxis.FontName = 'Times New Roman';
        ax.YAxis.FontName = 'Times New Roman';
        ax.ZAxis.FontName = 'Times New Roman';
        ax.XAxis.TickDirection = 'out';
        ax.YAxis.TickDirection = 'out';
        ax.ZAxis.TickDirection = 'out';
        % Add a legend
        legend(fake_handles, legend_entries, 'Location', 'northeastoutside');
        hold off
        % Plot the right subfigure
        subplot(1,2,2)
        which_index = 1 : 4;
        h_index = [0.1, 1, 2.5, 5, 7.5, 10];
        x = which_index;
        y = flip(10 - h_index);
        [X,Y] = meshgrid(x,y);
        % Set different color mappings and transparency
        C1 = [250,127,121]/256;
        C2 = [255,190,122]/256;
        C3 = [130,176,210]/256;
        data = {m_subtrans_scmean, s_subtrans_scmean, fusion_scmean};
        transparencies = [0.8, 0.5, 0.3];
        handles = gobjects(1, 3);
        handles(1) = surf(X,Y,data{3},'EdgeColor', [C3],'LineWidth',2);
        handles(1).CData = repmat(reshape(repmat(C3, [size(data{1}, 1),1]), [], 1, 3), [1 size(data{1}, 2) 1]);
        alpha(handles(1), transparencies(3));
        hold on;
        handles(2) = surf(X,Y,data{2},'EdgeColor', [C2],'LineWidth',2);
        handles(2).CData = repmat(reshape(repmat(C2, [size(data{2}, 1),1]), [], 1, 3), [1 size(data{2}, 2) 1]);
        alpha(handles(2), transparencies(2));
        hold on;
        handles(3) = surf(X,Y,data{1},'EdgeColor', [C1],'LineWidth',2);
        handles(3).CData = repmat(reshape(repmat(C1, [size(data{3}, 1),1]), [], 1, 3), [1 size(data{3}, 2) 1]);
        alpha(handles(3), transparencies(1));
        hold on;
        lighting phong;
        view(46, 8); %56 20
        axis vis3d;
        grid on;
        % Create dummy legend handles for custom legend color blocks
        fake_handles = gobjects(1, 3);
        fake_handles(1) = plot(nan, nan, 'Color', C1, 'LineWidth', 5); % 使用颜色映射的最后一个颜色
        fake_handles(2) = plot(nan, nan, 'Color', C2, 'LineWidth', 5); % 使用颜色映射的最后一个颜色
        fake_handles(3) = plot(nan, nan, 'Color', C3, 'LineWidth', 5); % 使用颜色映射的最后一个颜色
        legend_entries = {'M-SUBTRANS', 'S-SUBTRANS', 'FUSION',};
        % Add axis labels and title
        xlabel('Number of auxiliary studies(K)', 'FontSize', 14, 'FontName', 'Times New Roman');
        ylabel('Signal(h)', 'FontSize', 14, 'FontName', 'Times New Roman');
        zlabel('SC', 'FontSize', 14, 'FontName', 'Times New Roman');
        title('Trend of SC with the signal h and the number of auxiliary studies K', 'FontSize', 16, 'FontName', 'Times New Roman');
        % Set the axis ticks, font, and size
        ax = gca;
        ax.XAxis.FontSize = 14;
        ax.YAxis.FontSize = 14;
        ax.ZAxis.FontSize = 14;
        ax.XAxis.FontName = 'Times New Roman';
        ax.YAxis.FontName = 'Times New Roman';
        ax.ZAxis.FontName = 'Times New Roman';
        ax.XAxis.TickDirection = 'out';
        ax.YAxis.TickDirection = 'out';
        ax.ZAxis.TickDirection = 'out';
        % Add a legend
        legend(fake_handles, legend_entries, 'Location', 'northeastoutside');
        hold off;
        set(gca, 'Position', pos1);
        pos2 = pos1;
        pos2(1) = pos2(1) + pos1(3) + 0.1;
        set(gca, 'Position', pos2);
    end
end
% ----------------------------------------------------------------------------------------------


















% ------------------------------------Produce  Table 1----------------------------------------
if strcmp(OBJ,'T1') == 1
    if load_if ~= 1
        h_index = [0.1, 1, 10];
        for t = 1:length(h_index)
            h = h_index(1, t);
            % Data generation for the target domain
            [simulation_size, X_target, Y_target, sample_size, p, beta_target_real]...
                = Data_generation_setting_1(0, h);
            % Get the estimation results for target domain using RESPCLUST
            % when the signal equal to h and the number of auxilary domains is which_index
            results_dir = fullfile(NewAddress,...
                ['h=' num2str(h)]);
            [Result_table, Result_table_summary] = RESPCLUST(simulation_size, X_target, Y_target,...
                sample_size, p, beta_target_real, [], [], []);
            save(fullfile(results_dir, 'results_resp.mat'), 'Result_table', 'Result_table_summary');
        end
        h_index = [0.1, 1, 10];
        for t = 1:length(h_index)
            h = h_index(1, t);
            % Data generation for the target domain
            [simulation_size, X_target, Y_target, sample_size, p, beta_target_real]...
                = Data_generation_setting_1(0, h);
            % Get the estimation results for target domain using RESIDUAL
            % when the signal equal to h and the number of auxilary domains is which_index
            results_dir = fullfile(NewAddress,...
                ['h=' num2str(h)]);
            [Result_table, Result_table_summary] = RESIDUAL(simulation_size, X_target, Y_target, sample_size, p, beta_target_real,...
                [], [], []);
            save(fullfile(results_dir, 'results_resi.mat'), 'Result_table', 'Result_table_summary');
        end
        for which_index = 0:4
            h_index = [0.1, 1, 10];
            for i = 1:length(h_index)
                h = h_index(1, i);
                % Data generation for 1 target domain and 4 auxilary domains in a loop
                [simulation_size, X_target, Y_target, sample_size, p, beta_target_real]...
                    = Data_generation_setting_1(which_index, h);
                % Get the estimation results for 1 target domain and 4 auxilary domains
                % by using FUSION method
                Aux = Auxiliary_heterogeneity_regression();
                [Results,results,Results_opt,results_opt,...
                    initial,Results_list_opt,Results_single_opt,Class_summary]...
                    =  Aux.Auxiliary_heterogeneous_regression(simulation_size,[],...
                    X_target,Y_target,sample_size,p,beta_target_real,[],[],[],...
                    [],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],...
                    [],[],[],[],[],[],[],[],[],[],[]);
                % Save the estimation results for transfer method
                results_dir = fullfile(NewAddress,...
                    ['h=' num2str(h)]);
                save(fullfile(results_dir, ['results_' num2str(which_index) '.mat']), 'Results',...
                    'results', 'Results_opt', 'results_opt', 'initial', 'Results_list_opt',...
                    'Results_single_opt', 'Class_summary');
            end
        end
        multi_if = 0;
        for which_index = 1:4
            h_index = [0.1, 1, 10];
            for t = 1:length(h_index)
                h = h_index(1, t);
                % Data generation for the target domain
                [simulation_size, X_target, Y_target, sample_size, p, beta_target_real]...
                    = Data_generation_setting_1(0, h);
                % Get the estimation results for target domain using S-TRANS
                % when the signal equal to h and the number of auxilary domains is which_index
                Trans = Transfer_heterogeneity_regression_pro();
                Auxiliary_dateset_number = which_index;
                Auxiliary = struct;
                Auxiliary(Auxiliary_dateset_number).position_provide = [];
                results_dir = fullfile(NewAddress,...
                    ['h=' num2str(h)]);
                for i = 1 : which_index
                    Auxiliary(i).position_provide = results_dir;
                end
                [Results,results,Results_opt,results_opt,...
                    initial,Results_list_opt,Results_single_opt,Class_summary]...
                    = Trans.Transfer_heterogeneous_regression_pro(Auxiliary_dateset_number,...
                    Auxiliary,multi_if,simulation_size,[],X_target,Y_target,sample_size,p,...
                    beta_target_real,[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],...
                    [],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[]);
                save(fullfile(results_dir, ['results_' strjoin(arrayfun(@num2str, 1:which_index, 'UniformOutput', false), '') '_multi_if=0.mat']),...
                    'Results', 'results', 'Results_opt', 'results_opt', 'initial',...
                    'Results_list_opt', 'Results_single_opt', 'Class_summary');
            end
        end
        multi_if = 1;
        for which_index = 1:4
            h_index = [0.1, 1, 10];
            for t = 1:length(h_index)
                h = h_index(1, t);
                % Data generation for the target domain
                [simulation_size, X_target, Y_target, sample_size, p, beta_target_real]...
                    = Data_generation_setting_1(0, h);
                % Get the estimation results for target domain using M-TRANS
                % when the signal equal to h and the number of auxilary domains is which_index
                Trans = Transfer_heterogeneity_regression_pro();
                Auxiliary_dateset_number = which_index;
                Auxiliary = struct;
                Auxiliary(Auxiliary_dateset_number).position_provide = [];
                results_dir = fullfile(NewAddress,...
                    ['h=' num2str(h)]);
                for i = 1 : which_index
                    Auxiliary(i).position_provide = results_dir;
                end
                [Results,results,Results_opt,results_opt,...
                    initial,Results_list_opt,Results_single_opt,Class_summary]...
                    = Trans.Transfer_heterogeneous_regression_pro(Auxiliary_dateset_number,...
                    Auxiliary,multi_if,simulation_size,[],X_target,Y_target,sample_size,p,...
                    beta_target_real,[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],...
                    [],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[]);
                save(fullfile(results_dir, ['results_' strjoin(arrayfun(@num2str, 1:which_index, 'UniformOutput', false), '') '_multi_if=1.mat']),...
                    'Results', 'results', 'Results_opt', 'results_opt', 'initial',...
                    'Results_list_opt', 'Results_single_opt', 'Class_summary');
            end
        end
    else
        % Load the date for Table 1
        simulation_size = 100;
        max_num_aux = 4;
        h_index = [0.1, 1, 10];
        for t = 1:length(h_index)
            h = h_index(1, t);
            results_dir = fullfile(NewAddress,...
                ['h=' num2str(h)]);
            load(fullfile(results_dir, 'results_0.mat'))
            mse_beta = [];
            for l = 1:simulation_size
                mse_beta = [mse_beta;sqrt(results_opt(l).data.mse_beta)];
            end
            eval([ ['fusion' strrep(num2str(h), '.', '') '_msemean'] ' =  mean(mse_beta);']);
            eval([ ['fusion' strrep(num2str(h), '.', '') '_msesd'] ' =  std(mse_beta);']);
            eval([ ['fusion' strrep(num2str(h), '.', '') '_scmean'] ' =  Results_single_opt.SC_mean;']);
            eval([ ['fusion' strrep(num2str(h), '.', '') '_scsd'] ' =  Results_single_opt.SC_sd;']);
            eval([['fusion' strrep(num2str(h), '.', '') '_Mmean'] ' =  Results_single_opt.K_mean;']);
            eval([['fusion' strrep(num2str(h), '.', '') '_Msd'] ' =  Results_single_opt.K_sd;']);
            eval([['fusion' strrep(num2str(h), '.', '') '_Per'] ' =  Results_single_opt.per;']);
            results_dir = fullfile(NewAddress,...
                ['h=' num2str(h)]);
            load(fullfile(results_dir, 'results_resp.mat'))
            eval([ ['resp' strrep(num2str(h), '.', '') '_msemean'] ' =  mean(sqrt(Result_table.mse_beta));']);
            eval([ ['resp' strrep(num2str(h), '.', '') '_msesd'] ' =  std(sqrt(Result_table.mse_beta));']);
            eval([ ['resp' strrep(num2str(h), '.', '') '_scmean'] ' =  mean(Result_table.sc);']);
            eval([ ['resp' strrep(num2str(h), '.', '') '_scsd'] ' =  std(Result_table.sc);']);
            eval([['resp' strrep(num2str(h), '.', '') '_Mmean'] ' =  mean(Result_table.k);']);
            eval([['resp' strrep(num2str(h), '.', '') '_Msd'] ' =  std(Result_table.k);']);
            eval([['resp' strrep(num2str(h), '.', '') '_Per'] ' =  mean(Result_table.per_if);']);
            results_dir = fullfile(NewAddress,...
                ['h=' num2str(h)]);
            load(fullfile(results_dir, 'results_resi.mat'))
            eval([ ['resi' strrep(num2str(h), '.', '') '_msemean'] ' =  mean(sqrt(Result_table.mse_beta));']);
            eval([ ['resi' strrep(num2str(h), '.', '') '_msesd'] ' =  std(sqrt(Result_table.mse_beta));']);
            eval([ ['resi' strrep(num2str(h), '.', '') '_scmean'] ' =  mean(Result_table.sc);']);
            eval([ ['resi' strrep(num2str(h), '.', '') '_scsd'] ' =  std(Result_table.sc);']);
            eval([['resi' strrep(num2str(h), '.', '') '_Mmean'] ' =  mean(Result_table.k);']);
            eval([['resi' strrep(num2str(h), '.', '') '_Msd'] ' =  std(Result_table.k);']);
            eval([['resi' strrep(num2str(h), '.', '') '_Per'] ' =  mean(Result_table.per_if);']);

            eval([['s_subtrans' strrep(num2str(h), '.', '') '_Mmean'] ' =  zeros(max_num_aux,1);']);
            eval([['s_subtrans' strrep(num2str(h), '.', '') '_Msd'] ' =  zeros(max_num_aux,1);']);
            eval([['s_subtrans' strrep(num2str(h), '.', '') '_Per'] ' =  zeros(max_num_aux,1);']);
            eval([['s_subtrans' strrep(num2str(h), '.', '') '_msemean'] ' =  zeros(max_num_aux,1);']);
            eval([['s_subtrans' strrep(num2str(h), '.', '') '_msesd'] ' =  zeros(max_num_aux,1);']);
            eval([['s_subtrans' strrep(num2str(h), '.', '') '_scmean'] ' =  zeros(max_num_aux,1);']);
            eval([['s_subtrans' strrep(num2str(h), '.', '') '_scsd'] ' =  zeros(max_num_aux,1);']);
            for which_index = 1:max_num_aux
                results_dir = fullfile(NewAddress,...
                    ['h=' num2str(h)]);
                load(fullfile(results_dir, ['results_' strjoin(arrayfun(@num2str,...
                    1:which_index, 'UniformOutput', false), '') '_multi_if=0.mat']));
                mse_beta = [];
                for l = 1:simulation_size
                    mse_beta = [mse_beta;sqrt(results_opt(l).data.mse_beta)];
                end
                eval([ ['s_subtrans' strrep(num2str(h), '.', '') '_msemean(which_index,1)'] ' =  mean(mse_beta);']);
                eval([ ['s_subtrans' strrep(num2str(h), '.', '') '_msesd(which_index,1)'] ' =  std(mse_beta);']);
                eval([ ['s_subtrans' strrep(num2str(h), '.', '') '_scmean(which_index,1)'] ' =  Results_single_opt.SC_mean;']);
                eval([ ['s_subtrans' strrep(num2str(h), '.', '') '_scsd(which_index,1)'] ' =  Results_single_opt.SC_sd;']);
                eval([ ['s_subtrans' strrep(num2str(h), '.', '') '_Mmean(which_index,1)'] ' =  Results_single_opt.K_mean;']);
                eval([ ['s_subtrans' strrep(num2str(h), '.', '') '_Msd(which_index,1)'] ' =  Results_single_opt.K_sd;']);
                eval([ ['s_subtrans' strrep(num2str(h), '.', '') '_Per(which_index,1)'] ' =  Results_single_opt.per;']);
            end
            eval([['m_subtrans' strrep(num2str(h), '.', '') '_Mmean'] ' =  zeros(max_num_aux,1);']);
            eval([['m_subtrans' strrep(num2str(h), '.', '') '_Msd'] ' =  zeros(max_num_aux,1);']);
            eval([['m_subtrans' strrep(num2str(h), '.', '') '_Per'] ' =  zeros(max_num_aux,1);']);
            eval([['m_subtrans' strrep(num2str(h), '.', '') '_msemean'] ' =  zeros(max_num_aux,1);']);
            eval([['m_subtrans' strrep(num2str(h), '.', '') '_msesd'] ' =  zeros(max_num_aux,1);']);
            eval([['m_subtrans' strrep(num2str(h), '.', '') '_scmean'] ' =  zeros(max_num_aux,1);']);
            eval([['m_subtrans' strrep(num2str(h), '.', '') '_scsd'] ' =  zeros(max_num_aux,1);']);
            for which_index = 1:max_num_aux
                results_dir = fullfile(NewAddress,...
                    ['h=' num2str(h)]);
                load(fullfile(results_dir, ['results_' strjoin(arrayfun(@num2str,...
                    1:which_index, 'UniformOutput', false), '') '_multi_if=1.mat']));
                mse_beta = [];
                for l = 1:simulation_size
                    mse_beta = [mse_beta;sqrt(results_opt(l).data.mse_beta)];
                end
                eval([ ['m_subtrans' strrep(num2str(h), '.', '') '_msemean(which_index,1)'] ' =  mean(mse_beta);']);
                eval([ ['m_subtrans' strrep(num2str(h), '.', '') '_msesd(which_index,1)'] ' =  std(mse_beta);']);
                eval([ ['m_subtrans' strrep(num2str(h), '.', '') '_scmean(which_index,1)'] ' =  Results_single_opt.SC_mean;']);
                eval([ ['m_subtrans' strrep(num2str(h), '.', '') '_scsd(which_index,1)'] ' =  Results_single_opt.SC_sd;']);
                eval([ ['m_subtrans' strrep(num2str(h), '.', '') '_Mmean(which_index,1)'] ' =  Results_single_opt.K_mean;']);
                eval([ ['m_subtrans' strrep(num2str(h), '.', '') '_Msd(which_index,1)'] ' =  Results_single_opt.K_sd;']);
                eval([ ['m_subtrans' strrep(num2str(h), '.', '') '_Per(which_index,1)'] ' =  Results_single_opt.per;']);
            end
        end
        % Load the date for Table 1
        max_num_aux = 4;
        h_index = [0.1, 1, 10];
        for t = 1:length(h_index)
            h = h_index(1, t);
            M_mean = [eval(['resp' strrep(num2str(h), '.', '') '_Mmean']);eval(['resi' strrep(num2str(h), '.', '') '_Mmean']);...
                eval(['fusion' strrep(num2str(h), '.', '') '_Mmean']);];
            for i = 1 : max_num_aux
                M_mean = [M_mean;eval(['s_subtrans' strrep(num2str(h), '.', '') '_Mmean(i,1)'])];
                M_mean = [M_mean;eval(['m_subtrans' strrep(num2str(h), '.', '') '_Mmean(i,1)'])];
            end
            eval([ ['M_mean' strrep(num2str(h), '.', '')] ' =  M_mean;']);

            M_sd = [eval(['resp' strrep(num2str(h), '.', '') '_Msd']);eval(['resi' strrep(num2str(h), '.', '') '_Msd']);...
                eval(['fusion' strrep(num2str(h), '.', '') '_Msd']);];
            for i = 1 : max_num_aux
                M_sd = [M_sd;eval(['s_subtrans' strrep(num2str(h), '.', '') '_Msd(i,1)'])];
                M_sd = [M_sd;eval(['m_subtrans' strrep(num2str(h), '.', '') '_Msd(i,1)'])];
            end
            eval([ ['M_sd' strrep(num2str(h), '.', '')] ' =  M_sd;']);

            Per = [eval(['resp' strrep(num2str(h), '.', '') '_Per']);eval(['resi' strrep(num2str(h), '.', '') '_Per']);...
                eval(['fusion' strrep(num2str(h), '.', '') '_Per']);];
            for i = 1 : max_num_aux
                Per = [Per;eval(['s_subtrans' strrep(num2str(h), '.', '') '_Per(i,1)'])];
                Per = [Per;eval(['m_subtrans' strrep(num2str(h), '.', '') '_Per(i,1)'])];
            end
            eval([ ['Per' strrep(num2str(h), '.', '')] ' =  Per;']);

            Mse_mean = [eval(['resp' strrep(num2str(h), '.', '') '_msemean']);eval(['resi' strrep(num2str(h), '.', '') '_msemean']);...
                eval(['fusion' strrep(num2str(h), '.', '') '_msemean']);];
            for i = 1 : max_num_aux
                Mse_mean = [Mse_mean;eval(['s_subtrans' strrep(num2str(h), '.', '') '_msemean(i,1)'])];
                Mse_mean = [Mse_mean;eval(['m_subtrans' strrep(num2str(h), '.', '') '_msemean(i,1)'])];
            end
            eval([ ['Mse_mean' strrep(num2str(h), '.', '')] ' =  Mse_mean;']);

            Mse_sd = [eval(['resp' strrep(num2str(h), '.', '') '_msesd']);eval(['resi' strrep(num2str(h), '.', '') '_msesd']);...
                eval(['fusion' strrep(num2str(h), '.', '') '_msesd']);];
            for i = 1 : max_num_aux
                Mse_sd = [Mse_sd;eval(['s_subtrans' strrep(num2str(h), '.', '') '_msesd(i,1)'])];
                Mse_sd = [Mse_sd;eval(['m_subtrans' strrep(num2str(h), '.', '') '_msesd(i,1)'])];
            end
            eval([ ['Mse_sd' strrep(num2str(h), '.', '')] ' =  Mse_sd;']);

            Sc_mean = [eval(['resp' strrep(num2str(h), '.', '') '_scmean']);eval(['resi' strrep(num2str(h), '.', '') '_scmean']);...
                eval(['fusion' strrep(num2str(h), '.', '') '_scmean']);];
            for i = 1 : max_num_aux
                Sc_mean = [Sc_mean;eval(['s_subtrans' strrep(num2str(h), '.', '') '_scmean(i,1)'])];
                Sc_mean = [Sc_mean;eval(['m_subtrans' strrep(num2str(h), '.', '') '_scmean(i,1)'])];
            end
            eval([ ['Sc_mean' strrep(num2str(h), '.', '')] ' =  Sc_mean;']);

            Sc_sd = [eval(['resp' strrep(num2str(h), '.', '') '_scsd']);eval(['resi' strrep(num2str(h), '.', '') '_scsd']);...
                eval(['fusion' strrep(num2str(h), '.', '') '_scsd']);];
            for i = 1 : max_num_aux
                Sc_sd = [Sc_sd;eval(['s_subtrans' strrep(num2str(h), '.', '') '_scsd(i,1)'])];
                Sc_sd = [Sc_sd;eval(['m_subtrans' strrep(num2str(h), '.', '') '_scsd(i,1)'])];
            end
            eval([ ['Sc_sd' strrep(num2str(h), '.', '')] ' =  Sc_sd;']);
        end
        Table_1 = [];
        h_index = [0.1, 1, 10];
        for t = 1:length(h_index)
            h = h_index(1, t);
            Table_1 = [Table_1, [eval(['M_mean' strrep(num2str(h), '.', '')]), eval(['M_sd' strrep(num2str(h), '.', '')]),...
                eval(['Per' strrep(num2str(h), '.', '')])]];
        end
        disp(Table_1)
    end
end
% ----------------------------------------------------------------------------------------------




















% ------------------------------------Produce  Table S1----------------------------------------
if strcmp(OBJ,'TS1') == 1
    if load_if ~= 1
        h_index = [0.1, 1, 2.5];
        for t = 1:length(h_index)
            h = h_index(1, t);
            % Data generation for the target domain
            [simulation_size, X_target, Y_target, sample_size, p, beta_target_real]...
                = Data_generation_setting_1(0, h);
            % Get the estimation results for target domain using M-TRANS
            % when the signal equal to h and the number of auxilary domains is which_index
            results_dir = fullfile(NewAddress,...
                ['h=' num2str(h)]);
            [Result_table, Result_table_summary] = RESPCLUST(simulation_size, X_target, Y_target,...
                sample_size, p, beta_target_real, [], [], []);
            save(fullfile(results_dir, 'results_resp.mat'), 'Result_table', 'Result_table_summary');
        end
        h_index = [0.1, 1, 2.5];
        for t = 1:length(h_index)
            h = h_index(1, t);
            % Data generation for the target domain
            [simulation_size, X_target, Y_target, sample_size, p, beta_target_real]...
                = Data_generation_setting_1(0, h);
            % Get the estimation results for target domain using M-TRANS
            % when the signal equal to h and the number of auxilary domains is which_index
            results_dir = fullfile(NewAddress,...
                ['h=' num2str(h)]);
            [Result_table, Result_table_summary] = RESIDUAL(simulation_size, X_target, Y_target, sample_size, p, beta_target_real,...
                [], [], []);
            save(fullfile(results_dir, 'results_resi.mat'), 'Result_table', 'Result_table_summary');
        end
        for which_index = 0:4
            h_index = [0.1, 1, 2.5];
            for i = 1:length(h_index)
                h = h_index(1, i);
                % Data generation for 1 target domain and 4 auxilary domains in a loop
                [simulation_size, X_target, Y_target, sample_size, p, beta_target_real]...
                    = Data_generation_setting_1(which_index, h);
                % Get the estimation results for 1 target domain and 4 auxilary domains
                % by using FUSION method
                Aux = Auxiliary_heterogeneity_regression();
                [Results,results,Results_opt,results_opt,...
                    initial,Results_list_opt,Results_single_opt,Class_summary]...
                    =  Aux.Auxiliary_heterogeneous_regression(simulation_size,[],...
                    X_target,Y_target,sample_size,p,beta_target_real,[],[],[],...
                    [],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],...
                    [],[],[],[],[],[],[],[],[],[],[]);
                % Save the estimation results for transfer method
                results_dir = fullfile(NewAddress,...
                    ['h=' num2str(h)]);
                save(fullfile(results_dir, ['results_' num2str(which_index) '.mat']), 'Results',...
                    'results', 'Results_opt', 'results_opt', 'initial', 'Results_list_opt',...
                    'Results_single_opt', 'Class_summary');
            end
        end
        multi_if = 0;
        for which_index = 1:4
            h_index = [0.1, 1, 2.5];
            for t = 1:length(h_index)
                h = h_index(1, t);
                % Data generation for the target domain
                [simulation_size, X_target, Y_target, sample_size, p, beta_target_real]...
                    = Data_generation_setting_1(0, h);
                % Get the estimation results for target domain using S-TRANS
                % when the signal equal to h and the number of auxilary domains is which_index
                Trans = Transfer_heterogeneity_regression_pro();
                Auxiliary_dateset_number = which_index;
                Auxiliary = struct;
                Auxiliary(Auxiliary_dateset_number).position_provide = [];
                results_dir = fullfile(NewAddress,...
                    ['h=' num2str(h)]);
                for i = 1 : which_index
                    Auxiliary(i).position_provide = results_dir;
                end
                [Results,results,Results_opt,results_opt,...
                    initial,Results_list_opt,Results_single_opt,Class_summary]...
                    = Trans.Transfer_heterogeneous_regression_pro(Auxiliary_dateset_number,...
                    Auxiliary,multi_if,simulation_size,[],X_target,Y_target,sample_size,p,...
                    beta_target_real,[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],...
                    [],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[]);
                save(fullfile(results_dir, ['results_' strjoin(arrayfun(@num2str, 1:which_index, 'UniformOutput', false), '') '_multi_if=0.mat']),...
                    'Results', 'results', 'Results_opt', 'results_opt', 'initial',...
                    'Results_list_opt', 'Results_single_opt', 'Class_summary');
            end
        end
        multi_if = 1;
        for which_index = 1:4
            h_index = [0.1, 1, 2.5];
            for t = 1:length(h_index)
                h = h_index(1, t);
                % Data generation for the target domain
                [simulation_size, X_target, Y_target, sample_size, p, beta_target_real]...
                    = Data_generation_setting_1(0, h);
                % Get the estimation results for target domain using M-TRANS
                % when the signal equal to h and the number of auxilary domains is which_index
                Trans = Transfer_heterogeneity_regression_pro();
                Auxiliary_dateset_number = which_index;
                Auxiliary = struct;
                Auxiliary(Auxiliary_dateset_number).position_provide = [];
                results_dir = fullfile(NewAddress,...
                    ['h=' num2str(h)]);
                for i = 1 : which_index
                    Auxiliary(i).position_provide = results_dir;
                end
                [Results,results,Results_opt,results_opt,...
                    initial,Results_list_opt,Results_single_opt,Class_summary]...
                    = Trans.Transfer_heterogeneous_regression_pro(Auxiliary_dateset_number,...
                    Auxiliary,multi_if,simulation_size,[],X_target,Y_target,sample_size,p,...
                    beta_target_real,[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],...
                    [],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[]);
                save(fullfile(results_dir, ['results_' strjoin(arrayfun(@num2str, 1:which_index, 'UniformOutput', false), '') '_multi_if=1.mat']),...
                    'Results', 'results', 'Results_opt', 'results_opt', 'initial',...
                    'Results_list_opt', 'Results_single_opt', 'Class_summary');
            end
        end
    else
        % Load the date for Table S1
        simulation_size = 100;
        max_num_aux = 4;
        h_index = [0.1, 1, 2.5];
        for t = 1:length(h_index)
            h = h_index(1, t);
            results_dir = fullfile(NewAddress,...
                ['h=' num2str(h)]);
            load(fullfile(results_dir, 'results_0.mat'))
            mse_beta = [];
            for l = 1:simulation_size
                mse_beta = [mse_beta;sqrt(results_opt(l).data.mse_beta)];
            end
            eval([ ['fusion' strrep(num2str(h), '.', '') '_msemean'] ' =  mean(mse_beta);']);
            eval([ ['fusion' strrep(num2str(h), '.', '') '_msesd'] ' =  std(mse_beta);']);
            eval([ ['fusion' strrep(num2str(h), '.', '') '_scmean'] ' =  Results_single_opt.SC_mean;']);
            eval([ ['fusion' strrep(num2str(h), '.', '') '_scsd'] ' =  Results_single_opt.SC_sd;']);
            eval([['fusion' strrep(num2str(h), '.', '') '_Mmean'] ' =  Results_single_opt.K_mean;']);
            eval([['fusion' strrep(num2str(h), '.', '') '_Msd'] ' =  Results_single_opt.K_sd;']);
            eval([['fusion' strrep(num2str(h), '.', '') '_Per'] ' =  Results_single_opt.per;']);
            results_dir = fullfile(NewAddress,...
                ['h=' num2str(h)]);
            load(fullfile(results_dir, 'results_resp.mat'))
            eval([ ['resp' strrep(num2str(h), '.', '') '_msemean'] ' =  mean(sqrt(Result_table.mse_beta));']);
            eval([ ['resp' strrep(num2str(h), '.', '') '_msesd'] ' =  std(sqrt(Result_table.mse_beta));']);
            eval([ ['resp' strrep(num2str(h), '.', '') '_scmean'] ' =  mean(Result_table.sc);']);
            eval([ ['resp' strrep(num2str(h), '.', '') '_scsd'] ' =  std(Result_table.sc);']);
            eval([['resp' strrep(num2str(h), '.', '') '_Mmean'] ' =  mean(Result_table.k);']);
            eval([['resp' strrep(num2str(h), '.', '') '_Msd'] ' =  std(Result_table.k);']);
            eval([['resp' strrep(num2str(h), '.', '') '_Per'] ' =  mean(Result_table.per_if);']);
            results_dir = fullfile(NewAddress,...
                ['h=' num2str(h)]);
            load(fullfile(results_dir, 'results_resi.mat'))
            eval([ ['resi' strrep(num2str(h), '.', '') '_msemean'] ' =  mean(sqrt(Result_table.mse_beta));']);
            eval([ ['resi' strrep(num2str(h), '.', '') '_msesd'] ' =  std(sqrt(Result_table.mse_beta));']);
            eval([ ['resi' strrep(num2str(h), '.', '') '_scmean'] ' =  mean(Result_table.sc);']);
            eval([ ['resi' strrep(num2str(h), '.', '') '_scsd'] ' =  std(Result_table.sc);']);
            eval([['resi' strrep(num2str(h), '.', '') '_Mmean'] ' =  mean(Result_table.k);']);
            eval([['resi' strrep(num2str(h), '.', '') '_Msd'] ' =  std(Result_table.k);']);
            eval([['resi' strrep(num2str(h), '.', '') '_Per'] ' =  mean(Result_table.per_if);']);

            eval([['s_subtrans' strrep(num2str(h), '.', '') '_Mmean'] ' =  zeros(max_num_aux,1);']);
            eval([['s_subtrans' strrep(num2str(h), '.', '') '_Msd'] ' =  zeros(max_num_aux,1);']);
            eval([['s_subtrans' strrep(num2str(h), '.', '') '_Per'] ' =  zeros(max_num_aux,1);']);
            eval([['s_subtrans' strrep(num2str(h), '.', '') '_msemean'] ' =  zeros(max_num_aux,1);']);
            eval([['s_subtrans' strrep(num2str(h), '.', '') '_msesd'] ' =  zeros(max_num_aux,1);']);
            eval([['s_subtrans' strrep(num2str(h), '.', '') '_scmean'] ' =  zeros(max_num_aux,1);']);
            eval([['s_subtrans' strrep(num2str(h), '.', '') '_scsd'] ' =  zeros(max_num_aux,1);']);
            for which_index = 1:max_num_aux
                results_dir = fullfile(NewAddress,...
                    ['h=' num2str(h)]);
                load(fullfile(results_dir, ['results_' strjoin(arrayfun(@num2str,...
                    1:which_index, 'UniformOutput', false), '') '_multi_if=0.mat']));
                mse_beta = [];
                for l = 1:simulation_size
                    mse_beta = [mse_beta;sqrt(results_opt(l).data.mse_beta)];
                end
                eval([ ['s_subtrans' strrep(num2str(h), '.', '') '_msemean(which_index,1)'] ' =  mean(mse_beta);']);
                eval([ ['s_subtrans' strrep(num2str(h), '.', '') '_msesd(which_index,1)'] ' =  std(mse_beta);']);
                eval([ ['s_subtrans' strrep(num2str(h), '.', '') '_scmean(which_index,1)'] ' =  Results_single_opt.SC_mean;']);
                eval([ ['s_subtrans' strrep(num2str(h), '.', '') '_scsd(which_index,1)'] ' =  Results_single_opt.SC_sd;']);
                eval([ ['s_subtrans' strrep(num2str(h), '.', '') '_Mmean(which_index,1)'] ' =  Results_single_opt.K_mean;']);
                eval([ ['s_subtrans' strrep(num2str(h), '.', '') '_Msd(which_index,1)'] ' =  Results_single_opt.K_sd;']);
                eval([ ['s_subtrans' strrep(num2str(h), '.', '') '_Per(which_index,1)'] ' =  Results_single_opt.per;']);
            end
            eval([['m_subtrans' strrep(num2str(h), '.', '') '_Mmean'] ' =  zeros(max_num_aux,1);']);
            eval([['m_subtrans' strrep(num2str(h), '.', '') '_Msd'] ' =  zeros(max_num_aux,1);']);
            eval([['m_subtrans' strrep(num2str(h), '.', '') '_Per'] ' =  zeros(max_num_aux,1);']);
            eval([['m_subtrans' strrep(num2str(h), '.', '') '_msemean'] ' =  zeros(max_num_aux,1);']);
            eval([['m_subtrans' strrep(num2str(h), '.', '') '_msesd'] ' =  zeros(max_num_aux,1);']);
            eval([['m_subtrans' strrep(num2str(h), '.', '') '_scmean'] ' =  zeros(max_num_aux,1);']);
            eval([['m_subtrans' strrep(num2str(h), '.', '') '_scsd'] ' =  zeros(max_num_aux,1);']);
            for which_index = 1:max_num_aux
                results_dir = fullfile(NewAddress,...
                    ['h=' num2str(h)]);
                load(fullfile(results_dir, ['results_' strjoin(arrayfun(@num2str,...
                    1:which_index, 'UniformOutput', false), '') '_multi_if=1.mat']));
                mse_beta = [];
                for l = 1:simulation_size
                    mse_beta = [mse_beta;sqrt(results_opt(l).data.mse_beta)];
                end
                eval([ ['m_subtrans' strrep(num2str(h), '.', '') '_msemean(which_index,1)'] ' =  mean(mse_beta);']);
                eval([ ['m_subtrans' strrep(num2str(h), '.', '') '_msesd(which_index,1)'] ' =  std(mse_beta);']);
                eval([ ['m_subtrans' strrep(num2str(h), '.', '') '_scmean(which_index,1)'] ' =  Results_single_opt.SC_mean;']);
                eval([ ['m_subtrans' strrep(num2str(h), '.', '') '_scsd(which_index,1)'] ' =  Results_single_opt.SC_sd;']);
                eval([ ['m_subtrans' strrep(num2str(h), '.', '') '_Mmean(which_index,1)'] ' =  Results_single_opt.K_mean;']);
                eval([ ['m_subtrans' strrep(num2str(h), '.', '') '_Msd(which_index,1)'] ' =  Results_single_opt.K_sd;']);
                eval([ ['m_subtrans' strrep(num2str(h), '.', '') '_Per(which_index,1)'] ' =  Results_single_opt.per;']);
            end
        end
        % Load the date for Table S1
        max_num_aux = 4;
        h_index = [0.1, 1, 2.5];
        for t = 1:length(h_index)
            h = h_index(1, t);
            M_mean = [eval(['resp' strrep(num2str(h), '.', '') '_Mmean']);eval(['resi' strrep(num2str(h), '.', '') '_Mmean']);...
                eval(['fusion' strrep(num2str(h), '.', '') '_Mmean']);];
            for i = 1 : max_num_aux
                M_mean = [M_mean;eval(['s_subtrans' strrep(num2str(h), '.', '') '_Mmean(i,1)'])];
                M_mean = [M_mean;eval(['m_subtrans' strrep(num2str(h), '.', '') '_Mmean(i,1)'])];
            end
            eval([ ['M_mean' strrep(num2str(h), '.', '')] ' =  M_mean;']);

            M_sd = [eval(['resp' strrep(num2str(h), '.', '') '_Msd']);eval(['resi' strrep(num2str(h), '.', '') '_Msd']);...
                eval(['fusion' strrep(num2str(h), '.', '') '_Msd']);];
            for i = 1 : max_num_aux
                M_sd = [M_sd;eval(['s_subtrans' strrep(num2str(h), '.', '') '_Msd(i,1)'])];
                M_sd = [M_sd;eval(['m_subtrans' strrep(num2str(h), '.', '') '_Msd(i,1)'])];
            end
            eval([ ['M_sd' strrep(num2str(h), '.', '')] ' =  M_sd;']);

            Per = [eval(['resp' strrep(num2str(h), '.', '') '_Per']);eval(['resi' strrep(num2str(h), '.', '') '_Per']);...
                eval(['fusion' strrep(num2str(h), '.', '') '_Per']);];
            for i = 1 : max_num_aux
                Per = [Per;eval(['s_subtrans' strrep(num2str(h), '.', '') '_Per(i,1)'])];
                Per = [Per;eval(['m_subtrans' strrep(num2str(h), '.', '') '_Per(i,1)'])];
            end
            eval([ ['Per' strrep(num2str(h), '.', '')] ' =  Per;']);

            Mse_mean = [eval(['resp' strrep(num2str(h), '.', '') '_msemean']);eval(['resi' strrep(num2str(h), '.', '') '_msemean']);...
                eval(['fusion' strrep(num2str(h), '.', '') '_msemean']);];
            for i = 1 : max_num_aux
                Mse_mean = [Mse_mean;eval(['s_subtrans' strrep(num2str(h), '.', '') '_msemean(i,1)'])];
                Mse_mean = [Mse_mean;eval(['m_subtrans' strrep(num2str(h), '.', '') '_msemean(i,1)'])];
            end
            eval([ ['Mse_mean' strrep(num2str(h), '.', '')] ' =  Mse_mean;']);

            Mse_sd = [eval(['resp' strrep(num2str(h), '.', '') '_msesd']);eval(['resi' strrep(num2str(h), '.', '') '_msesd']);...
                eval(['fusion' strrep(num2str(h), '.', '') '_msesd']);];
            for i = 1 : max_num_aux
                Mse_sd = [Mse_sd;eval(['s_subtrans' strrep(num2str(h), '.', '') '_msesd(i,1)'])];
                Mse_sd = [Mse_sd;eval(['m_subtrans' strrep(num2str(h), '.', '') '_msesd(i,1)'])];
            end
            eval([ ['Mse_sd' strrep(num2str(h), '.', '')] ' =  Mse_sd;']);

            Sc_mean = [eval(['resp' strrep(num2str(h), '.', '') '_scmean']);eval(['resi' strrep(num2str(h), '.', '') '_scmean']);...
                eval(['fusion' strrep(num2str(h), '.', '') '_scmean']);];
            for i = 1 : max_num_aux
                Sc_mean = [Sc_mean;eval(['s_subtrans' strrep(num2str(h), '.', '') '_scmean(i,1)'])];
                Sc_mean = [Sc_mean;eval(['m_subtrans' strrep(num2str(h), '.', '') '_scmean(i,1)'])];
            end
            eval([ ['Sc_mean' strrep(num2str(h), '.', '')] ' =  Sc_mean;']);

            Sc_sd = [eval(['resp' strrep(num2str(h), '.', '') '_scsd']);eval(['resi' strrep(num2str(h), '.', '') '_scsd']);...
                eval(['fusion' strrep(num2str(h), '.', '') '_scsd']);];
            for i = 1 : max_num_aux
                Sc_sd = [Sc_sd;eval(['s_subtrans' strrep(num2str(h), '.', '') '_scsd(i,1)'])];
                Sc_sd = [Sc_sd;eval(['m_subtrans' strrep(num2str(h), '.', '') '_scsd(i,1)'])];
            end
            eval([ ['Sc_sd' strrep(num2str(h), '.', '')] ' =  Sc_sd;']);
        end
        Table_S1_1 = [];
        h_index = [0.1, 1, 2.5];
        for t = 1:length(h_index)
            h = h_index(1, t);
            Table_S1_1 = [Table_S1_1, [eval(['M_mean' strrep(num2str(h), '.', '')]), eval(['M_sd' strrep(num2str(h), '.', '')]),...
                eval(['Per' strrep(num2str(h), '.', '')]), nan*ones(length(eval(['Per' strrep(num2str(h), '.', '')])),1)]];
        end
        Table_S1_2 = [];
        h_index = [0.1, 1, 2.5];
        for t = 1:length(h_index)
            h = h_index(1, t);
            Table_S1_2 = [Table_S1_2, [eval(['Sc_mean' strrep(num2str(h), '.', '')]), eval(['Sc_sd' strrep(num2str(h), '.', '')]),...
                eval(['Mse_mean' strrep(num2str(h), '.', '')]), eval(['Mse_sd' strrep(num2str(h), '.', '')])]];
        end
        Table_S1 = [Table_S1_1; Table_S1_2];
        disp(Table_S1)
    end
end
% ----------------------------------------------------------------------------------------------
















% ------------------------------------Produce  Table S2----------------------------------------
if strcmp(OBJ,'TS2') == 1
    if load_if ~= 1
        h_index = [5, 7.5, 10];
        for t = 1:length(h_index)
            h = h_index(1, t);
            % Data generation for the target domain
            [simulation_size, X_target, Y_target, sample_size, p, beta_target_real]...
                = Data_generation_setting_1(0, h);
            % Get the estimation results for target domain using M-TRANS
            % when the signal equal to h and the number of auxilary domains is which_index
            results_dir = fullfile(NewAddress,...
                ['h=' num2str(h)]);
            [Result_table, Result_table_summary] = RESPCLUST(simulation_size, X_target, Y_target,...
                sample_size, p, beta_target_real, [], [], []);
            save(fullfile(results_dir, 'results_resp.mat'), 'Result_table', 'Result_table_summary');
        end
        h_index = [5, 7.5, 10];
        for t = 1:length(h_index)
            h = h_index(1, t);
            % Data generation for the target domain
            [simulation_size, X_target, Y_target, sample_size, p, beta_target_real]...
                = Data_generation_setting_1(0, h);
            % Get the estimation results for target domain using M-TRANS
            % when the signal equal to h and the number of auxilary domains is which_index
            results_dir = fullfile(NewAddress,...
                ['h=' num2str(h)]);
            [Result_table, Result_table_summary] = RESIDUAL(simulation_size, X_target, Y_target, sample_size, p, beta_target_real,...
                [], [], []);
            save(fullfile(results_dir, 'results_resi.mat'), 'Result_table', 'Result_table_summary');
        end
        for which_index = 0:4
            h_index = [5, 7.5, 10];
            for i = 1:length(h_index)
                h = h_index(1, i);
                % Data generation for 1 target domain and 4 auxilary domains in a loop
                [simulation_size, X_target, Y_target, sample_size, p, beta_target_real]...
                    = Data_generation_setting_1(which_index, h);
                % Get the estimation results for 1 target domain and 4 auxilary domains
                % by using FUSION method
                Aux = Auxiliary_heterogeneity_regression();
                [Results,results,Results_opt,results_opt,...
                    initial,Results_list_opt,Results_single_opt,Class_summary]...
                    =  Aux.Auxiliary_heterogeneous_regression(simulation_size,[],...
                    X_target,Y_target,sample_size,p,beta_target_real,[],[],[],...
                    [],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],...
                    [],[],[],[],[],[],[],[],[],[],[]);
                % Save the estimation results for transfer method
                results_dir = fullfile(NewAddress,...
                    ['h=' num2str(h)]);
                save(fullfile(results_dir, ['results_' num2str(which_index) '.mat']), 'Results',...
                    'results', 'Results_opt', 'results_opt', 'initial', 'Results_list_opt',...
                    'Results_single_opt', 'Class_summary');
            end
        end
        multi_if = 0;
        for which_index = 1:4
            h_index = [5, 7.5, 10];
            for t = 1:length(h_index)
                h = h_index(1, t);
                % Data generation for the target domain
                [simulation_size, X_target, Y_target, sample_size, p, beta_target_real]...
                    = Data_generation_setting_1(0, h);
                % Get the estimation results for target domain using S-TRANS
                % when the signal equal to h and the number of auxilary domains is which_index
                Trans = Transfer_heterogeneity_regression_pro();
                Auxiliary_dateset_number = which_index;
                Auxiliary = struct;
                Auxiliary(Auxiliary_dateset_number).position_provide = [];
                results_dir = fullfile(NewAddress,...
                    ['h=' num2str(h)]);
                for i = 1 : which_index
                    Auxiliary(i).position_provide = results_dir;
                end
                [Results,results,Results_opt,results_opt,...
                    initial,Results_list_opt,Results_single_opt,Class_summary]...
                    = Trans.Transfer_heterogeneous_regression_pro(Auxiliary_dateset_number,...
                    Auxiliary,multi_if,simulation_size,[],X_target,Y_target,sample_size,p,...
                    beta_target_real,[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],...
                    [],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[]);
                save(fullfile(results_dir, ['results_' strjoin(arrayfun(@num2str, 1:which_index, 'UniformOutput', false), '') '_multi_if=0.mat']),...
                    'Results', 'results', 'Results_opt', 'results_opt', 'initial',...
                    'Results_list_opt', 'Results_single_opt', 'Class_summary');
            end
        end
        multi_if = 1;
        for which_index = 1:4
            h_index = [5, 7.5, 10];
            for t = 1:length(h_index)
                h = h_index(1, t);
                % Data generation for the target domain
                [simulation_size, X_target, Y_target, sample_size, p, beta_target_real]...
                    = Data_generation_setting_1(0, h);
                % Get the estimation results for target domain using M-TRANS
                % when the signal equal to h and the number of auxilary domains is which_index
                Trans = Transfer_heterogeneity_regression_pro();
                Auxiliary_dateset_number = which_index;
                Auxiliary = struct;
                Auxiliary(Auxiliary_dateset_number).position_provide = [];
                results_dir = fullfile(NewAddress,...
                    ['h=' num2str(h)]);
                for i = 1 : which_index
                    Auxiliary(i).position_provide = results_dir;
                end
                [Results,results,Results_opt,results_opt,...
                    initial,Results_list_opt,Results_single_opt,Class_summary]...
                    = Trans.Transfer_heterogeneous_regression_pro(Auxiliary_dateset_number,...
                    Auxiliary,multi_if,simulation_size,[],X_target,Y_target,sample_size,p,...
                    beta_target_real,[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],...
                    [],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[]);
                save(fullfile(results_dir, ['results_' strjoin(arrayfun(@num2str, 1:which_index, 'UniformOutput', false), '') '_multi_if=1.mat']),...
                    'Results', 'results', 'Results_opt', 'results_opt', 'initial',...
                    'Results_list_opt', 'Results_single_opt', 'Class_summary');
            end
        end
    else
        % Load the date for Table S2
        simulation_size = 100;
        max_num_aux = 4;
        h_index = [5, 7.5, 10];
        for t = 1:length(h_index)
            h = h_index(1, t);
            results_dir = fullfile(NewAddress,...
                ['h=' num2str(h)]);
            load(fullfile(results_dir, 'results_0.mat'))
            mse_beta = [];
            for l = 1:simulation_size
                mse_beta = [mse_beta;sqrt(results_opt(l).data.mse_beta)];
            end
            eval([ ['fusion' strrep(num2str(h), '.', '') '_msemean'] ' =  mean(mse_beta);']);
            eval([ ['fusion' strrep(num2str(h), '.', '') '_msesd'] ' =  std(mse_beta);']);
            eval([ ['fusion' strrep(num2str(h), '.', '') '_scmean'] ' =  Results_single_opt.SC_mean;']);
            eval([ ['fusion' strrep(num2str(h), '.', '') '_scsd'] ' =  Results_single_opt.SC_sd;']);
            eval([['fusion' strrep(num2str(h), '.', '') '_Mmean'] ' =  Results_single_opt.K_mean;']);
            eval([['fusion' strrep(num2str(h), '.', '') '_Msd'] ' =  Results_single_opt.K_sd;']);
            eval([['fusion' strrep(num2str(h), '.', '') '_Per'] ' =  Results_single_opt.per;']);
            results_dir = fullfile(NewAddress,...
                ['h=' num2str(h)]);
            load(fullfile(results_dir, 'results_resp.mat'))
            eval([ ['resp' strrep(num2str(h), '.', '') '_msemean'] ' =  mean(sqrt(Result_table.mse_beta));']);
            eval([ ['resp' strrep(num2str(h), '.', '') '_msesd'] ' =  std(sqrt(Result_table.mse_beta));']);
            eval([ ['resp' strrep(num2str(h), '.', '') '_scmean'] ' =  mean(Result_table.sc);']);
            eval([ ['resp' strrep(num2str(h), '.', '') '_scsd'] ' =  std(Result_table.sc);']);
            eval([['resp' strrep(num2str(h), '.', '') '_Mmean'] ' =  mean(Result_table.k);']);
            eval([['resp' strrep(num2str(h), '.', '') '_Msd'] ' =  std(Result_table.k);']);
            eval([['resp' strrep(num2str(h), '.', '') '_Per'] ' =  mean(Result_table.per_if);']);
            results_dir = fullfile(NewAddress,...
                ['h=' num2str(h)]);
            load(fullfile(results_dir, 'results_resi.mat'))
            eval([ ['resi' strrep(num2str(h), '.', '') '_msemean'] ' =  mean(sqrt(Result_table.mse_beta));']);
            eval([ ['resi' strrep(num2str(h), '.', '') '_msesd'] ' =  std(sqrt(Result_table.mse_beta));']);
            eval([ ['resi' strrep(num2str(h), '.', '') '_scmean'] ' =  mean(Result_table.sc);']);
            eval([ ['resi' strrep(num2str(h), '.', '') '_scsd'] ' =  std(Result_table.sc);']);
            eval([['resi' strrep(num2str(h), '.', '') '_Mmean'] ' =  mean(Result_table.k);']);
            eval([['resi' strrep(num2str(h), '.', '') '_Msd'] ' =  std(Result_table.k);']);
            eval([['resi' strrep(num2str(h), '.', '') '_Per'] ' =  mean(Result_table.per_if);']);

            eval([['s_subtrans' strrep(num2str(h), '.', '') '_Mmean'] ' =  zeros(max_num_aux,1);']);
            eval([['s_subtrans' strrep(num2str(h), '.', '') '_Msd'] ' =  zeros(max_num_aux,1);']);
            eval([['s_subtrans' strrep(num2str(h), '.', '') '_Per'] ' =  zeros(max_num_aux,1);']);
            eval([['s_subtrans' strrep(num2str(h), '.', '') '_msemean'] ' =  zeros(max_num_aux,1);']);
            eval([['s_subtrans' strrep(num2str(h), '.', '') '_msesd'] ' =  zeros(max_num_aux,1);']);
            eval([['s_subtrans' strrep(num2str(h), '.', '') '_scmean'] ' =  zeros(max_num_aux,1);']);
            eval([['s_subtrans' strrep(num2str(h), '.', '') '_scsd'] ' =  zeros(max_num_aux,1);']);
            for which_index = 1:max_num_aux
                results_dir = fullfile(NewAddress,...
                    ['h=' num2str(h)]);
                load(fullfile(results_dir, ['results_' strjoin(arrayfun(@num2str,...
                    1:which_index, 'UniformOutput', false), '') '_multi_if=0.mat']));
                mse_beta = [];
                for l = 1:simulation_size
                    mse_beta = [mse_beta;sqrt(results_opt(l).data.mse_beta)];
                end
                eval([ ['s_subtrans' strrep(num2str(h), '.', '') '_msemean(which_index,1)'] ' =  mean(mse_beta);']);
                eval([ ['s_subtrans' strrep(num2str(h), '.', '') '_msesd(which_index,1)'] ' =  std(mse_beta);']);
                eval([ ['s_subtrans' strrep(num2str(h), '.', '') '_scmean(which_index,1)'] ' =  Results_single_opt.SC_mean;']);
                eval([ ['s_subtrans' strrep(num2str(h), '.', '') '_scsd(which_index,1)'] ' =  Results_single_opt.SC_sd;']);
                eval([ ['s_subtrans' strrep(num2str(h), '.', '') '_Mmean(which_index,1)'] ' =  Results_single_opt.K_mean;']);
                eval([ ['s_subtrans' strrep(num2str(h), '.', '') '_Msd(which_index,1)'] ' =  Results_single_opt.K_sd;']);
                eval([ ['s_subtrans' strrep(num2str(h), '.', '') '_Per(which_index,1)'] ' =  Results_single_opt.per;']);
            end
            eval([['m_subtrans' strrep(num2str(h), '.', '') '_Mmean'] ' =  zeros(max_num_aux,1);']);
            eval([['m_subtrans' strrep(num2str(h), '.', '') '_Msd'] ' =  zeros(max_num_aux,1);']);
            eval([['m_subtrans' strrep(num2str(h), '.', '') '_Per'] ' =  zeros(max_num_aux,1);']);
            eval([['m_subtrans' strrep(num2str(h), '.', '') '_msemean'] ' =  zeros(max_num_aux,1);']);
            eval([['m_subtrans' strrep(num2str(h), '.', '') '_msesd'] ' =  zeros(max_num_aux,1);']);
            eval([['m_subtrans' strrep(num2str(h), '.', '') '_scmean'] ' =  zeros(max_num_aux,1);']);
            eval([['m_subtrans' strrep(num2str(h), '.', '') '_scsd'] ' =  zeros(max_num_aux,1);']);
            for which_index = 1:max_num_aux
                results_dir = fullfile(NewAddress,...
                    ['h=' num2str(h)]);
                load(fullfile(results_dir, ['results_' strjoin(arrayfun(@num2str,...
                    1:which_index, 'UniformOutput', false), '') '_multi_if=1.mat']));
                mse_beta = [];
                for l = 1:simulation_size
                    mse_beta = [mse_beta;sqrt(results_opt(l).data.mse_beta)];
                end
                eval([ ['m_subtrans' strrep(num2str(h), '.', '') '_msemean(which_index,1)'] ' =  mean(mse_beta);']);
                eval([ ['m_subtrans' strrep(num2str(h), '.', '') '_msesd(which_index,1)'] ' =  std(mse_beta);']);
                eval([ ['m_subtrans' strrep(num2str(h), '.', '') '_scmean(which_index,1)'] ' =  Results_single_opt.SC_mean;']);
                eval([ ['m_subtrans' strrep(num2str(h), '.', '') '_scsd(which_index,1)'] ' =  Results_single_opt.SC_sd;']);
                eval([ ['m_subtrans' strrep(num2str(h), '.', '') '_Mmean(which_index,1)'] ' =  Results_single_opt.K_mean;']);
                eval([ ['m_subtrans' strrep(num2str(h), '.', '') '_Msd(which_index,1)'] ' =  Results_single_opt.K_sd;']);
                eval([ ['m_subtrans' strrep(num2str(h), '.', '') '_Per(which_index,1)'] ' =  Results_single_opt.per;']);
            end
        end
        % Load the date for Table S2
        max_num_aux = 4;
        h_index = [5, 7.5, 10];
        for t = 1:length(h_index)
            h = h_index(1, t);
            M_mean = [eval(['resp' strrep(num2str(h), '.', '') '_Mmean']);eval(['resi' strrep(num2str(h), '.', '') '_Mmean']);...
                eval(['fusion' strrep(num2str(h), '.', '') '_Mmean']);];
            for i = 1 : max_num_aux
                M_mean = [M_mean;eval(['s_subtrans' strrep(num2str(h), '.', '') '_Mmean(i,1)'])];
                M_mean = [M_mean;eval(['m_subtrans' strrep(num2str(h), '.', '') '_Mmean(i,1)'])];
            end
            eval([ ['M_mean' strrep(num2str(h), '.', '')] ' =  M_mean;']);

            M_sd = [eval(['resp' strrep(num2str(h), '.', '') '_Msd']);eval(['resi' strrep(num2str(h), '.', '') '_Msd']);...
                eval(['fusion' strrep(num2str(h), '.', '') '_Msd']);];
            for i = 1 : max_num_aux
                M_sd = [M_sd;eval(['s_subtrans' strrep(num2str(h), '.', '') '_Msd(i,1)'])];
                M_sd = [M_sd;eval(['m_subtrans' strrep(num2str(h), '.', '') '_Msd(i,1)'])];
            end
            eval([ ['M_sd' strrep(num2str(h), '.', '')] ' =  M_sd;']);

            Per = [eval(['resp' strrep(num2str(h), '.', '') '_Per']);eval(['resi' strrep(num2str(h), '.', '') '_Per']);...
                eval(['fusion' strrep(num2str(h), '.', '') '_Per']);];
            for i = 1 : max_num_aux
                Per = [Per;eval(['s_subtrans' strrep(num2str(h), '.', '') '_Per(i,1)'])];
                Per = [Per;eval(['m_subtrans' strrep(num2str(h), '.', '') '_Per(i,1)'])];
            end
            eval([ ['Per' strrep(num2str(h), '.', '')] ' =  Per;']);

            Mse_mean = [eval(['resp' strrep(num2str(h), '.', '') '_msemean']);eval(['resi' strrep(num2str(h), '.', '') '_msemean']);...
                eval(['fusion' strrep(num2str(h), '.', '') '_msemean']);];
            for i = 1 : max_num_aux
                Mse_mean = [Mse_mean;eval(['s_subtrans' strrep(num2str(h), '.', '') '_msemean(i,1)'])];
                Mse_mean = [Mse_mean;eval(['m_subtrans' strrep(num2str(h), '.', '') '_msemean(i,1)'])];
            end
            eval([ ['Mse_mean' strrep(num2str(h), '.', '')] ' =  Mse_mean;']);

            Mse_sd = [eval(['resp' strrep(num2str(h), '.', '') '_msesd']);eval(['resi' strrep(num2str(h), '.', '') '_msesd']);...
                eval(['fusion' strrep(num2str(h), '.', '') '_msesd']);];
            for i = 1 : max_num_aux
                Mse_sd = [Mse_sd;eval(['s_subtrans' strrep(num2str(h), '.', '') '_msesd(i,1)'])];
                Mse_sd = [Mse_sd;eval(['m_subtrans' strrep(num2str(h), '.', '') '_msesd(i,1)'])];
            end
            eval([ ['Mse_sd' strrep(num2str(h), '.', '')] ' =  Mse_sd;']);

            Sc_mean = [eval(['resp' strrep(num2str(h), '.', '') '_scmean']);eval(['resi' strrep(num2str(h), '.', '') '_scmean']);...
                eval(['fusion' strrep(num2str(h), '.', '') '_scmean']);];
            for i = 1 : max_num_aux
                Sc_mean = [Sc_mean;eval(['s_subtrans' strrep(num2str(h), '.', '') '_scmean(i,1)'])];
                Sc_mean = [Sc_mean;eval(['m_subtrans' strrep(num2str(h), '.', '') '_scmean(i,1)'])];
            end
            eval([ ['Sc_mean' strrep(num2str(h), '.', '')] ' =  Sc_mean;']);

            Sc_sd = [eval(['resp' strrep(num2str(h), '.', '') '_scsd']);eval(['resi' strrep(num2str(h), '.', '') '_scsd']);...
                eval(['fusion' strrep(num2str(h), '.', '') '_scsd']);];
            for i = 1 : max_num_aux
                Sc_sd = [Sc_sd;eval(['s_subtrans' strrep(num2str(h), '.', '') '_scsd(i,1)'])];
                Sc_sd = [Sc_sd;eval(['m_subtrans' strrep(num2str(h), '.', '') '_scsd(i,1)'])];
            end
            eval([ ['Sc_sd' strrep(num2str(h), '.', '')] ' =  Sc_sd;']);
        end
        Table_S2_1 = [];
        h_index = [5, 7.5, 10];
        for t = 1:length(h_index)
            h = h_index(1, t);
            Table_S2_1 = [Table_S2_1, [eval(['M_mean' strrep(num2str(h), '.', '')]), eval(['M_sd' strrep(num2str(h), '.', '')]),...
                eval(['Per' strrep(num2str(h), '.', '')]), nan*ones(length(eval(['Per' strrep(num2str(h), '.', '')])),1)]];
        end
        Table_S2_2 = [];
        h_index = [5, 7.5, 10];
        for t = 1:length(h_index)
            h = h_index(1, t);
            Table_S2_2 = [Table_S2_2, [eval(['Sc_mean' strrep(num2str(h), '.', '')]), eval(['Sc_sd' strrep(num2str(h), '.', '')]),...
                eval(['Mse_mean' strrep(num2str(h), '.', '')]), eval(['Mse_sd' strrep(num2str(h), '.', '')])]];
        end
        Table_S2 = [Table_S2_1; Table_S2_2];
        disp(Table_S2)
    end
end
% ----------------------------------------------------------------------------------------------













% ------------------------------------Produce  Figure 3----------------------------------------
if strcmp(OBJ,'F3') == 1
    if load_if ~= 1
        h_index = [0.1, 0.25, 0.5, 0.75, 1, 1.25,...
            7.5, 10, 12.5, 15, 17.5, 20];
        for target_if = 0 : 1
            if target_if == 0
                for i = 1:length(h_index)
                    h = h_index(1, i);
                    % Data generation for 12 auxilary domains in a loop
                    [simulation_size, X_target, Y_target, sample_size, p, beta_target_real]...
                        = Data_generation_setting_2(target_if, h);
                    % Get the estimation results for 12 auxilary domains
                    % by using FUSION method
                    Aux = Auxiliary_heterogeneity_regression();
                    [Results,results,Results_opt,results_opt,...
                        initial,Results_list_opt,Results_single_opt,Class_summary]...
                        =  Aux.Auxiliary_heterogeneous_regression(simulation_size,[],...
                        X_target,Y_target,sample_size,p,beta_target_real,[],[],[],...
                        [],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],...
                        [],[],[],[],[],[],[],[],[],[],[]);
                    % Save the estimation results for transfer method
                    results_dir = fullfile(NewAddress2,...
                        ['h=' num2str(h)]);
                    save(fullfile(results_dir, ['results_1' '.mat']), 'Results',...
                        'results', 'Results_opt', 'results_opt', 'initial', 'Results_list_opt',...
                        'Results_single_opt', 'Class_summary');
                end
            else
                % Data generation for the target domain
                [simulation_size, X_target, Y_target, sample_size, p, beta_target_real]...
                    = Data_generation_setting_2(1, 0);
                % Get the estimation results for the target domain
                % by using FUSION method
                Aux = Auxiliary_heterogeneity_regression();
                [Results,results,Results_opt,results_opt,...
                    initial,Results_list_opt,Results_single_opt,Class_summary]...
                    =  Aux.Auxiliary_heterogeneous_regression(simulation_size,[],...
                    X_target,Y_target,sample_size,p,beta_target_real,[],[],[],...
                    [],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],...
                    [],[],[],[],[],[],[],[],[],[],[]);
                % Save the estimation results for FUSION method
                for i = 1:length(h_index)
                    h = h_index(1, i);
                    results_dir = fullfile(NewAddress2,...
                        ['h=' num2str(h)]);
                    save(fullfile(results_dir, ['results_fusion' '.mat']), 'Results',...
                        'results', 'Results_opt', 'results_opt', 'initial', 'Results_list_opt',...
                        'Results_single_opt', 'Class_summary');
                end
            end
        end
        multi_if = 0;
        for t = 1:length(h_index)
            h = h_index(1, t);
            % Data generation for the target domain
            [simulation_size, X_target, Y_target, sample_size, p, beta_target_real]...
                = Data_generation_setting_2(1, 0);
            % Get the estimation results for target domain using S-TRANS
            % when the signal equal to h
            Trans = Transfer_heterogeneity_regression_pro();
            Auxiliary_dateset_number = 1;
            Auxiliary = struct;
            Auxiliary(Auxiliary_dateset_number).position_provide = [];
            results_dir = fullfile(NewAddress2,...
                ['h=' num2str(h)]);
            for i = 1 : 1
                Auxiliary(i).position_provide = results_dir;
            end
            [Results,results,Results_opt,results_opt,...
                initial,Results_list_opt,Results_single_opt,Class_summary]...
                = Trans.Transfer_heterogeneous_regression_pro(Auxiliary_dateset_number,...
                Auxiliary,multi_if,simulation_size,[],X_target,Y_target,sample_size,p,...
                beta_target_real,[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],...
                [],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[]);
            save(fullfile(results_dir, ['results_multi_if=0', '.mat']),...
                'Results', 'results', 'Results_opt', 'results_opt', 'initial',...
                'Results_list_opt', 'Results_single_opt', 'Class_summary');
        end
        multi_if = 1;
        for t = 1:length(h_index)
            h = h_index(1, t);
            % Data generation for the target domain
            [simulation_size, X_target, Y_target, sample_size, p, beta_target_real]...
                = Data_generation_setting_2(1, 0);
            % Get the estimation results for target domain using M-TRANS
            % when the signal equal to h
            Trans = Transfer_heterogeneity_regression_pro();
            Auxiliary_dateset_number = 1;
            Auxiliary = struct;
            Auxiliary(Auxiliary_dateset_number).position_provide = [];
            results_dir = fullfile(NewAddress2,...
                ['h=' num2str(h)]);
            for i = 1 : 1
                Auxiliary(i).position_provide = results_dir;
            end
            [Results,results,Results_opt,results_opt,...
                initial,Results_list_opt,Results_single_opt,Class_summary]...
                = Trans.Transfer_heterogeneous_regression_pro(Auxiliary_dateset_number,...
                Auxiliary,multi_if,simulation_size,[],X_target,Y_target,sample_size,p,...
                beta_target_real,[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],...
                [],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[]);
            save(fullfile(results_dir, ['results_multi_if=1', '.mat']),...
                'Results', 'results', 'Results_opt', 'results_opt', 'initial',...
                'Results_list_opt', 'Results_single_opt', 'Class_summary');
        end
        % Data generation for the target domain
        [simulation_size, X_target, Y_target, sample_size, p, beta_target_real]...
            = Data_generation_setting_2(1, 0);
        % Get the estimation results for target domain using RESPCLUST
        % when the signal equal to h
        [Result_table, Result_table_summary] = RESPCLUST(simulation_size, X_target, Y_target, sample_size, p, beta_target_real,...
            [], [], []);
        for t = 1:length(h_index)
            h = h_index(1, t);
            results_dir = fullfile(NewAddress2,...
                ['h=' num2str(h)]);
            save(fullfile(results_dir, 'results_resp.mat'), 'Result_table', 'Result_table_summary');
        end
        % Data generation for the target domain
        [simulation_size, X_target, Y_target, sample_size, p, beta_target_real]...
            = Data_generation_setting_2(1, 0);
        % Get the estimation results for target domain using RESIDUAL
        % when the signal equal to h
        [Result_table, Result_table_summary] = RESIDUAL(simulation_size, X_target, Y_target, sample_size, p, beta_target_real,...
            [], [], []);
        for t = 1:length(h_index)
            h = h_index(1, t);
            results_dir = fullfile(NewAddress2,...
                ['h=' num2str(h)]);
            save(fullfile(results_dir, 'results_resi.mat'), 'Result_table', 'Result_table_summary');
        end
    else
        h_index = [0.1, 0.25, 0.5, 0.75, 1, 1.25];
        for t = 1: length(h_index)
            h = h_index(t);
            results_dir = fullfile(NewAddress2,...
                ['h=' num2str(h)]);
            load(fullfile(results_dir, ['results_multi_if=0', '.mat']));
            eval([ ['psc_' strrep(num2str(h), '.', '') ] ' =  Results_list_opt.select_real_if;']);
            med = eval(['psc_' strrep(num2str(h), '.', '') ]);
            eval([ ['psc_median_' strrep(num2str(h), '.', '') ] ' =  median(med);']);
        end
        y = [];
        PSC = [];
        for t = 1: length(h_index)
            h = h_index(t);
            y = [y, eval(['psc_median_' strrep(num2str(h), '.', '') ])];
            PSC = [PSC, eval(['psc_' strrep(num2str(h), '.', '') ])];
        end
        figure;
        set(gcf, 'Position', [100, 100, 800, 600]);
        colors = [250,127,121]/256;
        subplot(1, 2, 1);
        handle = boxplot(PSC, 'Colors', colors, 'Widths', 0.6, 'BoxStyle',...
            'outline','MedianStyle','line','PlotStyle','traditional');
        set(handle, 'LineWidth', 1.5);
        xlim([0 7]);
        ylim([0.75, 1]);
        set(gca, 'XTickLabel', {'Y values'}, 'FontSize', 12);
        title('Boxplot', 'FontSize', 14);
        grid on;
        boxes = findobj(handle, 'Tag', 'Box');
        for j = 1:length(boxes)
            patch(get(boxes(j), 'XData'), get(boxes(j), 'YData'), colors, 'FaceAlpha', 0.6);
        end
        set(gcf, 'Color', [1 1 1]);
        hold on;
        plot(1:length(y), y, 'r--*', 'LineWidth', 2);
        h_index = [7.5, 10, 12.5, 15, 17.5, 20];
        for t = 1: length(h_index)
            h = h_index(t);
            results_dir = fullfile(NewAddress2,...
                ['h=' num2str(h)]);
            load(fullfile(results_dir, ['results_multi_if=0', '.mat']));
            eval([ ['psc_' strrep(num2str(h), '.', '') ] ' =  Results_list_opt.select_real_if;']);
            med = eval(['psc_' strrep(num2str(h), '.', '') ]);
            eval([ ['psc_median_' strrep(num2str(h), '.', '') ] ' =  median(med);']);
        end
        y = [];
        PSC = [];
        for t = 1: length(h_index)
            h = h_index(t);
            y = [y, eval(['psc_median_' strrep(num2str(h), '.', '') ])];
            PSC = [PSC, eval(['psc_' strrep(num2str(h), '.', '') ])];
        end
        colors = [250,127,121]/256;
        subplot(1, 2, 2);
        handle = boxplot(PSC, 'Colors', colors, 'Widths', 0.5, 'BoxStyle',...
            'outline','MedianStyle','line','PlotStyle','traditional');
        xlim([0 7]);
        ylim([-0.5, 0.5]);
        set(gca, 'XTickLabel', {'Y values'}, 'FontSize', 12);
        title('Boxplot', 'FontSize', 14);
        grid on;
        boxes = findobj(handle, 'Tag', 'Box');
        for j = 1:length(boxes)
            patch(get(boxes(j), 'XData'), get(boxes(j), 'YData'), colors, 'FaceAlpha', 0.6);
        end
        set(gcf, 'Color', [1 1 1]);
        hold on;
        plot(1:length(y), y, 'r--*', 'LineWidth', 2);
    end
end
% ----------------------------------------------------------------------------------------------














% ------------------------------------Produce  Table S3----------------------------------------
if strcmp(OBJ,'TS3') == 1
    if load_if ~= 1
        h_index = [0.1, 0.25, 0.5, 0.75, 1, 1.25,...
            7.5, 10, 12.5, 15, 17.5, 20];
        for target_if = 0 : 1
            if target_if == 0
                for i = 1:length(h_index)
                    h = h_index(1, i);
                    % Data generation for 12 auxilary domains in a loop
                    [simulation_size, X_target, Y_target, sample_size, p, beta_target_real]...
                        = Data_generation_setting_2(target_if, h);
                    % Get the estimation results for 12 auxilary domains
                    % by using FUSION method
                    Aux = Auxiliary_heterogeneity_regression();
                    [Results,results,Results_opt,results_opt,...
                        initial,Results_list_opt,Results_single_opt,Class_summary]...
                        =  Aux.Auxiliary_heterogeneous_regression(simulation_size,[],...
                        X_target,Y_target,sample_size,p,beta_target_real,[],[],[],...
                        [],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],...
                        [],[],[],[],[],[],[],[],[],[],[]);
                    % Save the estimation results for transfer method
                    results_dir = fullfile(NewAddress2,...
                        ['h=' num2str(h)]);
                    save(fullfile(results_dir, ['results_1' '.mat']), 'Results',...
                        'results', 'Results_opt', 'results_opt', 'initial', 'Results_list_opt',...
                        'Results_single_opt', 'Class_summary');
                end
            else
                % Data generation for the target domain
                [simulation_size, X_target, Y_target, sample_size, p, beta_target_real]...
                    = Data_generation_setting_2(1, 0);
                % Get the estimation results for the target domain
                % by using FUSION method
                Aux = Auxiliary_heterogeneity_regression();
                [Results,results,Results_opt,results_opt,...
                    initial,Results_list_opt,Results_single_opt,Class_summary]...
                    =  Aux.Auxiliary_heterogeneous_regression(simulation_size,[],...
                    X_target,Y_target,sample_size,p,beta_target_real,[],[],[],...
                    [],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],...
                    [],[],[],[],[],[],[],[],[],[],[]);
                % Save the estimation results for FUSION method
                for i = 1:length(h_index)
                    h = h_index(1, i);
                    results_dir = fullfile(NewAddress2,...
                        ['h=' num2str(h)]);
                    save(fullfile(results_dir, ['results_fusion' '.mat']), 'Results',...
                        'results', 'Results_opt', 'results_opt', 'initial', 'Results_list_opt',...
                        'Results_single_opt', 'Class_summary');
                end
            end
        end
        multi_if = 0;
        for t = 1:length(h_index)
            h = h_index(1, t);
            % Data generation for the target domain
            [simulation_size, X_target, Y_target, sample_size, p, beta_target_real]...
                = Data_generation_setting_2(1, 0);
            % Get the estimation results for target domain using S-TRANS
            % when the signal equal to h
            Trans = Transfer_heterogeneity_regression_pro();
            Auxiliary_dateset_number = 1;
            Auxiliary = struct;
            Auxiliary(Auxiliary_dateset_number).position_provide = [];
            results_dir = fullfile(NewAddress2,...
                ['h=' num2str(h)]);
            for i = 1 : 1
                Auxiliary(i).position_provide = results_dir;
            end
            [Results,results,Results_opt,results_opt,...
                initial,Results_list_opt,Results_single_opt,Class_summary]...
                = Trans.Transfer_heterogeneous_regression_pro(Auxiliary_dateset_number,...
                Auxiliary,multi_if,simulation_size,[],X_target,Y_target,sample_size,p,...
                beta_target_real,[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],...
                [],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[]);
            save(fullfile(results_dir, ['results_multi_if=0', '.mat']),...
                'Results', 'results', 'Results_opt', 'results_opt', 'initial',...
                'Results_list_opt', 'Results_single_opt', 'Class_summary');
        end
        multi_if = 1;
        for t = 1:length(h_index)
            h = h_index(1, t);
            % Data generation for the target domain
            [simulation_size, X_target, Y_target, sample_size, p, beta_target_real]...
                = Data_generation_setting_2(1, 0);
            % Get the estimation results for target domain using M-TRANS
            % when the signal equal to h
            Trans = Transfer_heterogeneity_regression_pro();
            Auxiliary_dateset_number = 1;
            Auxiliary = struct;
            Auxiliary(Auxiliary_dateset_number).position_provide = [];
            results_dir = fullfile(NewAddress2,...
                ['h=' num2str(h)]);
            for i = 1 : 1
                Auxiliary(i).position_provide = results_dir;
            end
            [Results,results,Results_opt,results_opt,...
                initial,Results_list_opt,Results_single_opt,Class_summary]...
                = Trans.Transfer_heterogeneous_regression_pro(Auxiliary_dateset_number,...
                Auxiliary,multi_if,simulation_size,[],X_target,Y_target,sample_size,p,...
                beta_target_real,[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],...
                [],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[]);
            save(fullfile(results_dir, ['results_multi_if=1', '.mat']),...
                'Results', 'results', 'Results_opt', 'results_opt', 'initial',...
                'Results_list_opt', 'Results_single_opt', 'Class_summary');
        end
        % Data generation for the target domain
        [simulation_size, X_target, Y_target, sample_size, p, beta_target_real]...
            = Data_generation_setting_2(1, 0);
        % Get the estimation results for target domain using RESPCLUST
        % when the signal equal to h
        [Result_table, Result_table_summary] = RESPCLUST(simulation_size, X_target, Y_target, sample_size, p, beta_target_real,...
            [], [], []);
        for t = 1:length(h_index)
            h = h_index(1, t);
            results_dir = fullfile(NewAddress2,...
                ['h=' num2str(h)]);
            save(fullfile(results_dir, 'results_resp.mat'), 'Result_table', 'Result_table_summary');
        end
        % Data generation for the target domain
        [simulation_size, X_target, Y_target, sample_size, p, beta_target_real]...
            = Data_generation_setting_2(1, 0);
        % Get the estimation results for target domain using RESIDUAL
        % when the signal equal to h
        [Result_table, Result_table_summary] = RESIDUAL(simulation_size, X_target, Y_target, sample_size, p, beta_target_real,...
            [], [], []);
        for t = 1:length(h_index)
            h = h_index(1, t);
            results_dir = fullfile(NewAddress2,...
                ['h=' num2str(h)]);
            save(fullfile(results_dir, 'results_resi.mat'), 'Result_table', 'Result_table_summary');
        end
    else
        h = 0.1;
        results_dir = fullfile(NewAddress2,...
            ['h=' num2str(h)]);
        load(fullfile(results_dir, 'results_resp.mat'));
        M_mean = [];M_sd = [];Per = [];SC_mean = [];SC_sd = [];MSE_mean = [];MSE_sd = [];Pcs_mean = [];Pcs_sd = [];
        M_mean = [M_mean;Result_table_summary.K_mean];
        M_sd = [M_sd;Result_table_summary.K_sd];
        Per = [Per;Result_table_summary.per];
        SC_mean = [SC_mean;Result_table_summary.SC_mean];
        SC_sd = [SC_sd;Result_table_summary.SC_sd];
        MSE_mean = [MSE_mean;mean(sqrt(Result_table.mse_beta))];
        MSE_sd = [MSE_sd;std(sqrt(Result_table.mse_beta))];
        Pcs_mean = [Pcs_mean;nan];
        Pcs_sd = [Pcs_sd;nan];
        results_dir = fullfile(NewAddress2,...
            ['h=' num2str(h)]);
        load(fullfile(results_dir, 'results_resi.mat'));
        M_mean = [M_mean;Result_table_summary.K_mean];
        M_sd = [M_sd;Result_table_summary.K_sd];
        Per = [Per;Result_table_summary.per];
        SC_mean = [SC_mean;Result_table_summary.SC_mean];
        SC_sd = [SC_sd;Result_table_summary.SC_sd];
        MSE_mean = [MSE_mean;mean(sqrt(Result_table.mse_beta))];
        MSE_sd = [MSE_sd;std(sqrt(Result_table.mse_beta))];
        Pcs_mean = [Pcs_mean;nan];
        Pcs_sd = [Pcs_sd;nan];
        results_dir = fullfile(NewAddress2,...
            ['h=' num2str(h)]);
        load(fullfile(results_dir, 'results_fusion.mat'));
        M_mean = [M_mean;Results_single_opt.K_mean];
        M_sd = [M_sd;Results_single_opt.K_sd];
        Per = [Per;Results_single_opt.per];
        SC_mean = [SC_mean;Results_single_opt.SC_mean];
        SC_sd = [SC_sd;Results_single_opt.SC_sd];
        MSE_mean = [MSE_mean;mean(sqrt(Results_list_opt.mse_beta))];
        MSE_sd = [MSE_sd;std(sqrt(Results_list_opt.mse_beta))];
        Pcs_mean = [Pcs_mean;nan];
        Pcs_sd = [Pcs_sd;nan];
        h_index = [0.1, 0.25, 0.5, 0.75, 1, 1.25,...
            7.5, 10, 12.5, 15, 17.5, 20];
        for t = 1:length(h_index)
            h = h_index(t);
            results_dir = fullfile(NewAddress2,...
                ['h=' num2str(h)]);
            load(fullfile(results_dir, 'results_multi_if=0.mat'));
            M_mean = [M_mean;Results_single_opt.K_mean];
            M_sd = [M_sd;Results_single_opt.K_sd];
            Per = [Per;Results_single_opt.per];
            SC_mean = [SC_mean;Results_single_opt.SC_mean];
            SC_sd = [SC_sd;Results_single_opt.SC_sd];
            MSE_mean = [MSE_mean;mean(sqrt(Results_list_opt.mse_beta))];
            MSE_sd = [MSE_sd;std(sqrt(Results_list_opt.mse_beta))];
            Pcs_mean = [Pcs_mean;Results_single_opt.select_real_if_mean];
            Pcs_sd = [Pcs_sd;Results_single_opt.select_real_if_sd];
        end
        Table_S3 = [M_mean, M_sd, Per, SC_mean, SC_sd, MSE_mean, MSE_sd, Pcs_mean, Pcs_sd];
        disp(Table_S3);
    end
end