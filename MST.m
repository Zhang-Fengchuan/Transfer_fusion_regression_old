%%
rng(123)
p = 13;
n = 112;
index = [1,3];
load('C:\Users\张丰川\Documents\MATLAB\Transfer_fusion_regression\real_data\Result\auxiliary\results_1.mat')
point1 = Class_summary.coef{1, 1};
point2 = Class_summary.coef{1, 2};
point3 = Class_summary.coef{1, 3};
point4 = Class_summary.coef{1, 4};
load('C:\Users\张丰川\Documents\MATLAB\Transfer_fusion_regression\real_data\Result\auxiliary\results_2.mat')
point5 = Class_summary.coef{1, 1};
point6 = Class_summary.coef{1, 2};
point7 = Class_summary.coef{1, 3};
load('C:\Users\张丰川\Documents\MATLAB\Transfer_fusion_regression\real_data\Result\auxiliary\results_0.mat')
pointb1 = Class_summary.coef{1, 1};
pointb2 = Class_summary.coef{1, 2};
load('C:\Users\张丰川\Documents\MATLAB\Transfer_fusion_regression\real_data\Result\transfer\results_0.mat')
pointl1 = Class_summary.coef{1, 1};
pointl2 = Class_summary.coef{1, 2};
pointl3 = Class_summary.coef{1, 3};
points_whole = [point1'; point2'; point3'; point4'; point5'; point6'; point7'];
distances_whole = pdist(points_whole);
dist_matrix_whole = squareform(distances_whole);
G_whole = graph(dist_matrix_whole);
T_whole = minspantree(G_whole);
points = [point1(index ,1)'; point2(index ,1)'; point3(index ,1)'; point4(index ,1)';...
            point5(index ,1)'; point6(index ,1)'; point7(index ,1)'];
figure;
h = plot(T_whole, 'XData', points(:,1), 'YData', points(:,2),...
    'LineWidth', 3, 'LineStyle','-');
hold on;
scatter(points(:,1), points(:,2), 80, 'o', 'filled');
hold on
cut = 1*sqrt(p*log(p)/n);
%highlight_edge_idx = T_whole.Edges{:,1}(find(T_whole.Edges{:,2}>cut),:); % 需要改变的边的索引
highlight(h, 'Edges', find(T_whole.Edges{:,2}>cut), 'LineStyle', '--', 'LineWidth', 2);
hold on
pointb = [pointb1(index ,1)'; pointb2(index ,1)'];
scatter(pointb(:,1), pointb(:,2), 80, 'g', 'filled');
pointl = [pointl1(index ,1)'; pointl2(index ,1)'; pointl3(index ,1)'];
scatter(pointl(:,1), pointl(:,2), 80, 'b', 'filled');
hold on
xlabel('是否手术');
ylabel('是否化疗');









%%
% 生成三个簇中心
centers = [1.5 1; 4.8 6.5; 8 2];

% 每个簇生成10个随机点
num_points = 8;
points1 = centers(1,:) + 0.6*randn(num_points, 2);
points2 = centers(2,:) + 0.6*randn(num_points, 2);
points3 = centers(3,:) + 0.6*randn(num_points, 2);

% 合并所有点
points = [points1; points2; points3];
points(5,1) = points(5,1) + 0.5;
points = [points; [0.6-1, 5-0.5]; [0.5-1,4.3-0.5]];
% 计算距离矩阵
distances = pdist(points);

% 生成距离矩阵的方阵形式
dist_matrix = squareform(distances);

% 创建图结构
G = graph(dist_matrix);

% 计算最小生成树
T = minspantree(G);

% 绘制最小生成树
figure;
h = plot(T, 'XData', points(:,1), 'YData', points(:,2), 'LineWidth', 3, 'LineStyle','-');
%h.NodeLabel = {}; % 不显示点的标号
hold on;
%line_color = [130,176,210]/256;
%h.EdgeColor = line_color;
scatter(points(:,1), points(:,2), 20, 'o', 'filled');
% hold on
% scatter(0.8, 4.5, 50,'r', 'filled');
title('Minimum Spanning Tree of Points');
xlabel('X');
ylabel('Y');
hold on;
scatter(centers(:,1), centers(:,2), 100,'g','filled','o')

row_indices1 = find(ismember(minspantree(G).Edges.EndNodes, [9,25], 'rows'));
row_indices2 = find(ismember(minspantree(G).Edges.EndNodes, [8,26], 'rows'));
row_indices3 = find(ismember(minspantree(G).Edges.EndNodes, [15,21], 'rows'));

highlight_edge_idx = [row_indices1;row_indices2;row_indices3]; % 需要改变的边的索引
highlight(h, 'Edges', highlight_edge_idx, 'LineStyle', '--', 'LineWidth', 2);
hold on


index_1 = [1,2,3,4,5,6,7,8];
points_clust_1 = points(index_1,:);

index_2 = [9,10,11,12,13,14,15,16];
points_clust_2 = points(index_2,:);

index_3 = [17,18,19,20,21,22,23,24];
points_clust_3 = points(index_3,:);

index_4 = [25,26];
points_clust_4 = points(index_4,:);

center_clust_1 = mean(points_clust_1);
center_clust_2 = mean(points_clust_2);
center_clust_3 = mean(points_clust_3);
center_clust_4 = mean(points_clust_4);
center_clust = [center_clust_1;center_clust_2;center_clust_3;center_clust_4];

scatter(center_clust(:,1), center_clust(:,2), 200,'r','filled','p')


ols = [2.2,1;4.6,2;3.2,4.95];
color_ols = [20, 184, 53]/256;
scatter(ols(:,1), ols(:,2), 100, color_ols,'filled','o')





% 设置特定点的坐标和圆的半径
center_point = centers(1,:); % 圆心坐标
radius = 1.3; % 圆的半径

% 绘制圆
theta = linspace(0, 2*pi, 100);
x_circle = center_point(1) + radius * cos(theta);
y_circle = center_point(2) + radius * sin(theta);
plot(x_circle, y_circle, 'g--', 'LineWidth', 1); % 'b-'表示蓝色实线
hold on





% 设置特定点的坐标和圆的半径
center_point = centers(2,:); % 圆心坐标
radius = 1.3; % 圆的半径

% 绘制圆
theta = linspace(0, 2*pi, 100);
x_circle = center_point(1) + radius * cos(theta);
y_circle = center_point(2) + radius * sin(theta);
plot(x_circle, y_circle, 'g--', 'LineWidth', 1); % 'b-'表示蓝色实线
hold on






% 设置特定点的坐标和圆的半径
center_point = centers(3,:); % 圆心坐标
radius = 1.3; % 圆的半径

% 绘制圆
theta = linspace(0, 2*pi, 100);
x_circle = center_point(1) + radius * cos(theta);
y_circle = center_point(2) + radius * sin(theta);
plot(x_circle, y_circle, 'g--', 'LineWidth', 1); % 'b-'表示蓝色实线
hold on

axis equal


