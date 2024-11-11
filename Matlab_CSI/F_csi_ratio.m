function [csi_ratio,level] = F_csi_ratio(ant1,ant2,ant3)
%计算csi_ratio值，并筛选天线
%   此处显示详细说明
%去除三根天线中无穷大值和0值点
ant1 = removecomplexzero( ant1 );
ant2 = removecomplexzero( ant2 );
ant3 = removecomplexzero( ant3 );

% 天线选择，其中用均值/标准差更关注信噪比，均值更关注信号的强度，方差更关注信号的稳定性
%求均值
mean_ant1 = mean(mean(abs(ant1)));
mean_ant2 = mean(mean(abs(ant2)));
mean_ant3 = mean(mean(abs(ant3)));
%求标准差
var_ant1 = mean(sqrt(var(abs(ant1))));
var_ant2 = mean(sqrt(var(abs(ant2))));
var_ant3 = mean(sqrt(var(abs(ant3))));
%求均值/标准差
ant1_mean_var_ratio = mean_ant1./var_ant1;
ant2_mean_var_ratio = mean_ant2./var_ant2;
ant3_mean_var_ratio = mean_ant3./var_ant3;
level = [ant1_mean_var_ratio ant2_mean_var_ratio ant3_mean_var_ratio];

% Create Gramm object
g = gramm('x', {'n1', 'n2', 'n3'}, 'y', [ant1_mean_var_ratio, ant2_mean_var_ratio, ant3_mean_var_ratio]);

% Plot bar graph with blue color
g.geom_bar('facecolor', [0.173 0.102 0.75]);

% Add axis labels and title, increase label font size
g.set_names('x', ' ', 'y', 'Mean-to-Variance Ratio');
g.axe_property('FontSize', 25); % Increase font size

% Add grid lines
g.axe_property('GridColor', [0.5 0.5 0.5]); % Gray color for grid lines
g.axe_property('GridAlpha', 0.5); % Set transparency for grid lines

% Bolden axis labels
g.axe_property('FontWeight', 'bold'); % Set font weight to bold for axis labels


% Display the plot
figure;
g.draw();

%求均值/标准差最大的天线
[~,max_index]=max (level);
[~,min_index]=min (level);

csi = [ant1,ant2,ant3];
%性能最好的和最差的天线的csi值
csi_max = csi(:,(max_index-1)*30+1:(max_index)*30);
csi_min = csi(:,(min_index-1)*30+1:(min_index)*30);
%csi_ratio
csi_ratio = csi_min ./ csi_max;
end

