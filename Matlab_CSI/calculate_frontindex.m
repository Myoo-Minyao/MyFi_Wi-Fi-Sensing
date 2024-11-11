function [min_index] = calculate_frontindex(csi_phase_butt_unwrapped)
%CALCULATE_FRONTINDEX 计算第运动开始值
%   此处显示详细说明
% 计算前半部分的长度
half_length = floor(length(csi_phase_butt_unwrapped) / 2);

% 获取前半部分数据
csi_phase_first_half = csi_phase_butt_unwrapped(1:half_length);

% 找到前半部分的最小值及其索引
[min_value, min_index] = min(csi_phase_first_half);

% 根据最小值索引获取后半部分数据
csi_phase_second_half = csi_phase_butt_unwrapped(min_index:half_length);

% 在后半部分找到最大值及其索引
[max_value, max_index] = max(csi_phase_second_half);

% 将最小值索引映射到原始数组中
min_index_global = min_index;

% 将最大值索引映射到全局索引中
max_index_global = min_index_global + max_index - 1;

% % 打印结果
% disp(['前半部分最小值：', num2str(min_value)]);
% disp(['前半部分最小值索引：', num2str(min_index_global)]);
% disp(['从最小值开始的后半部分最大值：', num2str(max_value)]);
% disp(['后半部分最大值索引：', num2str(max_index_global)]);
end

