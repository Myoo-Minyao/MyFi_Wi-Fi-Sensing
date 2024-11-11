function [min_index_global_1] = calculate_backindex(csi_phase_butt_unwrapped)
%CALCULATE_BACKINDEX 计算运动结束时
%   此处显示详细说明
% 获取后半部分数据
all_length = floor(length(csi_phase_butt_unwrapped));
half_length = floor(length(csi_phase_butt_unwrapped) / 2);
csi_phase_second_half = csi_phase_butt_unwrapped(half_length:all_length);

% 找到后半部分的最小值及其索引
[min_value_1, min_index_1] = min(csi_phase_second_half);
min_index_1=min_index_1+half_length;
% 根据最小值索引获取后半部分数据
csi_phase_second_half = csi_phase_butt_unwrapped(half_length:min_index_1);

% 在后半部分找到最大值及其索引
[max_value_1, max_index_1] = max(csi_phase_second_half);

% 将最小值索引映射到原始数组中
min_index_global_1 = min_index_1;

% 将最大值索引映射到全局索引中
max_index_global_1 =  half_length + max_index_1 ;

% % 打印结果
% disp(['后半部分最小值：', num2str(min_value_1)]);
% disp(['后半部分最小值索引：', num2str(min_index_global_1)]);
% disp(['从中间值到最小值后半部分最大值：', num2str(max_value_1)]);
% disp(['后半部分最大值索引：', num2str(max_index_global_1)]);
end

