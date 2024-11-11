function [csi_ratio_linear_filtering] = bandpass_butter(csi_ratio_linear,half_sample_rate,upper_f,upper_order,lower_order,lower_f)
%butterworth滤波器设计，过滤csi_ratio 数据
% 滤波器设计
[b,a] = butter(upper_order,upper_f/half_sample_rate,"low" );%低通滤波器设计
[d,c] = butter(lower_order,lower_f/half_sample_rate,'high');% 高通滤波器设计

% 遍历数据的每一列（每个信号通道），并对每列应用滤波器
for jj = 1:size(csi_ratio_linear,2)
    csi_ratio_linear(:,jj) = filter(b,a,csi_ratio_linear(:,jj));% 使用上通滤波器
    csi_ratio_linear(:,jj) = filter(d,c,csi_ratio_linear(:,jj));% 使用下通滤波器
end
csi_ratio_linear_filtering=csi_ratio_linear;
end

