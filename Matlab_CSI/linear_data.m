function [csi_ratio_linear] = linear_data(ant,csi_ratio)
%线性插值，对csi_ratio数据进行线性插值
%   此处显示详细说明
x = 1:length(ant);% 使用 1 到 ant1 的整数作为横坐标 x
%删除了NaN值
csi_ratio_linear = interp1(x(~isnan(x)),csi_ratio(~isnan(x), :), x);
% 附加 NaN 结尾，以确保数据流的末尾没有 NaN 值
idx = length(x);
% 计算每行中包含 NaN 值的数量
nanflag = sum(~isnan(csi_ratio_linear), 2);
while ~nanflag(idx, 1)
    idx = idx - 1;
end
csi_ratio_linear(idx, :) = csi_ratio_linear(idx-1, :);
% 寻找最后一个不包含 NaN 值的行
while idx ~= length(x)
    csi_ratio_linear(idx+1, :) = csi_ratio_linear(idx, :);%用前一个值去填充NaN值
    idx = idx + 1;
end
end

