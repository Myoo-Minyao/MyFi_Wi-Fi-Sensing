function [unbias_time_seq,csi_seq] = time(time_duration,cfr_array)
% 获取时序并筛选把错误点删除
%time1 = (1:size(time_duration))./200;
 time_span = [0, inf];
 unbias_time_seq = (double(time_duration) - double(time_duration(1)))./10^6; %和time1 是转置关系
 indicator_vec = unbias_time_seq <= time_span(2) & unbias_time_seq >= time_span(1);%找到时序错误的时间点
 indicator_vec =  ~indicator_vec;%筛选出错误时间点
 unbias_time_seq(indicator_vec) = [];%将时序错误点进行删除
 cfr_array(:,indicator_vec) = [];%将csi中时序错误点删除
 unbias_time_seq = unbias_time_seq - unbias_time_seq(1); %删除后的时序长度
 csi_seq = cfr_array;
end

