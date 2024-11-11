function [period, locs, diffs] = autocorrelation(detect_period)
    remove_average = detect_period - mean(detect_period);
    x = 1:length(remove_average);
    [~,locs] = findallpeaksdistance(x,remove_average,min(remove_average)-1,400);
% 计算信号的自相关函数，找到周期性并返回周期、峰值位置和峰值之间的差异。
% 参数说明:
% - detect_period: 输入信号
%
% 返回值:
% - period: 信号的周期
% - locs: 自相关函数中的峰值位置
% - diffs: 峰值之间的差异
    %     % calculate autocorrelation
    % 从输入信号中去除均值
%     remove_average = detect_period - mean(detect_period);
%     time_len = length(remove_average);
%     % element-wise multiplication
%     R = zeros(1, time_len);
%     for pau=1:time_len
%         R(1, pau) = sum(remove_average(1:time_len-pau+1) .* remove_average(pau:time_len));
%     end
% %     plot(R)
%    
%     % [~, locs] = min(R); %     % [~, locs] = min(R); % 其实我觉得按论文里，这个应该是locs后面那个再上去的峰值
%     % 按照论文中，半周期应该是第二个峰值点
% %     [~, locs] = findpeaks(R);
% 计算峰值位置
%     x = [1:length(R)];
%     P = findvalleys(x, R, 0, min(R)-1, 1, 10, 1); % 250 is the half width of period
%     period = P(2, 2) - P(1, 2);
% 如果找到的峰值数量不足2个，则无法计算周期
    if length(locs) <= 1
        period = 0;
        diffs = [0, 0];
    else
         % 计算峰值之间的差异
        diffs = diff(locs);
        % 计算周期，取峰值之间差异的平均值
        period = mean(diffs);
    end
end

