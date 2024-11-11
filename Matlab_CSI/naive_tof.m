function [tof_mat] = naive_tof(csi_data)
    % naive_tof
    % Input:
    %   - csi_data is the CSI used for ranging; [T S A L]
     %   - csi_data 是用于测距的CSI数据; [T S A L]
    %     T：时间实例数
    %     S：子载波数
    %     A：天线数
    %     L：数据包数
    % Output:
    %   - tof_mat is the rough time-of-flight estimation result; [T A]
     %   - tof_mat 是粗略的飞行时间估计结果; [T A]
    %     T：时间实例数
    %     A：天线数

    % The bandwidth parameter.% 带宽参数。
    global bw;
    [~, subcarrier_num, ~, ~] = size(csi_data);
    % Exponential powers of 2, based on the rounded up subcarrier number.
     % 计算向上取整后的子载波数的2的幂。
    ifft_point = power(2, ceil(log2(subcarrier_num)));
    % Get CIR from each packet and each antenna by ifft(CFR);
     % 通过ifft(CFR)获取每个数据包和每个天线的CIR。
    cir_sequence = ifft(csi_data, ifft_point, 2); % [T ifft_point A L]
    cir_sequence = squeeze(mean(cir_sequence, 4)); % [T ifft_point A]
    % Only consider half of the ifft points.
    % 仅考虑ifft点的一半
    half_point = ifft_point / 2;
    half_sequence = cir_sequence(:, 1:half_point, :); % [T half_point A]
    % Find the peak of the CIR sequence.
     % 找到CIR序列的峰值。
    [~, peak_indices] = max(half_sequence, [], 2); % [T 1 A]
    peak_indices = squeeze(peak_indices); % [T A]
    % Calculate ToF for each packet and each antenna, based on the CIR peak.
    % 基于CIR峰值，计算每个数据包和每个天线的ToF。
    tof_mat = peak_indices .* subcarrier_num ./ (ifft_point .* bw); % [T A]
end