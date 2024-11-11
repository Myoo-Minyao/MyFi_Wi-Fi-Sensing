function [fc_contour] = movement_speed(doppler_spectrum_squeezed,time1,doppler_frequencies,freq_axis_bin)
%移动速度
timebin_len=size(time1,2);
freqbin_len=size(doppler_frequencies,2);
% 存储移动轮廓频率的集合
fc_contour = [];fc_contour = [];
%找到满足功率大于总功率3%的频点
for i=1:timebin_len
    % 满足要求的某一个时间点的freq points集合
    contour_f_points = [];
    power_sum = sum(doppler_spectrum_squeezed(:, i));% 计算总功率
    for k=1:freqbin_len
        power_point = doppler_spectrum_squeezed(k, i);
        if power_point / power_sum > 0.03% 当前频率点的功率
            % 将当前频率点添加到频率轮廓集合中
            contour_f_points = [contour_f_points, freq_axis_bin(k)];
        end

    end

    if size(contour_f_points) == 0
        % 将0添加到人体移动频率轮廓集合中
        fc_contour = [fc_contour, 0];
    else
        % % 提取符合条件的频率最大值,轮廓
        fc_contour = [fc_contour, max(contour_f_points)];
    end

end


end

