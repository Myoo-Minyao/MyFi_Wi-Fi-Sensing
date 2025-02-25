function [csi_amplitude,csi_phase_unwrapped,csi_phase] = caculate_rawcsi_message(cfr_array,time1,start_time,end_time,j,i,k)
%RAW_CSI_MESSAGE 原始数据的振幅和相位信息
% start_time,end_time为切片时间段，i为子载波序列
% 计算振幅信息
csi_amplitude = abs(cfr_array);
% 计算相位信息
csi_phase = angle(cfr_array);
%csi_phase_unwrapped = csi_phase;
csi_phase_unwrapped = unwrap(csi_phase);

% 找到对应的索引
start_index = find(time1 >= start_time, 1);
end_index = find(time1 <= end_time, 1, 'last');

if j==0
    csi_amplitude_subset = csi_amplitude;
    csi_phase_unwrapped_subset = csi_phase_unwrapped;
    if k==1
        figure;
        plot(time1, csi_amplitude_subset(:, i))
        ylabel('CSI Amplitude(dB)', 'FontWeight', 'bold');
        xlabel('Time(s)', 'FontWeight', 'bold');
    elseif k==2
        figure;
        plot(time1, csi_phase_unwrapped_subset(:, i))
        ylabel('Phase(radian)', 'FontWeight', 'bold');
        xlabel('Time(s)', 'FontWeight', 'bold');
    elseif k==3
        figure;
        plot(time1, csi_amplitude_subset(:, i))
        ylabel('CSI Amplitude(dB)', 'FontWeight', 'bold');
        xlabel('Time(s)', 'FontWeight', 'bold');
        figure;
        plot(time1, csi_phase_unwrapped_subset(:, i))
        ylabel('Phase(radian)', 'FontWeight', 'bold');
        xlabel('Time(s)', 'FontWeight', 'bold');
    end
elseif j==1
    % 截取CSI数据
    csi_amplitude_subset = csi_amplitude(start_index:end_index, :);
    csi_phase_unwrapped_subset = csi_phase_unwrapped(start_index:end_index, :);
    if k==1
        figure;
        plot(time1(start_index:end_index), csi_amplitude_subset(:, i))
        ylabel('CSI Amplitude(dB)', 'FontWeight', 'bold');
        xlabel('Time(s)', 'FontWeight', 'bold');
    elseif k==2
        figure;
        plot(time1(start_index:end_index), csi_phase_unwrapped_subset(:, i))
        ylabel('Phase(radian)', 'FontWeight', 'bold');
        xlabel('Time(s)', 'FontWeight', 'bold');
    elseif k==3
        figure;
        plot(time1(start_index:end_index), csi_amplitude_subset(:, i))
        ylabel('CSI Amplitude(dB)', 'FontWeight', 'bold');
        xlabel('Time(s)', 'FontWeight', 'bold');
        figure;
        plot(time1(start_index:end_index), csi_phase_unwrapped_subset(:, i))
        ylabel('Phase(radian)', 'FontWeight', 'bold');
        xlabel('Time(s)', 'FontWeight', 'bold');
    end

end
end

