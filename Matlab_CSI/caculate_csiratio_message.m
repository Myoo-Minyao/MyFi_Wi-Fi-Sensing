function [csi_ratio_amplitude,csi_ratio_phase_unwrapped,csi_ratio_phase] = caculate_csiratio_message(csi_ratio,time1,start_time,end_time,j,i,k)
%CACULATE_CSIRATIO_MESSAGE 计算csi_ratio值的振幅比和相位差。

csi_ratio_amplitude = abs(csi_ratio);
csi_ratio_phase = angle(csi_ratio);
csi_ratio_phase_unwrapped = unwrap(csi_ratio_phase, pi, 2);
% 找到对应的索引
start_index = find(time1 >= start_time, 1);
end_index = find(time1 <= end_time, 1, 'last');

if j==0
    csi_amplitude_csiratio_subset = csi_ratio_amplitude;
    csi_phase_unwrapped_csiratio_subset =csi_ratio_phase_unwrapped;
    if k==1
        figure;
        plot(time1, csi_amplitude_csiratio_subset(:, i))
        ylabel('CSI Amplitude(dB)', 'FontWeight', 'bold');
        xlabel('Time(s)', 'FontWeight', 'bold');
    elseif k==2
        figure;
        plot(time1, csi_phase_unwrapped_csiratio_subset(:, i))
        ylabel('Phase(radian)', 'FontWeight', 'bold');
        xlabel('Time(s)', 'FontWeight', 'bold');
    elseif k==3
        figure;
        plot(time1, csi_amplitude_csiratio_subset(:, i))
        ylabel('CSI Amplitude(dB)', 'FontWeight', 'bold');
        xlabel('Time(s)', 'FontWeight', 'bold');
        figure;
        plot(time1, csi_phase_unwrapped_csiratio_subset(:, i))
        ylabel('Phase(radian)', 'FontWeight', 'bold');
        xlabel('Time(s)', 'FontWeight', 'bold');
    end
elseif j==1
    % 截取CSI数据
    csi_amplitude_csiratio_subset = csi_ratio_amplitude(start_index:end_index, :);
    csi_phase_unwrapped_csiratio_subset =csi_ratio_phase_unwrapped(start_index:end_index, :);
    if k==1
        figure;
        plot(time1(start_index:end_index), csi_amplitude_csiratio_subset(:, i))
        ylabel('CSI Amplitude(dB)', 'FontWeight', 'bold');
        xlabel('Time(s)', 'FontWeight', 'bold');
    elseif k==2
        figure;
        plot(time1(start_index:end_index), csi_phase_unwrapped_csiratio_subset(:, i))
        ylabel('Phase(radian)', 'FontWeight', 'bold');
        xlabel('Time(s)', 'FontWeight', 'bold');
    elseif k==3
        figure;
        plot(time1(start_index:end_index), csi_amplitude_csiratio_subset(:, i))
        ylabel('CSI Amplitude(dB)', 'FontWeight', 'bold');
        xlabel('Time(s)', 'FontWeight', 'bold');
        figure;
        plot(time1(start_index:end_index), csi_phase_unwrapped_csiratio_subset(:, i))
        ylabel('Phase(radian)', 'FontWeight', 'bold');
        xlabel('Time(s)', 'FontWeight', 'bold');
    end
end
