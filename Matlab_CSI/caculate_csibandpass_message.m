function [csi_amplitude_butt,csi_phase_butt_unwrapped] = caculate_csibandpass_message(csi_ratio_linear_filtering,time1,start_time,end_time,j,i,k)
%caculate_csibandpass_message 计算滤波后的振幅和相位信息
% 计算滤波后振幅信息
csi_amplitude_butt = abs(csi_ratio_linear_filtering);
% 计算滤波后相位信息
csi_phase_butt = angle(csi_ratio_linear_filtering);
csi_phase_butt_unwrapped = unwrap(csi_phase_butt);
% 找到对应的索引
start_index = find(time1 >= start_time, 1);
end_index = find(time1 <= end_time, 1, 'last');

if j==0
    csi_amplitude_csiratio_butt_subset = csi_amplitude_butt;
    csi_phase_unwrapped_csiratio_butt_subset =csi_phase_butt_unwrapped;
    if k==1
        figure;
        plot(time1, csi_amplitude_csiratio_butt_subset(:, i))
        ylabel('CSI Amplitude(dB)', 'FontWeight', 'bold');
        xlabel('Time(s)', 'FontWeight', 'bold');
    elseif k==2
        figure;
        plot(time1, csi_phase_unwrapped_csiratio_butt_subset(:, i))
        ylabel('Phase(radian)', 'FontWeight', 'bold');
        xlabel('Time(s)', 'FontWeight', 'bold');
    elseif k==3
        figure;
        plot(time1, csi_amplitude_csiratio_butt_subset(:, i))
        ylabel('CSI Amplitude(dB)', 'FontWeight', 'bold');
        xlabel('Time(s)', 'FontWeight', 'bold');
        figure;
        plot(time1, csi_phase_unwrapped_csiratio_butt_subset(:, i))
        ylabel('Phase(radian)', 'FontWeight', 'bold');
        xlabel('Time(s)', 'FontWeight', 'bold');
    end
elseif j==1
    % 截取CSI数据
    csi_amplitude_csiratio_butt_subset = csi_amplitude_butt(start_index:end_index, :);
    csi_phase_unwrapped_csiratio_butt_subset =csi_phase_butt_unwrapped(start_index:end_index, :);
    if k==1
        figure;
        plot(time1(start_index:end_index), csi_amplitude_csiratio_butt_subset(:, i))
        ylabel('CSI Amplitude(dB)', 'FontWeight', 'bold');
        xlabel('Time(s)', 'FontWeight', 'bold');
    elseif k==2
        figure;
        plot(time1(start_index:end_index), csi_phase_unwrapped_csiratio_butt_subset(:, i))
        ylabel('Phase(radian)', 'FontWeight', 'bold');
        xlabel('Time(s)', 'FontWeight', 'bold');
    elseif k==3
        figure;
        plot(time1(start_index:end_index), csi_amplitude_csiratio_butt_subset(:, i))
        ylabel('CSI Amplitude(dB)', 'FontWeight', 'bold');
        xlabel('Time(s)', 'FontWeight', 'bold');
        figure;
        plot(time1(start_index:end_index), csi_phase_unwrapped_csiratio_butt_subset(:, i))
        ylabel('Phase(radian)', 'FontWeight', 'bold');
        xlabel('Time(s)', 'FontWeight', 'bold');
    end
end

