function [estimated_steady_footstep1,estimated_steady_footstep2,time1,doppler_frequencies,rearranged_array] = Mult_data(path)

%% 加载数据文件
[cfr_array,ant1,ant2,ant3,Subcarriers,Tx,Rx,time_duration] = csi_loading(path);
%% 计算CSI信息
[csi_amplitude,csi_phase_unwrapped] = caculate_rawcsi_message(cfr_array);
time1 = (1:size(time_duration))./1000;

%% csi ratio 去除硬件噪声
[csi_ratio,level] = F_csi_ratio(ant1,ant2,ant3);
%% 计算CSI_ratio值
[csi_ratio_amplitude,csi_ratio_phase_unwrapped] = caculate_csiratio_message(csi_ratio);
%% 线性插值
[csi_ratio_linear] = linear_data(ant1,csi_ratio);
% 此时的csi ratio是去除了硬件噪声和空值的数据了
%% 画插值后的图

csi_ratio_linear_amplitude = abs(csi_ratio_linear);
csi_ratio_linear_phase = angle(csi_ratio_linear);
csi_ratio_linear_phase_unwrapped = unwrap(csi_ratio_linear_phase);

%% 设计带通滤波器(Butterworth 滤波器) 去除带外噪声
sample_rate = 1000;% 采样率
half_sample_rate = 0.5*sample_rate;% 分析频率
upper_order = 6;% 上通滤波器阶数
upper_f = 100; % 上截止频率
lower_order = 3;% 下通滤波器阶数
lower_f = 2;% 下截止频率
[csi_ratio_linear_filtering] = bandpass_butter(csi_ratio_linear,half_sample_rate,upper_f,upper_order,lower_order,lower_f);
%% 计算滤波后振幅相位
[csi_amplitude_butt,csi_phase_butt_unwrapped] = caculate_csibandpass_message(csi_ratio_linear_filtering);
%% 多普勒谱
WaveLength = 299792458 / 5.825e9;%波长
pca_num = 20;% PCA分析的主成分数量
window_size = 128;% 窗口大小

%频率轴对应的正0:half_sample_rate-1负频率序列
freq_bins_unwrap = [0:half_sample_rate-1 -half_sample_rate:-1]'/sample_rate;

% 设置频率轴为带通滤波后的频率轴-100到100HZ
freq_lpf_sele = freq_bins_unwrap <= upper_f / sample_rate & freq_bins_unwrap >= -upper_f / sample_rate;
%正频率数量
freq_lpf_positive_max = sum(freq_lpf_sele(2:length(freq_lpf_sele)/2));
%负频率
freq_lpf_negative_min = sum(freq_lpf_sele(length(freq_lpf_sele)/2:end));

% 创建频率和速度轴
freq_axis_bin = [0: freq_lpf_positive_max -1*freq_lpf_negative_min:-1];%频率轴对应的频率序列
velocity_axis_bin = freq_axis_bin * WaveLength/2; %  速度轴对应的速度序列

% 初始化多普勒谱矩阵
test = floor(size(csi_ratio_linear_filtering, 1));
doppler_spectrum = zeros(Rx, pca_num,1+freq_lpf_positive_max + freq_lpf_negative_min,...
    floor(size(csi_ratio_linear_filtering, 1)));
%% PCA降维
pca_coef = pca(csi_ratio_linear_filtering);
signal_pca= csi_ratio_linear_filtering * pca_coef(:,1:pca_num);

%% STFT时频分析
% 循环遍历主成分（pca_index 可能是主成分的索引）
for pca_index=1:pca_num
    % 创建一个时间实例数组，长度为 signal_pca 中的数据列数
    time_instance = 1:length(signal_pca(:,pca_index));
    % 如果 window_size 是偶数，将其增加 1，确保它是奇数
    if(~mod(window_size,2))
        window_size = window_size + 1;
    end
    % 使用 tfrsp 函数计算信号的时频表示，进行时频分析，采用高斯窗
    freq_time_prof_allfreq = tfrsp(signal_pca(:,pca_index), time_instance,...
        sample_rate, tftb_window(window_size, 'gauss'));
    % 选择带通滤波完频率内 2——100HZ
    freq_time_prof = freq_time_prof_allfreq(freq_lpf_sele, :);
    % 谱归一化
    freq_time_prof = abs(freq_time_prof) ./ repmat(sum(abs(freq_time_prof),1),...
        size(freq_time_prof,1), 1);
    % 存储多普勒谱
    if(size(freq_time_prof,3) >= size(doppler_spectrum,4))
        doppler_spectrum(1,pca_index,:,:) = freq_time_prof(:,1:size(doppler_spectrum,4));
    else
        doppler_spectrum(1,pca_index,:,:) = [freq_time_prof zeros(size(doppler_spectrum,3),...
            size(doppler_spectrum,4) - size(freq_time_prof,2))];
    end
end
%% 处理时间计算

% 创建一个与多普勒谱相同大小的全零矩阵 after_doppler_spectrum 用于存储处理后的图像
after_doppler_spectrum = zeros(size(doppler_spectrum));
% 循环，一次迭代（i=1），处理多普勒谱的主成分
for i=1:pca_num % pca_num
    % 去除服从高斯分布的带限高斯噪声
    % 创建一个高斯滤波器，5 表示滤波器大小，0.8 是滤波器的标准差
    mdlSpecLpf = fspecial('gaussian', 5, 0.8);
    % 使用滤波器对多普勒谱进行平滑处理
    outSpec = imfilter(doppler_spectrum(1,i,:,:),mdlSpecLpf);
    % 将处理后的数据添加到 after_doppler_spectrum 中
    after_doppler_spectrum (1,i,:,:) = after_doppler_spectrum (1,i,:,:) + outSpec;
end
% 将多余的维度进行压缩，得到 doppler_spectrum_squeezed 数据
doppler_spectrum_squeezed = squeeze(after_doppler_spectrum (1,1,:,:));
% 假设你的数组为 data_array
% 例如，创建一个5x3的示例数组
data_array = doppler_spectrum_squeezed;

% 找到数组的最中间的行数
middle_row = ceil(size(data_array, 1) / 2);

% 对数组进行重新排列
rearranged_array = [data_array(middle_row-1:end, :); data_array(1:middle_row, :)];


%doppler_spectrum_squeezed=squeeze(doppler_spectrum(1, 4, :, :));%
% 画多谱勒频谱图
% doppler_spectrum 的维度应为 (Rx, pca_num, 频率数, 时间点)
% 创建多普勒频率坐标
doppler_frequencies = sort(freq_axis_bin);
% 绘制多普勒谱的时频图
% figure(1);
% imagesc(time1,doppler_frequencies,abs(rearranged_array));
% shading interp; % 平滑颜色
% colormap('jet'); % 选择颜色映射
% colorbar; % 添加颜色刻度条
% xlabel('Time(s)', 'FontWeight', 'bold');
% ylabel('Frequency(Hz)', 'FontWeight', 'bold');
% axis xy;
%% 人体移动速度提取
%人体移动速度对应频率提取
freqbin_len = size(doppler_spectrum_squeezed, 1);% 获取频率维度的长度
timebin_len = size(doppler_spectrum_squeezed, 2);% 获取时间维度的长度
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

% 计算人体移动速度
fc = 300;
fs = 1000;
[b,a] = butter(10, fc/(fs/2));
% 使用滤波器对轮廓频率进行滤波
fc_filter = filter(b,a,fc_contour);
% 使用高斯平滑对滤波后的数据进行平滑处理，时间窗口大小为 250
fc_contour_freq_smooth = smoothdata(fc_filter,'gaussian',250); % 250大概是500ms的一半，也就是近似步态周期的一半
% 速度计算
movement_speed = fc_contour_freq_smooth*WaveLength/2;
%{
    %% 躯干速度提取
%躯干移动速度对应频率提取
freqbin_len = size(doppler_spectrum_squeezed, 1);% 获取频率维度的长度
timebin_len = size(doppler_spectrum_squeezed, 2);% 获取时间维度的长度
% 存储移动轮廓频率的集合
torso_freq = [];
%找到满足功率大于总功率50%的频点
 for i=1:timebin_len
     % 满足要求的某一个时间点的freq points集合
     contour_f_points = [];
        power_sum = sum(doppler_spectrum_squeezed(:, i));% 计算总功率
        for k=1:freqbin_len
            cumulated_sum = doppler_spectrum_squeezed(k, i);
            if cumulated_sum / power_sum > 0.2% 当前频率点的功率
                % 将当前频率点添加到频率轮廓集合中
                torso_freq = [contour_f_points, freq_axis_bin(k)];
            end
        end
        if size(contour_f_points) == 0
            % 将0添加到人体移动频率轮廓集合中
            torso_freq = [torso_freq, 0];
        else
            % % 提取符合条件的频率最大值,轮廓
            torso_freq = [torso_freq, max(contour_f_points)];
        end
 end
 % 计算躯体移动速度
    % 使用相同的滤波器对躯干频率进行滤波
    torso_filter = filter(b,a,torso_freq);
    % 使用高斯平滑对滤波后的数据进行平滑处理，时间窗口大小为 250   
    torso_freq_smooth = smoothdata(torso_filter,'gaussian',250); % 250大概是500ms的一半，也就是近似步态周期的一半
    % 将躯干频率转换为躯干速度
    torso_speed = torso_freq_smooth*WaveLength/2;
%% 腿部速度提取
%腿部移动速度对应频率提取
freqbin_len = size(doppler_spectrum_squeezed, 1);% 获取频率维度的长度
timebin_len = size(doppler_spectrum_squeezed, 2);% 获取时间维度的长度
% 存储移动轮廓频率的集合
leg_freq = []; 
%找到满足功率大于总功率95%的频点
 for i=1:timebin_len
     contour_f_points = [];
     % 满足要求的某一个时间点的freq points集合
        power_sum = sum(doppler_spectrum_squeezed(:, i));% 计算总功率
        for k=1:freqbin_len
            cumulated_sum = doppler_spectrum_squeezed(k, i);
            if cumulated_sum / power_sum > 0.95% 当前频率点的功率
                % 将当前频率点添加到频率轮廓集合中
                leg_freq = [contour_f_points, freq_axis_bin(k)];
            end
        end
         
 end
 % 计算腿部移动速度
leg_filter = filter(b,a,leg_freq);
leg_freq_smooth = smoothdata(leg_filter,'gaussian',250); % 250大概是500ms的一半，也就是近似步态周期的一半
leg_speed = leg_freq_smooth*WaveLength/2;
%x = squeeze(1:size(doppler_spectrum, 4));
%}
% figure(2);
% plot(time1,movement_speed,'LineWidth', 1.5,'Color', 'b')
% xlabel('Time(s)', 'FontWeight', 'bold');
% ylabel('Velocity(m/s)', 'FontWeight', 'bold');
% grid on ;
%{
%saveas(gcf,  ['v_' num2str(i) '.png']);
%{
figure;
time = 1:size(doppler_spectrum_squeezed, 2);
doppler_frequencies = freq_axis_bin;

% 绘制多普勒谱的时频图
imagesc(time1,doppler_frequencies,abs(rearranged_array));
shading interp; % 平滑颜色
colormap('jet'); % 选择颜色映射
colorbar; % 添加颜色刻度条
xlabel('Time(s)', 'FontWeight', 'bold');
ylabel('Frequency(Hz)', 'FontWeight', 'bold');
axis xy;

% 获取速度曲线的最大值和最小值
max_speed = max(movement_speed);
min_speed = min(movement_speed);

% 将速度曲线数值缩放到更大的范围以便显示
scaled_speed = (movement_speed - min_speed) / (max_speed - min_speed) * (max(doppler_frequencies) - min(doppler_frequencies)) + min(doppler_frequencies);

% 使用 hold on 保持当前坐标轴状态
hold on;

% 绘制速度曲线
plot(time1, scaled_speed, 'cyan', 'LineWidth', 2); % 使用白色线条

% 如果需要保存图像
% saveas(gcf, 'doppler_spectrum_with_scaled_velocity.png');


%{
figure;
plot(torso_speed)

figure;
plot(leg_speed);

%}
%}
%}
%% autocorrelation
% extract period where speed is no less than 80% maxmium speed
% 定义速度为0的阈值
velocity_threshold = 0.1; % 根据实际情况调整阈值

% 寻找速度开始从0增加的点
start_increasing_velocity_index = find(diff(movement_speed > velocity_threshold) > 0, 1, 'first');

% 判断在窗口内速度是否一直为0
if ~isempty(start_increasing_velocity_index)
    fprintf('速度开始从0增加的点是：%d\n', start_increasing_velocity_index);
else
    fprintf('未找到速度开始从0增加的点\n');
end
%% 寻找速度持续小于0.1的最后一个点
threshold_value = 0;

% 找到速度持续小于0.1的点
low_speed_points = find(movement_speed < threshold_value);

if ~isempty(low_speed_points)
    % 找到最后一个速度持续小于0.1的点
    end_decreasing_velocity_index = low_speed_points(end);
    fprintf('速度持续小于0.1的最后一个点是：%d\n', end_decreasing_velocity_index);
else
    fprintf('未找到速度持续小于0.1的点\n');
end

movement_speed1=movement_speed(:,start_increasing_velocity_index:end_decreasing_velocity_index);
max_speed_80 = max(movement_speed1) * 0.8;
max_speed_75 = max(movement_speed1) * 0.75;
left_index = 0;
right_index = 0;
for i=1:length(movement_speed1)
    if movement_speed1(1, i) > max_speed_80
        left_index = i;
        break;
    end
end
for i=length(movement_speed1):-1:1
    if movement_speed1(1, i) > max_speed_75
        right_index = i;
        break;
    end
end
   detect_period = movement_speed(left_index:right_index);
%     period = autocorrelation(detect_period);
    remove_average = detect_period - mean(detect_period);
    [period, locs, diffs] = autocorrelation(remove_average);
    
    average_speed = mean(movement_speed);
    steady_average_speed = mean(detect_period);
%% extract walking period
% 存储四个周期的均值
fisrt_magnitude_array = [];
second_magnitude_array = [];
third_magnitude_array = [];
fourth_magnitude_array = [];
period_start = round(locs(1) - period/2);
for i=1:length(locs)
    %         period_start = round(locs(i) - period/2);
    if period_start < 1
        period_start = 1;
    end
    period_end = locs(i);
    slice_len = period_end - period_start;
    one_period = doppler_spectrum_squeezed(:, period_start:period_end);
    fisrt_stage = one_period(:, 1:round(slice_len/4)); % 半周期的1/4
    second_stage = one_period(:, round(slice_len/4):round(slice_len/2));
    third_stage = one_period(:, round(slice_len/2):round(3*slice_len/4));
    fourth_stage = one_period(:, round(3*slice_len/4):end);

    fisrt_magnitude_array = [fisrt_magnitude_array; extract_spec_signature(fisrt_stage)];
    second_magnitude_array = [second_magnitude_array; extract_spec_signature(second_stage)];
    third_magnitude_array = [third_magnitude_array; extract_spec_signature(third_stage)];
    fourth_magnitude_array = [fourth_magnitude_array; extract_spec_signature(fourth_stage)];

    % 更新周期
    % 到达边界，停止检测
    if i < length(locs) && locs(i+1) < length(doppler_spectrum_squeezed)
        period_between = movement_speed(locs(i):locs(i+1));
        valley_loc = findvalleys([1:length(period_between)], period_between, 0, min(period_between)-1, 1, 1, 1);
        if size(valley_loc, 1) > 1
            period_start = round(locs(i) + period / 2);
        else
            period_start = locs(i) + valley_loc(1, 2);
        end
    else
        break;
    end
end
gait_cycle = period / 1000; % 步态周期
estimated_steady_footstep = steady_average_speed *2* period / 1000; % 模拟步长
estimated_steady_footstep1 = movement_speed * period / 1000; % 模拟步长
length(estimated_steady_footstep1)
%{
% 绘制图形
figure;
plot(time1,abs(estimated_steady_footstep1),'LineWidth', 1.5,'Color', 'b')
title('模拟步长随时间的变化');
xlabel('时间');
ylabel('模拟步长');
grid on;
%}

%%
% stride Length 步态周期
gait_cycle1 = 2 * period / 1000;
% 计算模拟步长
estimated_steady_footstep1 = fc_contour_freq_smooth * WaveLength / 2 * gait_cycle1;
estimated_steady_footstep1 = zeros(size(movement_speed));

% 划分窗口，计算步长
window_size = round(gait_cycle1 * 1000); % 将步态周期转换为毫秒
num_windows = floor(length(movement_speed) / window_size);

for i = 1:num_windows
    window_start = (i - 1) * window_size + 1;
    window_end = i * window_size;

    % 在窗口内计算平均速度
    window_average_speed = mean(movement_speed(window_start:window_end));

    % 计算模拟步长并赋值到相应位置
    estimated_steady_footstep1(window_start:window_end) = window_average_speed * gait_cycle1;
end

% step length
gait_cycle2 = period / 1000;
% 计算模拟步长
estimated_steady_footstep2 = fc_contour_freq_smooth * WaveLength / 2 * gait_cycle2;
estimated_steady_footstep2 = zeros(size(movement_speed));

% 划分窗口，计算步长
window_size = round(gait_cycle2 * 1000); % 将步态周期转换为毫秒
num_windows = floor(length(movement_speed) / window_size);

for i = 1:num_windows
    window_start = (i - 1) * window_size + 1;
    window_end = i * window_size;

    % 在窗口内计算平均速度
    window_average_speed = mean(movement_speed(window_start:window_end));

    % 计算模拟步长并赋值到相应位置
    estimated_steady_footstep2(window_start:window_end) = window_average_speed * gait_cycle2;
end

end

