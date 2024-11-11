clear all;
addpath(genpath(pwd));%把子文件夹加入到路径中
%% 加载数据文件
[cfr_array, ant1, ant2, ant3, Subcarriers, Tx, Rx, time_duration,csi_trace] = csi_loading('sampledata.dat');

%% 获取时间序列
[time1,cfr_array] = time(time_duration,cfr_array);

%% 子载波选择 
power = abs(cfr_array).^2 ;
sum_power = sum(power, 1);
[ ~ ,subcarrier_idx] = maxk(sum_power, 1);%找到子载波中功率最大的

%% 数据切片
start_time=1;
end_time=281;
% 计算CSI信息
% (cfr_array,time1,start_time,end_time,j,i,k) j是切片开关,i是子载波序列，k是1振幅，2相位，3振幅相位
[csi_amplitude, csi_phase_unwrapped,csi_phase] = caculate_rawcsi_message(cfr_array,time1,start_time,end_time,2,subcarrier_idx,3);

%% csi ratio 去除硬件噪声
[csi_ratio,level] = F_csi_ratio(ant1,ant2,ant3);
% 计算CSI_ratio值
[csi_ratio_amplitude,csi_ratio_phase_unwrapped,csi_ratio_phase] = caculate_csiratio_message(csi_ratio,time1,start_time,end_time,2,1,3);

%% 线性相位
m=[-28,-26,-24,-22,-20,-18,-16,-14,-12,-10,-8,-6,-4,-2,-1,1,3,5,7,9,11,13,15,17,19,21,23,25,27,28];%30子载波
for t=1:30
    afterphase(:,t)=csi_ratio_phase_unwrapped(:,t)-(csi_ratio_phase_unwrapped(:,30)-csi_ratio_phase_unwrapped(:,1))/56*m(t)-1/30*sum(csi_ratio_phase_unwrapped,2);%linear transformation
end

figure;
plot(time1, afterphase(:,1),'LineWidth', 2,'Color',[0.173 0.102 0.75])
ylabel('Phase', 'FontWeight', 'bold','FontSize', 50);
xlabel('Time(s)', 'FontWeight', 'bold','FontSize', 50);
% 设置图形窗口位置和大小
x = 100; % 图形窗口左上角 x 坐标
y = 100; % 图形窗口左上角 y 坐标
width = 800; % 图形窗口宽度
height = 600; % 图形窗口高度
set(gca,'fontsize',25);
set(gcf, 'Position', [x y width height]); % 设置图形窗口位置和大小
grid on;

%% 天线选择图
figure;
categories = {'Antenna 1', 'Antenna 2', 'Antenna 3'};
bar(categories, level,'BarWidth',0.5);
ylabel('Average / Standard Deviation', 'FontWeight', 'bold');
text(1:length(categories), level, num2str(level'), 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');
set(gca, 'FontSize', 14);  % 增大坐标轴标签的字体大小
% 增大数据标签的字体大小
text_objects = findall(gcf, 'Type', 'text');
set(text_objects, 'FontSize', 12); % 增大数据标签的字体大小

%% 线性插值
[csi_ratio_linear] = linear_data(ant1,csi_ratio);
% 此时的csi ratio是去除了硬件噪声和空值的数据了
%% 设计带通滤波器(Butterworth 滤波器) 去除带外噪声
sample_rate = 1000;% 采样率
half_sample_rate = 0.5*sample_rate;% 分析频率
upper_order = 6;% 上通滤波器阶数
upper_f = 100; % 上截止频率
lower_order = 3;% 下通滤波器阶数
lower_f = 2;% 下截止频率
[csi_ratio_linear_filtering] = bandpass_butter(csi_ratio_linear,half_sample_rate,upper_f,upper_order,lower_order,lower_f);

%% 计算滤波后振幅相位
[csi_amplitude_butt,csi_phase_butt_unwrapped] = caculate_csibandpass_message(csi_ratio_linear_filtering,time1,start_time,end_time,2,1,3);

figure;
plot(time1(start_time:end_time), csi_amplitude_butt(start_time:end_time,1),'LineWidth', 2,'Color',[0.173 0.102 0.75])
ylabel('CSI', 'FontWeight', 'bold','FontSize', 50);
xlabel('Time(s)', 'FontWeight', 'bold','FontSize', 50);
% 设置图形窗口位置和大小
x = 100; % 图形窗口左上角 x 坐标
y = 100; % 图形窗口左上角 y 坐标
width = 800; % 图形窗口宽度
height = 600; % 图形窗口高度
set(gca,'fontsize',25);
set(gcf, 'Position', [x y width height]); % 设置图形窗口位置和大小
grid on;

% min_index = calculate_frontindex(csi_phase_butt_unwrapped);
% max_index = calculate_backindex(csi_phase_butt_unwrapped);

% % 绘制从最小值到最大值之间的相位数据
% figure;
% plot(time1(min_index:max_index), csi_phase_butt_unwrapped(min_index:max_index));
% xlabel('时间');
% ylabel('相位');
% title('从最小值到最大值之间的相位数据');
WaveLength = 299792458 / 5.825e9;%波长
d =( csi_phase_butt_unwrapped / WaveLength )*0.0525 ;
figure;
plot(d)

% 计算CSI数据的能量
csi_energy = abs(csi_ratio_linear_filtering).^2;

% 计算时间轴上的能量平均值
average_energy = mean(csi_energy, 2); % 对每一行（时间点）求平均值

% 指定移动平均窗口大小
window_size = 10;

% 对能量数据进行移动平均处理
smoothed_energy = movmean(average_energy, window_size);

% 绘制平滑后的时间-能量图
figure;
plot(time1, abs(log10(smoothed_energy))); % 对平滑后的能量取对数，以便更好地展示
xlabel('Time(s)', 'FontWeight', 'bold');
ylabel('Smoothed Energy (dB)', 'FontWeight', 'bold');
title('Smoothed Time-Energy Plot');

figure;
plot(time1, csi_amplitude); % 对平滑后的能量取对数，以便更好地展示
xlabel('Time(s)', 'FontWeight', 'bold');
ylabel('Smoothed Energy (dB)', 'FontWeight', 'bold');

%%

sequence =30:-1:1;
figure;
%imagesc(time1,sequence,power');
imagesc(time1,sequence,csi_amplitude_butt');
shading interp; % 平滑颜色
colormap('jet'); % 选择颜色映射
xlabel('Time(s)', 'FontWeight', 'bold');
ylabel('Subcarriers', 'FontWeight', 'bold');
x = 100; % 图形窗口左上角 x 坐标
y = 100; % 图形窗口左上角 y 坐标
width = 800; % 图形窗口宽度
height = 600; % 图形窗口高度
set(gca,'fontsize',25);
set(gcf, 'Position', [x y width height]); % 设置图形窗口位置和大小
axis xy;

%% 多普勒谱
WaveLength = 299792458 / 5.825e9;%波长
pca_num = 20;% PCA分析的主成分数量
window_size = 128;% 窗口大小

[doppler_spectrum,doppler_frequencies,doppler_rearranged_array,doppler_spectrum_squeezed,freq_axis_bin] = doppler(csi_ratio_linear_filtering,Rx);

draw_doppler(time1,doppler_frequencies,doppler_rearranged_array);

%% 提取出对象运动时间，减去环境影响部分。
% 计算频域中的功率分布
power_distribution = abs(fft(csi_ratio_linear_filtering)).^2;%将每个频率分量的振幅平方

% 平滑功率分布
smoothed_power_distribution = smoothdata(power_distribution, 'gaussian', 2500); % 使用高斯平滑，窗口大小为20
soot_1 =abs(log10(smoothed_power_distribution(:,1)));

% 计算 soot_1 数组的一阶导数
soot_1_derivative = diff(soot_1);

% 绘制yi阶导数曲线
figure;
plot(time1(1:end-1), soot_1_derivative,'LineWidth', 2,'Color',[0.173 0.102 0.75]);
xlabel('time');
ylabel('Power Distribution Derivative');
% 设置图形窗口位置和大小
x = 100; % 图形窗口左上角 x 坐标
y = 100; % 图形窗口左上角 y 坐标
width = 800; % 图形窗口宽度
height = 600; % 图形窗口高度
set(gca,'fontsize',25);
set(gcf, 'Position', [x y width height]); % 设置图形窗口位置和大小
grid on;

% 找到曲线的负极值点
[~, min_indices] = findpeaks(-soot_1);

% 找到大于4的索引
greater_than_4_indices = find(soot_1 > 3.5);

% 找到连续段的起始和结束点
start_end_points = [greater_than_4_indices(1)];
for i = 2:length(greater_than_4_indices)
    if greater_than_4_indices(i) ~= greater_than_4_indices(i-1) + 1
        start_end_points = [start_end_points, greater_than_4_indices(i-1), greater_than_4_indices(i)];
    end
end
start_end_points = [start_end_points, greater_than_4_indices(end)];

% 将起始和结束点显示在图上
figure;
plot(time1, soot_1,'LineWidth', 5,'Color',[0.173 0.102 0.75]);
% hold on;
% scatter(time1(min_indices), soot_1(min_indices), 'ro');
% scatter(time1(start_end_points), soot_1(start_end_points), 'g*');
% hold off;
xlabel('Time(s)');
ylabel('Dynamic Power Distribution');
% legend('Smoothed Power distribution', 'Minima', 'Points > o');
% 设置图形窗口位置和大小
x = 100; % 图形窗口左上角 x 坐标
y = 100; % 图形窗口左上角 y 坐标
width = 800; % 图形窗口宽度
height = 600; % 图形窗口高度
set(gca,'fontsize',25);
set(gcf, 'Position', [x y width height]); % 设置图形窗口位置和大小
grid on;

% 获取开始指针和结束指针
start_index = start_end_points(3);
end_index = start_end_points(end);

% 显示结果
disp(['Start Index: ', num2str(start_index)]);
disp(['End Index: ', num2str(end_index)]);

smoothed_csi_phase_butt_unwrapped_distribution = smoothdata(csi_phase_butt_unwrapped(:,1), 'gaussian', 450); % 使用高斯平滑，窗口大小为20
% 计算 soot_1 数组的一阶导数
csi_phase_butt_unwrapped_derivative = diff(smoothed_csi_phase_butt_unwrapped_distribution);

% 绘制yi阶导数曲线
figure;
plot(time1(start_index:end_index), csi_phase_butt_unwrapped_derivative(start_index:end_index));



%{
holdmath = 0.1;
[contour_t_points,contour_f_points] = gen_speed(time1,doppler_frequencies,doppler_rearranged_array,holdmath);
% figure;
% plot(contour_f_points);
mindata=min(contour_t_points(:,2));
mindata_indices = find(contour_t_points(:, 2) == mindata)
std_data=std(contour_t_points(:,2))
range_data=range(contour_t_points(:,2))
mean_Data = mean(contour_t_points(15:20,2))
fc_good = []; % 初始化 fc_good
for i = 1:size(contour_t_points, 1) % 遍历 contour_t_points 的每一行
    if contour_t_points(i, 2) - mindata < 0.000001 % 判断每一行的第二列是否小于阈值
        fc_good = [fc_good; contour_t_points(i, :)]; % 将符合条件的行添加到 fc_good 中
    end
end

figure;
plot(fc_good(:,1), fc_good(:,2)); % 绘制符合条件的数据点

figure;
plot(contour_t_points(:, 1), contour_t_points(:, 2)-mindata);
fc = 300;
fs = 1000;
[b,a] = butter(10, fc/(fs/2));
% 使用滤波器对轮廓频率进行滤波
fc_filter = filter(b,a,contour_f_points);
%平滑
contour_f_points_smooth = smoothdata(fc_filter,'gaussian',250);
% figure;
% plot(contour_f_points_smooth);
% 频率谱转换速度谱

wave_length = 299792458 / 5.825e9;
speed = contour_f_points_smooth * wave_length / 2;

figure;
plot(speed);
plot(time1(start_index:end_index), speed(start_index:end_index));
%}
%% 人体移动速度提取
%人体移动速度对应频率提取
% 计算人体移动速度
fc_contour = movement_speed(doppler_spectrum_squeezed,time1,doppler_frequencies,freq_axis_bin);


fc = 300;
fs = 1000;
[b,a] = butter(10, fc/(fs/2));
% 使用滤波器对轮廓频率进行滤波
fc_filter = filter(b,a,fc_contour);
% 使用高斯平滑对滤波后的数据进行平滑处理，时间窗口大小为 250
fc_contour_freq_smooth = smoothdata(fc_filter,'gaussian',250); % 250大概是500ms的一半，也就是近似步态周期的一半
% 速度计算

movement_speed1 = fc_contour_freq_smooth*WaveLength/2;
smoothed_movement_speed1 = smoothdata(movement_speed1, 'gaussian', 150); % 使用高斯平滑，窗口大小为20
figure;
plot(time1(start_index:end_index),smoothed_movement_speed1(start_index:end_index),'LineWidth', 3,'Color',[0.173 0.102 0.75]);
xlabel('Time(s)');
ylabel('Velocity(m/s)');
% 设置图形窗口位置和大小
x = 100; % 图形窗口左上角 x 坐标
y = 100; % 图形窗口左上角 y 坐标
width = 800; % 图形窗口宽度
height = 600; % 图形窗口高度
set(gca,'fontsize',25);
set(gcf, 'Position', [x y width height]); % 设置图形窗口位置和大小
grid on;
figure;
imagesc(time1,movement_speed1,csi_amplitude);
%导数

soot_2 =abs(log10(smoothed_movement_speed1));
% 计算 soot_1 数组的一阶导数
soot_2_derivative = diff(soot_2);
figure;
plot(time1(start_index:end_index),soot_2(start_index:end_index),'LineWidth', 3,'Color',[0.173 0.102 0.75]);
xlabel('Time(s)');
ylabel('Velocity Derivative');
% 设置图形窗口位置和大小
x = 100; % 图形窗口左上角 x 坐标
y = 100; % 图形窗口左上角 y 坐标
width = 800; % 图形窗口宽度
height = 600; % 图形窗口高度
set(gca,'fontsize',25);
set(gcf, 'Position', [x y width height]); % 设置图形窗口位置和大小
grid on;



%{
%% 自相关
% 计算最大移动速度的80%和75%
max_speed_80 = max(movement_speed) * 0.8;
max_speed_75 = max(movement_speed) * 0.75;

% 初始化检测周期的左右索引
left_index = 0;
right_index = 0;

% 寻找检测周期的左索引
for i = 1:length(movement_speed)
    if movement_speed(1, i) > max_speed_80
        left_index = i;
        break;
    end
end

% 寻找检测周期的右索引
for i = length(movement_speed):-1:1
    if movement_speed(1, i) > max_speed_75
        right_index = i;
        break;
    end
end

% 根据左右索引提取检测周期
detect_period = movement_speed(left_index:right_index);

% 从检测周期中去除平均值
remove_average = detect_period - mean(detect_period);

% 计算自相关并获取周期、位置和差异
[period, locs, diffs] = autocorrelation(remove_average);

% 计算整个移动的平均速度
average_speed = mean(movement_speed);

% 计算检测周期内的稳定平均速度
steady_average_speed = mean(detect_period);

figure;
plot(steady_average_speed);
%% 自相关提取步态周期

%test = floor(size(cfr_array, 1));
%doppler_spectrum = zeros(Rx, pca_num,1+freq_lpf_positive_max + freq_lpf_negative_min,...
%    floor(size(cfr_array, 1)));

%% 运行时间计算 tic为秒表计时器，toc是计算现在到tic所开始的计时器的时间间隔。
%{
loading_data_time=0;%初始化数据加载时间
signal_processing_time=0;%初始化信号处理时间
loading_data_time=loading_data_time+toc;%计算数据加载时间
tic;
%}



%% 处理时间计算
%{
signal_processing_time=signal_processing_time+toc;
loading_data_time=loading_data_time/Rx;
signal_processing_time=signal_processing_time/Rx;
%}
%{
% 找到频率轴上的最大值的索引
    [~,idx] = max(freq_axis_bin);
% 计算圆形的长度，这里是基于频率轴上的索引位置
    circle_length = length(freq_axis_bin) - idx;
% doppler_spectrum_noflip = circshift(doppler_spectrum, [0 0 circle_length 0]);
% 将多普勒谱进行循环移位，以便处理数据
    doppler_spectrum = circshift(doppler_spectrum, [0 0 circle_length 0]);
% 对频率轴进行循环移位
    freq_axis_bin=circshift(freq_axis_bin, [0 circle_length]);
% 对速度轴进行循环移位
    velocity_axis_bin=circshift(velocity_axis_bin, [0 circle_length]);
%}

% % 绘制多普勒谱的时频图

% % 保存图像
%saveas(gcf,  ['d_' num2str(i) '.png']);

% % 创建频率轴和速度轴
% freq_axis_bin = [0:freq_lpf_positive_max -1*freq_lpf_negative_min:-1];
% [S1, F1, T1] = spectrogram (user1, 256, 128, 512, fs);
% [S2, F2, T2] = spectrogram (user2, 256, 128, 512, fs);
% 
% % 创建一个图形窗口，使用 subplot 函数将其分为两个子图
% figure;
% subplot (2, 1, 1); % 第一个子图，显示用户 1 的频率-时间图
% 
% % 使用 imagesc 函数绘制用户 1 的频率-时间图，注意矩阵 S1 需要转置
% % 您可以根据您的需要修改颜色映射（这里是 jet）、坐标轴范围（这里是 [-60, 60] Hz）和坐标轴标签
% imagesc (T1, F1, abs (S1.'));
% colormap (jet);
% axis xy;
% ylim ([-60, 60]);
% xlabel ('Time (s)');
% ylabel ('Frequency (Hz)');
% 
% 
% % 使用 colorbar 函数添加一个色条，表示频率强度的范围
% % 您可以根据您的需要修改色条的位置和标签
% colorbar ('eastoutside');
% caxis ([0, 0.18]);
% clabel = 'Frequency Intensity';
% h = get (gca, 'ylabel');
% p = get (h, 'position');
% set (gca, 'ylabel', text (p(1), p(2), clabel, 'rotation', -90, 'verticalalignment', 'bottom'));
% 
% % 使用 text 函数添加一些注释，指示特定的事件
% % 您可以根据您的需要修改注释的位置、内容和样式
% text (0.5, -50, 'Limb movement', 'color', 'w', 'fontsize', 12, 'horizontalalignment', 'center');
% text (2.5, -50, 'Action interval', 'color', 'w', 'fontsize', 12, 'horizontalalignment', 'center');
% text (4.5, -50, 'Motion change', 'color', 'w', 'fontsize', 12, 'horizontalalignment', 'center');
% 
% % 重复上述步骤，绘制用户 2 的频率-时间图
% subplot (2, 1, 2); % 第二个子图，显示用户 2 的频率-时间图
% imagesc (T2, F2, abs (S2.'));
% colormap (jet);
% axis xy;
% ylim ([-60, 60]);
% xlabel ('Time (s)');
% ylabel ('Frequency (Hz)');
% title ('User 2');
% colorbar ('eastoutside');
% caxis ([0, 0.18]);
% clabel = 'Frequency Intensity';
% h = get (gca, 'ylabel');
% p = get (h, 'position');
% set (gca, 'ylabel', text (p(1), p(2), clabel, 'rotation', -90, 'verticalalignment', 'bottom'));
% text (0.5, -50, 'Limb movement', 'color', 'w', 'fontsize', 12, 'horizontalalignment', 'center');
% text (2.5, -50, 'Action interval', 'color', 'w', 'fontsize', 12, 'horizontalalignment', 'center');
% text (4.5, -50, 'Motion change', 'color', 'w', 'fontsize', 12, 'horizontalalignment', 'center');


% %{
%     %% 躯干速度提取
% %躯干移动速度对应频率提取
% freqbin_len = size(doppler_spectrum_squeezed, 1);% 获取频率维度的长度
% timebin_len = size(doppler_spectrum_squeezed, 2);% 获取时间维度的长度
% % 存储移动轮廓频率的集合
% torso_freq = [];
% %找到满足功率大于总功率50%的频点
%  for i=1:timebin_len
%      % 满足要求的某一个时间点的freq points集合
%      contour_f_points = [];
%         power_sum = sum(doppler_spectrum_squeezed(:, i));% 计算总功率
%         for k=1:freqbin_len
%             cumulated_sum = doppler_spectrum_squeezed(k, i);
%             if cumulated_sum / power_sum > 0.2% 当前频率点的功率
%                 % 将当前频率点添加到频率轮廓集合中
%                 torso_freq = [contour_f_points, freq_axis_bin(k)];
%             end
%         end
%         if size(contour_f_points) == 0
%             % 将0添加到人体移动频率轮廓集合中
%             torso_freq = [torso_freq, 0];
%         else
%             % % 提取符合条件的频率最大值,轮廓
%             torso_freq = [torso_freq, max(contour_f_points)];
%         end
%  end
%  % 计算躯体移动速度
%     % 使用相同的滤波器对躯干频率进行滤波
%     torso_filter = filter(b,a,torso_freq);
%     % 使用高斯平滑对滤波后的数据进行平滑处理，时间窗口大小为 250   
%     torso_freq_smooth = smoothdata(torso_filter,'gaussian',250); % 250大概是500ms的一半，也就是近似步态周期的一半
%     % 将躯干频率转换为躯干速度
%     torso_speed = torso_freq_smooth*WaveLength/2;
% %% 腿部速度提取
% %腿部移动速度对应频率提取
% freqbin_len = size(doppler_spectrum_squeezed, 1);% 获取频率维度的长度
% timebin_len = size(doppler_spectrum_squeezed, 2);% 获取时间维度的长度
% % 存储移动轮廓频率的集合
% leg_freq = []; 
% %找到满足功率大于总功率95%的频点
%  for i=1:timebin_len
%      contour_f_points = [];
%      % 满足要求的某一个时间点的freq points集合
%         power_sum = sum(doppler_spectrum_squeezed(:, i));% 计算总功率
%         for k=1:freqbin_len
%             cumulated_sum = doppler_spectrum_squeezed(k, i);
%             if cumulated_sum / power_sum > 0.95% 当前频率点的功率
%                 % 将当前频率点添加到频率轮廓集合中
%                 leg_freq = [contour_f_points, freq_axis_bin(k)];
%             end
%         end
% 
%  end
%  % 计算腿部移动速度
% leg_filter = filter(b,a,leg_freq);
% leg_freq_smooth = smoothdata(leg_filter,'gaussian',250); % 250大概是500ms的一半，也就是近似步态周期的一半
% leg_speed = leg_freq_smooth*WaveLength/2;
% %x = squeeze(1:size(doppler_spectrum, 4));
% %}
figure;
plot(time1,movement_speed,'LineWidth', 1.5,'Color', 'b')
xlabel('Time(s)', 'FontWeight', 'bold');
ylabel('Velocity(m/s)', 'FontWeight', 'bold');
grid on ;
% 
%  % 计算速度谱的自相关
% autocorr_velocity = xcorr(movement_speed, 'coeff'); % 'coeff'用于归一化自相关
% 
% % 绘制自相关函数
% figure;
% plot(autocorr_velocity);
% xlabel('时间滞后');
% ylabel('自相关');
% title('速度谱的自相关函数');
% % 假设您已经计算了速度谱的自相关函数 autocorr_velocity
% 
% % 假设您已经计算了速度谱的自相关函数 autocorr_velocity
% 
% % 取自相关函数的负数来找到谷值
% neg_autocorr_velocity = -autocorr_velocity;
% 
% % 找到自相关函数的谷值
% [valleys, valley_locations] = findpeaks(neg_autocorr_velocity);
% 
% % 找到最小谷值的位置
% [min_valley, min_index] = min(valleys);
% min_valley_location = valley_locations(min_index);
% 
% % 提取相应的速度曲线
% window_size = min(length(movement_speed), length(autocorr_velocity) - min_valley_location + 1);
% min_correlated_speed = movement_speed(1:window_size);
% 
% % 绘制最小自相关程度的速度曲线
% figure;
% plot(min_correlated_speed, 'LineWidth', 1.5);
% xlabel('时间');
% ylabel('速度');
% title('最小自相关程度的速度曲线');
% grid on;
% 
% 
% 
% 
% %{
% %saveas(gcf,  ['v_' num2str(i) '.png']);
% %{
% figure;
% time = 1:size(doppler_spectrum_squeezed, 2);
% doppler_frequencies = freq_axis_bin;
% 
% % 绘制多普勒谱的时频图
% imagesc(time1,doppler_frequencies,abs(rearranged_array));
% shading interp; % 平滑颜色
% colormap('jet'); % 选择颜色映射
% colorbar; % 添加颜色刻度条
% xlabel('Time(s)', 'FontWeight', 'bold');
% ylabel('Frequency(Hz)', 'FontWeight', 'bold');
% axis xy;
% 
% % 获取速度曲线的最大值和最小值
% max_speed = max(movement_speed);
% min_speed = min(movement_speed);
% 
% % 将速度曲线数值缩放到更大的范围以便显示
% scaled_speed = (movement_speed - min_speed) / (max_speed - min_speed) * (max(doppler_frequencies) - min(doppler_frequencies)) + min(doppler_frequencies);
% 
% % 使用 hold on 保持当前坐标轴状态
% hold on;
% 
% % 绘制速度曲线
% plot(time1, scaled_speed, 'cyan', 'LineWidth', 2); % 使用白色线条
% 
% % 如果需要保存图像
% % saveas(gcf, 'doppler_spectrum_with_scaled_velocity.png');
% 
% 
% %{
% figure;
% plot(torso_speed)
% 
% figure;
% plot(leg_speed);
% 
% %}
% %}
% %}
% 
% %% autocorrelation
% % extract period where speed is no less than 80% maxmium speed
% velocity_threshold = 0.1; % 阈值
% 
% % 寻找速度开始从0增加的点
% start_increasing_velocity_index = find(diff(movement_speed > velocity_threshold) > 0, 1, 'first');
% 
% % 判断在窗口内速度是否一直为0
% if ~isempty(start_increasing_velocity_index)
%     fprintf('速度开始从0增加的点是：%d\n', start_increasing_velocity_index);
% else
%     fprintf('未找到速度开始从0增加的点\n');
% end
% %% 寻找速度持续小于0.1的最后一个点
% threshold_value = 0;
% 
% % 找到速度持续小于0.1的点
% low_speed_points = find(movement_speed < threshold_value);
% 
% if ~isempty(low_speed_points)
%     % 找到最后一个速度持续小于0.1的点
%     end_decreasing_velocity_index = low_speed_points(end);
%     fprintf('速度持续小于0.1的最后一个点是：%d\n', end_decreasing_velocity_index);
% else
%     fprintf('未找到速度持续小于0.1的点\n');
% end
% 
% movement_speed1=movement_speed(:,start_increasing_velocity_index:end_decreasing_velocity_index);
% max_speed_80 = max(movement_speed1) * 0.8;
% max_speed_75 = max(movement_speed1) * 0.75;
% 
% left_index = 0;
% right_index = 0;
% for i=1:length(movement_speed1)
%     if movement_speed1(1, i) > max_speed_80
%         left_index = i;
%         break;
%     end
% end
% for i=length(movement_speed1):-1:1
%     if movement_speed1(1, i) > max_speed_75
%         right_index = i;
%         break;
%     end
% end
%    detect_period = movement_speed(left_index:right_index);
%      period1 = autocorrelation(detect_period);
%     remove_average = detect_period - mean(detect_period);
%     [period, locs, diffs] = autocorrelation(remove_average);
% 
% 
% %%
% % 计算加速度
% acceleration1 = diff(movement_speed) ./ diff(time1);
% time_acceleration = time1(1:end-1);
% x=exp(acceleration1);
% 
% figure;
% plot(time_acceleration, x);
% xlabel('时间');
% ylabel('加速度');
% title('加速度随时间的变化');
% grid on;
% 
% %%
% 
% % 时间范围
% start_time =8; 
% end_time =20;   
% % 找到对应时间范围内的索引
% start_index = find(time1 >= start_time, 1, 'first');
% end_index = find(time1 <= end_time, 1, 'last');
% 
% % 截取指定时间范围内的速度数据
% selected_time_range = time1(start_index:end_index);
% selected_velocity_range = movement_speed(start_index:end_index);
% 
% % 映射速度数据到正半轴
% mapped_velocity = abs(selected_velocity_range);
% 
% % 使用 findpeaks 函数提取全部极大值点的位置
% [~, maxima_locations] = findpeaks(mapped_velocity);
% 
% % 根据位置获取极大值点的速度值
% maxima_values = mapped_velocity(maxima_locations);
% 
% % 计算极大值点的速度均值
% average_speed_at_maxima = mean(maxima_values);
% 
% disp(['速度均值: ' num2str(average_speed_at_maxima)]);
% 
% figure;
% plot(selected_time_range, selected_velocity_range, 'LineWidth', 1.5, 'Color', 'b');
% hold on;
% scatter(selected_time_range(maxima_locations), selected_velocity_range(maxima_locations), 100, 'r', 'filled');
% xlabel('Time(s)', 'FontWeight', 'bold');
% ylabel('Velocity(m/s)', 'FontWeight', 'bold');
% title('Velocity Data and Maxima');
% legend('Velocity', 'Maxima');
% grid on;
% estimated_steady_footstep1 = average_speed_at_maxima * period / 1000; % 步长
% disp(['步长均值: ' num2str(estimated_steady_footstep1)]);

%{

    average_speed = mean(abs(movement_speed));
    steady_average_speed = mean(abs(detect_period));
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
path = 'Env3_13_test_b2.dat';
[estimated_steady_footstep3,estimated_steady_footstep4,time2] = Mult_data(path);

% % 绘制图形
% figure;
% 
% % 创建第一张图 (subplot 1x2，第一张图)
% subplot(2, 1, 1);
% 
% % 绘制第一条曲线（蓝色实线）
% plot(time1, abs(estimated_steady_footstep2), '-o', 'Color', 'b');
% hold on;
% 
% % 绘制第二条曲线（红色实线）
% plot(time1, abs(estimated_steady_footstep1), '-o', 'Color', 'r');
% 
% 
% ylabel('Estimated Length(m)', 'FontWeight', 'bold');
% grid on;
% 
% legend('Cross LOS Step Length', 'Cross LOS Stride Length');
% 
% % 创建第二张图 (subplot 1x2，第二张图)
% subplot(2, 1, 2);
% 
% % 绘制第三条曲线（蓝色虚线带星号）
% plot(time2, abs(estimated_steady_footstep4), '--*b');
% hold on;
% 
% % 绘制第四条曲线（红色虚线带星号）
% plot(time2, abs(estimated_steady_footstep3), '--*r');
% 
% xlabel('Time(s)', 'FontWeight', 'bold');
% ylabel('Estimated Length(m)', 'FontWeight', 'bold');
% grid on;
% 
% legend('None Cross LOS Step Length', 'None Cross LOS Stride Length');

%}
%}
%}