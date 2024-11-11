function [doppler_spectrum,doppler_frequencies,doppler_rearranged_array,doppler_spectrum_squeezed,freq_axis_bin] = doppler(csi_ratio_linear_filtering,Rx)
%DOPPLER 此处显示有关此函数的摘要
%   此处显示详细说明
WaveLength = 299792458 / 5.825e9;%波长
pca_num = 20;% PCA分析的主成分数量
window_size = 128;% 窗口大小
sample_rate = 1000;% 采样率
half_sample_rate = 0.5*sample_rate;% 分析频率
upper_order = 6;% 上通滤波器阶数
upper_f = 100; % 上截止频率
lower_order = 3;% 下通滤波器阶数
lower_f = 2;% 下截止频率
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
data_array = doppler_spectrum_squeezed;
% 找到数组的最中间的行数
middle_row = ceil(size(data_array, 1) / 2);

% 对数组进行重新排列
doppler_rearranged_array = [data_array(middle_row-1:end, :); data_array(1:middle_row, :)];
doppler_frequencies = sort(freq_axis_bin);
end