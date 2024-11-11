clear;
path1 = 'Env3_13_test_c2.dat' ;
path2 = 'Env3_13_test_b2.dat' ;
[estimated_steady_footstep2,estimated_steady_footstep2,time2,doppler_frequencies2,rearranged_array2] =Mult_data(path1);
[estimated_steady_footstep,estimated_steady_footstep,time1,doppler_frequencies,rearranged_array] =Mult_data(path2);
%doppler_spectrum_squeezed=squeeze(doppler_spectrum(1, 4, :, :));%
% 画多谱勒频谱图
% doppler_spectrum 的维度应为 (Rx, pca_num, 频率数, 时间点)
% 创建多普勒频率坐标
% % 绘制多普勒谱的时频图
figure();
subplot(2,1,1)
imagesc(time1,doppler_frequencies,abs(rearranged_array));
shading interp; % 平滑颜色
colormap('jet'); % 选择颜色映射

ylabel('Frequency(Hz)', 'FontWeight', 'bold');
axis xy;
subplot(2,1,2)
imagesc(time2,doppler_frequencies2,abs(rearranged_array2));
shading interp; % 平滑颜色
colormap('jet'); % 选择颜色映射
xlabel('Time(s)', 'FontWeight', 'bold');
ylabel('Frequency(Hz)', 'FontWeight', 'bold');
axis xy;
% 为所有子图添加一个共享的颜色栏
h = colorbar; % 创建一个颜色栏对象
set(h, 'Position', [.91 .11 .0320 .8150]); % 调整颜色栏的位置和大小
set(gcf, 'Position'); % 调整图窗的位置和大小

