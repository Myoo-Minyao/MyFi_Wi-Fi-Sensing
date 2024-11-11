function [] = draw_doppler(time1,doppler_frequencies,rearranged_array)
figure();
imagesc(time1,doppler_frequencies,abs(rearranged_array));
shading interp; % 平滑颜色
colormap('jet'); % 选择颜色映射
xlabel('Time(s)', 'FontWeight', 'bold');
ylabel('Frequency(Hz)', 'FontWeight', 'bold');
x = 100; % 图形窗口左上角 x 坐标
y = 100; % 图形窗口左上角 y 坐标
width = 800; % 图形窗口宽度
height = 600; % 图形窗口高度
set(gca,'fontsize',25);
set(gcf, 'Position', [x y width height]); % 设置图形窗口位置和大小
axis xy;
end

