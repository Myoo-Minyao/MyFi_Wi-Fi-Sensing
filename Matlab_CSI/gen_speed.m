function [contour_t_points,contour_f_points] = gen_speed(time1,doppler_frequencies,doppler_rearranged_array,holdmath)
%躯干速度提取函数
timebin_len=size(time1,2);
freqbin_len=size(doppler_frequencies,2);
fc_contour = [];
contour_t_points = zeros(timebin_len, 2);
contour_f_points = zeros(freqbin_len, 2);
contour_f_points = [];
for i=1:timebin_len
    power_sum = sum(doppler_rearranged_array(:,i));
    contour_t_points(i, :) = [time1(i), power_sum];
    % end
    % size(contour_t_points);
    % figure;
    % plot(contour_t_points(:, 1), contour_t_points(:, 2));
    % xlabel('Time');
    % ylabel('Power Sum');
    for k=1:freqbin_len
        sum_power_point = sum(doppler_rearranged_array(1:k,i));
        % contour_f_points(:,k) = [sum_power_point,doppler_frequencies(k)];
        % end
        % figure;
        % plot(contour_f_points(:, 2), contour_f_points(:, 1));
        % xlabel('Time');
        % ylabel('Power Sum');
        if sum_power_point / power_sum >holdmath
            contour_f_points =[contour_f_points,doppler_frequencies(k)];
            break;
        end
    end
end

end

