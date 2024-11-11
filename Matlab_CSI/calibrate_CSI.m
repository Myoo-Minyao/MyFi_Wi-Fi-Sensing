function output = calibrate_CSI(Rx_angle_diff_x, Rx_angle_diff_y, middle_line)
% 校准 CSI 数据的函数

    % 如果没有传入中线值，则计算数据的平均值作为中线值
    if isnan(middle_line)
        Rx_angle_diff_y_mean = mean(Rx_angle_diff_y);
    else
        Rx_angle_diff_y_mean = middle_line;
    end
    if Rx_angle_diff_y_mean > 0
        downer_Rx_angle_range_1 = Rx_angle_diff_y_mean;
        upper_Rx_angle_range_1 = -180 + Rx_angle_diff_y_mean;

        condition_index = Rx_angle_diff_y>downer_Rx_angle_range_1 | Rx_angle_diff_y< upper_Rx_angle_range_1;
        line_1_x = Rx_angle_diff_x(condition_index);
        line_1_y = Rx_angle_diff_y(condition_index);
        line_1_y(line_1_y < upper_Rx_angle_range_1) = line_1_y(line_1_y < upper_Rx_angle_range_1) + 180 + 180;
        changed_line_1_y = line_1_y - mean(line_1_y);

        condition_index = ~condition_index;
        line_2_x = Rx_angle_diff_x(condition_index);
        line_2_y = Rx_angle_diff_y(condition_index);
        changed_line_2_y = line_2_y - mean(line_2_y);
        
        combined_line = [line_1_x, line_2_x; changed_line_1_y, changed_line_2_y].';
        sort_combined_line = sortrows(combined_line, 1);
    else
        downer_Rx_angle_range_1 = Rx_angle_diff_y_mean;
        upper_Rx_angle_range_1 = 180 + Rx_angle_diff_y_mean;

        condition_index = Rx_angle_diff_y>downer_Rx_angle_range_1 & Rx_angle_diff_y< upper_Rx_angle_range_1;
        line_1_x = Rx_angle_diff_x(condition_index);
        line_1_y = Rx_angle_diff_y(condition_index);
        changed_line_1_y = line_1_y - mean(line_1_y);

        condition_index = ~condition_index;
        line_2_x = Rx_angle_diff_x(condition_index);
        line_2_y = Rx_angle_diff_y(condition_index);
        line_2_y(line_2_y > upper_Rx_angle_range_1) = - (180 - line_2_y(line_2_y > upper_Rx_angle_range_1))- 180;
        changed_line_2_y = line_2_y - mean(line_2_y);
        
        combined_line = [line_1_x, line_2_x; changed_line_1_y, changed_line_2_y].';
        sort_combined_line = sortrows(combined_line, 1);
    end
    output = sort_combined_line(:,2).';
     % 根据中线将数据分成两段，分别校准
    % 注意：角度的变化是周期性的，因此需要特殊处理
    % 在这里，将角度限制在 (-180, 180] 区间内
    % 以及处理周期性的角度变化
    % 这部分代码相当于对数据进行了周期性处理，使得处理后的数据更加平滑
end