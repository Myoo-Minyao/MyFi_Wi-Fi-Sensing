function [ox_vec, oy_vec] = removeNaN(x_vec, y_vec) 
% 移除 NaN 值的函数

    % 找到 NaN 值的索引，并将对应的数据从原始数据中删除
    nanIdx = isnan(y_vec);
    ox_vec = x_vec;
    oy_vec = y_vec;
    ox_vec(nanIdx) = [];
    oy_vec(nanIdx) = [];
end

