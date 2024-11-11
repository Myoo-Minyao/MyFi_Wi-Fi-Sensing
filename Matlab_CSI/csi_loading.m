function [cfr_array,ant1,ant2,ant3,Subcarriers,Tx,Rx,time_duration,csi_trace] = csi_loading(file)
% 读取CSI文件信息
% 输入：CSI文件名
% 输出：CSI数组，三根接收天线数据，子载波个数，发送天线，接收天线，发射天线，时间间隔。
csi_trace = read_bf_file(file) ;
time_duration = zeros(length(csi_trace),1);%数据时长（行，列） 该数组长度为数据包个数
%size(time_duration)( 数据包， 1)
% cfr_array为所有采集数据，格式为（数据包，RX*TX*子载波）
% ant1 ant2 ant3 分别对应3根天线数据 格式为（数据包，子载波）
%每个数据包都包含了一整个struct，csi_trace{1}
%{
 包含以下字段的 struct:
    timestamp_low: 1.1805e+09
       bfee_count: 526
              Nrx: 3
              Ntx: 1
           rssi_a: 39
           rssi_b: 36
           rssi_c: 32
            noise: -127
              agc: 45
             perm: [1 2 3]
             rate: 257
              csi: [1×3×30 double]
%}
scaled_first_csi = get_scaled_csi(csi_trace{1});%获取struct里的csi信息
% csi: [1×3×30 double]
Tx = size(scaled_first_csi, 1);%获取发射天线
Rx = size(scaled_first_csi, 2);%获取接收天线
Subcarriers = size(scaled_first_csi, 3);%获取子载波个数
cfr_array = zeros(length(csi_trace), Rx*Tx*Subcarriers);% 为存储 CSI 数据创建零矩阵
ant1 = zeros(length(csi_trace),Subcarriers);%（数据包长度，子载波数）单个天线数据
ant2 = zeros(length(csi_trace),Subcarriers);
ant3 = zeros(length(csi_trace),Subcarriers);
for k = 1 : length(csi_trace)%遍历每个数据包
    csi_entry = csi_trace{k};%获取第K个数据包的CSI结构体
    %csi_test = get_scaled_csi(csi_entry);%获取第K个数据包的CSI数据,30行3列
    csi_sque = squeeze(get_scaled_csi(csi_entry)).';%转置数组 30行3列 .'复数共轭转置，squeeze是删去维度为1的维数，降维
    csi = [csi_sque(:,1);csi_sque(:,2);csi_sque(:,3)].';%把天线对应加入每一列中
    ant1(k,:)=[csi_sque(:,1)];%第一根接收天线
    ant2(k,:)=[csi_sque(:,2)];%第二根接收天线
    ant3(k,:)=[csi_sque(:,3)];%第三根接收天线
    time_duration(k) = csi_entry.timestamp_low;%获取时间戳
    cfr_array(k,:) = csi;
end
end

