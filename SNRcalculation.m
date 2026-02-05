function SNR = SNRcalculation( data, Nmax, flag0, flag1, flag2 )
% 计算频域的信噪比
% data 传入函数的频谱(需要转化为20倍logdB)
% Nmax 剔除的极大值
% flag0 是否需要计算成归一化的幅值，flag0=0为不需要，flag0=1需要
% flag1 计算模式，flag1=0为全部计算，flag1=0为选取中位数以上的数计算 flag2为选取峰值
% SNR 函数计算的信噪比
N = length(data);
if flag0 == 0
    data =data;
else flag0 == 1
    data = 20*log10(abs(data)/max(abs(data)));
end

% 是否指定信号的谱
if flag2 == 0
    datas = max(data);
else
    datas = data(flag2);
end

% 计算
if flag1 == 0
    data1 = sort(data,'descend') ; % 该波束幅值排序
    dataNmax = data1(1:Nmax);      % 该波束幅值的前5大
    sumdata = sum(data1)- sum(dataNmax);    % 除了前5个大值的总和
    Background = sumdata/( N - Nmax ) ;  % 背景（除了前5个大值的平均）
    SNR = datas - Background;
elseif flag1 == 1
    data1 = sort(data,'descend') ; % 该波束幅值排序
    data1 = data1(Nmax+1:floor(N/3) );    % 前一半中位数
    dataNmax = data1(1:Nmax);      % 该波束幅值的前5大
    sumdata = sum(data1)- sum(dataNmax);    % 除了前5个大值的总和
    Background = sumdata/( floor(N/3) - Nmax ) ;  % 背景（除了前5个大值的平均）
    SNR = datas - Background;
elseif flag1 == 2 
    datap = findpeaks(data);
    datap1 = sort(datap,'descend') ; % 该波束幅值排序
    dataNmax = datap1(1:Nmax);      % 该波束幅值的前5大
    sumdata = sum(datap1)- sum(dataNmax);    % 除了前5个大值的总和
    Background = sumdata/( N - Nmax ) ;  % 背景（除了前5个大值的平均）
    SNR = datas - Background;
    
end

end
    