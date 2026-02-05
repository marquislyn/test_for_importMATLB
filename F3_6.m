clc;
clear;
close all;
%% 信号参数设置
f = 180;              % 信号频率
fs = 2000;            % 采样频率
Ts = 1/fs;             % 采样间隔

T = 50;                % 时间

Z = fs*T;                           % 数字信号的长度
t = 0:Ts:T-Ts;                      % 时间信号时域向量
sin_Amp = 1;                        % 正弦信号的幅度
randn_Amp = 2;                      % 随机噪声的标准差
NFFT = 2*fs;

or = 0.5;                             % overlap rate
R = NFFT*(1-or);

L2 = 1+floor((Z-NFFT)/R);           % 分段数

S = sin(2*pi*f*t);                  % 信号

f_p = fs*(0:(NFFT-1))/NFFT;         % 或者f=Fs/2*linspace(0,1,L1/2+1)

SNR = 0:10:20;
MC = 1;

CC = length(SNR)*MC;            % 需要循环的总次数
count = 0;                      % 开始

SNR_AVGPR_MC = zeros(MC,length(SNR));
SNR_PCAVG_MC = zeros(MC,length(SNR));
SNR_RDA_MC = zeros(MC,length(SNR));
SNR_AVGPR = zeros(1,length(SNR));
SNR_PCAVG = zeros(1,length(SNR));
SNR_RDA = zeros(1,length(SNR));

fL1 = 10;
fU1 = 500;
h = fir1(512,[fL1,fU1]/fs*2);   % 滤波器参数

for ii = 1:length(SNR)
    for qq = 1:MC
        
        %% 加高斯白噪声
        
%                 [X,NOISE]= noisegen(S,SNR(ii));
        
        % 高斯白噪声
        Noise_C0 = randn(1,Z+1*fs);
        Noise_C = (filter(h,1,Noise_C0));
        Noise_C = Noise_C/sqrt(var(Noise_C(1,:)));
        
        Sigma_Sig1 = 1*10^(SNR(ii)/20);
        Sigma_Noi = 1;
        
        BW = fU1 - fL1;
        Noise = Sigma_Noi*sqrt(BW)*Noise_C(1,fs+1:end);
        X = Sigma_Sig1 * S + Noise;
        
        %% 对时域信号做FFT
        YY = zeros(L2,NFFT);
        YY_abs = zeros(L2,NFFT);
        Y_FFT_1 = zeros(L2,NFFT);
        
        for i = 1:L2
            xl = (i-1)*R+1;
            xr = (i-1)*R+NFFT;
            Y = fft(X(xl:xr),NFFT);        % 对有1000个离散点的时域信号做1000点的FFT
            
            Y_FFT_1(i,:) = abs(Y);
            YY(i,:) = Y;
            YY_abs(i,:) = abs(Y);
        end
        %% RDA法 Rotational Differential Average, RDA
        YY_angle = angle(YY);
        
        Fhi2 = zeros(L2,NFFT);
        for n = 1:NFFT
            for l = 3:L2
                Fhi2(l,n) = YY_angle(l,n)-2*YY_angle(l-1,n)+YY_angle(l-2,n);
            end
        end
        
        R1 = YY_abs.*cos(Fhi2);
        R2 = YY_abs.*sin(Fhi2);
        
        R1_sum = (1/L2)*sum(R1);
        R2_sum = (1/L2)*sum(R2);
        
        R1_sum = R1_sum.^2;
        R2_sum = R2_sum.^2;
        
        T_RDA = R1_sum + R2_sum;
        T_RDA_dB = 10*log10(T_RDA);
        T_RDA_dB1 = 10*log10(T_RDA./max(T_RDA));
        
        %% PCAVG
        Y_PCAVG = YY;
        [T_PCAVG] = PCAVG(Y_PCAVG,or);
        T_PCAVG = abs(T_PCAVG);
        
        T_PCAVG_dB1 = 10*log10((T_PCAVG)/max((T_PCAVG)));
        %% AVGPR平均周期图法
        
        T_AVGPR = (1/(L2))*sum(YY_abs.*YY_abs);
        
        T_AVGPR_dB = 10*log10(T_AVGPR);
        T_AVGPR_dB1 = 10*log10(T_AVGPR./max(T_AVGPR));
        
        %% 计算信噪比 
        xianpu_xuanqu = (f*NFFT/fs+1)-fL1*NFFT/fs;
        
        SNR_AVGPR_MC(qq,ii) = SNRcalculation(T_AVGPR_dB1(fL1*NFFT/fs+1:fU1*NFFT/fs+1), 5, 0, 1 ,xianpu_xuanqu);
        SNR_PCAVG_MC(qq,ii) = SNRcalculation(T_PCAVG_dB1(fL1*NFFT/fs+1:fU1*NFFT/fs+1), 5, 0, 1 ,xianpu_xuanqu);
        SNR_RDA_MC(qq,ii) = SNRcalculation(T_RDA_dB1(fL1*NFFT/fs+1:fU1*NFFT/fs+1), 5, 0, 1 ,xianpu_xuanqu);
        count = count + 1;
        fprintf('目前的次数%.0f，需要的总次数%.0f\n',count,CC);
         
    end
    SNR_AVGPR(ii) = sum(SNR_AVGPR_MC(:,ii))/MC;
    SNR_PCAVG(ii) = sum(SNR_PCAVG_MC(:,ii))/MC;
    SNR_RDA(ii) = sum(SNR_RDA_MC(:,ii))/MC;
end
figure
plot(SNR,SNR_AVGPR,'b--*','linewidth',2.5);
hold on;
plot(SNR,SNR_PCAVG,'k--x','linewidth',2.5);
hold on;
plot(SNR,SNR_RDA,'r--s','linewidth',2.5);
hold on;
grid on;
legend('\fontname{times new roman}AVGPR','\fontname{times new roman}PCAVG','\fontname{times new roman}PAV');
xlabel('\fontname\{宋体}输入信噪比\fontname{times new roman}(dB)');
ylabel('\fontname\{宋体}输出信噪比\fontname{times new roman}(dB)');
xlim([SNR(1) SNR(end)]);
set(gca,'FontSize',16,'FontName','Times New Roman');

Gain_AVGPR = SNR_AVGPR-SNR;
Gain_PCAVG = SNR_PCAVG-SNR;
Gain_RDA = SNR_RDA-SNR;

Gain_AVGPR2 = SNR_AVGPR-SNR_AVGPR;
Gain_PCAVG2 = SNR_PCAVG-SNR_AVGPR;
Gain_RDA2 = SNR_RDA-SNR_AVGPR;

figure
plot(SNR,Gain_AVGPR,'b--*','linewidth',2.5);
hold on;
plot(SNR,Gain_PCAVG,'k--x','linewidth',2.5);
hold on;
plot(SNR,Gain_RDA,'r--s','linewidth',2.5);
hold on;
grid on;
legend('\fontname{times new roman}AVGPR','\fontname{times new roman}PCAVG','\fontname{times new roman}PAV');
xlabel('\fontname\{宋体}输入信噪比\fontname{times new roman}(dB)');
ylabel('\fontname\{宋体}增益\fontname{times new roman}(dB)');
xlim([SNR(1) SNR(end)]);
set(gca,'FontSize',16,'FontName','Times New Roman');

figure
plot(SNR,Gain_AVGPR2,'b--*','linewidth',2.5);
hold on;
plot(SNR,Gain_PCAVG2,'k--x','linewidth',2.5);
hold on;
plot(SNR,Gain_RDA2,'r--s','linewidth',2.5);
hold on;
grid on;
legend('\fontname{times new roman}AVGPR','\fontname{times new roman}PCAVG','\fontname{times new roman}PAV');
xlabel('\fontname\{宋体}输入信噪比\fontname{times new roman}(dB)');
ylabel('\fontname\{宋体}增益\fontname{times new roman}(dB)');
xlim([SNR(1) SNR(end)]);
set(gca,'FontSize',16,'FontName','Times New Roman');