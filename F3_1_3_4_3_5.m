
clc;
clear;
close all;

f = 180;              
fs = 2000;             % 采样频率
Ts = 1/fs;             % 采样间隔
SNR_in  = 0;
T = 50;                % 时间

Z = fs*T;               % 数字信号的长度
t = 0:Ts:T-Ts;          % 时间信号时域向量
sin_Amp = 1;            % 正弦信号的幅度
randn_Amp = 2;          % 随机噪声的标准差
NFFT = 2*fs;

or = 0.5;                 % overlap rate
R = NFFT*(1-or);

L2 = 1+floor((Z-NFFT)/R);           % 分段数

S = sin_Amp*sin(2*pi*f*t);          % 信号
%% 带噪信号的生成
fL1 = 10;
fU1 = 500;
h = fir1(512,[fL1,fU1]/fs*2);   % 滤波器参数

% [X,NOISE]= noisegen(S,10);
% 高斯白噪声
Noise_C0 = randn(1,Z+1*fs);
Noise_C = (filter(h,1,Noise_C0));
Noise_C = Noise_C/sqrt(var(Noise_C(1,:)));

Sigma_Sig1 = 1*10^(SNR_in/20);
Sigma_Noi = 1;

BW = fU1 - fL1;
Noise = Sigma_Noi*sqrt(BW)*Noise_C(1,fs+1:end);
X = Sigma_Sig1 * S + Noise;
%% 对时域信号做FFT
f_p = fs*(0:(NFFT-1))/NFFT;         

YY = zeros(L2,NFFT);
YY_abs = zeros(L2,NFFT);

for i = 1:L2
    xl = (i-1)*R+1;
    xr = (i-1)*R+NFFT;
    Y = fft(X(xl:xr),NFFT);        

    YY(i,:) = Y;
    YY_abs(i,:) = abs(Y);
end

% SNR_AVGPR = SNRcalculation(YY_abs(1,fL1*NFFT/fs+1:fU1*NFFT/fs+1), 5, 1, 1 ,xianpu_xuanqu);
% YY_abs_sum  = sum(YY_abs);
% figure;
% plot(f_p,20*log10(YY_abs_sum/max(YY_abs_sum)))
% SNR1 = SNRcalculation(YY_abs_sum(fL1*NFFT/fs+1:fU1*NFFT/fs+1), 5, 1, 1 ,xianpu_xuanqu);

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

T_RDA = R1_sum+R2_sum;
T_RDA_dB = 10*log10(T_RDA);
T_RDA_dB1 = 10*log10(T_RDA./max(T_RDA));

%% PCAVG
Y_PCAVG = YY;
[T_PCAVG] = PCAVG(Y_PCAVG,or);
T_PCAVG = abs(T_PCAVG);

T_PCAVG_dB1 = 10*log10((T_PCAVG)/max((T_PCAVG)));
%% AVGPR平均周期图法
T_AVGPR = (1/(L2))*sum(YY_abs.*YY_abs);
T_AWSUM = ((1/(L2))*sum(YY_abs.^(-4))).^(-1/4);
T_AVGPR_dB = 10*log10(T_AVGPR);
T_AVGPR_dB1 = 10*log10(T_AVGPR./max(T_AVGPR));
T_AWSUM_dB1 = 10*log10(T_AWSUM./max(T_AWSUM));
%% AVGPRO法
T_AVGPRO = (1/(L2))*sum(YY_abs.^2.4);

%% 计算信噪比
xianpu_xuanqu = (f*NFFT/fs+1)-fL1*NFFT/fs;

xianpu_xuanqu2 = (f*NFFT/fs+1);

SNR_AVGPR = SNRcalculation(T_AVGPR_dB1(fL1*NFFT/fs+1:fU1*NFFT/fs+1), 5, 0, 1 ,xianpu_xuanqu);
SNR_PCAVG = SNRcalculation(T_PCAVG_dB1(fL1*NFFT/fs+1:fU1*NFFT/fs+1), 5, 0, 1 ,xianpu_xuanqu);
SNR_RDA = SNRcalculation(T_RDA_dB1(fL1*NFFT/fs+1:fU1*NFFT/fs+1), 5, 0, 1 ,xianpu_xuanqu);
SNR_AWSUM = SNRcalculation(T_AWSUM_dB1(fL1*NFFT/fs+1:fU1*NFFT/fs+1), 5, 0, 1 ,xianpu_xuanqu);
%% 画图
figure 
plot(f_p,T_AVGPR_dB1,'Linewidth',1.8);
title('\fontname{times new roman}AVGPR');
xlabel('\fontname\{宋体}频率\fontname\{times new roman}(Hz)');
ylabel('\fontname\{宋体}幅度\fontname\{times new roman}(dB)');
set(gca,'FontSize',16,'FontName','Times New Roman');
axis([fL1,fU1,-70,0]);
grid on;

figure 
plot(f_p,10*log10(T_RDA./max(T_RDA)),'Linewidth',1.8);
title('\fontname{times new roman}PAV');
xlabel('\fontname\{宋体}频率\fontname\{times new roman}(Hz)');
ylabel('\fontname\{宋体}幅度\fontname\{times new roman}(dB)');
set(gca,'FontSize',16,'FontName','Times New Roman');
axis([fL1,fU1,-70,0]);
grid on;

figure 
plot(f_p,T_PCAVG_dB1,'Linewidth',1.8);
title('\fontname{times new roman}PCAVG');
xlabel('\fontname\{宋体}频率\fontname\{times new roman}(Hz)');
ylabel('\fontname\{宋体}幅度\fontname\{times new roman}(dB)');
set(gca,'FontSize',16,'FontName','Times New Roman');
axis([fL1,fU1,-70,0]);
grid on;

figure 
plot(f_p,10*log10(T_AWSUM./max(T_AWSUM)),'Linewidth',1.8);
title('\fontname{times new roman}T_AWSUM');
xlabel('\fontname\{宋体}频率\fontname\{times new roman}(Hz)');
ylabel('\fontname\{宋体}幅度\fontname\{times new roman}(dB)');
set(gca,'FontSize',16,'FontName','Times New Roman');
axis([fL1,fU1,-70,0]);
grid on;
figure 
plot(f_p,10*log10(T_AVGPR/max(T_AVGPR)),'Linewidth',1.8);
hold on;
plot(f_p,T_PCAVG_dB1,'k-.','Linewidth',1.8);
hold on;
plot(f_p,10*log10(T_RDA/max(T_RDA)),'r:','Linewidth',1.8)
xlabel('\fontname\{宋体}频率\fontname\{times new roman}(Hz)');
ylabel('\fontname\{宋体}幅度\fontname\{times new roman}(dB)');
legend('\fontname{times new roman}AVGPR','\fontname{times new roman}PCAVG','\fontname{times new roman}PAV','Location','south');
set(gca,'FontSize',16,'FontName','Times New Roman');
axis([fL1,fU1,-50,0]);
grid on;


%% 画相位的
signal_phase_yuan = YY_angle(:,xianpu_xuanqu2)*180/pi;
signal_phase_duiqi = Fhi2(:,xianpu_xuanqu2)*180/pi;

noise_phase_yuan = YY_angle(:,xianpu_xuanqu2+153)*180/pi;
noise_phase_duiqi = Fhi2(:,xianpu_xuanqu2+153)*180/pi;

figure
subplot 211
scatter((1:L2),signal_phase_yuan,'Linewidth',2);
legend('\fontname\{宋体}原相位')
set(gca,'FontSize',15,'FontName','Times New Roman');
xlim([1 L2]);
xlabel('\fontname\{宋体}分段号');
ylabel('\fontname\{宋体}角度\fontname\{times new roman}(\circ)');

subplot 212
scatter((1:L2),signal_phase_duiqi,'r*','Linewidth',2);
legend('\fontname\{宋体}对齐相位')
axis([1 L2 -500 500]);
% xlim([1 L2]);
set(gca,'FontSize',15,'FontName','Times New Roman');
xlabel('\fontname\{宋体}分段号');
ylabel('\fontname\{宋体}角度\fontname\{times new roman}(\circ)');

figure
subplot 211
scatter((1:L2),noise_phase_yuan,'Linewidth',2);
legend('\fontname\{宋体}原相位')
set(gca,'FontSize',15,'FontName','Times New Roman');
xlim([1 L2]);
xlabel('\fontname\{宋体}分段号');
ylabel('\fontname\{宋体}角度\fontname\{times new roman}(\circ)');

subplot 212
scatter((1:L2),noise_phase_duiqi,'r*','Linewidth',2);
legend('\fontname\{宋体}对齐相位')
axis([1 L2 -500 500]);
set(gca,'FontSize',15,'FontName','Times New Roman');
xlabel('\fontname\{宋体}分段号');
ylabel('\fontname\{宋体}角度\fontname\{times new roman}(\circ)');