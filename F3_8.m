clc;
clear;
close all;
%% 阵列设置
M = 12;
r = 1;

phi_s0 = (0:M-1)*2*pi/M;
theta_s = pi/2;
posi = r*[sin(theta_s)*cos(phi_s0'), sin(theta_s)*sin(phi_s0'), cos(theta_s)*ones(M,1)];
posi(:,3) = 0;
SensorPosition = [posi'];

%% 参数设置
c = 1500;
N = 20000;                                
K = 2;                                   

f = [160 162 151 418 197]; 
k0 = 2*pi*f(1)/c;
k1 = 2*pi*f(2)/c;

ff = f/1000;

fs = 5000;                            
Ts = 1/fs;
T = 4;
tt = 0:Ts:T;
t = tt(1:N).';

% t1=(0:M-1)*Ts;
n = 1:N;

%% 产生声压通道信号
% Generate sensor signals
params.c = c;
params.fs = fs;
params.N_phi = 64;

phi0 = 70/180*pi;
theta0 =90/180*pi;

k0_vector = -1*k0*[sin(theta0)*sin(phi0); sin(theta0)*cos(phi0); cos(theta0)];
p0_0 = exp(-1*1i*posi*k0_vector);

k1_vector = -1*k1*[sin(theta0)*sin(phi0); sin(theta0)*cos(phi0); cos(theta0)];
p0_1 = exp(-1*1i*posi*k1_vector);
%% 循环
SNR = 0:10:20;
MC = 1;
f_Es1_MUSIC = zeros(length(SNR),MC);
f_Es2_MUSIC = zeros(length(SNR),MC);

f_Es1_Capon = zeros(length(SNR),MC);
f_Es2_Capon = zeros(length(SNR),MC);

f_Es1_FC = zeros(length(SNR),MC);
f_Es2_FC = zeros(length(SNR),MC);

f_Es1_TRCIS_FC = zeros(length(SNR),MC);
f_Es2_TRCIS_FC = zeros(length(SNR),MC);

%% 加噪
fL1 = 100;
fU1 = 200;

h = fir1(512,[fL1,fU1]/fs*2);           % fir1 设计fir滤波器 带通和低通ftype无需输入 高通输入high 带阻输入stop 512为滤波器阶数 [fL1,fU1]/fs*2为带通范围

for tt = 1:length(SNR)
    C1 = 0;
    C2 = 0;
    C3 = 0;
    C4 = 0;
    for e = 1:MC
        Sigma_Sig1 = 1*10^(SNR(tt)/20);
        Sigma_Sig2 = 1*10^(SNR(tt)/20);
        Sigma_Sig3 = 1*10^(-3/20);
        Sigma_Noi = 1;
        
        Noise_C01 = sinf_3D(posi.',N,params);
        Noise_C1 = (filter(h,1,Noise_C01.')).';
        Noise_C1 = Noise_C1/sqrt(var(Noise_C1(1,:)));
        
        BW = fU1 - fL1;
        
        X = exp(1i*2*pi*f(1)*t);    % 信号
        data1 = hilbert(X.').';
        
        X2 = exp(1i*2*pi*f(2)*t);    % 信号
        data2 = hilbert(X2.').';
        
        Signals1 = p0_0 * data1.';
        Signals1 = Signals1/sqrt(var(Signals1(1,:)));
        
        Signals2 = p0_1 * data2.';
        Signals2 = Signals2/sqrt(var(Signals2(1,:)));
        
        data = Sigma_Sig1*Signals1 + Sigma_Sig2*Signals2 + Sigma_Noi*sqrt(BW)*Noise_C1 ;
        
        Noise_12 = Sigma_Noi*sqrt(BW)*Noise_C1;
        Noise = Noise_12(6,:).';
        X1 = data(6,:).';
        
        NFFT = fs;
        y_fft = fft(X1,NFFT);
        y_f = fs*(0:NFFT-1)/NFFT;
        
        NFFT2 = length(X1);
        y_fft2 = fft(X1,NFFT2);
        y_f2 = fs*(0:NFFT2-1)/NFFT2;
        %% 协方差矩阵的求解
        M = 1024;                   
        for i = 1:110:N-M+1
            yyy(:,i) = X1(i:1:M+i-1);
        end
        
        L = (N-M+1);                             
        RR = yyy*yyy'/L;
        
        %% TRCIS
        t_l = (1:N)*Ts;        
        
        t_r1 = -fliplr(t_l);
        t_r1(1) = [];
        t_l1 = t_l;
        t_l1(end) = [];
        t_l1 = [0 t_l1];
        t_2 = [t_r1 t_l1];
        
        R_X = xcorr(X1);                 
        R_S = xcorr(X1);
        R_N = xcorr(Noise);
        
        E = ones(1,2*N-1);
        for k = 1:(2*N-1)
            if k<=N-1
                E(k) = 2./(T+t_r1(k));
            else
                E(k) = 2./(T-t_l1(k-N+1));
            end
        end
        
        range_Interference_rej = 10;
        
        I = ones(1,2*N-1);
        I(N-range_Interference_rej:N+range_Interference_rej)=0;
        
        y_sfjj = R_X.'.*E.*I;
        
        
        R_sfjj = zeros(1,N);
        for i = 1:N
            R_sfjj(i) = y_sfjj(N+i-1);
        end
        
        for i = 1:110:(N)-M+1
            yyy_sfjj(:,i) = R_sfjj(i:1:M+i-1);
        end
        
        RR_sfjj = yyy_sfjj*yyy_sfjj'/L;         
        
        %% FIM-Capon法 (将空域DOA估计中的FIM-Capon方法引入到谱估计中)
        M1 = 900;        
        
        RR_FC = my_FIM_Capon(RR,M1);
        RR_sfjj_FC = my_FIM_Capon(RR_sfjj,M1);
        %% 各种方法频谱计算 CBF Capon FC TRCIS
        KK = 1001;
        f_scan_pinlv = linspace(fL1,fU1,KK);  
        
        f_scan = f_scan_pinlv/fs;             
        
        AA = exp(1i*(0:M-1).'*2*pi*f_scan);
        
        P_CBF = abs(diag(AA'*RR*AA));
        P_sfjj = abs(diag(AA'*RR_sfjj*AA));
        
        P_Capon = abs(1./diag(AA'*(RR+0.001*eye(M))^-1*AA));
        
        AA_FC = exp(-1i*(0:M1-1).'*2*pi*f_scan);
        P_FC = abs(1./diag(AA_FC'*RR_FC^-1*AA_FC));
        P_TRCIS_FC = abs(1./diag(AA_FC'*RR_sfjj_FC^-1*AA_FC));
        
        P_CBF_dB = 10*log10(abs(P_CBF)/max(abs(P_CBF)));
        P_Capon_dB = 10*log10(abs(P_Capon)/max(abs(P_Capon)));
        P_sfjj_dB = 10*log10(abs(P_sfjj)/max(abs(P_sfjj)));
        P_FC_dB = 10*log10(abs(P_FC)/max(abs(P_FC)));
        P_TRCIS_FC_dB = 10*log10(abs(P_TRCIS_FC)/max(abs(P_TRCIS_FC)));
        
        % MUSIC
        [EV,D] = eig(RR);                     
        
        N_num_En = K;                          
        En = EV(:,M-N_num_En:-1:1);
        
        P_MUSIC = 1./diag(AA'*(En*En')*AA);
        P_MUSIC_dB = 10*log10(abs(P_MUSIC)/max(abs(P_MUSIC)));
        %% 寻找峰值点
      
        [Max_peak_Capon,T_peak_Capon] = max_peak(P_Capon_dB,f_scan_pinlv,2);                % Capon
        [Max_peak_MUSIC,T_peak_MUSIC] = max_peak(P_MUSIC_dB,f_scan_pinlv,2);                % MUSIC
        [Max_peak_FC,T_peak_FC] = max_peak(P_FC_dB,f_scan_pinlv,2);                         % FIM-Capon
        [Max_peak_TRCIS_FC,T_peak_TRCIS_FC] = max_peak(P_TRCIS_FC_dB,f_scan_pinlv,2);       % TRCIS-FC
                
        f_Es1_Capon(tt,e) = T_peak_Capon(1);
        f_Es2_Capon(tt,e) = T_peak_Capon(2);
        if abs(f_Es1_Capon(tt,e)-f(1)) + abs(f_Es2_Capon(tt,e)-f(2))<abs(f(1)-f(2))
            C1 = C1+1;
        end
        
        f_Es1_MUSIC(tt,e) = T_peak_MUSIC(1);
        f_Es2_MUSIC(tt,e) = T_peak_MUSIC(2);
        if abs(f_Es1_MUSIC(tt,e)-f(1)) + abs(f_Es2_MUSIC(tt,e)-f(2))<abs(f(1)-f(2))
            C2 = C2+1;
        end
        
        f_Es1_FC(tt,e) = T_peak_FC(1);
        f_Es2_FC(tt,e) = T_peak_FC(2);
        if abs(f_Es1_FC(tt,e)-f(1)) + abs(f_Es2_FC(tt,e)-f(2))<abs(f(1)-f(2))
            C3 = C3+1;
        end
        
        f_Es1_TRCIS_FC(tt,e) = T_peak_TRCIS_FC(1);
        f_Es2_TRCIS_FC(tt,e) = T_peak_TRCIS_FC(2);
        if abs(f_Es1_TRCIS_FC(tt,e)-f(1)) + abs(f_Es2_TRCIS_FC(tt,e)-f(2))<abs(f(1)-f(2))
            C4 = C4+1;
        end
       
    end
    
    RP1(tt) = C1/MC;
    RP2(tt) = C2/MC;
    
    RP3(tt) = C3/MC;
    RP4(tt) = C4/MC;
    
end
%% 画图

t1_RMSE = f(1)*ones(length(SNR),1);
t2_RMSE = f(2)*ones(length(SNR),1);

RMSE_Capon = RMSE(f_Es1_Capon,f_Es2_Capon,t1_RMSE,t2_RMSE,MC);
RMSE_MUSIC = RMSE(f_Es1_MUSIC,f_Es2_MUSIC,t1_RMSE,t2_RMSE,MC);
RMSE_FC = RMSE(f_Es1_FC,f_Es2_FC,t1_RMSE,t2_RMSE,MC);
RMSE_TRCIS_FC = RMSE(f_Es1_TRCIS_FC,f_Es2_TRCIS_FC,t1_RMSE,t2_RMSE,MC);

figure 
plot(SNR,RMSE_Capon,'b--*','linewidth',1.5);
hold on;
plot(SNR,RMSE_MUSIC,'r--s','linewidth',1.5);
hold on;
plot(SNR,RMSE_FC,'k--d','linewidth',1.5);
hold on;
plot(SNR,RMSE_TRCIS_FC,'m--x','linewidth',1.5);
hold on;
grid on;
axis([SNR(1) SNR(end) 0 2]);
xlabel('\fontname{times new roman}SNR(dB)');
ylabel('\fontname{times new roman}RMSE(\circ)');
legend('\fontname{Times new roman}Capon','\fontname{Times new roman}MUSIC','\fontname{Times new roman}TA','\fontname{Times new roman}NSG-TA','FontSize',15);
set(gca,'FontSize',15);

figure 
plot(SNR,RP1,'b--*','linewidth',1.5);
hold on;
plot(SNR,RP2,'r--s','linewidth',1.5);
hold on;
plot(SNR,RP3,'k--d','linewidth',1.5);
hold on;
plot(SNR,RP4,'m--x','linewidth',1.5);
hold on;
grid on;
axis([SNR(1) SNR(end) 0 1]);
xlabel('\fontname{times new roman}SNR(dB)');
ylabel('\fontname{times new roman}RP');
legend('\fontname{Times new roman}Capon','\fontname{Times new roman}MUSIC','\fontname{Times new roman}TA','\fontname{Times new roman}NSG-TA','FontSize',15);
% legend('\fontname{Times new roman}Capon','\fontname{Times new roman}FIM-Capon','\fontname{Times new roman}CMR-APC','FontSize',13);
set(gca,'FontSize',15);

RMSE_0402(1,:) = RMSE_Capon;
RMSE_0402(2,:) = RMSE_MUSIC;
RMSE_0402(3,:) = RMSE_FC;
RMSE_0402(4,:) = RMSE_TRCIS_FC;
RMSE_0402(5,:) = SNR;

RP_0402(1,:) = RP1;
RP_0402(2,:) = RP2;
RP_0402(3,:) = RP3;
RP_0402(4,:) = RP4;
RP_0402(5,:) = SNR;
