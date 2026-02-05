
clc;
clear all;
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

n = 1:N;
%% 产生声压通道信号

params.c = c;
params.fs = fs;
params.N_phi = 64;

phi0 = 70/180*pi;
theta0 = 90/180*pi;

k0_vector = -1*k0*[sin(theta0)*sin(phi0); sin(theta0)*cos(phi0); cos(theta0)];
p0_0 = exp(-1*1i*posi*k0_vector);

k1_vector = -1*k1*[sin(theta0)*sin(phi0); sin(theta0)*cos(phi0); cos(theta0)];
p0_1 = exp(-1*1i*posi*k1_vector);
%% 加噪
fL1 = 100;
fU1 = 200;

h = fir1(512,[fL1,fU1]/fs*2);         

Sigma_Sig1 = 1*10^(10/20);
Sigma_Sig2 = 1*10^(10/20);
Sigma_Sig3 = 1*10^(-3/20);
Sigma_Noi = 1;

Noise_C01 = sinf_3D(posi.',N,params);
Noise_C1 = (filter(h,1,Noise_C01.')).';
Noise_C1 = Noise_C1/sqrt(var(Noise_C1(1,:)));

BW = fU1 - fL1;
X = exp(1i*2*pi*f(1)*t);    % 信号
data1 = X;

X2 = exp(1i*2*pi*f(2)*t);    % 信号
data2 = X2;

Signals1 = p0_0 * data1.';
Signals1 = Signals1/sqrt(var(Signals1(1,:)));

Signals2 = p0_1 * data2.';
Signals2 = Signals2/sqrt(var(Signals2(1,:)));

S_12 = Sigma_Sig1*Signals1 + Sigma_Sig2*Signals2;
Noise_12 = Sigma_Noi*sqrt(BW)*Noise_C1;

data = Sigma_Sig1*Signals1 + Sigma_Sig2*Signals2 + Noise_12 ;


%%
S = S_12(6,:).';
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
for i=1:110:N-M+1   
    yyy(:,i) = X1(i:1:M+i-1);
end
                      
L = (N-M+1);                            

RR = yyy*yyy'/L;

autocorr = RR(:,1);
autocorr_FFT = fft(autocorr,length(autocorr));
f_p = (0:length(autocorr)-1)/length(autocorr)*fs;
%% TRCIS
t_l = (1:N)*Ts;        

t_r1 = -fliplr(t_l);
t_r1(1) = [];
t_l1 = t_l;
t_l1(end)=[];
t_l1 = [0 t_l1];
t_2 = [t_r1 t_l1];

R_X = xcorr(X1);                  
R_S = xcorr(S);
R_N = xcorr(Noise);

E = ones(1,2*N-1);
for k=1:(2*N-1)
    if k<=N-1
        E(k) = 2./(T+t_r1(k));
    else
        E(k) = 2./(T-t_l1(k-N+1));
    end
end

range_Interference_rej = 100;

I = ones(1,2*N-1);
I(N-range_Interference_rej:N+range_Interference_rej) = 0;

y_sfjj = R_X.'.*E.*I;

s_sfjj_junheng = R_S.'.*E;    
n_sfjj_junheng = R_N.'.*E;   

s_sfjj = s_sfjj_junheng.*I;
n_sfjj = n_sfjj_junheng.*I;

R_sfjj=zeros(1,N);
for i=1:N
   R_sfjj(i) = y_sfjj(N+i-1); 
end

for i=1:110:(N)-M+1   
    yyy_sfjj(:,i) = R_sfjj(i:1:M+i-1);
end

RR_sfjj = yyy_sfjj*yyy_sfjj'/L;         
%% FIM-Capon法 (将空域DOA估计中的FIM-Capon方法引入到谱估计中)
M1 = 900;         

RR_FC = my_FIM_Capon(RR,M1);
RR_sfjj_FC = my_FIM_Capon(RR_sfjj,M1);
%% 各种方法频谱计算 CBF Capon FC TRCIS
KK = 5001;
f_scan_pinlv = linspace(fL1,fU1,KK);  

f_scan = f_scan_pinlv/fs;             

AA = exp(1i*(0:M-1).'*2*pi*f_scan); 

P_CBF = abs(diag(AA'*RR*AA));
% P_sfjj=abs(diag(AA'*RR_sfjj*AA));

P_sfjj = abs(1./diag(AA'*(RR_sfjj+0.01*eye(M))^-1*AA));
P_Capon = abs(1./diag(AA'*(RR+0.01*eye(M))^-1*AA));

AA_FC = exp(-1i*(0:M1-1).'*2*pi*f_scan);
P_FC = abs(1./diag(AA_FC'*RR_FC^-1*AA_FC));
P_TRCIS_FC = abs(1./diag(AA_FC'*(RR_sfjj_FC+0.01*eye(M1))^-1*AA_FC));

P_CBF_dB = 10*log10(abs(P_CBF)/max(abs(P_CBF)));
P_Capon_dB = 10*log10(abs(P_Capon)/max(abs(P_Capon)));
P_sfjj_dB = 10*log10(abs(P_sfjj)/max(abs(P_sfjj)));
P_FC_dB = 10*log10(abs(P_FC)/max(abs(P_FC)));
P_TRCIS_FC_dB = 10*log10(abs(P_TRCIS_FC)/max(abs(P_TRCIS_FC)));

% MUSIC
[EV,D] = eig(RR);                     % 特征值分解

N_num_En = K;                          % 选取前几个特征向量为噪声子空间
En = EV(:,M-N_num_En:-1:1);

P_MUSIC = 1./diag(AA'*(En*En')*AA);
P_MUSIC_dB = 10*log10(abs(P_MUSIC)/max(abs(P_MUSIC)));

% MUSIC
[EV2,D2] = eig(RR_sfjj);                     % 特征值分解

N_num_En2 = K;                          % 选取前几个特征向量为噪声子空间
En2 = EV2(:,M-N_num_En2:-1:1);

P_NSG_MUSIC = 1./diag(AA'*(En2*En2')*AA);
P_NSG_MUSIC_dB = 10*log10(abs(P_NSG_MUSIC)/max(abs(P_NSG_MUSIC)));

%% 寻找峰值点
[Max_peak_CBF,T_peak_CBF] = max_peak(P_CBF_dB,f_scan_pinlv,2);                             % CBF
[Max_peak_Capon,T_peak_Capon] = max_peak(P_Capon_dB,f_scan_pinlv,2);                       % Capon
[Max_peak_FC,T_peak_FC] = max_peak(P_FC_dB,f_scan_pinlv,2);                                % FIM-Capon
[Max_peak_TRCIS_FC,T_peak_TRCIS_FC] = max_peak(P_TRCIS_FC_dB,f_scan_pinlv,2);              % TRCIS-FC

[Max_peak_MUSIC,T_peak_MUSIC] = max_peak(P_MUSIC_dB,f_scan_pinlv,2);                       % MUSIC
%% 画图

figure
plot(f_scan_pinlv ,P_Capon_dB,'b','linewidth',2.5);
hold on;
plot(f_scan_pinlv ,P_MUSIC_dB,'k-.','linewidth',2.5);
hold on;
plot(f_scan_pinlv ,P_FC_dB,'r:','linewidth',2.5);
hold on;
plot(f_scan_pinlv ,P_TRCIS_FC_dB,'m--','linewidth',2.5);
hold on;
for i=1:K
    plot([f(i) f(i)],[-70 0],'g','linewidth',1.5);
    hold on;
end
legend('\fontname{times new roman}Capon','\fontname{times new roman}MUSIC','\fontname{times new roman}TA','\fontname{times new roman}NSG-TA','\fontname{Times new roman}Real-Fre','Location','south');
xlabel('\fontname\{宋体}频率\fontname\{times new roman}(Hz)');
ylabel('\fontname\{宋体}幅度\fontname\{times new roman}(dB)');
grid on;
set(gca,'FontSize',15);
axis([140 180 -50 0]);

figure
plot(t_2,real(R_S)./max(real(R_S)));
xlabel('\fontname\{宋体}时间\fontname\{times new roman}(s)');
ylabel('\fontname\{宋体}归一化幅度');
set(gca,'FontSize',15);

figure
plot(t_2,R_N./max(R_N));
xlabel('\fontname\{宋体}时间\fontname\{times new roman}(s)');
ylabel('\fontname\{宋体}归一化幅度');
set(gca,'FontSize',16);

figure
plot(t_2,I,'linewidth',1.5);
xlabel('\fontname\{宋体}时间\fontname\{times new roman}(s)');
ylabel('\fontname\{宋体}归一化幅度');
set(gca,'FontSize',16);
ylim([-0.5 2]);

figure
plot(t_2(200:end-200),E(200:end-200),'linewidth',1.5);
xlabel('\fontname\{宋体}时间\fontname\{times new roman}(s)');
ylabel('\fontname\{宋体}归一化幅度');
set(gca,'FontSize',16);
