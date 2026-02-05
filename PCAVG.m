function [T_PCAVG] = PCAVG(Y_PCAVG,ol)

% ff -- 预先假定的线谱的频率（对于PCAVG和CGLRT关键的数据）
% ol -- 重叠率
% 此PCAVG算法将所有频点补偿
%%
L2 = length(Y_PCAVG(:,1));       % L2 -- 分的段数
N = length(Y_PCAVG(1,:));        % N  -- 每段的数据个数

Elenum = 1;
R = N*(1-ol);

% 第一次补偿
for ii = 0:N-1
    BC1 = zeros(L2,1);
    for i = 1:L2
        BC1(i) = exp(-1i*2*pi*(ii/N)*R*(i-1));
    end
    for n = 1:Elenum
        Y_PCAVG(:,ii+1,n) = Y_PCAVG(:,ii+1,n).*BC1;
    end
    
    % 第二次补偿
    for n = 1:Elenum
        A_FFT1 = fft(Y_PCAVG(:,ii+1,n),L2);
        A_FFT1 = abs(A_FFT1);
        [~,II1] = max(A_FFT1);
        
        if II1>(L2+1)/2
            kmax1 = L2+1-II1;
            alpha1g = -(N*kmax1)/(R*L2);
        else
            kmax1 = II1-1;
            alpha1g = (N*kmax1)/(R*L2);
        end
        
        BAS1 = zeros(L2,1);
        
        for i = 1:(L2)
            BAS1(i) = exp(-1i*2*pi*(alpha1g/N)*R*(i-1));
        end
        
        Y_PCAVG(:,ii+1,n) = Y_PCAVG(:,ii+1,n).*BAS1;
    end
    
    YY_PCAVG = sum(Y_PCAVG.*Y_PCAVG);
    T_PCAVG1 = 1/L2*YY_PCAVG;
    for n = 1:Elenum
        T_PCAVG1(1,1,n) = 0;
    end
    
    T_PCAVG = zeros(Elenum,N);
    for k = 1:Elenum
        T_PCAVG(k,:) = T_PCAVG1(:,:,k);
    end
end
% v0S = (1/L2)*sum(abs(Y_PCAVG).*abs(Y_PCAVG));
% uS = (1/L2)*sum(Y_PCAVG);
% umS = zeros(L2,N,Elenum);
% for n = 1:Elenum
%     for i = 1:L2
%         umS(i,:,n) = uS(:,:,n);
%     end
% end
% vmmS = Y_PCAVG-umS;
% v1S = (1/L2)*sum(abs(vmmS).*abs(vmmS));
%
% T_CGLRT1 = (-L2)*log(v1S./v0S);
% for n = 1:Elenum
%     T_CGLRT1(1,1,n) = 0;
% end
%
% T_CGLRT = zeros(Elenum,N);
% for k = 1:Elenum
%     T_CGLRT(k,:) = T_CGLRT1(:,:,k);
% end
end