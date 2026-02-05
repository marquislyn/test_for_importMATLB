function R_SA=my_FIM_Capon(RR,M1)
% RR 为需要做平均的协方差矩阵
% M1 为子阵的长度，即子阵平滑长度
%%
M = length(RR(1,:));

C1 = 0;
C = zeros(2*M-1,1);
for  i = -M+1:-1
    for n = abs(i)+1:M
        C1 = RR(n,n+i)+C1;
    end
    C(i+M) = 1/(M-abs(i))*C1;
    C1=0;
end

for i = 0:M-1
    for n = 1:M-abs(i)
        C1 = RR(n,n+i)+C1;
    end
    C(i+M) = 1/(M-abs(i))*C1;
    C1 = 0;
end

R_FIM = C*C';                      % R_FIM矩阵

M_subarray = 2*M-M1;               % 相同信息子阵数量

R_SA = zeros(M1,M1);

for m = 1:M_subarray
    R_SA = R_FIM(m:m+M1-1,m:m+M1-1)+R_SA;
end

R_SA = R_SA/M_subarray;            % R_SA矩阵
end