function Mse=RMSE(DOA1,DOA2,theta_DOA1,theta_DOA2,MC)
% DOA1:h*MC  为衡量的指标向量长度（如SNR、AI、NVAE、f）*Monte Carlo实验次数
%

MSE1_Capon=(DOA1-theta_DOA1).^2;
MSE2_Capon=(DOA2-theta_DOA2).^2;
Mse=sqrt((sum(MSE1_Capon,2)+sum(MSE2_Capon,2))/(2*MC));
% Mse(Mse>4)=4;
end