function [M,T]=max_peak(P,theta_angle,model)
% P=abs(P);
if model==1   % 找出第一个极大值点（最值点）
    P_fengzhi2=-1e30*ones(1,length(theta_angle));
    T_fengzhi2=-1e30*ones(1,length(theta_angle));
    for i=2:length(theta_angle)-1
        if P(i)>P(i-1)&&P(i)>P(i+1)
            P_fengzhi2(i)=(P(i));
            T_fengzhi2(i)=theta_angle(i);
        end
    end
    [MAX,I]=max(P_fengzhi2);
    M=MAX;
    T=theta_angle(I);
end

if model==2   % 找出前两个极大值点（最大的两个极值）
    P_fengzhi2=-1e30*ones(1,length(theta_angle));
    T_fengzhi2=-1e30*ones(1,length(theta_angle));
    for i=2:length(theta_angle)-1
        if P(i)>P(i-1)&&P(i)>P(i+1)
            P_fengzhi2(i)=P(i);
            T_fengzhi2(i)=theta_angle(i);
        end
    end
    [MAX1,I1]=max(P_fengzhi2);
    M1=MAX1;
    P_fengzhi2(I1)=-1e30;
    [MAX2,I2]=max(P_fengzhi2);
    M2=MAX2;
    M=[M1 M2];
    T=[theta_angle(I1) theta_angle(I2)];
    T=sort(T);
end

if model==3     % 找出前三个极大值点（最大的三个极值）
    P_fengzhi2=-1e30*ones(1,length(theta_angle));
    T_fengzhi2=-1e30*ones(1,length(theta_angle));
    for i=2:length(theta_angle)-1
        if P(i)>P(i-1)&&P(i)>P(i+1)
            P_fengzhi2(i)=P(i);
            T_fengzhi2(i)=theta_angle(i);
        end
    end
    [MAX1,I1]=max(P_fengzhi2);
    M1=MAX1;
    P_fengzhi2(I1)=-1e30;
    [MAX2,I2]=max(P_fengzhi2);
    M2=MAX2;
    P_fengzhi2(I2)=-1e30;
    [MAX3,I3]=max(P_fengzhi2);
    M3=MAX3;
    
    M=[M1 M2 M3];
    T=[theta_angle(I1) theta_angle(I2) theta_angle(I3)];
    T=sort(T);
end

end