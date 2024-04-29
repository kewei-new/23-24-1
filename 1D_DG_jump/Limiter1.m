function uh_mod = Limiter1(uh,h,M,G,space_order)
% M：关于二阶导数的上界，一般取M=max(abs(u_0,xx))，不知道取什么
% h：剖分步长，必须是均匀剖分

N = length(uh)/space_order;
uh_mod = zeros(N,space_order);


for i = 1:N
    
    % 当前单元的uh系数
    C = uh((i-1)*space_order+1:(i)*space_order,1);
    temp = G*C;
    
    % 求带入minmod的第一个参数
    ubar = temp(1,1);%当前单元的均值
    utilde = temp(2,1);% 右减中 u_{j+1/2} - ubar
    uhat =temp(3,1);% 中减左 ubar - u_{j-1/2}
    
    % 求左右单元的uh参数
    % 使用周期边界条件
    if i==1
        l = N;
        r = i+1;
    elseif i==N
        l = i-1;
        r = 1;
    else
        l=i-1;
        r=i+1;
    end

    C_R = uh((r-1)*space_order+1:(r)*space_order,1);
    C_L = uh((l-1)*space_order+1:(l)*space_order,1);
    
    % 求均值的右差分和左差分，即输入minmod的第二个和第三个参数
    DeltaR = G(1,:)*(C_R-C);
    DeltaL = G(1,:)*(C-C_L);
    
    % 装载成向量
    X1 = [utilde,DeltaR,DeltaL];
    X2 = [uhat,DeltaR,DeltaL];
    
    M1 = max((2/3)*M,(2/9)*(3+10*M)*M*(h^2/(h^2+abs(DeltaL)+abs(DeltaR))));
    utildemod = mod_minmod(X1,M1,h);
    uhatmod = mod_minmod(X2,M1,h);
    
    %解方程即可
    b = [ubar;utildemod;uhatmod];
    uh_mod(i,:) = G\b;
    

end

uh_mod = reshape(uh_mod',N*space_order,1);

end

