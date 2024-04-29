function uh_mod = Limiter(uh,P,M,space_order,Gauss_reference_coefficient,Gauss_reference_point)
% M：关于二阶导数的上界，一般取M=max(abs(u_0,xx))
% h：剖分步长，必须是均匀剖分

h = P(2)-P(1);
N = length(P)-1;
uh_mod = zeros(N,space_order);

for i = 1:N
    
    mesh_point = P(1,[i,i+1]);
    G = generate_limiter_matrix(mesh_point,space_order,Gauss_reference_coefficient,Gauss_reference_point);
    % 当前单元的uh系数
    C = uh((i-1)*space_order+1:(i)*space_order,1);
    temp = G*C;
    
    % 求带入minmod的第一个参数
    ubar = (1/h)*temp(1,1);%当前单元的均值
    utilde = temp(2,1) - ubar;
    uhat = ubar - temp(2,1);
    
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
    
    M = max((2/3)*M,(2/9)*(3+10*M)*M*(h^2/(h^2+abs(DeltaL)+abs(DeltaR))));
    utildemod = ubar + mod_minmod(X1,M,h);
    uhatmod = ubar - mod_minmod(X2,M,h);
    
    %解方程即可
    b = [ubar;utildemod;uhatmod];
    uh_mod(i,:) = G\b;
end

uh_mod = reshape(uh_mod',N*space_order,1);

end

