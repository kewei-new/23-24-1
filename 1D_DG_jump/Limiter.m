function uh_mod = Limiter(uh,h,M,G,space_order)
% M：关于二阶导数的上界，一般取M=max(abs(u_0,xx))，不知道取什么
% h：剖分步长，必须是均匀剖分


N = length(uh)/space_order;
uh_mod = zeros(N,space_order);
for i = 1:N
    
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
    
    C = uh((i-1)*space_order+1:(i)*space_order,1);
    deltaUR = C(2) + (2/3)*C(3);
    deltaUL = C(2) - (2/3)*C(3);
    deltaURM = C_R(1) - C(1);
    deltaULM = C(1) - C_L(1);
    
    X1 = [deltaUR,deltaURM,deltaULM];
    X2 = [deltaUL,deltaURM,deltaULM];
    
    M1 = max((2/3)*M,(2/9)*(3+10*M)*M*(h^2/(h^2+abs(deltaURM)+abs(deltaULM))));
    deltaUR1 = mod_minmod(X1,M1,h);
    deltaUL1 = mod_minmod(X2,M1,h);
    
    uh_mod(i,2) = (deltaUR1 + deltaUL1)/2;
    uh_mod(i,3) = 3*(deltaUR1 - deltaUL1)/4;
    uh_mod(i,1) = C(1);
    

end

uh_mod = reshape(uh_mod',N*space_order,1);



end

