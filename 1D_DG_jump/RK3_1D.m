function Q = RK3_1D(function_name,h_N,CFL,T_last,space_order,T_partion,matrix_E,M,Gauss_reference_coefficient,Gauss_reference_point)

% 求初值，需要使用L2投影(uh0,v) = (u,v)
% 由于uh0可以表示成C1*P1 + C2*P2 + C3*P3，同时由于使用了标准正交基，所以可以得到Ci = (u,Pi)

Q=[];
for i = 1:h_N

    mesh_point = T_partion(1,[i,i+1]);
    [Gauss_local_coefficient,Gauss_local_point]=generate_local_Guass_1D(mesh_point,Gauss_reference_coefficient,Gauss_reference_point);
    
    for k = 1:space_order
        
        Q1 = 0;
        Q2 = 0;
        % 求(u,Pj)
        for j = 1:length(Gauss_local_coefficient)
    
            Q2 = Q2 + Gauss_local_coefficient(j)*local_basis(mesh_point,Gauss_local_point(j),k-1,0)*local_basis(mesh_point,Gauss_local_point(j),k-1,0);
            Q1 = Q1 + Gauss_local_coefficient(j)*local_basis(mesh_point,Gauss_local_point(j),k-1,0)*feval(function_name,Gauss_local_point(j));
    
        end

        Q=[Q;Q1/Q2];

    end

end

h = max(diff(T_partion));
T = 0;

% 求limiter的系数矩阵,如果不是均匀剖分需要注释掉
G = generate_limiter_matrix(T_partion(1,[1,2]),space_order,Gauss_reference_coefficient,Gauss_reference_point);

while T < T_last

    t = CFL*h^((space_order + 1)/3);
    if T + t >= T_last
        t = T_last - T;
    end
    T = T+t;
    

    % RK3
    K1 = Q + t*matrix_E*Q;
    % K1 = Limiter1(K1,h,M,G,space_order);

    K2 = 3/4*Q + 1/4*K1 + 1/4*t*matrix_E*K1;
    % K2 = Limiter1(K2,h,M,G,space_order);

    Q = 1/3*Q + 2/3*K2 + 2/3*t*matrix_E*K2;
    Q = Limiter1(Q,h,M,G,space_order);

end

