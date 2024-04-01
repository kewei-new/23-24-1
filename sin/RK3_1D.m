function C = RK3_1D(function_name,h_N,CFL,T_last,space_order,T_partion,matrix_E,Gauss_reference_coefficient,Gauss_reference_point)

% 求初值，需要使用L2投影(uh0,v) = (u,v)
% 由于uh0可以表示成C1*P1 + C2*P2 + C3*P3，同时由于使用了标准正交基，所以可以得到Ci = (u,Pi)

Q=[];
for i = 1:h_N

    mesh_point = T_partion(1,[i,i+1]);
    [Gauss_local_coefficient,Gauss_local_point]=generate_local_Guass_1D(mesh_point,Gauss_reference_coefficient,Gauss_reference_point);
    
    for k = 1:space_order
        
        Q1 = 0;

        % 求(u,Pj)
        for j = 1:length(Gauss_local_coefficient)
    
            Q1 = Q1 + Gauss_local_coefficient(j)*local_basis(mesh_point,Gauss_local_point(j),k-1,0)*feval(function_name,Gauss_local_point(j));
    
        end

        Q=[Q;Q1];

    end

end



C = Q;
h = max(diff(T_partion));
T = 0;
while T(end) < T_last

    t = CFL*h^((space_order + 1)/3);
    if T(end) + t >= T_last
        t = T_last - T(end);
    end
    T = [T;T(end) + t];

    K1 = Q + t*matrix_E*Q;
    K2 = 3/4*Q + 1/4*K1 + 1/4*t*matrix_E*K1;
    Q = 1/3*Q + 2/3*K2 + 2/3*t*matrix_E*K2;
    C = [C,Q];

end

