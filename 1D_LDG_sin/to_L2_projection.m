function u0 = to_L2_projection(function_name,h_N,T_partion,space_order,Gauss_reference_coefficient,Gauss_reference_point)

u0=[];
for i = 1:h_N

    mesh_point = T_partion(1,[i,i+1]);
    [Gauss_local_coefficient,Gauss_local_point]=generate_local_Guass_1D(mesh_point,Gauss_reference_coefficient,Gauss_reference_point);
    
    for k = 1:space_order
        
        Q1 = 0;

        % 求(u,Pj)
        for j = 1:length(Gauss_local_coefficient)
    
            Q1 = Q1 + Gauss_local_coefficient(j)*local_basis(mesh_point,Gauss_local_point(j),k-1,0)*feval(function_name,Gauss_local_point(j));
    
        end

        u0=[u0;Q1];

    end
end