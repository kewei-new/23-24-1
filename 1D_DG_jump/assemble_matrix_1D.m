function result = assemble_matrix_1D(N,T_partion,space_order,derivative_order,Gauss_reference_coefficient,Gauss_reference_point)

result = sparse(N*space_order);

if derivative_order == 1

    for i = 1:N
        
        mesh_point = T_partion(1,[i,i+1]);
        [Gauss_local_coefficient,Gauss_local_point]=generate_local_Guass_1D(mesh_point,Gauss_reference_coefficient,Gauss_reference_point);

        for alpha = 1:space_order % trial

            for beta = 1:space_order % test

                result((i-1)*space_order+beta,(i-1)*space_order+alpha)=Gauss_int_1D(mesh_point,alpha-1,beta-1,Gauss_local_coefficient,Gauss_local_point,0,0);
            end
        end
    end

elseif derivative_order == 0

    mesh_point_0 = T_partion(1,[N,N+1]);
    for i = 1:N
        
        mesh_point = T_partion(1,[i,i+1]);
        [Gauss_local_coefficient,Gauss_local_point]=generate_local_Guass_1D(mesh_point,Gauss_reference_coefficient,Gauss_reference_point);

        for beta = 1:space_order % test

            for alpha = 1:space_order % trial

                result((i-1)*space_order+beta,(i-1)*space_order+alpha)=-Gauss_int_1D(mesh_point,alpha-1,beta-1,Gauss_local_coefficient,Gauss_local_point,0,1)+...
                                                                        local_basis(mesh_point,mesh_point(2),alpha-1,0)*local_basis(mesh_point,mesh_point(2),beta-1,0);
                if i == 1%周期边界
                    j=N;
                else
                    j=i-1;
                end

                result((i-1)*space_order+beta,(j-1)*space_order+alpha)=-local_basis(mesh_point_0,mesh_point_0(2),alpha-1,0)*local_basis(mesh_point,mesh_point(1),beta-1,0);
                

            end
        end
        mesh_point_0 = mesh_point;
    end

else
    result =  "derivate_order error!"
end

end