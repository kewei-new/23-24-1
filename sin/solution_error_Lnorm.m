function result=solution_error_Lnorm(real_fun,Norm_order,h_N,uh,space_order,T_partion,Gauss_reference_coefficient,Gauss_reference_point)

result = 0;

for i = 1:h_N
    
    mesh_point = T_partion(1,[i,i+1]);
    [Gauss_local_coefficient,Gauss_local_point]=generate_local_Guass_1D(mesh_point,Gauss_reference_coefficient,Gauss_reference_point);
    j = (i-1)*space_order;
    uh_local=uh(j+1:j+space_order,:);

    result = result + Gauss_int_error_1D(real_fun,uh_local,Norm_order,mesh_point,space_order,Gauss_local_coefficient,Gauss_local_point);

end

result = result^(1/Norm_order);