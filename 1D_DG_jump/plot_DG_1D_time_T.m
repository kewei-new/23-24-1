function plot_DG_1D_time_T(T_partion,uh,space_order,N)

for i = 1:N

    mesh_point = T_partion(:,[i,i+1]);
    j = (i-1)*space_order;
    uh_local=uh(j+1:j+space_order,:);

    estimated_value_temp = 0;
    for k = 1:space_order

       estimated_value_temp = estimated_value_temp + uh_local(k)*local_basis(mesh_point,mesh_point(1),k-1,0);

    end

    if i == 1
        estimated_value = estimated_value_temp;
    else
        estimated_value = [estimated_value,estimated_value_temp];
    end


end
% estimated_value_0 = ENO(estimated_value,N);

plot(T_partion(1:N),estimated_value,'-');
% hold on
% plot(T_partion(1:N+1),estimated_value_0,'bo');
