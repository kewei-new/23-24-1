function u = min_limiter(u_pre)

[Gauss_coefficient_reference_1D,Gauss_point_reference_1D]=generate_Gauss_reference(4);

for i = 1:N

    if (i+1)>N
        next_mesh_point = T_partion(1,[1,2]);
    elseif i == 1
    % 初始pre的数据
        h_local = diff(T_partion(1,[N,N+1]);
        j = (N-1)*space_order;
        uh_local=u_pre(j+1:j+space_order,:);
        [Gauss_local_coefficient,Gauss_local_point]=generate_local_Guass_1D(T_partion(1,[N,N+1]),Gauss_coefficient_reference_1D,Gauss_point_reference_1D);    
        u_average_pre = 0;

        for k = 1:space_order   
            u_average_pre = u_average_pre + uh_local(k)*Gauss_int_1D(T_partion(1,[N,N+1]),k-1,0,Gauss_local_coefficient,Gauss_local_point,0,0)/h_local;
        end
    % 初始present的数据
        h_local = diff(T_partion(1,[i,i+1]);
        uh_local=u_pre(1:space_order,:);
        [Gauss_local_coefficient,Gauss_local_point]=generate_local_Guass_1D(T_partion(1,[i,i+1]),Gauss_coefficient_reference_1D,Gauss_point_reference_1D);    
        u_average = 0;   
        estimated_value_left = 0;
        estimated_value_right = 0;
        for k = 1:space_order   
            u_average_pre = u_average_pre + uh_local(k)*Gauss_int_1D(T_partion(1,[1,2]),k-1,0,Gauss_local_coefficient,Gauss_local_point,0,0)/h_local;
            estimated_value_left = estimated_value_left + uh_local(k)*local_basis(T_partion(1,[1,2]),mesh_point(1),k-1,0);
            estimated_value_right = estimated_value_right + uh_local(k)*local_basis(T_partion(1,[1,2]),mesh_point(2),k-1,0);
        end

        u_wan_1 = estimated_value_right - u_average;
        u_wan_2 = u_average - estimated_value_left;
    % 初始next网格信息
        next_mesh_point = T(1,[i+1,i+2]);
    else
        next_mesh_point = T(1,[i+1,i+2]);
    end
    
    h_local = diff(next_mesh_point);
    j = (i)*space_order;
    uh_local=u_pre(j+1:j+space_order,:);
    
    [Gauss_local_coefficient,Gauss_local_point]=generate_local_Guass_1D(next_mesh_point,Gauss_coefficient_reference_1D,Gauss_point_reference_1D);
    
    u_average_next = 0;
    estimated_value_left_next = 0;
    estimated_value_right_next = 0;
    
    for k = 1:space_order   
        u_average_next = u_average_next + uh_local(k)*Gauss_int_1D(next_mesh_point,k-1,0,Gauss_local_coefficient,Gauss_local_point,0,0)/h_local;
        estimated_value_left_next = estimated_value_left_next + uh_local(k)*local_basis(mesh_point_next,mesh_point(1),k-1,0);
        estimated_value_right_next = estimated_value_right_next + uh_local(k)*local_basis(mesh_point_next,mesh_point(2),k-1,0);
    end
    
    u_wan_1_next = estimated_value_right_next - u_average_next;
    u_wan_2_next = u_average_nex - estimated_value_left_next;
    
    u_wan_1_mod = min_mod(u_wan_1,u_average_next-u_average,u_average-u_average_pre);
    u_wan_2_mod = min_mod(u_wan_2,u_average_next-u_average,u_average-u_average_pre);
    
    % 需要更新系数
    
    u_average_pre = u_average;
    u_wan_1 = u_wan_1_next;
    u_wan_2 = u_wan_2_next;
    estimated_value_right = estimated_value_right_next;
    estimated_value_left = estimated_value_left_next;


end