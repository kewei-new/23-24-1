function result = generate_limiter_matrix(mesh_point,space_order,Gauss_reference_coefficient,Gauss_reference_point)

[Gauss_coefficient_local_1D,Gauss_point_local_1D]=generate_local_Guass_1D(mesh_point,Gauss_reference_coefficient,Gauss_reference_point);

% 矩阵大小是3*space_order，3是因为要求取均值、右端点值，左端点值这三个参数，space_order是多项式空间的多项式个数，即需要space_order个多项式的结果相加
result = sparse(3,space_order);
% 空间步长，用于求均值的时候使用
h = diff(mesh_point);

for k=1:space_order
    
    % 均值就是在区间I上求积分并除以
    result(1,k) = 1/h*Gauss_int_1D(mesh_point,k-1,0,Gauss_coefficient_local_1D,Gauss_point_local_1D,0,0);
    % 求右端点值
    result(2,k) = local_basis(mesh_point,mesh_point(2),k-1,0) - result(1,k);
    % 求左端点值
    result(3,k) = result(1,k) - local_basis(mesh_point,mesh_point(1),k-1,0);

end



