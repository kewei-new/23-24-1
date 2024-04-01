function result = Gauss_int_1D(mesh_point,trial_order,test_order,Gauss_local_coefficient,Gauss_local_point,trial_derivate_order,test_derivate_order)
% 求质量或者刚度矩阵
result = 0;
Gpn = length(Gauss_local_coefficient);

for i = 1:Gpn

    % 第一个basis是试探函数、第二个basis是测试函数
    result = result + Gauss_local_coefficient(i)*local_basis(mesh_point,Gauss_local_point(i),trial_order,trial_derivate_order)*local_basis(mesh_point,Gauss_local_point(i),test_order,test_derivate_order);

end