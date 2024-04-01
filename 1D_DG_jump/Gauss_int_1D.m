function result = Gauss_int_1D(mesh_point,trial_order,test_order,Gauss_local_coefficient,Gauss_local_point,trial_derivate_order,test_derivate_order)

result = 0;
Gpn = length(Gauss_local_coefficient);

for i = 1:Gpn

    result = result + Gauss_local_coefficient(i)*local_basis(mesh_point,Gauss_local_point(i),trial_order,trial_derivate_order)*local_basis(mesh_point,Gauss_local_point(i),test_order,test_derivate_order);

end