function result = Gauss_int_error_1D(real_fun,uh_local,Norm_order,mesh_point,space_order,Gauss_local_coefficient,Gauss_local_point)

result = 0;
Gpn = length(Gauss_local_coefficient);

for i = 1:Gpn

    estimated_value = 0;

   for j = 1:space_order

       estimated_value = estimated_value + uh_local(j)*local_basis(mesh_point,Gauss_local_point(i),j-1,0);

   end
   result = result + Gauss_local_coefficient(i)*abs(feval(real_fun,Gauss_local_point(i),2*pi)-estimated_value)^Norm_order;
end

