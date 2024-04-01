function result = solution_error_abs_max(exact_function_name,h_N,uh,space_order,T_partion)

result = 0;
temp=[];

for i = 1:h_N

   mesh_point = T_partion(1,[i,i+1]);
   estimated_value = 0;
   k = (i-1)*space_order;
   uh_local=uh(k+1:k+space_order,:);

   for j = 1:space_order

       estimated_value = estimated_value + uh_local(j)*local_basis(mesh_point,mesh_point(1),j-1,0);

   end
    
    result = abs(feval(exact_function_name,mesh_point(1),2*pi)-estimated_value);
    temp=[temp,result];
end

result = max(temp);