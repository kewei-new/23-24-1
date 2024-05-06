function result = local_basis(mesh_point,x,fun_order,der_order)
% 这里用的是勒让德多项式
% 和一般勒让德的区别是，Pn和Pn的积分结果为1，具体操作是在Pn上乘了系数

a = max(mesh_point);
b = min(mesh_point);
h = (a-b)/2;
center = (a+b)/2;
t = (x-center)/h;

if fun_order == 0

    if der_order == 0
        
        result = 1;

    elseif der_order == 1

        result = 0;

    end
    
elseif fun_order == 1
    
    if der_order == 0
        
        result = t;

    elseif der_order == 1

        result = 1;
        
    end

elseif fun_order == 2
   
    if der_order == 0
        
        result = 3/2*(t^2 - 1/3);

    elseif der_order == 1

        result = 3*t;
        
    end

else

    result = "fun_order Error!"

end

if der_order == 0

    result = 2*result;

elseif der_order == 1

    result = 2/h*result;

end


