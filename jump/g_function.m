function result = g_function(x)
% 初值函数

result = zeros(1,length(x));

for i = 1:length(x)

    if x(i)>=pi/2 && x(i)<=3*pi/2
    
        result(1,i) = 1;

    else

        result(1,i) = 0;

    end

end