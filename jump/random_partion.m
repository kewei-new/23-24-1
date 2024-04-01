function result = random_partion(T_partion_base,N)
% 先生成一系列均值为0，方差为1的随机数
% 再用最大值将其限制在（-1，1）之间，然后除以三限制在（-1/2,1/2）之间
% 然后再除以剖分数

num_random = randn(1,N-1);
num_random = num_random/(2*N*max(abs(num_random)));

result = T_partion_base;

for i = 2:N
    
    result(1,i) = result(1,i) + num_random(1,i-1);

end
