function result = random_partion(T_partion_base,N)

% 生成一串区间在（-1/2，1/2）的随机数
num_random = rand(1,N-1)-0.5;
% 为了防止加在相邻两个网格点上的随机数导致区间长度为非正数，需要把（-1/2，1/2）再放缩到（-h/2，h/2）
% h = (right - left)/N
% 10^(-1)是对随机数再次放缩
num_random = num_random/(10*N);

% 初始化剖分网格点（均匀的）
result = T_partion_base;

% 除了两侧端点值，其余网格点值都加上随机数
% 因此随机数只生成了N+1-2=N-1个
for i = 2:N
    
    result(1,i) = result(1,i) + num_random(1,i-1);

end
