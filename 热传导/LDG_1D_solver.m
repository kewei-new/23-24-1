function [u0,p0,T] = LDG_1D_solver(u0,h_N,T_partion,space_order,T_last,CFL,Gauss_reference_coefficient,Gauss_reference_point)

%% 先求第二个方程的系数 P-Ux=0
% 通过观察发现，第二个方程的系数矩阵和DG连续初值的那个程序的系数矩阵相同，只是取相反的符号
% 因此直接调用assemble_matrix函数

% A1是P的系数矩阵
A1 = assemble_matrix_1D(h_N,T_partion,space_order,1,Gauss_reference_coefficient,Gauss_reference_point);
% A2是u的系数矩阵
A2 = assemble_matrix_1D(h_N,T_partion,space_order,0,Gauss_reference_coefficient,Gauss_reference_point);
E = A1\A2;

%% 然后利用上述RK3求解
% 第一个方程是Ut = Px-Ux，左端系数矩阵是A1，右端项P的系数矩阵是A3，U的系数矩阵是A2
% 需要先求出Ph，即Ph=E*u0;

% 这里由于辅助变量P的数值通量是P+，因此需要写一个向后区间的系数矩阵
A3 = assemble_matrix_1D(h_N,T_partion,space_order,2,Gauss_reference_coefficient,Gauss_reference_point);
h = max(diff(T_partion));

% T是用来存储时间，方便随时暂停检查运算到哪了
T = 0;
% matrix_E是第一个方程的迭代矩阵，因为ut=[-A2,A3]*[u;p]
matrix_E = [-A1\A2,A1\A3];

% 求解T0时刻的ph值
p0=E*u0;

%用来储存各个时刻计算的结果，按列存储
u=u0;
p=p0;
while T(end) < T_last
    
    t = CFL*h^2;
    if T(end) + t >= T_last
        t = T_last - T(end);
    end
    T = [T;T(end) + t];
    
    % RK3求解下一时刻的u
    u0=RK3_1D_without_projection(u0,p0,E,matrix_E,t);
    % 把u保存用来检验数值解阶数
    u=[u,u0];
    % 由第二个方程求辅助变量p
    p0 = E*u0;
    % 把p保存用来检验辅助变量阶数
    p=[p,p0];

end

