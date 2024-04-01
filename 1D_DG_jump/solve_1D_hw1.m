function [result_T,T_partion] = solve_1D_hw1(left,right,T_last,h_N,CFL,space_order,partion_type)
% left,right 左右边界
% T_last最大时间
% h_N空间剖分数
% t_N时间剖分数

[Gauss_reference_coefficient,Gauss_reference_point]=generate_Gauss_reference(4);

%% 剖分

T_partion_base = linspace(left,right,h_N+1);

if partion_type == 101
    T_partion = T_partion_base;

elseif partion_type == 102
    T_partion = random_partion(T_partion_base,h_N);

else
    T_partion = "partion_type error!"
end

%% 组装系数矩阵
% A1是C导数的系数矩阵，A2是C的稀疏矩阵,A1*dC + A2*C = 0
A1 = assemble_matrix_1D(h_N,T_partion,space_order,1,Gauss_reference_coefficient,Gauss_reference_point);
A2 = assemble_matrix_1D(h_N,T_partion,space_order,0,Gauss_reference_coefficient,Gauss_reference_point);

%% 3阶龙格库塔
matrix_E = -A1\A2;
result = RK3_1D('g_function',h_N,CFL,T_last,space_order,T_partion,matrix_E,Gauss_reference_coefficient,Gauss_reference_point);
result_T = result(:,end);%输出最后时刻

end




