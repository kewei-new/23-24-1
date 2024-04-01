function [U,P,T_partion]=solve_1D_hw2(left,right,T_last,h_N,CFL,space_order,partion_type)
% ut + ux - uxx = 0
% u(x,0) = sin(x)
% 周期边界条件,求解区域[0,2*pi]

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

%% L2投影求初值
u0 = to_L2_projection('Initial_function',h_N,T_partion,space_order,Gauss_reference_coefficient,Gauss_reference_point);

%% LDG求解
[U,P,T] = LDG_1D_solver(u0,h_N,T_partion,space_order,T_last,CFL,Gauss_reference_coefficient,Gauss_reference_point);
