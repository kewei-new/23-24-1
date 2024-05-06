%% 各种参数
left = 0;
right = 2*pi;
T_last = 2*pi;
h_N = 20;
CFL = 0.1;
M = 1;
space_order = 2;
partion_type = 101;


%% 求误差阶数
[Gauss_reference_coefficient,Gauss_reference_point]=generate_Gauss_reference(4);
for i = 1:4
   
    N = h_N*2^(i-1);
    [uh,T_partion] = solve_1D_hw1(left,right,T_last,N,CFL,M,space_order,partion_type);
    L1_error = solution_error_Lnorm('exact_function',1,N,uh,space_order,T_partion,Gauss_reference_coefficient,Gauss_reference_point);
    L2_error = solution_error_Lnorm('exact_function',2,N,uh,space_order,T_partion,Gauss_reference_coefficient,Gauss_reference_point);
    L_infinite_error = solution_error_abs_max('exact_function',N,uh,space_order,T_partion);
    h=max(diff(T_partion));
    if i>1
        L1_order = log(L1_error_last/L1_error)/log(h_last/h);
        L2_order = log(L2_error_last/L2_error)/log(h_last/h);
        L_infinite_order = log(L_infinite_error_last/L_infinite_error)/log(h_last/h);
    
        L1_order_list = [L1_order_list;L1_order];
        L2_order_list = [L2_order_list;L2_order];
        L_infinite_order_list = [L_infinite_order_list;L_infinite_order];
        L1_error_list = [L1_error_list;L1_error];
        L2_error_list = [L2_error_list;L2_error];
        L_infinite_error_list = [L_infinite_error_list;L_infinite_error];

        h_last = h;
    else
        L1_order_list = 0;
        L2_order_list = 0;
        L_infinite_order_list = 0;
        L1_error_list = L1_error;
        L2_error_list = L2_error;
        L_infinite_error_list = L_infinite_error;
        h_last = h;

    end

    L1_error_last = L1_error;
    L2_error_last = L2_error;
    L_infinite_error_last = L_infinite_error;
    
    hold on
    plot_DG_1D_time_T(T_partion,uh,space_order,N)

end

hold on
j=1;
T_temp=[];
for i = 1:length(T_partion)

    if mod(j-1,4)==0
    T_temp = [T_temp,T_partion(j)];
    end

    j=j+1;
end
plot(T_partion,feval("exact_function",T_partion,T_last),'-')
% legend('N=20','N=40','N=80','N=160','真实解')