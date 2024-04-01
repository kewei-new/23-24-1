%% 各种参数
left = 0;
right = 2*pi;
T_last = 100*pi;
h_N = 40;
CFL = 0.1;
space_order = 2;
partion_type = 101;

%% 求误差阶数
[Gauss_reference_coefficient,Gauss_reference_point]=generate_Gauss_reference(4);
for i = 1:1
   
    N = h_N*2^(i-1);
    [uh,T_partion] = solve_1D_hw1(left,right,T_last,N,CFL,space_order,partion_type);
    L1_error = solution_error_Lnorm('exact_function',1,N,uh,space_order,T_partion,Gauss_reference_coefficient,Gauss_reference_point);
    L2_error = solution_error_Lnorm('exact_function',2,N,uh,space_order,T_partion,Gauss_reference_coefficient,Gauss_reference_point);

    if i>1
        L1_order = log(L1_error_last/L1_error)/log(2);
        L2_order = log(L2_error_last/L2_error)/log(2);
    
        L1_order_list = [L1_order_list;L1_order];
        L2_order_list = [L2_order_list;L2_order];
    else
        L1_order_list = 0;
        L2_order_list = 0;

    end

    L1_error_last = L1_error;
    L2_error_last = L2_error;

end

L1_order_list
L2_order_list
plot_DG_1D_time_T(T_partion,uh,space_order,N)

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
