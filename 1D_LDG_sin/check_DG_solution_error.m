%% 各种参数
left = 0;
right = 2*pi;
T_last = 2*pi; % 时间最后时刻，最好不要改，
               % 因为算误差的函数里我固定为了2pi,需要自己去里面找到相应参数改下
h_N = 20;
space_order = 3;% 实际空间维数是space_order-1，即==2时是一次多项式，==3时是二次多项式
partion_type = 101; %均匀剖分101，随机剖分102

%% 求误差阶数
if space_order==3
    CFL = 0.015;
elseif space_order==2
    CFL = 0.06;
end

[Gauss_reference_coefficient,Gauss_reference_point]=generate_Gauss_reference(4);
% i>3会很慢
for i = 1:3
   
    N = h_N*2^(i-1);
    % uh是数值解，ph是辅助变量
    [uh,ph,T_partion] = solve_1D_hw2(left,right,T_last,N,CFL,space_order,partion_type);
    % 上面输出的是所有时刻的计算结果，取最后时刻进行检查
    uh =uh(:,end);
    ph =ph(:,end);
    % 下面求的是辅助变量的阶，如果要求uh的阶
    % 把'exact_der_1_function'换成'exact_function'，把ph换成uh
    % L1_order_list\L2_order_list\L_infinite_order_list 存储着阶
    L1_error = solution_error_Lnorm('exact_function',T_last,1,N,uh,space_order,T_partion,Gauss_reference_coefficient,Gauss_reference_point);
    L2_error = solution_error_Lnorm('exact_function',T_last,2,N,uh,space_order,T_partion,Gauss_reference_coefficient,Gauss_reference_point);
    L_infinite_error = solution_error_abs_max('exact_function',T_last,N,uh,space_order,T_partion);
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
    
    % plot_DG_1D_time_T(T_partion,uh,space_order,N,i)
end

% plot_DG_1D_time_T(T_partion,uh,space_order,N)
% hold on
% j=1;
% T_temp=[];
% for i = 1:length(T_partion)
% 
%     if mod(j-1,4)==0
%     T_temp = [T_temp,T_partion(j)];
%     end
% 
%     j=j+1;
% end
% plot(T_partion,feval("exact_function",T_partion,T_last),'-o')
% legend('N=20','N=40','N=80','N=160','真实解')
% legend('K=1','K=2','真实解')