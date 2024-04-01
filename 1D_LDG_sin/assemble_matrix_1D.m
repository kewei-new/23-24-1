function result = assemble_matrix_1D(N,T_partion,space_order,derivative_order,Gauss_reference_coefficient,Gauss_reference_point)
% 该函数是生成u-ux=0的两个系数矩阵

% 矩阵大小为（剖分区间数*试探函数空间维数）^2
% 矩阵的行和测试函数相关，列和试探函数相关

% derivative_order==1是质量矩阵
% derivative_order==0是flux（u）=u-由ux生成的矩阵
% derivative_order==2是flux（u）=u+由ux生成的矩阵


result = zeros(N*space_order,N*space_order);

if derivative_order == 1

    % 遍历每个剖分区间
    for i = 1:N
        
        % 第i个剖分区间的网格点
        mesh_point = T_partion(1,[i,i+1]);
        % 根据reference_Gauss点和系数生成对应网格偶分区间的点和系数
        [Gauss_local_coefficient,Gauss_local_point]=generate_local_Guass_1D(mesh_point,Gauss_reference_coefficient,Gauss_reference_point);
        
        % 遍历试探函数
        for alpha = 1:space_order % trial
            % 遍历测试函数
            for beta = 1:space_order % test
                % 因为单纲矩阵是（函数空间维数）^2的方阵，所以位置就是对角线上第i个矩阵块
                % 这里求的是质量矩阵，因为basis函数设置是根据次数（0次多项式、1次多项式等），需要alpha-1,beta-1
                % 又因为质量矩阵basis函数不用求导，最后两个参数都是0。
                result((i-1)*space_order+beta,(i-1)*space_order+alpha)=Gauss_int_1D(mesh_point,alpha-1,beta-1,Gauss_local_coefficient,Gauss_local_point,0,0);

            end
        end
    end

elseif derivative_order == 0

    % 第一行的i-1（也就是1-1），根据周期条件是N
    mesh_point_0 = T_partion(1,[N,N+1]);
    for i = 1:N
        
        mesh_point = T_partion(1,[i,i+1]);
        [Gauss_local_coefficient,Gauss_local_point]=generate_local_Guass_1D(mesh_point,Gauss_reference_coefficient,Gauss_reference_point);

        for beta = 1:space_order % test

            for alpha = 1:space_order % trial
               
                % 这里第一个-Gauss_int_1D求的是trial为0阶导数，test为1阶导数的积分
                % 第二部分是间断处右端点的值
                % 位置和第一个矩阵相同是ui的系数
                result((i-1)*space_order+beta,(i-1)*space_order+alpha)=-Gauss_int_1D(mesh_point,alpha-1,beta-1,Gauss_local_coefficient,Gauss_local_point,0,1)+...
                                                                        local_basis(mesh_point,mesh_point(2),alpha-1,0)*local_basis(mesh_point,mesh_point(2),beta-1,0);
                if i == 1%周期边界
                    j=N;
                else
                    j=i-1;
                end
                % 这里求的是前一个网格区间的右端点值，和当前区间左端点值的乘积
                % 由于trial是前一个区间，所以需要i-1列上的矩阵块
                result((i-1)*space_order+beta,(j-1)*space_order+alpha)=-local_basis(mesh_point_0,mesh_point_0(2),alpha-1,0)*local_basis(mesh_point,mesh_point(1),beta-1,0);
                

            end
        end
        mesh_point_0 = mesh_point;
    end

elseif derivative_order==2
    % 第N行的i+1（也就是N+1），根据周期条件是N,实际上是第1个有限元区间
    mesh_point_1 = T_partion(1,[1,2]);
    for i = N:-1:1
        
        mesh_point = T_partion(1,[i,i+1]);
        [Gauss_local_coefficient,Gauss_local_point]=generate_local_Guass_1D(mesh_point,Gauss_reference_coefficient,Gauss_reference_point);

        for beta = 1:space_order % test

            for alpha = 1:space_order % trial
               
                % 这里第一个-Gauss_int_1D求的是trial为0阶导数，test为1阶导数的积分
                % 第二部分是间断处左端点的值
                % 位置和第一个矩阵相同是ui的系数
                result((i-1)*space_order+beta,(i-1)*space_order+alpha)=-Gauss_int_1D(mesh_point,alpha-1,beta-1,Gauss_local_coefficient,Gauss_local_point,0,1)...
                                                                       -local_basis(mesh_point,mesh_point(1),alpha-1,0)*local_basis(mesh_point,mesh_point(1),beta-1,0);
                if i == N%周期边界
                    j=1;
                else
                    j=i+1;
                end
                % 这里求的是后一个网格区间的右端点值，和当前区间左端点值的乘积
                % 由于trial是后一个区间，所以需要i+1列上的矩阵块
                result((i-1)*space_order+beta,(j-1)*space_order+alpha)=local_basis(mesh_point_1,mesh_point_1(1),alpha-1,0)*local_basis(mesh_point,mesh_point(2),beta-1,0);
                

            end
        end
        mesh_point_1 = mesh_point;
    end
   

else
    result =  "derivate_order error!"
end

result(abs(result)<10^(-10))=0;
end