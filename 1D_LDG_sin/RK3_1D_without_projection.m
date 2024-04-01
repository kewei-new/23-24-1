function Q=RK3_1D_without_projection(Q,Q1,E,matrix_E,t)
    
    % 主体部分

    K1 = Q + t*matrix_E*[Q;Q1];
    K2 = 3/4*Q + 1/4*K1 + 1/4*t*matrix_E*[K1;E*K1];
    Q = 1/3*Q + 2/3*K2 + 2/3*t*matrix_E*[K2;E*K2];

end
