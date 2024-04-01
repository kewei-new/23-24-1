function YH = ENO(Y,N)

YH = zeros(1,N + 1); % 重构的解

YL = [Y(end),Y(1:end - 1)]; % 左边
YR = [Y(2:end),Y(1)]; % 右边
YL2 = [Y(end - 1:end),Y(1:end - 2)]; % 左边2个
YR2 = [Y(3:end),Y(1:2)]; % 右边2个

for i = 1:N
    % 判断添加哪个节点
    if abs(Y(i) - YL(i)) < abs(YR(i) - Y(i))
        % 模板为 i - 1,i
        if abs(Y(i) - 2*YL(i) + YL2(i)) < abs(YR(i) - 2*Y(i) + YL(i))
            % 模板为 i - 2,i - 1,i
            YH(i + 1) = (1/3)*YL2(i) - (7/6)*YL(i) + (11/6)*Y(i);
        else
            % 模板为 i - 1,i,i + 1
            YH(i + 1) = -(1/6)*YL(i) + (5/6)*Y(i) + (1/3)*YR(i);
        end
    else
        % 模板为 i,i + 1
        if abs(YR(i) - 2*Y(i) + YL(i)) < abs(YR2(i) - 2*YR(i) + Y(i))
            % 模板为 i - 1,i,i + 1
            YH(i + 1) = -(1/6)*YL(i) + (5/6)*Y(i) + (1/3)*YR(i);
        else
            % 模板为 i,i + 1,i + 2
            YH(i + 1) = (1/3)*Y(i) + (5/6)*YR(i) - (1/6)*YR2(i);
        end
    end
end