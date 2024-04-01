function YH = ENO(Y,N)

YH = zeros(1,N + 1); % �ع��Ľ�

YL = [Y(end),Y(1:end - 1)]; % ���
YR = [Y(2:end),Y(1)]; % �ұ�
YL2 = [Y(end - 1:end),Y(1:end - 2)]; % ���2��
YR2 = [Y(3:end),Y(1:2)]; % �ұ�2��

for i = 1:N
    % �ж�����ĸ��ڵ�
    if abs(Y(i) - YL(i)) < abs(YR(i) - Y(i))
        % ģ��Ϊ i - 1,i
        if abs(Y(i) - 2*YL(i) + YL2(i)) < abs(YR(i) - 2*Y(i) + YL(i))
            % ģ��Ϊ i - 2,i - 1,i
            YH(i + 1) = (1/3)*YL2(i) - (7/6)*YL(i) + (11/6)*Y(i);
        else
            % ģ��Ϊ i - 1,i,i + 1
            YH(i + 1) = -(1/6)*YL(i) + (5/6)*Y(i) + (1/3)*YR(i);
        end
    else
        % ģ��Ϊ i,i + 1
        if abs(YR(i) - 2*Y(i) + YL(i)) < abs(YR2(i) - 2*YR(i) + Y(i))
            % ģ��Ϊ i - 1,i,i + 1
            YH(i + 1) = -(1/6)*YL(i) + (5/6)*Y(i) + (1/3)*YR(i);
        else
            % ģ��Ϊ i,i + 1,i + 2
            YH(i + 1) = (1/3)*Y(i) + (5/6)*YR(i) - (1/6)*YR2(i);
        end
    end
end