function result = minmod(F_pre)
% 最小值限制子，符号相同取绝对值最小，符号不同取0

temp = sum(sign(F_pre));
s = sign(temp);

if abs(temp)/length(F_pre)==1
    
    result = s*min(abs(F_pre));

else

    result = 0;

end