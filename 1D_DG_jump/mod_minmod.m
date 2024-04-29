function result=mod_minmod(F_pre,M,h)
% 用来处理存在critical points的限制子，在最小值限制子的基础上进行了一点改良

if abs(F_pre(1))<=M*h^2

    result = F_pre(1);

else

    result = minmod(F_pre);


end