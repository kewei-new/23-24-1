function result = min_mod(u_wan,u_add_ave,u_plus_ave)

if sign(u_wan)==sign(u_add_ave)==sign(u_plus_ave)
    result = sign(u_wan)*min([abs(u_wan),abs(u_add_ave),abs(u_plus_ave)]);
else
    result = 0;
end

