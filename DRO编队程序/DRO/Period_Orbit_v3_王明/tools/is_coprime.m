function [flag, gcd] = is_coprime(x , y)
% 判断两正整数是否互质
% 输出
% flag: 是否互质； gcd : x 和 y 的最大公约数

if x == y 
    gcd = x;
    flag = true;
elseif x<y
    temp = x;
    x = y;
    y = temp;
end
% 辗转相除发求最大公约数
while (true)
    z = mod(x,y);
    if z ==0
        break;
    else
        x = y;
        y = z;
    end
end
if y==1
    flag = true;
    gcd = 1;
else
    flag = false;
    gcd = y;
end

