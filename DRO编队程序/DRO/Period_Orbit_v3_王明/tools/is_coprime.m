function [flag, gcd] = is_coprime(x , y)
% �ж����������Ƿ���
% ���
% flag: �Ƿ��ʣ� gcd : x �� y �����Լ��

if x == y 
    gcd = x;
    flag = true;
elseif x<y
    temp = x;
    x = y;
    y = temp;
end
% շת����������Լ��
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

