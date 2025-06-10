function vec_I = inertial2synodic(vec_S,t)
vec_S = vec_S(:);
theta = t;
if length(vec_S) == 2
    M = [cos(theta),-sin(theta);
        sin(theta),cos(theta)];
else
    M = [cos(theta),-sin(theta),0;
        sin(theta),cos(theta),0;
        0,0,1];
end
vec_I = M*vec_S;
vec_I = reshape(vec_I,size(vec_S));
