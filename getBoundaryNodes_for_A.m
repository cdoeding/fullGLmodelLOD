function [Bx,By] = getBoundaryNodes_for_A(p,x_a,x_b,condition)
%GETBOUNDARYNODES Summary of this function goes here
%   Detailed explanation goes here

N = size(p,1);
Bx = ones(N,1);
By = ones(N,1);

if strcmp(condition,'non-natural')
    for i = 1:N
        if p(i,1) == x_a || p(i,1) == x_b
            Bx(i) = 0;
        end
        
        if p(i,2) == x_a || p(i,2) == x_b
            By(i) = 0;
        end
    end
end

Bx = spdiags(Bx,0,N,N);
By = spdiags(By,0,N,N);

end

