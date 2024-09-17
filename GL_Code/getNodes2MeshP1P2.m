function [bn,bn_x,bn_y] = getNodes2MeshP1P2(Nd,Nd_P2,x_a,x_b,boundary_u,boundary_A)
%GETNODES2MESHP1P2 Summary of this function goes here
%   Detailed explanation goes here

N = size(Nd,1);
N_P2 = size(Nd_P2,1);

bn = zeros(N,1); % indicator if node is boundary node
bn_x = zeros(N_P2,1); % indicator if node is boundary node
bn_y = zeros(N_P2,1); % indicator if node is boundary node

% generate nodes and boundary indicator
%% u-equation P1
counter_mesh = 1;
for i = 1:N
    if (Nd(i,1) ~= x_a && Nd(i,2) ~= x_a && Nd(i,1) ~= x_b && Nd(i,2) ~= x_b) || strcmp(boundary_u,'Neumann')
        bn(i) = counter_mesh;
        counter_mesh = counter_mesh + 1;
    end
end

%% A-equation P2
counter_mesh_x = 1;
counter_mesh_y = 1;
for i = 1:N_P2
    % is node boundary node?
    if (Nd_P2(i,1) ~= x_a && Nd_P2(i,1) ~= x_b) || strcmp(boundary_A,'natural')
        bn_x(i) = counter_mesh_x;
        counter_mesh_x = counter_mesh_x + 1;
    end
    
    if (Nd_P2(i,2) ~= x_a && Nd_P2(i,2) ~= x_b) || strcmp(boundary_A,'natural')
        bn_y(i) = counter_mesh_y;
        counter_mesh_y = counter_mesh_y + 1;
    end
end

end

