function F = assemble_RHS_curlH_P2(H,T,Nd,nodes2mesh_x,nodes2mesh_y)
%ASSEMBLE_RHS_CURLH_P2 Summary of this function goes here
%   Detailed explanation goes here

Nx = sum(logical(nodes2mesh_x));
Ny = sum(logical(nodes2mesh_y));
F = zeros(Nx + Ny,1);

grad = { @(x) [4*x(1) + 4*x(2) - 3; 4*x(1) + 4*x(2) - 3];
    @(x) [-4*(2*x(1) + x(2) - 1); -4*x(1)];
    @(x) [4*x(1) - 1; 0];
    @(x) [4*x(2); 4*x(1)];
    @(x) [0; 4*x(2) - 1];
    @(x) [-4*x(2); -4*(x(1) + 2*x(2) - 1)]};

[quad,w] = getQuadrature(3);
no_of_basis = size(grad,1);
no_of_quad_points = length(w);
grad_in_quad = zeros(2,no_of_quad_points,no_of_basis);
for i = 1:no_of_basis
    for j = 1:no_of_quad_points
        grad_in_quad(:,j,i) = grad{i}(quad(:,j));
    end
end

H_in_quad = zeros(1,no_of_quad_points);

for k = 1:size(T,1)
    tri = T(k,:); %node index of triangle
    z1 = Nd(T(k,1),:); %coordinates of 1st triangle node
    z2 = Nd(T(k,3),:); %coordinates of 2nd triangle node
    z3 = Nd(T(k,5),:); %coordinates of 3rd triangle node
    
    %% transformation to refenrence triangle
    BT = [z2(1)-z1(1), z3(1)-z1(1); ...
        z2(2)-z1(2), z3(2)-z1(2)];
    
    detBT = BT(1,1)*BT(2,2)-BT(1,2)*BT(2,1);
    b = [z1(1); z1(2)];
    
    BTinv = inv(BT)';
    grad_in_BT = pagemtimes(BTinv,grad_in_quad);

    for i = 1:no_of_quad_points
        H_in_quad(i) = H(BT*quad(:,i) + b);
    end

    for i = 1:6
        index_i = nodes2mesh_x(tri(i));
        if index_i ~= 0
            e = -detBT*evaluateQuadrature(H_in_quad.*grad_in_BT(2,:,i),w);
            F(index_i) = F(index_i) + e;
        end
        
        index_i = nodes2mesh_y(tri(i));
        if index_i ~= 0
            e = detBT*evaluateQuadrature(H_in_quad.*grad_in_BT(1,:,i),w);
            F(Nx + index_i) = F(Nx + index_i) + e;
        end
    end
end