function F = assemble_RHS_u_nabla_u_P2(u,T,T_P1,Nd,nodes2mesh,nodes2mesh_x,nodes2mesh_y)
%ASSEMBLE_RHS_U_NABLA_U_P2 Summary of this function goes here
%   Detailed explanation goes here

Nx = sum(logical(nodes2mesh_x));
Ny = sum(logical(nodes2mesh_y));
F = zeros(Nx + Ny,1);

grad = {[-1; -1]; [1; 0]; [0; 1]};

phi = { @(x) -x(1) - x(2) + 1;
    @(x) x(1);
    @(x) x(2)};

phi_P2 = {@(x) (1 - x(1) - x(2))*(1-2*x(1)-2*x(2));
    @(x) 4*x(1)*(1-x(1)-x(2));
    @(x) x(1)*(2*x(1)-1);
    @(x) 4*x(1)*x(2);
    @(x) x(2)*(2*x(2)-1);
    @(x) 4*x(2)*(1-x(1)-x(2))};

[quad,w] = getQuadrature(3);
no_of_basis = size(grad,1);
no_of_basis_P2 = size(phi_P2,1);
no_of_quad_points = length(w);

grad_in_quad = zeros(2,no_of_quad_points,no_of_basis);
phi_in_quad = zeros(no_of_quad_points,no_of_basis);
phi_P2_in_quad = zeros(no_of_quad_points,no_of_basis_P2);

for j = 1:no_of_quad_points
    for i = 1:no_of_basis
        grad_in_quad(:,j,i) = grad{i};
        phi_in_quad(j,i) = phi{i}(quad(:,j));
    end

    for i = 1:no_of_basis_P2
        phi_P2_in_quad(j,i) = phi_P2{i}(quad(:,j));
    end
end

for k = 1:size(T,1)
    tri = T(k,:); %node index of triangle
    z1 = Nd(T(k,1),:); %coordinates of 1st triangle node
    z2 = Nd(T(k,3),:); %coordinates of 2nd triangle node
    z3 = Nd(T(k,5),:); %coordinates of 3rd triangle node

    %% transformation to refenrence triangle
    BT = [z2(1)-z1(1), z3(1)-z1(1); ...
        z2(2)-z1(2), z3(2)-z1(2)];

    detBT = BT(1,1)*BT(2,2)-BT(1,2)*BT(2,1);

    BTinv = inv(BT)';
    grad_in_BT = pagemtimes(BTinv,grad_in_quad);

    %% construct u on triangle
    tri_P1 = T_P1(k,:);
    u_in_quad = zeros(no_of_quad_points,1);
    gradu_in_quad = zeros(no_of_quad_points,2);

    for j = 1:no_of_basis
        index = nodes2mesh(tri_P1(j));
        if index ~= 0
            u_in_quad = u_in_quad + u(index)*phi_in_quad(:,j);
            gradu_in_quad = gradu_in_quad + u(index)*grad_in_BT(:,:,j)';
        end
    end

    u_in_quad = abs(u_in_quad).^2;
    f1_in_quad = conj(u_in_quad).*gradu_in_quad(:,1) - u_in_quad.*conj(gradu_in_quad(:,1));
    f2_in_quad = conj(u_in_quad).*gradu_in_quad(:,2) - u_in_quad.*conj(gradu_in_quad(:,2));

    %% assemble
    for i = 1:6
        index_i = nodes2mesh_x(tri(i));
        if index_i ~= 0
            F(index_i) = F(index_i) + detBT*evaluateQuadrature(transpose(f1_in_quad.*phi_P2_in_quad(:,i)),w);
        end

        index_i = nodes2mesh_y(tri(i));
        if index_i ~= 0
            F(Nx + index_i) = F(Nx + index_i) + detBT*evaluateQuadrature(transpose(f2_in_quad.*phi_P2_in_quad(:,i)),w);
        end
    end
end

end