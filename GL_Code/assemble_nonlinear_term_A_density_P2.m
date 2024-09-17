function F = assemble_nonlinear_term_A_density_P2(u,T,T_P1,Nd,nodes2mesh,nodes2mesh_x,nodes2mesh_y)
%ASSEMBLE_NONLINEAR_TERM_A_DENSITY_P2 Summary of this function goes here
%   Detailed explanation goes here

Nx = sum(logical(nodes2mesh_x));

F_i = zeros(2*36*size(T,1),1);
F_j = zeros(2*36*size(T,1),1);
FF = zeros(2*36*size(T,1),1);
ind = 1;

phi = { @(x) -x(1) - x(2) + 1;
    @(x) x(1);
    @(x) x(2)};

phi_P2 = {@(x) (1 - x(1) - x(2))*(1-2*x(1)-2*x(2));
    @(x) 4*x(1)*(1-x(1)-x(2));
    @(x) x(1)*(2*x(1)-1);
    @(x) 4*x(1)*x(2);
    @(x) x(2)*(2*x(2)-1);
    @(x) 4*x(2)*(1-x(1)-x(2))};

[quad,w] = getQuadrature(7);
no_of_basis = size(phi,1);
no_of_basis_P2 = size(phi_P2,1);
no_of_quad_points = length(w);

phi_in_quad = zeros(no_of_quad_points,no_of_basis);
phi_P2_in_quad = zeros(no_of_quad_points,no_of_basis_P2);

for j = 1:no_of_quad_points
    for i = 1:no_of_basis
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
    
    %% construct f = |u|^2 on triangle
    tri_P1 = T_P1(k,:);
    f_in_quad = zeros(no_of_quad_points,1);

    for j = 1:no_of_basis
        index = nodes2mesh(tri_P1(j));
        if index ~= 0
            f_in_quad = f_in_quad + u(index)*phi_in_quad(:,j);
        end
    end

    f_in_quad = abs(f_in_quad).^2;
    
    %% assemble
    for i = 1:6
        for j = 1:6
            index_i = nodes2mesh_x(tri(i));
            index_j = nodes2mesh_x(tri(j));
            if index_i ~= 0 && index_j ~= 0
                e = detBT*evaluateQuadrature(transpose(f_in_quad.*phi_P2_in_quad(:,i).*phi_P2_in_quad(:,j)),w);
                F_i(ind) = index_i;
                F_j(ind) = index_j;
                FF(ind) = e;
                ind = ind + 1;
            end
            
            index_i = nodes2mesh_y(tri(i));
            index_j = nodes2mesh_y(tri(j));
            if index_i ~= 0 && index_j ~= 0
                e = detBT*evaluateQuadrature(transpose(f_in_quad.*phi_P2_in_quad(:,i).*phi_P2_in_quad(:,j)),w);
                F_i(ind) = Nx + index_i;
                F_j(ind) = Nx + index_j;
                FF(ind) = e;
                ind = ind + 1;
            end
        end
    end
end

del_index = find(F_i==0,1);
if isempty(del_index)
    del_index = length(F_i)+1;
end
F = sparse(F_i(1:del_index-1),F_j(1:del_index-1),FF(1:del_index-1));
end

