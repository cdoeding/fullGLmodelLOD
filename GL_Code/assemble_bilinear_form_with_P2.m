function S = assemble_bilinear_form_with_P2(A,kappa,T,T_P2,Nd,nodes2mesh_x,nodes2mesh_y)
%ASSEMBLEGLOBALBILINEARFORM Summary of this function goes here
%   Detailed explanation goes here

Nx = sum(logical(nodes2mesh_x));

s_i = zeros(9*size(T,1),1);
s_j = zeros(9*size(T,1),1);
s_val = zeros(9*size(T,1),1);
ind = 1;

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
    z2 = Nd(T(k,2),:); %coordinates of 2nd triangle node
    z3 = Nd(T(k,3),:); %coordinates of 3rd triangle node

    %% transformation to refenrence triangle
    BT = [z2(1)-z1(1), z3(1)-z1(1); ...
        z2(2)-z1(2), z3(2)-z1(2)];

    detBT = BT(1,1)*BT(2,2)-BT(1,2)*BT(2,1);

    BTinv = inv(BT)';
    grad_in_BT = pagemtimes(BTinv,grad_in_quad);

    %% construct A on triangle
    tri_P2 = T_P2(k,:);
    A_in_quad = zeros(2,no_of_quad_points);

    for i = 1:no_of_basis_P2
        index = nodes2mesh_x(tri_P2(i));
        if index ~= 0
            A_in_quad(1,:) = A_in_quad(1,:) + A(index)*phi_P2_in_quad(:,i)';
        end

        index = nodes2mesh_y(tri_P2(i));
        if index ~= 0
            A_in_quad(2,:) = A_in_quad(2,:) + A(Nx + index)*phi_P2_in_quad(:,i)';
        end
    end

    %% assemble
    for i = 1:3
        for j = i:3
            int = diag((1i/kappa * grad_in_BT(:,:,i)+ phi_in_quad(:,i)'.*A_in_quad)'*(1i/kappa * grad_in_BT(:,:,j)+ phi_in_quad(:,j)'.*A_in_quad));
            e = detBT*evaluateQuadrature(transpose(int),w);
            if i == j
                s_i(ind) = tri(i);
                s_j(ind) = tri(i);
                s_val(ind) = e;
                ind = ind + 1;
            else
                s_i(ind) = tri(i);
                s_j(ind) = tri(j);
                s_val(ind) = conj(e);
                ind = ind + 1;

                s_i(ind) = tri(j);
                s_j(ind) = tri(i);
                s_val(ind) = e;
                ind = ind + 1;
            end
        end
    end
end

S = sparse(s_i,s_j,s_val);

end

