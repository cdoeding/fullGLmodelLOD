function [S,D,G,GA] = assembleLinearMatrices_An_Dn(A,D,kappa,T,T_P2,Nd,nodes2mesh_x,nodes2mesh_y,parallel)
%ASSEMBLEGLOBALBILINEARFORM Summary of this function goes here
%   Detailed explanation goes here

if parallel
    [S,D,G,GA] = assembleLinearMatrices_An_Dn_parallel(A,D,kappa,T,T_P2,Nd,nodes2mesh_x,nodes2mesh_y);
else
    [S,D,G,GA] = assembleLinearMatrices_An_Dn_sequential(A,D,kappa,T,T_P2,Nd,nodes2mesh_x,nodes2mesh_y);
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [S,D,G,GA] = assembleLinearMatrices_An_Dn_sequential(A,D,kappa,T,T_P2,Nd,nodes2mesh_x,nodes2mesh_y)
Nx = sum(logical(nodes2mesh_x));

s_i = zeros(9*size(T,1),1);
s_j = zeros(9*size(T,1),1);
s_val = zeros(9*size(T,1),1);
g_i = zeros(9*size(T,1),1);
g_j = zeros(9*size(T,1),1);
g_val = zeros(9*size(T,1),1);
a_i = zeros(9*size(T,1),1);
a_j = zeros(9*size(T,1),1);
a_val = zeros(9*size(T,1),1);
d_i = zeros(9*size(T,1),1);
d_j = zeros(9*size(T,1),1);
d_val = zeros(9*size(T,1),1);
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

[quad,w] = getQuadrature(7);
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
    D_in_quad = zeros(2,no_of_quad_points);

    for i = 1:no_of_basis_P2
        index = nodes2mesh_x(tri_P2(i));
        if index ~= 0
            A_in_quad(1,:) = A_in_quad(1,:) + A(index)*phi_P2_in_quad(:,i)';
            D_in_quad(1,:) = D_in_quad(1,:) + D(index)*phi_P2_in_quad(:,i)';
        end

        index = nodes2mesh_y(tri_P2(i));
        if index ~= 0
            A_in_quad(2,:) = A_in_quad(2,:) + A(Nx + index)*phi_P2_in_quad(:,i)';
            D_in_quad(2,:) = D_in_quad(2,:) + D(Nx + index)*phi_P2_in_quad(:,i)';
        end
    end

    %% assemble
    for i = 1:3
        for j = i:3
            int = diag((1i/kappa * grad_in_BT(:,:,i)+ phi_in_quad(:,i)'.*A_in_quad)'*(phi_in_quad(:,j)'.*D_in_quad));
            e = detBT*evaluateQuadrature(transpose(int),w);

            int_g = diag((phi_in_quad(:,i)'.*D_in_quad)'*(phi_in_quad(:,j)'.*D_in_quad));
            e_g = detBT*evaluateQuadrature(transpose(int_g),w);

            int_a = diag((phi_in_quad(:,i)'.*A_in_quad)'*(phi_in_quad(:,j)'.*D_in_quad));
            e_a = detBT*evaluateQuadrature(transpose(int_a),w);

            int_d = diag((phi_in_quad(:,i)'.*D_in_quad)'*(1i/kappa * grad_in_BT(:,:,j)+ phi_in_quad(:,j)'.*A_in_quad));
            e_d = detBT*evaluateQuadrature(transpose(int_d),w);


            if i == j
                s_i(ind) = tri(i);
                s_j(ind) = tri(i);
                s_val(ind) = e;
                g_i(ind) = tri(i);
                g_j(ind) = tri(i);
                g_val(ind) = e_g;
                a_i(ind) = tri(i);
                a_j(ind) = tri(i);
                a_val(ind) = e_a;
                d_i(ind) = tri(i);
                d_j(ind) = tri(i);
                d_val(ind) = e_d;
                ind = ind + 1;
            else
                s_i(ind) = tri(i);
                s_j(ind) = tri(j);
                s_val(ind) = conj(e);
                g_i(ind) = tri(i);
                g_j(ind) = tri(j);
                g_val(ind) = conj(e_g);
                a_i(ind) = tri(i);
                a_j(ind) = tri(j);
                a_val(ind) = conj(e_a);
                d_i(ind) = tri(i);
                d_j(ind) = tri(j);
                d_val(ind) = conj(e_d);
                ind = ind + 1;

                s_i(ind) = tri(j);
                s_j(ind) = tri(i);
                s_val(ind) = e;
                g_i(ind) = tri(j);
                g_j(ind) = tri(i);
                g_val(ind) = e_g;
                a_i(ind) = tri(j);
                a_j(ind) = tri(i);
                a_val(ind) = e_a;
                d_i(ind) = tri(j);
                d_j(ind) = tri(i);
                d_val(ind) = e_d;
                ind = ind + 1;
            end
        end
    end
end

S = sparse(s_i,s_j,s_val);
G = sparse(g_i,g_j,g_val);
GA = sparse(a_i,a_j,a_val);
D = sparse(d_i,d_j,d_val);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [S,D,G,GA] = assembleLinearMatrices_An_Dn_parallel(A,D,kappa,T,T_P2,Nd,nodes2mesh_x,nodes2mesh_y)
Nx = sum(logical(nodes2mesh_x));
dim = size(Nd,1);

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

[quad,w] = getQuadrature(7);
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

spmd
    my_index_start = floor(size(T,1)*(spmdIndex - 1)/spmdSize + 1);
    my_index_end = floor(size(T,1)*(spmdIndex)/spmdSize);

    s_i = zeros(9*size(T,1),1);
    s_j = zeros(9*size(T,1),1);
    s_val = zeros(9*size(T,1),1);
    g_i = zeros(9*size(T,1),1);
    g_j = zeros(9*size(T,1),1);
    g_val = zeros(9*size(T,1),1);
    a_i = zeros(9*size(T,1),1);
    a_j = zeros(9*size(T,1),1);
    a_val = zeros(9*size(T,1),1);
    d_i = zeros(9*size(T,1),1);
    d_j = zeros(9*size(T,1),1);
    d_val = zeros(9*size(T,1),1);
    ind = 1;

    for k = my_index_start:my_index_end
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
        D_in_quad = zeros(2,no_of_quad_points);

        for i = 1:no_of_basis_P2
            index = nodes2mesh_x(tri_P2(i));
            if index ~= 0
                A_in_quad(1,:) = A_in_quad(1,:) + A(index)*phi_P2_in_quad(:,i)';
                D_in_quad(1,:) = D_in_quad(1,:) + D(index)*phi_P2_in_quad(:,i)';
            end

            index = nodes2mesh_y(tri_P2(i));
            if index ~= 0
                A_in_quad(2,:) = A_in_quad(2,:) + A(Nx + index)*phi_P2_in_quad(:,i)';
                D_in_quad(2,:) = D_in_quad(2,:) + D(Nx + index)*phi_P2_in_quad(:,i)';
            end
        end

        %% assemble
        for i = 1:3
            for j = i:3
                int = diag((1i/kappa * grad_in_BT(:,:,i)+ phi_in_quad(:,i)'.*A_in_quad)'*(phi_in_quad(:,j)'.*D_in_quad));
                e = detBT*evaluateQuadrature(transpose(int),w);

                int_g = diag((phi_in_quad(:,i)'.*D_in_quad)'*(phi_in_quad(:,j)'.*D_in_quad));
                e_g = detBT*evaluateQuadrature(transpose(int_g),w);

                int_a = diag((phi_in_quad(:,i)'.*A_in_quad)'*(phi_in_quad(:,j)'.*D_in_quad));
                e_a = detBT*evaluateQuadrature(transpose(int_a),w);

                int_d = diag((phi_in_quad(:,i)'.*D_in_quad)'*(1i/kappa * grad_in_BT(:,:,j)+ phi_in_quad(:,j)'.*A_in_quad));
                e_d = detBT*evaluateQuadrature(transpose(int_d),w);

                if i == j
                    s_i(ind) = tri(i);
                    s_j(ind) = tri(i);
                    s_val(ind) = e;
                    g_i(ind) = tri(i);
                    g_j(ind) = tri(i);
                    g_val(ind) = e_g;
                    a_i(ind) = tri(i);
                    a_j(ind) = tri(i);
                    a_val(ind) = e_a;
                    d_i(ind) = tri(i);
                    d_j(ind) = tri(i);
                    d_val(ind) = e_d;
                    ind = ind + 1;
                else
                    s_i(ind) = tri(i);
                    s_j(ind) = tri(j);
                    s_val(ind) = conj(e);
                    g_i(ind) = tri(i);
                    g_j(ind) = tri(j);
                    g_val(ind) = conj(e_g);
                    a_i(ind) = tri(i);
                    a_j(ind) = tri(j);
                    a_val(ind) = conj(e_a);
                    d_i(ind) = tri(i);
                    d_j(ind) = tri(j);
                    d_val(ind) = conj(e_d);
                    ind = ind + 1;

                    s_i(ind) = tri(j);
                    s_j(ind) = tri(i);
                    s_val(ind) = e;
                    g_i(ind) = tri(j);
                    g_j(ind) = tri(i);
                    g_val(ind) = e_g;
                    a_i(ind) = tri(j);
                    a_j(ind) = tri(i);
                    a_val(ind) = e_a;
                    d_i(ind) = tri(j);
                    d_j(ind) = tri(i);
                    d_val(ind) = e_d;
                    ind = ind + 1;
                end
            end
        end
    end

    del_index = find(s_i==0,1);
    if isempty(del_index)
        del_index = length(s_i)+1;
    end
    S_worker = sparse(s_i(1:del_index-1),s_j(1:del_index-1),s_val(1:del_index-1),dim,dim);
    G_worker = sparse(g_i(1:del_index-1),g_j(1:del_index-1),g_val(1:del_index-1),dim,dim);
    GA_worker = sparse(a_i(1:del_index-1),a_j(1:del_index-1),a_val(1:del_index-1),dim,dim);
    D_worker = sparse(d_i(1:del_index-1),d_j(1:del_index-1),d_val(1:del_index-1),dim,dim);
end

% collect results from workers
S = sparse(dim,dim);
G = sparse(dim,dim);
GA = sparse(dim,dim);
D = sparse(dim,dim);
for j = 1:length(my_index_start)
    S = S + S_worker{j};
    G = G + G_worker{j};
    GA = GA + GA_worker{j};
    D = D + D_worker{j};
end

end

