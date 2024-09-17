function M = assemble_mass_A_P2(T,Nd,nodes2mesh_x,nodes2mesh_y)
%ASSEMBLE_MASS_A_P2 Summary of this function goes here
%   Detailed explanation goes here

Nx = sum(logical(nodes2mesh_x));

m_i = zeros(2*36*size(T,1),1);
m_j = zeros(2*36*size(T,1),1);
m = zeros(2*36*size(T,1),1);
ind = 1;

phi = { @(x) (1 - x(1) - x(2))*(1-2*x(1)-2*x(2));
    @(x) 4*x(1)*(1-x(1)-x(2));
    @(x) x(1)*(2*x(1)-1);
    @(x) 4*x(1)*x(2);
    @(x) x(2)*(2*x(2)-1);
    @(x) 4*x(2)*(1-x(1)-x(2)) };

[quad,w] = getQuadrature(7);
phi_in_quad = zeros(6,7);
for i = 1:6
    for j = 1:7
        phi_in_quad(i,j) = phi{i}(quad(:,j));
    end
end

M_loc = zeros(6,6);
for i = 1:6
    for j = 1:6
        M_loc(i,j) = evaluateQuadrature(phi_in_quad(i,:).*phi_in_quad(j,:),w);
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

    %% assemble
    for i = 1:6
        for j = 1:6
            index_i = nodes2mesh_x(tri(i));
            index_j = nodes2mesh_x(tri(j));
            if index_i ~= 0 && index_j ~= 0
                m_i(ind) = index_i;
                m_j(ind) = index_j;
                m(ind) = detBT*M_loc(i,j);
                ind = ind + 1;
            end

            index_i = nodes2mesh_y(tri(i));
            index_j = nodes2mesh_y(tri(j));
            if index_i ~= 0 && index_j ~= 0
                m_i(ind) = Nx + index_i;
                m_j(ind) = Nx + index_j;
                m(ind) = detBT*M_loc(i,j);
                ind = ind + 1;
            end
        end
    end
end

%A = sparse(A);
del_index = find(m_i==0,1);
if isempty(del_index)
    del_index = length(m_i)+1;
end
M = sparse(m_i(1:del_index-1),m_j(1:del_index-1),m(1:del_index-1));
end