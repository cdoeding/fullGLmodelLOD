function M = assemble_mass_matrix(T,Nd,nodes2mesh)
%ASSEMBLE_MASS_MATRIX Summary of this function goes here
%   Detailed explanation goes here

m_i = zeros(9*size(T,1),1);
m_j = zeros(9*size(T,1),1);
m = zeros(9*size(T,1),1);
ind = 1;

M_loc = (1/24)*[2 1 1; 1 2 1; 1 1 2];

for k = 1:size(T,1)
    tri = T(k,:); %node index of triangle
    z1 = Nd(T(k,1),:); %coordinates of 1st triangle node
    z2 = Nd(T(k,2),:); %coordinates of 2nd triangle node
    z3 = Nd(T(k,3),:); %coordinates of 3rd triangle node

    %% transformation to refenrence triangle
    BT = [z2(1)-z1(1), z3(1)-z1(1); ...
        z2(2)-z1(2), z3(2)-z1(2)];
    detBT = BT(1,1)*BT(2,2)-BT(1,2)*BT(2,1);

    for i = 1:3
        for j = 1:3
            m_i(ind) = nodes2mesh(tri(i),1);
            m_j(ind) = nodes2mesh(tri(j),1);
            m(ind) = detBT*M_loc(i,j);
            ind = ind + 1;
        end
    end
end
del_index = find(m_i==0,1);
if isempty(del_index)
    del_index = length(m_i)+1;
end
M = sparse(m_i(1:del_index-1),m_j(1:del_index-1),m(1:del_index-1));
end