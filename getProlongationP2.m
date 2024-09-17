function P = getProlongationP2(T_H,T_h,P0)
%GETPROLONGATION Summary of this function goes here
%   Detailed explanation goes here

tri_H = T_H.t;
tri_h = T_h.t;

Nd_H = T_H.p;
Nd_h = T_h.p;

p_i = zeros(100*size(tri_H,1),1);
p_j = zeros(100*size(tri_H,1),1);
p_val = zeros(100*size(tri_H,1),1);
ind = 1;

phi1 = @(x) (1 - x(1) - x(2))*(1-2*x(1)-2*x(2));
phi2 = @(x) 4*x(1)*(1-x(1)-x(2));
phi3 = @(x) x(1)*(2*x(1)-1);
phi4 = @(x) 4*x(1)*x(2);
phi5 = @(x) x(2)*(2*x(2)-1);
phi6 = @(x) 4*x(2)*(1-x(1)-x(2));

for k = 1:size(tri_H,1) % k = index of coarse triangle
    tri = tri_H(k,:); %node index of triangle
    z1 = Nd_H(tri(1),:); %coordinates of 1st triangle node
    z2 = Nd_H(tri(3),:); %coordinates of 2nd triangle node
    z3 = Nd_H(tri(5),:); %coordinates of 3rd triangle node

    %% transformation to refenrence triangle
    BT = [z2(1)-z1(1), z3(1)-z1(1); ...
        z2(2)-z1(2), z3(2)-z1(2)];

    BTinv = BT\eye(2);

    b = [z1(1); z1(2)];

    % get indecies of fine triangles on coarse triangle
    active_fine_triangles = find(P0(:,k));

    for l = 1:size(active_fine_triangles,1)
        for i = 1:6
            index_i = tri_h(active_fine_triangles(l),i);
            % j = 1
            index_j = tri_H(k,1);
            p_i(ind) = index_i;
            p_j(ind) = index_j;
            p_val(ind) = phi1(BTinv*(Nd_h(index_i,:)' - b));
            ind = ind + 1;

            % j = 2
            index_j = tri_H(k,2);
            p_i(ind) = index_i;
            p_j(ind) = index_j;
            p_val(ind) = phi2(BTinv*(Nd_h(index_i,:)' - b));
            ind = ind + 1;

            % j = 3
            index_j = tri_H(k,3);
            p_i(ind) = index_i;
            p_j(ind) = index_j;
            p_val(ind) = phi3(BTinv*(Nd_h(index_i,:)' - b));
            ind = ind + 1;

            % j = 4
            index_j = tri_H(k,4);
            p_i(ind) = index_i;
            p_j(ind) = index_j;
            p_val(ind) = phi4(BTinv*(Nd_h(index_i,:)' - b));
            ind = ind + 1;

            % j = 5
            index_j = tri_H(k,5);
            p_i(ind) = index_i;
            p_j(ind) = index_j;
            p_val(ind) = phi5(BTinv*(Nd_h(index_i,:)' - b));
            ind = ind + 1;

            % j = 6
            index_j = tri_H(k,6);
            p_i(ind) = index_i;
            p_j(ind) = index_j;
            p_val(ind) = phi6(BTinv*(Nd_h(index_i,:)' - b));
            ind = ind + 1;

        end
    end
end

[p_unique, p_ind] = unique([p_i, p_j], 'rows');
p_i_new = p_unique(:,1);
p_j_new = p_unique(:,2);
p_val_new = p_val(p_ind);

del_index = find(p_i_new==0,1);
if ~isempty(del_index)
    p_i_new(del_index) = [];
    p_j_new(del_index) = [];
    p_val_new(del_index) = [];
end
P = sparse(p_i_new,p_j_new,p_val_new);
P = P';
end



