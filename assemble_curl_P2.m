function A = assemble_curl_P2(T,Nd,nodes2mesh_x,nodes2mesh_y)
%ASSEMBLE_CURL Summary of this function goes here
%   Detailed explanation goes here

Nx = sum(logical(nodes2mesh_x));
a_i = zeros(2*36*size(T,1),1);
a_j = zeros(2*36*size(T,1),1);
a = zeros(2*36*size(T,1),1);
ind = 1;

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
    
    %% assemble
    for i = 1:6
        for j = 1:6
            index_i = nodes2mesh_x(tri(i));
            index_j = nodes2mesh_x(tri(j));
            if index_i ~= 0 && index_j ~= 0
                a_i(ind) = index_i;
                a_j(ind) = index_j;
                a(ind) = detBT*evaluateQuadrature(grad_in_BT(2,:,i).*grad_in_BT(2,:,j),w);
                ind = ind + 1;
            end
            
            index_i = nodes2mesh_y(tri(i));
            index_j = nodes2mesh_y(tri(j));
            if index_i ~= 0 && index_j ~= 0
                a_i(ind) = Nx + index_i;
                a_j(ind) = Nx + index_j;
                a(ind) = detBT*evaluateQuadrature(grad_in_BT(1,:,i).*grad_in_BT(1,:,j),w);
                ind = ind + 1;
            end
            
            index_i = nodes2mesh_x(tri(i));
            index_j = nodes2mesh_y(tri(j));
            if index_i ~= 0 && index_j ~= 0
                a_i(ind) = index_i;
                a_j(ind) = Nx + index_j;
                a(ind) = -detBT*evaluateQuadrature(grad_in_BT(2,:,i).*grad_in_BT(1,:,j),w);
                ind = ind + 1;
            end
            
            index_i = nodes2mesh_y(tri(i));
            index_j = nodes2mesh_x(tri(j));
            if index_i ~= 0 && index_j ~= 0
                a_i(ind) = Nx + index_i;
                a_j(ind) = index_j;
                a(ind) = -detBT*evaluateQuadrature(grad_in_BT(2,:,i).*grad_in_BT(1,:,j),w);
                ind = ind + 1;
            end     
        end
    end
end

del_index = find(a_i==0,1);
if isempty(del_index)
    del_index = length(a_i)+1;
end
A = sparse(a_i(1:del_index-1),a_j(1:del_index-1),a(1:del_index-1));
end

