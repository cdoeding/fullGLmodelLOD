function [curlA,E_curl] = get_curl_of_A_P2(H,A,T,Nd,nodes2mesh_x,nodes2mesh_y)

E_curl = 0;

Nx = sum(logical(nodes2mesh_x));
curlA = zeros(size(T,1),1);

grad = { @(x) [4*x(1) + 4*x(2) - 3; 4*x(1) + 4*x(2) - 3];
    @(x) [-4*(2*x(1) + x(2) - 1); -4*x(1)];
    @(x) [4*x(1) - 1; 0];
    @(x) [4*x(2); 4*x(1)];
    @(x) [0; 4*x(2) - 1];
    @(x) [-4*x(2); -4*(x(1) + 2*x(2) - 1)]};

[quad,w] = getQuadrature(7);
no_of_basis = size(grad,1);
no_of_quad_points = length(w);
grad_in_quad = zeros(2,no_of_quad_points,no_of_basis);
for i = 1:no_of_basis
    for j = 1:no_of_quad_points
        grad_in_quad(:,j,i) = grad{i}(quad(:,j));
    end
end

curlA_in_quad = zeros(1,no_of_quad_points);
H_in_quad = zeros(1,no_of_quad_points);

for k = 1:size(T,1)
    tri = T(k,:); %node index of triangle
    z1 = Nd(T(k,1),:); %coordinates of 1st triangle node
    z2 = Nd(T(k,3),:); %coordinates of 2nd triangle node
    z3 = Nd(T(k,5),:); %coordinates of 3rd triangle node
    
    % transformation to refenrence triangle
    BT = [z2(1)-z1(1), z3(1)-z1(1); ...
        z2(2)-z1(2), z3(2)-z1(2)];
    
    BTinv = inv(BT)';
    grad_in_BT = pagemtimes(BTinv,grad_in_quad);
    
    b = [z1(1); z1(2)];
    
    detBT = BT(1,1)*BT(2,2)-BT(1,2)*BT(2,1);

    for i = 1:no_of_quad_points
        H_in_quad(i) = H(BT*quad(:,i) + b);
    end
    
    %% assemble
    % construct A contributions on triangle
    index11 = nodes2mesh_x(tri(1));
    if index11 ~= 0
        %define contribution of phi1
        A11 = A(index11);
    else
        A11 = 0;
    end
    
    index12 = nodes2mesh_x(tri(2));
    if index12 ~= 0
        %define contribution of phi1
        A12 = A(index12);
    else
        A12 = 0;
    end
    index13 = nodes2mesh_x(tri(3));
    if index13 ~= 0
        %define contribution of phi3
        A13 = A(index13);
    else
        A13 = 0;
    end
    
    index14 = nodes2mesh_x(tri(4));
    if index14 ~= 0
        %define contribution of phi1
        A14 = A(index14);
    else
        A14 = 0;
    end
    
    index15 = nodes2mesh_x(tri(5));
    if index15 ~= 0
        % define contribution of phi1
        A15 = A(index15);
    else
        A15 = 0;
    end
    index16 = nodes2mesh_x(tri(6));
    if index16 ~= 0
        %define contribution of phi3
        A16 = A(index16);
    else
        A16 = 0;
    end
    
    index21 = nodes2mesh_y(tri(1));
    if index21 ~= 0
        %define contribution of phi1
        A21 = A(Nx + index21);
    else
        A21 = 0;
    end
    
    index22 = nodes2mesh_y(tri(2));
    if index22 ~= 0
        %define contribution of phi1
        A22 = A(Nx + index22);
    else
        A22 = 0;
    end
    index23 = nodes2mesh_y(tri(3));
    if index23 ~= 0
        %define contribution of phi3
        A23 = A(Nx + index23);
    else
        A23 = 0;
    end
    
    index24 = nodes2mesh_y(tri(4));
    if index24 ~= 0
        %define contribution of phi1
        A24 = A(Nx + index24);
    else
        A24 = 0;
    end
    
    index25 = nodes2mesh_y(tri(5));
    if index25 ~= 0
        % define contribution of phi1
        A25 = A(Nx + index25);
    else
        A25 = 0;
    end
    index26 = nodes2mesh_y(tri(6));
    if index26 ~= 0
        %define contribution of phi3
        A26 = A(Nx + index26);
    else
        A26 = 0;
    end
    
    curlA_in_quad = [A21, -A11]*(grad_in_BT(:,:,1)) ...
        + [A22, -A12]*(grad_in_BT(:,:,2)) ...
        + [A23, -A13]*(grad_in_BT(:,:,3)) ...
        + [A24, -A14]*(grad_in_BT(:,:,4)) ...
        + [A25, -A15]*(grad_in_BT(:,:,5)) ...
        + [A26, -A16]*(grad_in_BT(:,:,6));
    
    curlA(k) = curlA_in_quad(1);
    
    E_curl = E_curl + detBT*evaluateQuadrature((curlA_in_quad - H_in_quad).^2,w);
    
end

end

