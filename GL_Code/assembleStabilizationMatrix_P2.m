function G = assembleStabilizationMatrix_P2(A,T,Nd)
%   Summary of this function goes here
%   Detailed explanation goes here

% the local part of the stabilization matrix is a 4x4 matrix. 
% we will store the values in a vector of length 4x4=36.
G_i = zeros(36*size(T,1),1); % stores the global i-index of the matrix
G_j = zeros(36*size(T,1),1); % stores the global j-index of the matrix
G_val = zeros(36*size(T,1),1); % stores the value of the element contribution
ind = 1; % index counter

% basis functions P2
phi = {@(x) (1 - x(1) - x(2))*(1-2*x(1)-2*x(2));
    @(x) 4*x(1)*(1-x(1)-x(2));
    @(x) x(1)*(2*x(1)-1);
    @(x) 4*x(1)*x(2);
    @(x) x(2)*(2*x(2)-1);
    @(x) 4*x(2)*(1-x(1)-x(2))};

no_of_basis = size(phi,1);

% get quadrature points (quad) and and correpsonding quadrature weights (w)
[x_quad,w] = getQuadrature(7); % 7 = number of quadrature points
no_of_quad_points = length(w);

% values of the 6 basis functions in the quadrature nodes
phi_in_quad_point = zeros(no_of_basis,no_of_quad_points); % 6 basis functions, 7 quad points
for i = 1:no_of_basis
     for q = 1:no_of_quad_points
          phi_in_quad_point(i,q) = phi{i}(x_quad(:,q));
     end
end

for k = 1:size(T,1)

    tri = T(k,:); % Lagrange node indices of the triangle
    z1 = Nd(T(k,1),:); %coordinates of 1st triangle node
    z2 = Nd(T(k,3),:); %coordinates of 2nd triangle node
    z3 = Nd(T(k,5),:); %coordinates of 3rd triangle node
    
    %% transformation to refenrence triangle
    BT = [z2(1)-z1(1), z3(1)-z1(1); ...
        z2(2)-z1(2), z3(2)-z1(2)];
    
    detBT = BT(1,1)*BT(2,2)-BT(1,2)*BT(2,1);
    
    % transformation F : T_0 -> T (reference triangle T_0 to current triangle T)
    % F(x) = b + BT x
    b = [z1(1); z1(2)]; % b=z1'
    
    A_in_quad_point = zeros(2,7);
    A_in_quad_point(:,1) = A(BT*x_quad(:,1)+b);
    A_in_quad_point(:,2) = A(BT*x_quad(:,2)+b);
    A_in_quad_point(:,3) = A(BT*x_quad(:,3)+b);
    A_in_quad_point(:,4) = A(BT*x_quad(:,4)+b);
    A_in_quad_point(:,5) = A(BT*x_quad(:,5)+b);
    A_in_quad_point(:,6) = A(BT*x_quad(:,6)+b);
    A_in_quad_point(:,7) = A(BT*x_quad(:,7)+b);

    % weight * |A|^2 on reference element (unit triangle)
    abs_A_in_quad_point = zeros(7);
    abs_A_in_quad_point(1) = w(1) * norm( A_in_quad_point(:,1) ).^2;
    abs_A_in_quad_point(2) = w(2) * norm( A_in_quad_point(:,2) ).^2;
    abs_A_in_quad_point(3) = w(3) * norm( A_in_quad_point(:,3) ).^2;
    abs_A_in_quad_point(4) = w(4) * norm( A_in_quad_point(:,4) ).^2;
    abs_A_in_quad_point(5) = w(5) * norm( A_in_quad_point(:,5) ).^2;
    abs_A_in_quad_point(6) = w(6) * norm( A_in_quad_point(:,6) ).^2;       
    abs_A_in_quad_point(7) = w(7) * norm( A_in_quad_point(:,7) ).^2;    
        
    % f = @(x) norm( A(BT*x+b) ).^2;

    %% assemble
    for i = 1:no_of_basis
        for j = 1:no_of_basis
        
             % exploit symmetry of matrix 
            if ( i <= j )       
            
               index_i = tri(i); % get global index of local node i
               index_j = tri(j); % get global index of local node j

               G_i(ind) = index_i;
               G_j(ind) = index_j;
            
               % integrate on reference element
               e =      abs_A_in_quad_point(1) * phi_in_quad_point(i,1) * phi_in_quad_point(j,1) ... 
                       + abs_A_in_quad_point(2) * phi_in_quad_point(i,2) * phi_in_quad_point(j,2) ...
                       + abs_A_in_quad_point(3) * phi_in_quad_point(i,3) * phi_in_quad_point(j,3) ...
                       + abs_A_in_quad_point(4) * phi_in_quad_point(i,4) * phi_in_quad_point(j,4) ...
                       + abs_A_in_quad_point(5) * phi_in_quad_point(i,5) * phi_in_quad_point(j,5) ...
                       + abs_A_in_quad_point(6) * phi_in_quad_point(i,6) * phi_in_quad_point(j,6) ...
                       + abs_A_in_quad_point(7) * phi_in_quad_point(i,7) * phi_in_quad_point(j,7);                    
               e = detBT * e;
            
               %
               % e_TEST = detBT*integrate_unit_triangle(@(x) (f(x)*phi{i}(x)*phi{j}(x)) , 6 ); %
   
               G_val(ind) = e;
               ind = ind + 1;
               
               if ( i < j )
                      G_i(ind) = tri(j); % get global index of local node i
                      G_j(ind) = tri(i); % get global index of local node j
                      G_val(ind) = e;
                      ind = ind + 1;
               end
               
            end % end "if ( i <= j )"
               
        end % end i
    end % end j
    
end

G = sparse(G_i,G_j,G_val);

end

