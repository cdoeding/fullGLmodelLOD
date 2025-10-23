H_mag = @(x) 10*sin(pi*x(2))*sin(pi*x(1)); % external magnetic field

x_a = 0; % domain left/bottom end point
x_b = 1; % domain right/top end point
area = 1; % area of rectangle

H_level = 3; % coarse mesh size level u (LOD)
h_level = 9; % fine mesh size level u (LOD)
ht_level = 7; % mesh size level A (P2)

%% coarse, fine mesh and patches
boundary_u = 'Neumann';
boundary_A = 'non-natural';
[T_H,T_h,~,~] = getCoarseFineTriangulation(x_a,x_b,H_level,h_level);

% mesh for A equation (P1 & P2)
if h_level <= ht_level
    [~,T_ht,~,Pt0] = getCoarseFineTriangulation(x_a,x_b,h_level,ht_level);
    tic;
    [T_h_P2,T_ht_P2,P2] = getCoarseFineTriangulation_for_P2(T_h,T_ht,Pt0);
    toc;
else
    [T_ht,~,~,Pt0] = getCoarseFineTriangulation(x_a,x_b,ht_level,h_level);
    [T_ht_P2,T_h_P2,P2] = getCoarseFineTriangulation_for_P2(T_ht,T_h,Pt0);
end

% boundaries for H mesh (u equation)
B_H = getBoundaryNodes(T_H.p,x_a,x_b,boundary_u);
[nodes2mesh_H,nodes2mesh_Hx,nodes2mesh_Hy] = getNodes2Mesh(T_H.p,x_a,x_b,boundary_u,boundary_A);
Nx_H = sum(logical(nodes2mesh_Hx));

% boundaries for h mesh (u equation)
B_h = getBoundaryNodes(T_h.p,x_a,x_b,boundary_u);
[Bx_h,By_h] = getBoundaryNodes_for_A(T_h_P2.p,x_a,x_b,boundary_A);
Bdx_h = getBoundaryRestriction(Bx_h);
Bdy_h = getBoundaryRestriction(By_h);

% boundaries for ht mesh (A equation)
B_ht = getBoundaryNodes(T_ht.p,x_a,x_b,boundary_u);
[Bx_ht,By_ht] = getBoundaryNodes_for_A(T_ht_P2.p,x_a,x_b,boundary_A);
Bdx_ht = getBoundaryRestriction(Bx_ht);
Bdy_ht = getBoundaryRestriction(By_ht);

[nodes2mesh_h,nodes2mesh_hx,nodes2mesh_hy] = getNodes2MeshP1P2(T_h.p,T_h_P2.p,x_a,x_b,boundary_u,boundary_A);
[nodes2mesh_ht,nodes2mesh_htx,nodes2mesh_hty] = getNodes2MeshP1P2(T_ht.p,T_ht_P2.p,x_a,x_b,boundary_u,boundary_A);

Nx_ht = sum(logical(nodes2mesh_htx));

% assemble matricies for A equation
S_curl = assemble_curl_P2(T_ht_P2.t,T_ht_P2.p,nodes2mesh_htx,nodes2mesh_hty);
S_div = assemble_div_P2(T_ht_P2.t,T_ht_P2.p,nodes2mesh_htx,nodes2mesh_hty);

% assemble rhs
RHS_H_mag = assemble_RHS_curlH_P2(H_mag,T_ht_P2.t,T_ht_P2.p,nodes2mesh_htx,nodes2mesh_hty);

%% compute Astar via curl(A_star) = H and div(A_star) = 0
A_for_LOD_ht = (S_curl  + S_div)\RHS_H_mag;
A1_for_LOD_h = Bdx_h*P2'*Bdx_ht'*A_for_LOD_ht(1:Nx_ht); % prolongation
A2_for_LOD_h = Bdy_h*P2'*Bdy_ht'*A_for_LOD_ht(Nx_ht+1:end); % prolongation
A_for_LOD_h = [A1_for_LOD_h; A2_for_LOD_h];

clearvars -except A_for_LOD_h
save("A_star.mat",'-v7.3')


