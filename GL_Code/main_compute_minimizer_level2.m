% "A multiscale approach to the stationary Ginzburg-Landau equations of superconductivity"
% by Christian Doeding, Benjamin Doerich, and Patrick Henning
%
% Compute minimizers (u,A) of the Ginzburg-Landau energy in 2d
% E(u,A)=\int_\Omega |(i/kappa)*grad(u)+Au|^2+(1/2)*(1-|u|^2)^2+|curl(A)-H|^2+|div(A)|^2 dx
%
% with
% LOD discretization in order parameter u
% P2 Lagrange FEM discretization in vector potential A
% Sobolev gradient flow discretization for minimization
%
% on level 2 (H = 2^(-5), h = 2^(-4), l = 5).
%
% usage:
% - set model parameter and solver preferences in 'initialize parameters'
% - enable/disable parallel computing
% - minimizer from "main_compute_minimizer_level1.m" required
% - A_star from "main_compute_Astar.m" required

%% initialize parameters
parallel = true;

if parallel == true
    parpool('local','IdleTimeout', Inf)
end

kappa = 5; % GL parameter
H_mag = @(x) 10*sin(pi*x(2))*sin(pi*x(1)); % external magnetic field

x_a = 0; % domain left/bottom end point
x_b = 1; % domain right/top end point
area = 1; % area of rectangle

H_level = 5; % coarse mesh size level u (LOD)
h_level = 9; % fine mesh size level u (LOD)
ht_level = 4; % mesh size level A (P2)
ell = 5; % oversampling of ell-layers for u (LOD)

tol = 10^(-10); % tolerance for termination
i_max = 100000; % maximum number of iterations

save_path = "solution_level2_kappa"+kappa+".mat"; % path for save file

%% coarse, fine mesh and patches
boundary_u = 'Neumann';
boundary_A = 'non-natural';
[T_H,T_h,P1,P0] = getCoarseFineTriangulation(x_a,x_b,H_level,h_level);
T = T_H.t;
Nd = T_H.p;
P1 = P1';

% mesh for A equation (P1 & P2)
if h_level <= ht_level
    [~,T_ht,Pt1,Pt0] = getCoarseFineTriangulation(x_a,x_b,h_level,ht_level);
    tic;
    [T_h_P2,T_ht_P2,P2] = getCoarseFineTriangulation_for_P2(T_h,T_ht,Pt0);
    toc;
else
    [T_ht,~,Pt1,Pt0] = getCoarseFineTriangulation(x_a,x_b,ht_level,h_level);
    [T_ht_P2,T_h_P2,P2] = getCoarseFineTriangulation_for_P2(T_ht,T_h,Pt0);
end
Pt = Pt1';

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

% patches
patches = getPatches(T_H,ell); % patches_ij non-zero iff jth triangle is in patch of ith triangle

% assemble matricies for u equation
M_H = assemble_mass_matrix(T_H.t,T_H.p,nodes2mesh_H);
M_h = assemble_mass_matrix(T_h.t,T_h.p,nodes2mesh_h);
M_ht = assemble_mass_matrix(T_ht.t,T_ht.p,nodes2mesh_ht);
Stiff_h = assemble_stiffness_matrix(T_h.t,T_h.p,nodes2mesh_h);

% assemble matricies for A equation
S_curl = assemble_curl_P2(T_ht_P2.t,T_ht_P2.p,nodes2mesh_htx,nodes2mesh_hty);
S_div = assemble_div_P2(T_ht_P2.t,T_ht_P2.p,nodes2mesh_htx,nodes2mesh_hty);
M_A = assemble_mass_A_P2(T_ht_P2.t,T_ht_P2.p,nodes2mesh_htx,nodes2mesh_hty);

% assemble rhs
RHS_H_mag = assemble_RHS_curlH_P2(H_mag,T_ht_P2.t,T_ht_P2.p,nodes2mesh_htx,nodes2mesh_hty);

%% compute Astar via curl(A_star) = H and div(A_star) = 0
A_for_LOD_h = load("A_star.mat","A_for_LOD_h").A_for_LOD_h;

%% compute corrector (construct LOD space)
% assemble bilinear form
S_for_LOD_h = assemble_bilinear_form_with_P2(A_for_LOD_h,kappa,T_h.t,T_h_P2.t,T_h.p,nodes2mesh_hx,nodes2mesh_hy,parallel);

beta = 0;
% compute LOD-corrector (for parallization use 'getCorrectorMatrixParallel_with_P2' in l.132)
disp("compute corrector")
Q = getCorrectorMatrix_with_P2(T_H,T_h,T_h_P2,patches,A_for_LOD_h,kappa,beta,S_for_LOD_h,M_h,P1',P0,B_H,B_h,nodes2mesh_hx,nodes2mesh_hy,parallel);

%% get initial value
path_best = "solution_level1_kappa"+kappa+".mat";
z = load(path_best,"A_ht","u_h");

u_h = getBestApproximation_for_u(path_best,x_a,x_b,kappa,Stiff_h,M_h,h_level,9);
A_ht = getBestApproximation_for_A(path_best,x_a,x_b,ht_level,4,Bdx_ht,Bdy_ht,M_A,S_curl,S_div);

A1_h = Bdx_h*P2'*Bdx_ht'*A_ht(1:Nx_ht); % prolongation
A2_h = Bdy_h*P2'*Bdy_ht'*A_ht(Nx_ht+1:end); % prolongation
A_h = [A1_h; A2_h];

% initial value in LOD space (L2-projection) and representation
M_LOD = (P1 + Q)*M_h*(P1 + Q)';
u_LOD = M_LOD'\((P1+Q)*M_h'*u_h);
u_h = (P1+Q)'*u_LOD;

%% corrector computation and assembling
[S_h, G_h] = assemble_bilinear_form_and_A_density_with_P2(A_h,kappa,T_h.t,T_h_P2.t,T_h.p,nodes2mesh_hx,nodes2mesh_hy,parallel);
Fu_h = assemble_nonlinear_term_impl(u_h,T_h.t,T_h.p,nodes2mesh_h,parallel);

% assemble LOD matrices
S_LOD = (P1 + Q)*S_h*(P1 + Q)';
G_LOD = (P1 + Q)*G_h*(P1 + Q)';
Fu_LOD = (P1 + Q)*Fu_h*(P1 + Q)';

% compute initial matrices for A-equation
% get u on ht mesh
u_ht = M_ht'\(Pt*M_h'*u_h);
FA_ht = assemble_nonlinear_term_A_density_P2(u_ht,T_ht_P2.t,T_ht.t,T_ht_P2.p,nodes2mesh_ht,nodes2mesh_htx,nodes2mesh_hty,parallel);
RHS_u = assemble_RHS_u_nabla_u_P2_new(u_ht,T_ht_P2.t,T_ht.t,T_ht_P2.p,nodes2mesh_ht,nodes2mesh_htx,nodes2mesh_hty,parallel);

% compute energy
[curlA,E_curl] = get_curl_of_A_P2(H_mag,A_ht,T_ht_P2.t,T_ht_P2.p,nodes2mesh_htx,nodes2mesh_hty,parallel);
E = 0.5*real(u_LOD'*S_LOD*u_LOD + 0.5*area - u_LOD'*M_LOD*u_LOD + 0.5*u_LOD'*Fu_LOD*u_LOD + E_curl + A_ht'*S_div*A_ht);

%% computation of Sobolev gradient flow
delta = tol;
counter = 0;
stab = 0.1;

while abs(delta) >= tol && counter < i_max
    E_old = E;

    % compute Sobolev gradient
    g_u = (S_LOD + Fu_LOD + stab*M_LOD + G_LOD)\(((1+stab)*M_LOD + G_LOD)*u_LOD);
    g_A = (S_curl + S_div + FA_ht + stab*M_A)\(stab*M_A*A_ht - (0.5i/kappa)*RHS_u + RHS_H_mag);

    % descent direction
    d_u = g_u - u_LOD;
    d_A = g_A - A_ht;

    % find optimal step size tau
    d_uh = (P1 + Q)'*d_u;
    d1_h = Bdx_h*P2'*Bdx_ht'*d_A(1:Nx_ht); % prolongation
    d2_h = Bdy_h*P2'*Bdy_ht'*d_A(Nx_ht+1:end); % prolongation
    d_Ah = [d1_h; d2_h];

    [F_un_dn,F_dn_dn] = assembleNonlinearMatrices_un_dn_P1(u_h,d_uh,T_h.t,T_h.p,parallel);
    F_un_dn_LOD = (P1 + Q)*F_un_dn*(P1 + Q)';
    F_dn_dn_LOD = (P1 + Q)*F_dn_dn*(P1 + Q)';

    [SD_h, DS_h, DD_h, AD_h] = assembleLinearMatrices_An_Dn(A_h,d_Ah,kappa,T_h.t,T_h_P2.t,T_h.p,nodes2mesh_hx,nodes2mesh_hy,parallel);
    SD_LOD = (P1 + Q)*SD_h*(P1 + Q)';
    DS_LOD = (P1 + Q)*DS_h*(P1 + Q)';
    DD_LOD = (P1 + Q)*DD_h*(P1 + Q)';
    AD_LOD = (P1 + Q)*AD_h*(P1 + Q)';

    c0 = 0.5*real( u_LOD'*S_LOD*u_LOD ...
        + 0.5*area ...
        - u_LOD'*M_LOD*u_LOD  ...
        + 0.5*u_LOD'*Fu_LOD*u_LOD ...
        + E_curl + A_ht'*S_div*A_ht );

    c1 = 0.5*real( u_LOD'*S_LOD*d_u + d_u'*S_LOD*u_LOD...
        + u_LOD'*(SD_LOD + DS_LOD)*u_LOD ...
        - d_u'*M_LOD*u_LOD - u_LOD'*M_LOD*d_u ... 
        + 0.5*d_u'*Fu_LOD*u_LOD + 0.5*u_LOD'*Fu_LOD*d_u ...
        + u_LOD'*F_un_dn_LOD*u_LOD ...
        + 2*A_ht'*S_curl*d_A ...
        - 2*RHS_H_mag'*d_A ...
        + A_ht'*S_div*d_A + d_A'*S_div*A_ht );

    c2 = 0.5*real( d_u'*S_LOD*d_u ...
        + u_LOD'*(SD_LOD + DS_LOD)*d_u ...
        + d_u'*(SD_LOD + DS_LOD)*u_LOD ...
        + u_LOD'*DD_LOD*u_LOD ...
        - d_u'*M_LOD*d_u ...
        + 0.5*u_LOD'*F_dn_dn_LOD*u_LOD ...
        + u_LOD'*F_un_dn_LOD*d_u ...
        + d_u'*F_un_dn_LOD*u_LOD ...
        + 0.5*d_u'*Fu_LOD*d_u ...
        + d_A'*S_curl*d_A + d_A'*S_div*d_A);

    c3 = 0.5*real( d_u'*DD_LOD*u_LOD + u_LOD'*DD_LOD*d_u ...
        + d_u'*(SD_LOD + DS_LOD)*d_u ...
        + d_u'*F_un_dn_LOD*d_u ... 
        + 0.5*u_LOD'*F_dn_dn_LOD*d_u + 0.5*d_u'*F_dn_dn_LOD*u_LOD );

    c4 = 0.5*real( d_u'*DD_LOD*d_u ...
        + 0.5*d_u'*F_dn_dn_LOD*d_u );

    g_tau = @(t) ( c0 + c1*t + c2*t.^2 + c3*t.^3 + c4*t.^4 );

    upper_search_bound = 30;
    lower_search_bound = 0.1;
    [tau,g_tau_value] = golden_search_section(lower_search_bound,upper_search_bound,g_tau);

    % Sobolev gradient descent step
    u_LOD = u_LOD + tau*d_u;
    A_ht = A_ht + tau*d_A;

    u_h = (P1 + Q)'*u_LOD;

    % update matrices
    % matricies on fine mesh (u-equation)
    % get A on fine mesh
    A1_h = Bdx_h*P2'*Bdx_ht'*A_ht(1:Nx_ht); % prolongation
    A2_h = Bdy_h*P2'*Bdy_ht'*A_ht(Nx_ht+1:end); % prolongation
    A_h = [A1_h; A2_h];

    Fu_LOD = Fu_LOD + 2*tau*F_un_dn_LOD + (tau^2)*F_dn_dn_LOD;
    S_LOD = S_LOD + tau*(SD_LOD + DS_LOD) + (tau^2)*DD_LOD;
    G_LOD = G_LOD + 2*tau*AD_LOD + (tau^2)*DD_LOD;

    % compute energy
    E_curl = E_curl + 2*tau*(A_ht - tau*d_A)'*S_curl*d_A - 2*tau*RHS_H_mag'*d_A + (tau^2)*d_A'*S_curl*d_A;

    % compute energy
    E = 0.5*real(u_LOD'*S_LOD*u_LOD + 0.5*area - u_LOD'*M_LOD*u_LOD + 0.5*u_LOD'*Fu_LOD*u_LOD + E_curl + A_ht'*S_div*A_ht);

    % update matricies for A-equation
    % get u on ht mesh
    u_ht = M_ht'\(Pt*M_h'*u_h);
    FA_ht = assemble_nonlinear_term_A_density_P2(u_ht,T_ht_P2.t,T_ht.t,T_ht_P2.p,nodes2mesh_ht,nodes2mesh_htx,nodes2mesh_hty,parallel);
    RHS_u = assemble_RHS_u_nabla_u_P2_new(u_ht,T_ht_P2.t,T_ht.t,T_ht_P2.p,nodes2mesh_ht,nodes2mesh_htx,nodes2mesh_hty,parallel);

    delta = E_old - E;
    counter = counter + 1;
    disp(counter)
    disp('step size')
    disp(tau)
    disp('energy difference')
    disp(delta)
    disp('energy')
    disp(E)

end

%% save
clearvars -except save_path H_level ht_level delta counter tol E Q u_LOD u_h M_h T_h nodes2mesh_h A_ht S_curl S_div M_A nodes2mesh_htx T_ht T_ht_P2 parallel
save(save_path,'-v7.3')

if parallel == true
    delete(gcp())
end



