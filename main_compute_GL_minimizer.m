% "A multiscale approach to the stationary Ginzburg-Landau equations of superconductivity"
% by Christian Doeding, Benjamin Doerich, and Patrick Henning
%
% Compute minimizers (u,A) of the Ginzburg-Landau energy in 2d
% E(u,A)=\int_\Omega |(i/kappa)*grad(u)+Au|^2+(1/2)*(1-|u|^2)^2+|curl(A)-H|^2+|div(A)|^2 dx
%
% with
% LOD discretization in order parameter u
% P2 Lagrange FEM discretization in vector potential A
% L^2 gradient flow discretization for minimization
%
% usage:
% - set model parameter and solver preferences in 'initialize parameters' section (l.21-l.43)
% - initial values for gradient flow can be set in 'initial values' section (l.105-l.118)
% - parallel computing for corrector can be used in l.130 and l.194

clearvars
close all

%% initialize parameters
save_results = true; % save results to .mat-file
plot_results = true; % plot the solution after computation

kappa = 12; % GL parameter
H_mag = @(x) 10*sin(pi*x(2))*sin(pi*x(1)); % external magnetic field

x_a = 0; % domain left/bottom end point
x_b = 1; % domain right/top end point
area = 1; % area of rectangle

H_level = 4; % coarse mesh size level u (LOD)
h_level = 6; % fine mesh size level u (LOD)
ht_level = 4; % mesh size level A (P2)
ell = 3; % oversampling of ell-layers for u (LOD)

tau = 1; % step size for gradient flow
tol = 10^(-10); % tolerance for termination
i_max = 1000000; % maximum number of iterations
update_para1 = 10; % update corrector in first 'update_para1' steps
update_para2 = 100; % update corrector after each 'update_para2' steps

save_path = "./solution_kappa"+kappa+".mat"; % path for save file

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


%% initial values
% specify initial values here
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1. e.g. default values u = 0.8 + 0.6i and A = 0
u_h = (0.8 + 0.6i)*ones(size(M_h,1),1);
A_ht = zeros(size(S_curl,1),1);

% or

% 2. best approximation from file of previous calculation
%path_best = "./solution"+kappa+".mat";
%u_h = getBestApproximation_for_u(path_best,x_a,x_b,kappa,Stiff_h,M_h,h_level,h_ref_level);
%A_ht = getBestApproximation_for_A(path_best,x_a,x_b,ht_level,h_ref_level,Bdx_ht,Bdy_ht,M_A,S_curl,S_div);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

A1_h = Bdx_h*P2'*Bdx_ht'*A_ht(1:Nx_ht); % prolongation
A2_h = Bdy_h*P2'*Bdy_ht'*A_ht(Nx_ht+1:end); % prolongation
A_h = [A1_h; A2_h];


%% corrector computation and assembling
S_h = assemble_bilinear_form_with_P2(A_h,kappa,T_h.t,T_h_P2.t,T_h.p,nodes2mesh_hx,nodes2mesh_hy);
Fu_h = assemble_nonlinear_term_impl(u_h,T_h.t,T_h.p,nodes2mesh_h);

beta = 0;
% compute LOD-corrector (for parallization use 'getCorrectorMatrixParallel_with_P2' in l.132)
disp("compute corrector")
Q = getCorrectorMatrix_with_P2(T_H,T_h,T_h_P2,patches,A_h,kappa,beta,S_h,M_h,P1',P0,B_H,B_h,nodes2mesh_hx,nodes2mesh_hy);
%Q = getCorrectorMatrixParallel_with_P2(T_H,T_h,T_h_P2,patches,A_h,kappa,beta,S_h,M_h,P1',P0,B_H,B_h,nodes2mesh_hx,nodes2mesh_hy);

% assemble LOD matrices
M_LOD = (P1 + Q)*M_h*(P1 + Q)';
S_LOD = (P1 + Q)*S_h*(P1 + Q)';
Fu_LOD = (P1 + Q)*Fu_h*(P1 + Q)';

% initial value in LOD space
u_LOD = M_LOD'\((P1+Q)*M_h'*u_h);

% compute intial matrices for A-equation
% get u on ht mesh
u_ht = M_ht'\(Pt*M_h'*u_h);
FA_ht = assemble_nonlinear_term_A_density_P2(u_ht,T_ht_P2.t,T_ht.t,T_ht_P2.p,nodes2mesh_ht,nodes2mesh_htx,nodes2mesh_hty);
RHS_u = assemble_RHS_u_nabla_u_P2(u_ht,T_ht_P2.t,T_ht.t,T_ht_P2.p,nodes2mesh_ht,nodes2mesh_htx,nodes2mesh_hty);

% compute energy
[curlA,E_curl] = get_curl_of_A_P2(H_mag,A_ht,T_ht_P2.t,T_ht_P2.p,nodes2mesh_htx,nodes2mesh_hty);
%E_curl = get_energy_curl_contribution(H_mag,curlA,T,Nd);
E = real(u_LOD'*S_LOD*u_LOD + 0.5*area - u_LOD'*M_LOD*u_LOD + 0.5*u_LOD'*Fu_LOD*u_LOD) + E_curl;


%% computation of gradient flow
delta = tol;
counter = 0;

while delta >= tol && counter < i_max
    u_old = u_h;
    A_old = A_ht;
    E_old = E;

    % solve u equation
    u_LOD = (M_LOD + tau*(S_LOD + Fu_LOD - M_LOD))\(M_LOD*u_LOD);
    u_h = (P1 + Q)'*u_LOD;

    % solve A equation
    A_ht = (M_A + tau*(S_div + S_curl + FA_ht))\(M_A*A_ht + RHS_H_mag - (0.5i/kappa)*RHS_u);

    % compute energy (idenpendent on update)
    [curlA,E_curl] = get_curl_of_A_P2(H_mag,A_ht,T_ht_P2.t,T_ht_P2.p,nodes2mesh_htx,nodes2mesh_hty);
    E = real(u_LOD'*S_LOD*u_LOD + 0.5*area - u_LOD'*M_LOD*u_LOD + 0.5*u_LOD'*Fu_LOD*u_LOD) + E_curl;

    delta = abs(E_old - E);
    counter = counter + 1;
    disp(counter)
    disp(delta)
    disp(E)


    %% corrector update
    if delta >= tol && counter < i_max %only if necessary
        % matricies on fine mesh (u-equation)
        % get A on fine mesh
        A1_h = Bdx_h*P2'*Bdx_ht'*A_ht(1:Nx_ht); % prolongation
        A2_h = Bdy_h*P2'*Bdy_ht'*A_ht(Nx_ht+1:end); % prolongation
        A_h = [A1_h; A2_h];

        S_h = assemble_bilinear_form_with_P2(A_h,kappa,T_h.t,T_h_P2.t,T_h.p,nodes2mesh_hx,nodes2mesh_hy);
        Fu_h = assemble_nonlinear_term_impl(u_h,T_h.t,T_h.p,nodes2mesh_h);

        if mod(counter,update_para2)==0 || (counter <= update_para1)
            % compute LOD-corrector (for parallization use 'getCorrectorMatrixParallel_with_P2' in l.196)
            disp("compute corrector")
            Q = getCorrectorMatrix_with_P2(T_H,T_h,T_h_P2,patches,A_h,kappa,beta,S_h,M_h,P1',P0,B_H,B_h,nodes2mesh_hx,nodes2mesh_hy);
            %Q = getCorrectorMatrixParallel_with_P2(T_H,T_h,T_h_P2,patches,A_h,kappa,beta,S_h,M_h,P1',P0,B_H,B_h,nodes2mesh_hx,nodes2mesh_hy);

            % u in new LOD space
            u_LOD = M_LOD'\((P1+Q)*M_h'*u_h);
            M_LOD = (P1 + Q)*M_h*(P1 + Q)';
        end

        % assemble LOD-matricis
        S_LOD = (P1 + Q)*S_h*(P1 + Q)';
        Fu_LOD = (P1 + Q)*Fu_h*(P1 + Q)';

        % compute intial matricies for A-equation
        % get u on ht mesh
        u_ht = M_ht'\(Pt*M_h'*u_h);
        FA_ht = assemble_nonlinear_term_A_density_P2(u_ht,T_ht_P2.t,T_ht.t,T_ht_P2.p,nodes2mesh_ht,nodes2mesh_htx,nodes2mesh_hty);
        RHS_u = assemble_RHS_u_nabla_u_P2(u_ht,T_ht_P2.t,T_ht.t,T_ht_P2.p,nodes2mesh_ht,nodes2mesh_htx,nodes2mesh_hty);

    end
end


%% save
if save_results == true
    save(save_path,'-v7.3')
end


%% plot
if plot_results == true
    plot_solution(u_h,A_ht,curlA,T_h,T_ht,T_ht_P2,nodes2mesh_htx,nodes2mesh_hty)
end



