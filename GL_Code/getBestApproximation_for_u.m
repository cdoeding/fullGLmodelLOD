function u0 = getBestApproximation_for_u(path,x_a,x_b,kappa,Stiff_h,M_h,h_level,h_ref_level)
%GETBESTAPPROXIMATION Summary of this function goes here
%   Detailed explanation goes here

%kappa_str = num2str(kappa,'%i');
%path = strcat("./reference_solution_kappa",kappa_str,".mat");

%% get u
z_u = load(path,"u_h","M_h", "T_h", "nodes2mesh_h");

u_ref = z_u.u_h;
M_h_ref = z_u.M_h;
T_h_ref = z_u.T_h;
nodes2mesh_h_ref = z_u.nodes2mesh_h;
Stiff_h_ref = assemble_stiffness_matrix(T_h_ref.t,T_h_ref.p,nodes2mesh_h_ref);

[~,~,P1,~] = getCoarseFineTriangulation(x_a,x_b,h_level,h_ref_level);
P = P1';

u0 = ((1/kappa^2)*Stiff_h + M_h)'\(P*((1/kappa^2)*Stiff_h_ref + M_h_ref)'*u_ref);

end

