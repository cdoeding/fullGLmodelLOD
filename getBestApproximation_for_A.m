function A0 = getBestApproximation_for_A(path,x_a,x_b,ht_level,h_ref_level,Bdx_ht,Bdy_ht,M_A_ht,S_ht_curl,S_ht_div)
%GETBESTAPPROXIMATION Summary of this function goes here
%   Detailed explanation goes here

%kappa_str = num2str(kappa,'%i');
%path = strcat("./reference_solution_kappa",kappa_str,".mat");

%% get A
% load data
z_A = load(path, ...
    "A_ht", ...
    "S_curl", ...
    "S_div", ...
    "M_A", ...
    "nodes2mesh_htx", ...
    "T_ht", ...
    "T_ht_P2");

A_ref = z_A.A_ht;
S_curl_ref = z_A.S_curl;
S_div_ref = z_A.S_div;
M_A_ref = z_A.M_A;
T_h_ref = z_A.T_ht;
T_h_ref_P2 = z_A.T_ht_P2;

% Polongation operator from ref
[T_ht,~,~,Pt0] = getCoarseFineTriangulation(x_a,x_b,ht_level,h_ref_level);
[~,~,P2] = getCoarseFineTriangulation_for_P2(T_ht,T_h_ref,Pt0);
P = P2;

% mesh
Nx = sum(logical(z_A.nodes2mesh_htx));
[Bx_h_ref,By_h_ref] = getBoundaryNodes_for_A(T_h_ref_P2.p,x_a,x_b,'non-natural');
Bdx_h_ref = getBoundaryRestriction(Bx_h_ref);
Bdy_h_ref = getBoundaryRestriction(By_h_ref);

% project
A_elliptic = (S_curl_ref + S_div_ref + M_A_ref)'*A_ref;
A1 = Bdx_ht*P*Bdx_h_ref'*A_elliptic(1:Nx);
A2 = Bdy_ht*P*Bdy_h_ref'*A_elliptic(Nx+1:end);
A0 = (S_ht_curl + S_ht_div + M_A_ht)\[A1; A2];

end

