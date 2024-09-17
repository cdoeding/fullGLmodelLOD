function A_on_u = get_A_on_u_mesh(A,h_level,ht_level,P,M_A,M_A_h,Bdx_h,Bdy_h,Bdx_ht,Bdy_ht,nodes2mesh_htx)
%GET_U_ON_A_MESH Summary of this function goes here
%   Detailed explanation goes here

Nx = sum(logical(nodes2mesh_htx));

if h_level >= ht_level
    % polongation
    A1 = Bdx_h*P'*Bdx_ht'*A(1:Nx); % polongation
    A2 = Bdy_h*P'*Bdy_ht'*A(Nx+1:end); % polongation
    A_on_u = [A1; A2];

else
    % L2-projection
    A_mass = M_A'*A;
    A1 = Bdx_h*P*Bdx_ht'*A_mass(1:Nx);
    A2 = Bdy_h*P*Bdy_ht'*A_mass(Nx+1:end);
    A_on_u = M_A_h\[A1; A2];

end
