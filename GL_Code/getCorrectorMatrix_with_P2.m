function Q = getCorrectorMatrix_with_P2(T_H,T_h,T_P2,patches,A_h,kappa,beta,S_h,M_h,P1,P0,B_H,B_h,nodes2mesh_x,nodes2mesh_y)
%GETCORRECTORMA Summary of this function goes here
%   Detailed explanation goes here

N_h = size(T_h.p,1);
N_H = size(T_H.p,1);
NT_H = size(T_H.t,1);

C_h = P1'*M_h;
Q = sparse(N_H,N_h);

for l = 1:NT_H
    [Rl_H,Rl_h] = getRestriction(T_H,T_h,l,patches,P0,B_H,B_h);
    T = getNode2MeshMatrix(T_H,l);
    S_loc = assembleLocalBilinearForm_with_P2(A_h,kappa,T_h,T_P2,P0,l,nodes2mesh_x,nodes2mesh_y);
    M_loc = assembleLocalMassMatrix(T_h,P0,l);
    
    Nl_h = size(Rl_h,1);
    Nl_H = size(Rl_H,1);
    
    S = Rl_h*(S_h + beta^2*M_h)*Rl_h';
    C = Rl_H*C_h*Rl_h';
    
    rhs = -T*B_H*P1'*(S_loc + beta^2*M_loc)*Rl_h';
    
    Y = S\C';
    S_inv = (C*Y)\speye(Nl_H,Nl_H);
    
    w = zeros(3,Nl_h);
    for i = 1:3
        q = S\rhs(i,:)';
        lambda = S_inv*C*q;
        w(i,:) = (q - Y*lambda)';
    end
    
    Q = Q + T'*sparse(w)*Rl_h;
    
    if mod(l,500)==0
        fprintf(strjoin([{num2str(l,'%i')},{'patches of'},{num2str(NT_H,'%i')},{'computed \n'}]))
    end
end


