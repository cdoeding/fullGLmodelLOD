function [T_H_P2,T_h_P2,P_P2] = getCoarseFineTriangulation_for_P2(T_H,T_h,P0)
%GETCOARSEFINETRIANGULATION Summary of this function goes here
%   Detailed explanation goes here

%% Corse mesh
% triangles
tri_P1 = T_H.t;
tri_P2 = zeros(size(tri_P1,1),5);
tri_P2(:,[1 3 5]) = tri_P1;

% nodes
Nd_P1 = T_H.p;
Nd_P2 = T_H.p;
Nx = size(Nd_P1,1);

% generate edges mid-points
edges = Tri2Edge(T_H.p',T_H.t');
[xmid,ymid] = EdgeMidPoints(T_H.p',edges,T_H.t');
edges = edges(:,[3 1 2]) + Nx*ones(size(edges));
tri_P2(:,[2 4 6]) = edges;
Nd_P2(Nx+1:Nx+length(xmid),:) = [xmid', ymid'];

T_H_P2 = struct('t',tri_P2,'p',Nd_P2);

%% Fine mesh
% triangles
tri_P1 = T_h.t;
tri_P2 = zeros(size(tri_P1,1),5);
tri_P2(:,[1 3 5]) = tri_P1;

% nodes
Nd_P1 = T_h.p;
Nd_P2 = T_h.p;
Nx = size(Nd_P1,1);

% generate edges mid-points
edges = Tri2Edge(T_h.p',T_h.t');
[xmid,ymid] = EdgeMidPoints(T_h.p',edges,T_h.t');
edges = edges(:,[3 1 2]) + Nx*ones(size(edges));
tri_P2(:,[2 4 6]) = edges;
Nd_P2(Nx+1:Nx+length(xmid),:) = [xmid', ymid'];

T_h_P2 = struct('t',tri_P2,'p',Nd_P2);

%% get prolongation
P_P2 = getProlongationP2(T_H_P2,T_h_P2,P0);

end

