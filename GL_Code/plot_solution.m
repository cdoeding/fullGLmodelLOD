function plot_solution(u_h,A_ht,curlA,T_h,T_ht,T_ht_P2,nodes2mesh_htx,nodes2mesh_hty)
%PLOT_SOLUTION Summary of this function goes here
%   Detailed explanation goes here

model = createpde();
geometryFromMesh(model,T_h.p',T_h.t');

%% plot
font = 20;
cmap = flipud(hot);
figure(1)
h = gcf;
h.Position = [400 400 2*560 2*420];
hold on
clf;

%% plot u (1)
subplot(2,2,1)
pdeplot(model,"XYData",abs(u_h), "ZData",abs(u_h),'colormap',cmap,'Mesh','off','FaceAlpha', 1)
clim([0 1])
xlabel("$x$",'Interpreter', 'Latex', 'Fontsize', font)
ylabel("$y$",'Interpreter', 'Latex', 'Fontsize', font)
title("$|u_h|$", 'Interpreter', 'Latex', 'Fontsize', 20)
ax = gca;
ax.FontSize = font;

%% plot u (2)
subplot(2,2,2)
pdeplot(model,"XYData",abs(u_h), "ZData",abs(u_h),'colormap',cmap,'Mesh','off','FaceAlpha', 1)
view(0,90)
clim([0 1])
xlabel("$x$",'Interpreter', 'Latex', 'Fontsize', font)
ylabel("$y$",'Interpreter', 'Latex', 'Fontsize', font)
title("$|u_h|$", 'Interpreter', 'Latex', 'Fontsize', 20)
ax = gca;
ax.FontSize = font;

%% plot A
A1 = zeros(size(T_ht_P2.p,1),1);
A2 = zeros(size(T_ht_P2.p,1),1);

Nx = sum(logical(nodes2mesh_htx));
for j = 1:size(T_ht_P2.p,1)
    if nodes2mesh_htx(j) ~= 0
        A1(j) = A_ht(nodes2mesh_htx(j));
    end

    if nodes2mesh_hty(j) ~= 0
        A2(j) = A_ht(Nx + nodes2mesh_hty(j));
    end
end

subplot(2,2,3)
quiver(T_ht_P2.p(:,1),T_ht_P2.p(:,2),A1,A2)
xlabel("$x$",'Interpreter', 'Latex', 'Fontsize', font)
ylabel("$y$",'Interpreter', 'Latex', 'Fontsize', font)
title('$A$','Interpreter', 'Latex','Fontsize', font)
ax = gca;
ax.FontSize = font;

%% plot curlA
subplot(2,2,4)
plot_curl(curlA,T_ht.t,T_ht.p)
colormap(cmap);
xlabel("$x$",'Interpreter', 'Latex', 'Fontsize', font)
ylabel("$y$",'Interpreter', 'Latex', 'Fontsize', font)
title('$\mathrm{curl} A$','Interpreter', 'Latex','Fontsize', font)
ax = gca;
ax.FontSize = font;

end

