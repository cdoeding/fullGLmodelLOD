function plot_curl(curlA,T,Nd)
%PLOT_CURL Summary of this function goes here
%   Detailed explanation goes here

tri_x = zeros(size(T,1),1);
tri_y = zeros(size(T,1),1);

for i = 1:size(T,1)
    z1 = Nd(T(i,1),:); %coordinates of 1st triangle node
    z2 = Nd(T(i,2),:); %coordinates of 2nd triangle node
    z3 = Nd(T(i,3),:); %coordinates of 3rd triangle node
    
    tri_x(i) = (1/3)*(z1(1)+z2(1)+z3(1));
    tri_y(i) = (1/3)*(z1(2)+z2(2)+z3(2));
end

scatter3(tri_x,tri_y,curlA,50,curlA,'filled')

end

