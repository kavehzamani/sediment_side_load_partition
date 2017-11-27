
clc; clear all;clf

cbar_main = 1;
u = 0.01:0.05:2;
diam = 0.02; % in meter
h_water = 1; 
h_levee = 0.2:0.01:0.95; 
n = length (h_levee);
m = length (u);
alpha = zeros(n,m);
flag_shear=2;
ks = 2*max(diam);

for i = 1: n 
    for j = 1: m
        alpha(i,j) = sideload_coef_beta_shear(ks,diam,cbar_main,u(j),h_levee(i),h_water,flag_shear,2);
    end
end

surf (u*ks/0.001,h_levee,alpha)
title(['d_s :',num2str(diam)])
if flag_shear ==1
    xlabel('U_*')
elseif flag_shear ==2
    xlabel('$\frac{u_*K_s}{\nu}$','interpreter','latex')
end
ylabel('h_{levee}')
zlabel('\alpha')


