function [alpha] = sideload_coef_cohesive(u_bar,c_bar,h_water,h_sw,m,k,C_f,rho_w,rho_s, nu,C_internal,g,B,B_flag)
% alpha [dimensionless]: is the ratio of the cohesive sediment
% concentration in side channel to the main channel
% u_bar [m/s]: velocity in the main channel
% k:  cohesive sediment coefficent with appropriate dimension page 87 of
% "Dynamics of estuine mud" defult value is 0.001 in SI unit
% m:  cohesive sediment coefficent with appropriate dimension page 87 of
% "Dynamics of estuine mud" default value is 1 in SI unit
% cbar_main [kg/m3]: is the depth average concentration in the main channel
% h_sw [m]: is the height of side weir or levee
% h_water [m]: is the height of free water surface in the channel
% this subroutine only works for dilute mixtures of water and cohesive
% sediment when mass concentration of suspended mud is less than 3 kg/m^3 



%% Calculating W_50b
rho_e = rho_w + C_internal*(rho_s-rho_w);
length_scale = sqrt(19.8*rho_w*nu*(rho_s^m)*k/(g*(rho_e-rho_w)));
de = length_scale*c_bar^(0.5*m);
d_star = de * (g*(rho_e-rho_w)./(rho_w*nu^2))^(1/3);
w_50b = (nu/de)*(sqrt(10.36^2+1.049*((1-C_f)^4.7)*d_star^3)-10.36);
%% calculating B
if B_flag == 0
    B = m*w_50b/(0.0025*u_bar);
end
%% alpha
y = h_water - h_sw;
alpha = (1-B*y/2/m/h_water+B*B*(m+1)*y*y/(6*m*m*h_water*h_water))/(1-B/2/m+B*B*(m+1)/6/m/m);

return