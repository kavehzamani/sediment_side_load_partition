function [alpha] = sideload_coef_beta_shear(k_s,diam,cbar_main,u,h_levee,h_water,flag_shear,beta)

% alpha [dimensionless]: is the ratio of the sediment concentration in side channel to the
% k_s is the roughness height
% depth average sediment concentration in the cell next to the levee or
% weir (first cell of side channel) or alpha = cbar_side/cbar_main
% diam [m]: is the diameter of sediment particle
% cbar_main [m3/m3]: is the depth average concentration in the main channel
% ubar_main [m/s]: is the velocity in the cell next to side weir or shear
% velocity directly based on flag (Flag_shear = 1: main channel velocity; flag-shear=2: shear velocity directly) 
% h_sw [m]: is the height of side weir or levee
% h_water [m]: is the height of free water surface in the channel
% flag_shear (Flag_shear = 1: main channel velocity; flag-shear=2: shear velocity directly) 


%% constants
rho_s = 2650; % mass density of sediment grains (kg/m3)
rho_w = 1010;% mass density of water (kg/m3).
visco_kin = 10E-6; % kinematic viscosity of the water_sediment mixture (m2/s).
kappa = 0.4; %von Karman constant 

%% limitations 
%todo: put limitation on the minimum and maximum of side weir height

%% turbulent shmidt number

    

%% calculate rouse number
% settling velocity
[w_fall] = settlingvelocity_vanrijn(rho_s,diam,rho_w,visco_kin);
 
% shear velocity Ian King formula for RMA11 
if flag_shear ==1
    u_star = kappa*u/(log(12.27*h_water/k_s));
elseif flag_shear==2 
    u_star = u;
end

if nargin<6
    beta = 1; % + 2*(w_fall/u_star)^2; %van Rijn 1984
else
end

% Rouse number
ro_num = w_fall/(beta*kappa*u_star);
%% calculating the side weir coefficent

if (ro_num < 0.75) % washload are well-mixed
    alpha = 1;
elseif (ro_num > 2.5) % bedload particles would not raise more than 10 diameter of itself
    alpha = 0;
else
    fun1 = @(z)(((h_water-z).*k_s./(z.*(h_water-k_s))).^ro_num); % Rouse distribution
    int1 = integral(@(z)fun1(z),k_s,h_levee,'AbsTol',10e-4,'RelTol',10e-4);
    int2 = integral(@(z)fun1(z),h_levee,h_water,'AbsTol',10e-4,'RelTol',10e-4);
    c_sub_a = cbar_main/(int1+int2); % concentration at reference level in the rouse distribution
    cbar_side = c_sub_a*int2; % averge concentration on the top of the side spillway or levee
    %%%% 
    alpha = (h_water/(h_water-h_levee))* cbar_side/cbar_main; 
end

return 