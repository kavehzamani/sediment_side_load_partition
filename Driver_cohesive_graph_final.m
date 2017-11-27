clc; clear all; close

u = 0.5:0.1:5;
H = 1;
c_bar = 1;
m = 1.07;
k = 0.00099;
Cf = 0.032;
hw = (0.4:0.02:0.96)*H;
rho_w = 1027;
rho_s= 2650;
B = 123456789;
nu = 1.36e-6;
b_flag = 0;
g =9.81;
c_internal = 0.03;
 nu = length (u);
 nh = length(hw);
 
 alpha = zeros(nu,nh);
 
 for i=1:nu
     for j=1:nh
      alpha(i,j)=sideload_coef_cohesive(u(i),c_bar,H,hw(j),m,k,Cf,rho_w,rho_s,nu,c_internal,g,B,b_flag);
     end 
 end
 
 
surf (hw,u,alpha)
title(['k:',num2str(k),', Cf:',num2str(Cf),', m:',num2str(m)],'FontSize',14)
ylabel('Velocity (m/s)','FontSize',6)
xlabel('Side weir height ratio','FontSize',16)
zlabel('\alpha','FontSize',16)




