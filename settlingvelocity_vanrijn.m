function [wf] = settlingvelocity_vanrijn(rhom,diam,rho,viscocine)
%% This function calculates settling velocity of grains
% Based on the method suggestion by Leo Van Rijn:  
% Van Rijn, L.C. (1993). Principles of sediment transport in rivers,
% estuaries and coastal seas. Aqua Publications.
% Implemented by Kaveh Zamani at UNSW Sydney, School of  of Civil and Environmental Engg.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ___________________________________________________________________
% This routine finds the value of sediment fall velocity.
%
% 1. Inputs:
% rhom = mass density of sediment grains (kg/m3)
% diam = sediment grain diameter (m).
% rho = mass density of seawater (kg/m3).
% viscocine = kinematic viscosity (m2/s).
%
% 2. Outputs:
% wf = settling or fall velocity     (m/s)
% 
% 3. Example: 
%  [wf] = settlingvelocity_vanrijn(2650,0.2E-3,1025,1E-6)
% 
% wf = 
%   0.0238750897450597
%% 


dtr = ( ( rhom - rho ) / rho );
S =  ( ( diam / ( 4 * viscocine) ) * sqrt( ( (rhom / rho) - 1 ) * 9.81 * diam  ) );     % (dimensionless)

if S >= 300
            wf = 1.82 * sqrt( ( (rhom / rho) - 1 ) * 9.81 * diam );
elseif S  <= 0.8
            wf = ( 2 / 9 ) *sqrt( ( (rhom / rho) - 1 ) * 9.81 * diam ); 
else
            Dasterisk = ( ( 9.81 * ( (rhom / rho) -1 ) )/ viscocine ^2 ) ^ ( 1 / 3 ) *diam;
            Dcubic = Dasterisk ^3;
             Arquimids = ( dtr * 9.81 * diam ^ 3 ) / viscocine ^2;
             Cl = 0.055 * tanh( ( 12 * Arquimids ^-0.59 ) * exp( -0.0004 * Arquimids ) );
             Ct = 1.06 * tanh( ( 0.016 * Arquimids ^ 0.5 ) * exp( -120 / Arquimids ) );
             wf = ( ( Cl * dtr * 9.81 * diam ^2 ) / viscocine ) + ( Ct * sqrt( dtr * 9.81 * diam ) );

end
return 