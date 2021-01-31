function dydz=yderiv(z,y)
global supvel molmass ptot rhob pb0 Cp enthalpy U dt Tr R;
pA = y(1); 
T = y(2);
k=1.0e-10*exp(19.837-13636.0/T)/3600.0; 
r=k*pb0*pA; 
rhog=ptot*molmass/(R*T);
dTdz=1.0/(supvel*rhog*Cp)*(enthalpy*rhob*r-4.0*U/dt*(T-Tr));
dpdz=-1.0/supvel*(molmass*ptot*rhob/rhog*r)+pA/T*dTdz;
dydz=[dpdz;dTdz];
end


