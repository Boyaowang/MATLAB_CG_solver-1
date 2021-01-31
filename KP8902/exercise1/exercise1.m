clc;
clear all;

global supvel molmass ptot rhob pb0 Cp enthalpy U dt Tr R;

supvel = 1.0; %[m/s]
molmass = 29.48; %[kg/kmole] 
ptot = 1.0e5; %[Pa] 
rhob = 1300; %[kg/m^3]
pb0 = 0.211e5; %[Pa]
Cp = 0.992; %[kJ/kg*K] 
enthalpy = 1285409.0; %[kJ/kmole]
U = 0.096; %[kJ/m^2*s]
dt = 2.54e-2; %[m]
Tr = 625.0; %[K]
R = 8314.5; %[J/kmole*K]
zstart = 0; %[m]
zend = 3; %[m]
pA0 = 0.015e5; %[Pa]
T0 = 625; %[K]

%integration span
zspan=[zstart zend];

%initial conditions
y0=[pA0 T0];

[z,y]=ode15s(@yderiv,zspan,y0);

%plot the result

m = 2;
n = 1;
nr = 1;


subplot(m,n,nr);
plot(z,y(:,1))
title('Partial pressure profile')
xlabel('z [m]') 
ylabel('p [Pa]')

subplot(m,n,2);
plot(z,y(:,2))
title('Temperature profile')
xlabel('z [m]') 
ylabel('T [K]')



