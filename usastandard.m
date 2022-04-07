
 function [p,rho,T]=usastandard(z)

p0=1.01325*10^5 % Initial pressure
T0=288.16; % temperature in Kelvin
a=6.5*10^-3 % slope coefficient
g=9.81 % gravity acceleration m/s²
R=287.17 % Universal gas constant for atmospheric air
z1=0 % z1 altitude beginning of zone,
z2=11000%  z2 altitude to be computed
dz=z2-z1 

%Pressure calculation
p=p0*((T0-a*dz)/T0)^(g/(a*R))% N/m²
  
%Density calculation
rho0=p0/(R*T0)% Kg/m³
rho=rho0*((T0-a*dz)/T0)^((g/(a*R))-1)% Kg/m³

% Temperature calculation
T=T0-a*dz

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%      Flight conditions at 30000 m for the waverider 14-XB     %%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

a20=1*10^-3 % following US standard atmosphere
p11=p
p20=5474.875 % Pa
T11=T
T20=T11

z=30000
z20=20000
Dz=z30-z20

% Temperature calculation
T30=T20+a20*Dz

rho11=rho
rho20=p20/(R*T20)

%Pressure calculation
p30=p20*(T20/(T20+a20*Dz))^(g/(a20*R))% N/m²
rho30=rho20*((T20-a20*Dz)/T20)^((g/(a20*R))-1)% Kg/m³

%Thermo-physical calculations

mi= 1.458*10^-6*(T30)^(3/2)/(T30+110.4)  % kg/s.m
k= (1.99*10^-3)*((T30)^(3/2))/(T30+112) % W/mK

p=p30;
rho=rho30;
T=T30;

 end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%