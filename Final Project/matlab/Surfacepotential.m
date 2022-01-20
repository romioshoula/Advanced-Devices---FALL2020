%%Rami wail - Ahmed Alaa Mustafa Ali-Mohamed Hatem
%surface potential numerical calculation 
clear all;
close all;
clc;
%% parameters 

Tox=1.7*10^-7;
Tsi=1*10^-6;
Vg=0.9;
Vd=1.5;
L=3*10^-5;
Vs=0;
K=1.38*10^-23;
T=300;
vt=26*10^-3;
epsilon=8.85*10^-14;
es=11.7;
eox=3.9;
Ni=1.5*10^10;
Eg=1.12;
q=1.6*10^-19;
Nd=10^19;
%% calculations 
Cox=eox*epsilon/Tox
Phi_FP=vt*log(Nd/Ni)
Xdep=sqrt(4*es*epsilon*Phi_FP/(q*Nd))
Phi_MS=-(Eg/2+Phi_FP)
V_FB= -Phi_MS
Qsd_max = q*Nd*Xdep
Vth=Qsd_max/Cox + Phi_MS + 2*Phi_FP
m=1+3*Tox/Xdep
temp=(Vg-Vth)/m
Vys=0;
Vyd=temp-sqrt(temp^2-2*temp*Vd+Vd^2)

%syms Phi_O
Phi_O=0.1;
Phi_S=Phi_O-2*vt*log(   cos(    sqrt(q*Ni/(2*vt*es*epsilon))  *    Tsi/2*exp( abs(Phi_O-Vys)/vt)  ))

 %Phi_O=0.1;
%% numerical calculations  
f= -Vg + V_FB + Phi_S + Tox*es/eox*sqrt(    2*vt*Ni*q/(es*epsilon)*(   exp((Phi_S-Vys)/vt) -   exp((Phi_O-Vys)/vt) )   )
Phi_O=Phi_O+0.001
f_pos= -Vg + V_FB + Phi_S + Tox*es/eox*sqrt(    2*vt*Ni*q/(es*epsilon)*(   exp((Phi_S-Vys)/vt) -   exp((Phi_O-Vys)/vt) )   )
Phi_O=Phi_O-2*0.001
f_neg= -Vg + V_FB + Phi_S + Tox*es/eox*sqrt(    2*vt*Ni*q/(es*epsilon)*(   exp((Phi_S-Vys)/vt) -   exp((Phi_O-Vys)/vt) )   )
Phi_O=Phi_O+0.001

Phi_O=Phi_O-f/(f_pos+f_neg)/0.001;

f2dash=- 2*log(cos((Tsi*exp(abs(Phi_O - Vys)/vt)*((Ni*q)/(2*epsilon*es*vt))^(1/2))/2)) - (2*vt*sin((Tsi*exp(abs(Phi_O - Vys)/vt)*((Ni*q)/(2*epsilon*es*vt))^(1/2))/2)*((Tsi*exp(abs(Phi_O - Vys)/vt)*abs(Phi_O - Vys)*((Ni*q)/(2*epsilon*es*vt))^(1/2))/(2*vt^2) + (Ni*Tsi*q*exp(abs(Phi_O - Vys)/vt))/(8*epsilon*es*vt^2*((Ni*q)/(2*epsilon*es*vt))^(1/2))))/cos((Tsi*exp(abs(Phi_O - Vys)/vt)*((Ni*q)/(2*epsilon*es*vt))^(1/2))/2) - (2^(1/2)*Tox*es*((Ni*q*(exp((Phi_O - Vys)/vt) - exp(-(Vys - Phi_O + 2*vt*log(cos((Tsi*exp(abs(Phi_O - Vys)/vt)*((Ni*q)/(2*epsilon*es*vt))^(1/2))/2)))/vt)))/(epsilon*es) - (Ni*q*vt*(exp(-(Vys - Phi_O + 2*vt*log(cos((Tsi*exp(abs(Phi_O - Vys)/vt)*((Ni*q)/(2*epsilon*es*vt))^(1/2))/2)))/vt)*((Vys - Phi_O + 2*vt*log(cos((Tsi*exp(abs(Phi_O - Vys)/vt)*((Ni*q)/(2*epsilon*es*vt))^(1/2))/2)))/vt^2 - (2*log(cos((Tsi*exp(abs(Phi_O - Vys)/vt)*((Ni*q)/(2*epsilon*es*vt))^(1/2))/2)) + (2*vt*sin((Tsi*exp(abs(Phi_O - Vys)/vt)*((Ni*q)/(2*epsilon*es*vt))^(1/2))/2)*((Tsi*exp(abs(Phi_O - Vys)/vt)*abs(Phi_O - Vys)*((Ni*q)/(2*epsilon*es*vt))^(1/2))/(2*vt^2) + (Ni*Tsi*q*exp(abs(Phi_O - Vys)/vt))/(8*epsilon*es*vt^2*((Ni*q)/(2*epsilon*es*vt))^(1/2))))/cos((Tsi*exp(abs(Phi_O - Vys)/vt)*((Ni*q)/(2*epsilon*es*vt))^(1/2))/2))/vt) + (exp((Phi_O - Vys)/vt)*(Phi_O - Vys))/vt^2))/(epsilon*es)))/(2*eox*(-(Ni*q*vt*(exp((Phi_O - Vys)/vt) - exp(-(Vys - Phi_O + 2*vt*log(cos((Tsi*exp(abs(Phi_O - Vys)/vt)*((Ni*q)/(2*epsilon*es*vt))^(1/2))/2)))/vt)))/(epsilon*es))^(1/2));

 Phi_O=Phi_O-f/f2dash;
