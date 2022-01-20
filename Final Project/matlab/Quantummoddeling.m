%%Rami wail - Ahmed Alaa Mustafa Ali-Mohamed Hatem
%surface potential numerical calculation 
clear all;
close all;
clc;
%%parameters
q=1.6*10^-19;
h=6.62606896*10^-34;
mo=9.10938215*10^-31;
ml=0.98*mo;
mt=0.19*mo;
pi=3.14159;
h_p=h/(2*pi);
tsi=10^-8 ;
vt=0.026;

%descrite energy levels in the 1D 

E00=    h^2 /   (4*2*ml*tsi^2)/q
E01=(0+1)^2*pi^2*h_p^2/(2*mt*tsi^2)/q
E10=(1+1)^2*pi^2*h_p^2/(2*ml*tsi^2)/q
E11=(1+1)^2*pi^2*h_p^2/(2*mt*tsi^2)/q


%number of k states
N2d00=q*2*4*pi*vt*sqrt(mt*mt)/h^2
N2d01=q*4*4*pi*vt*sqrt(ml*mt)/h^2
N2d10=q*2*4*pi*vt*sqrt(mt*mt)/h^2
N2d11=q*4*4*pi*vt*sqrt(ml*mt)/h^2

%charge per unit area 
Q00=q*N2d00*log(1+exp((0.08-E00)/vt))
Q01=q*N2d00*log(1+exp((0.08-E01)/vt))
Q10=q*N2d00*log(1+exp((0.08-E10)/vt))
Q11=q*N2d00*log(1+exp((0.08-E11)/vt))

