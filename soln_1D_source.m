clear all
close all

%develop solver for 1D slab SN with isotropic scattering
%uniform mesh
%zero incoming flux BC
%use cross section without fission:
%%%%%CASE 1: USE 1-GROUP CROSS SECTIONS --> SOURCE ITERATION ONLY
%%%%%CASE 2: USE G-GROUP CROSS SECTIONS --> ADD THERMAL ITERATIONS

%source localized somewhere
H=100;

Nz=100;
dz=H/Nz;
z=[dz/2:dz:H-dz/2];
N=4; %sn order

%generation of S_N directions and weights
[mu,weig]=lgwt(N,-1.0,1.0);

%definition of source in the central part of the system
SS=zeros(Nz,1);SS(Nz/4+1:3*Nz/4)=1;

%source iteration
toll=1e-5;
err=1;
it=0;
while err>toll
    it=it+1;
    %transport sweep - positive directions

    %transport sweep - negative directions

    %DEFINITION OF ERROR???
    err

end