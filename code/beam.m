function [b_moms_def,fq] = beam(twists,pitchrad,pn,pt,stiffs_mass_cg,data,nu)
% credit for this source code goes to:
% Matias Sessarego, M.Sc., msessare@ucalgary.ca, msessare@gmail.com
%Organizing Turbine Data
x=nu*data.R;
beta=twists;
v=stiffs_mass_cg(:,3);
beta_v=(beta+v);

%Mass
m=stiffs_mass_cg(:,4);

%Stiffness about Principal Axes
EI1=stiffs_mass_cg(:,1);
EI2=stiffs_mass_cg(:,2);

%Loads about Tip Chord
py=(pt.*cos(pitchrad))+(pn.*sin(pitchrad));
pz=(pn.*cos(pitchrad))-(pt.*sin(pitchrad));
N=length(py);

%Calculate deflection----------------------------------
[beam_out] = deflect(py,pz,EI1,EI2,beta_v,N,x,pitchrad);

%Calculate centrifugal force Fx------------------------
Fx=zeros(N,1);
if data.loadcase==3
    omega=data.rotoroverspeed*data.rat_rpm*(pi/30);
    Fx=abs(trapz((x),m*(omega^2).*x)-cumtrapz((x),m*(omega^2).*x));
end

%Calculate Eigenmodes----------------------------------
pIni=zeros(N,1);
pIni(end)=sqrt(0.5);
%1st flapwise eigenmode
wSQf1=0;
pyf1=pIni; pzf1=pIni;
itermax=1000000;
for iter=1:1:itermax
    wSQf1old=wSQf1;
    def=deflect(pyf1,pzf1,EI1,EI2,beta_v,N,x,pitchrad);
    uyf1=def(:,1);  uzf1=def(:,2);
    wSQf1=(norm(pzf1(end))/norm(uzf1(end)*m(end)));
    deno=(((uyf1(end)^2)+(uzf1(end)^2)).^0.5);
    pyf1=(wSQf1*m.*uyf1)./deno;
    pzf1=(wSQf1*m.*uzf1)./deno;
    tol1=abs(wSQf1-wSQf1old);
    if (tol1<1e-8)
        break 
    end
end
wf1=wSQf1^0.5;  Uyf1=uyf1./deno;    Uzf1=uzf1./deno;

%Output results----------------------------------------
M1=beam_out(:,3);   M2=beam_out(:,4);   
un=beam_out(:,5);   ut=beam_out(:,6);  
K1=beam_out(:,7);   fq=wf1;

b_moms_def=cat(2,M1,M2,un,ut,-K1,Uyf1,Uzf1,Fx);
end

%Calculation of bending moments and deflections
function [beam_out] = deflect(py,pz,EI1,EI2,beta_v,N,x,pitchrad)
%Shear about Tip Chord
Ty=zeros(N,1);
Tz=zeros(N,1);
for i=N:-1:2
    Ty(i-1)=Ty(i)+(0.5*(py(i-1)+py(i))*(x(i)-x(i-1)));
    Tz(i-1)=Tz(i)+(0.5*(pz(i-1)+pz(i))*(x(i)-x(i-1)));
end

%Moment about Tip Chord
My=zeros(N,1);
Mz=zeros(N,1);
for i=N:-1:2  
    My(i-1)=My(i)-(Tz(i)*(x(i)-x(i-1)))-(((1/6*pz(i-1))+(1/3*pz(i)))*((x(i)-x(i-1))^2));
    Mz(i-1)=Mz(i)+(Ty(i)*(x(i)-x(i-1)))+(((1/6*py(i-1))+(1/3*py(i)))*((x(i)-x(i-1))^2));
end

%Moment about Principal Axes
M1=(My.*cos(beta_v))-(Mz.*sin(beta_v));
M2=(My.*sin(beta_v))+(Mz.*cos(beta_v));

%Curvatures
K1=(M1./EI1);
K2=(M2./EI2);

Kz=(-K1.*sin(beta_v))+(K2.*cos(beta_v));
Ky=(K1.*cos(beta_v))+(K2.*sin(beta_v));

%Angular Deformation
ad_y=zeros(N,1);
ad_z=zeros(N,1);
for i=1:1:(N-1)
    ad_y(i+1)=ad_y(i)+(0.5*(Ky(i+1)+Ky(i))*(x(i+1)-x(i)));
    ad_z(i+1)=ad_z(i)+(0.5*(Kz(i+1)+Kz(i))*(x(i+1)-x(i)));
end

%Deflection
uy=zeros(N,1);
uz=zeros(N,1);
for i=1:1:(N-1)
    uy(i+1)=uy(i)+(ad_z(i)*(x(i+1)-x(i)))+(((1/6*Kz(i+1))+(1/3*Kz(i)))*((x(i+1)-x(i))^2));
    uz(i+1)=uz(i)-(ad_y(i)*(x(i+1)-x(i)))-(((1/6*Ky(i+1))+(1/3*Ky(i)))*((x(i+1)-x(i))^2));
end

un=(-uy.*sin(pitchrad))+(uz.*cos(pitchrad));
ut=(uy.*cos(pitchrad))+(uz.*sin(pitchrad));

beam_out=cat(2,uy,uz,M1,M2,un,ut,K1,ad_y,ad_z,Ky,Kz);
end