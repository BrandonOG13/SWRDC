function [dcp_dcT] = bem(pitch,tsr,twists,chords,data,un_afdata,nu)
% credit for this source code goes to:
% Matias Sessarego
B=data.numBlades;
R=data.R;
relax=data.relax;

belements=length(nu);
lengthtsr=length(tsr);
thetap=pitch';
lengthp=length(thetap);
thetap=repmat(repmat(thetap,1,belements),lengthtsr,1);
thetap=sort(thetap,1);

tsr=repmat(repmat(tsr',1,belements),lengthp,1);
twists=repmat(repmat(twists',lengthtsr,1),lengthp,1);
chords=repmat(repmat(chords',lengthtsr,1),lengthp,1);
nu=repmat(repmat(nu',lengthtsr,1),lengthp,1);

cl_cd=zeros((lengthtsr*lengthp),belements,2);

iter=data.iter;
sigma=((chords*B)./(2*pi*nu*R));

%Step(1) Initialize a and a'
%a=0.25*(2+(pi*tsr.*nu.*sigma)-sqrt(4-(4*pi*tsr.*nu.*sigma)+(pi*((tsr.*nu).^2).*sigma.*(8*twists+pi*sigma))));
a=zeros((lengthtsr*lengthp),belements);
ap=zeros((lengthtsr*lengthp),belements);
anew=a;

for z=1:1:iter
%Step(2) Compute flow angle phi
v1=(1-a);
v2=(tsr.*nu.*(1+ap));
v_ratio=(v1./v2);

phi=atan(v_ratio);
sinphi=sin(phi);
cosphi=cos(phi);

%Step(3) Compute tip loss factor
f=abs((B./(2*sinphi)).*((1-nu)./nu));
pF=(2/pi)*acos(exp(-f));

%Step(4) Compute local angle of attack alpha
theta=(thetap+twists);
alpha=((180/pi)*(phi-theta));

%Check alpha for +-180 bounds
uplimit=(alpha>180);
lolimit=(alpha<-180);
alpha=(alpha+(uplimit*(-360)));
alpha=(alpha+(lolimit*(360)));
Vrel=((((data.vo_design*v1).^2)+((data.vo_design*v2).^2)).^(1/2));
Re=(data.rho*chords.*Vrel)/data.vis;
%Step(5) Read Cl(alpha) & Cd(alpha) from airfoil tables
%Reynolds number range
re_array=data.re_array;

%Step(5) Read Cl(alpha) & Cd(alpha) from airfoil tables
aoa=un_afdata(:,1,1);
re10_array=log10(re_array);
Re10=log10(Re);
for k=1:1:size(a,1);
    for i=1:1:belements
        [cli cdi]=Calc_Cl_Cd(Re10(k,i),re10_array,un_afdata,alpha(k,i),aoa);
        cl_cd(k,i,1)=cli;
        cl_cd(k,i,2)=cdi;
    end
end
cl=cl_cd(:,:,1);
cd=cl_cd(:,:,2);

%Step(6) Compute Cn and Ct
cn=(cl.*cosphi)+(cd.*sinphi);
ct=(cl.*sinphi)-(cd.*cosphi);

%Step(7) Compute anew and ap
sw=sigma.*(cn./(4.*(sinphi.^2)));
cT=4*sw.*((1-a).^2);
aind=(cT<(0.96*pF));
ath=(1-a).*sw./pF;
B0=(cT.*(50-(36*pF)))+((12*pF).*((3*pF)-4));
aco=((18*pF)-20-(3*(B0.^0.5)))./((36*pF)-50);
anew(aind)=ath(aind);
anew(~aind)=aco(~aind);
apnew=((((4*pF.*sinphi.*cosphi)./(sigma.*ct))- 1).^(-1));

%Step(8) a and ap tolerance
a_tol=abs(anew-a);
a=a+(relax*(anew-a));
ap=ap+(relax*(apnew-ap));

a_maxtol=max(max(a_tol));
if a_maxtol<1e-8 & z>20;    break;  end
end
%disp(z)

%Step(9) Compute the power & thrust coefficients
dcp=(B*(1-a).*(tsr.^2).*(1+ap).*chords.*ct.*(nu.^2))./(pi*R*sinphi.*cosphi);
dcT=((B*((1-a).^2).*chords.*cn)./(pi*R*(sinphi.^2)));
%Replace dcp=NaN & dcT=NaN with 0
indcp=isnan(dcp);
indcT=isnan(dcT);
dcp(indcp)=0;
dcT(indcT)=0;

dcp_dcT=cat(3,dcp,dcT);
end

function [cli, cdi] = Calc_Cl_Cd(Re,re_array,un_afdata,add,aoa)
if Re<min(re_array)
    cli=interpolate(aoa,un_afdata(:,2,1),add);
    cdi=interpolate(aoa,un_afdata(:,3,1),add);
    return
elseif Re>max(re_array)
    cli=interpolate(aoa,un_afdata(:,2,end),add);
    cdi=interpolate(aoa,un_afdata(:,3,end),add);
    return
end

I=find(re_array>Re,1);
K=(Re-re_array(I-1))/(re_array(I)-re_array(I-1));

cl_re=(un_afdata(:,2,I-1)*(1-K))+(un_afdata(:,2,I)*K);
cd_re=(un_afdata(:,3,I-1)*(1-K))+(un_afdata(:,3,I)*K);

cli=interpolate(aoa,cl_re,add);
cdi=interpolate(aoa,cd_re,add);
return
end

function [yi] = interpolate(X,Y,xi)
I=find(X > xi,1);
yi=Y(I-1)+( (Y(I)-Y(I-1))*((xi-X(I-1))/(X(I)-X(I-1))) ); 
end
