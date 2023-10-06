function [pn_worst,pt_worst,tempdata] = iecClass(pitch,twists,chords,un_afdata,tempdata,nu,rpm_max)
%IEC Regulations
%Safety Factor on Aerodynamic Loads:
SF=tempdata.loadSF;
R=tempdata.R;
Ve50=tempdata.Ve50;

loadcase=tempdata.loadcase;

if loadcase==1 || loadcase==2
    [pn_worst,pt_worst] = parked(pitch,twists,chords,un_afdata,tempdata,nu,SF,Ve50,R);
elseif loadcase==3
    tempdata.vo_design=Ve50;

    over_omega=tempdata.rotoroverspeed*tempdata.rat_rpm*(pi/30);
    tsr=(over_omega*R)/tempdata.vo_design;
    tempdata.lambda_design=tsr;

    [pn_worst,pt_worst] = overspeed(pitch,tsr,twists,chords,un_afdata,tempdata,nu,SF,rpm_max);
end

end

function [pn_worst,pt_worst] = parked(pitch,twists,chords,un_afdata,tempdata,nu,SF,Ve50,R)
belements=length(twists);

phi=(-180:45:180);
phi_rep=repmat(phi,belements,1);
phi_rep_deg=phi_rep*(pi/180);
lenphi=length(phi);

factor=repmat((0.5*tempdata.rho*(Ve50^2)*chords),1,lenphi);
twists_rep=repmat(twists,1,lenphi);
t_p=(twists_rep+pitch)*(180/pi);
ch_rep=repmat(chords,1,lenphi);

alpha=phi_rep-t_p;
%Check alpha for +-180 bounds
uplimit=(alpha>180);
lolimit=(alpha<-180);
alpha=(alpha+(uplimit*(-360)));
alpha=(alpha+(lolimit*(360)));

Re10=log10((tempdata.rho*ch_rep.*Ve50)/tempdata.vis);
re10_array=log10(tempdata.re_array);

%Step(5) Read Cl(alpha) & Cd(alpha) from airfoil tables
aoa=un_afdata(:,1,1);
cl_cd=zeros(belements,lenphi,2);
for k=1:1:lenphi;
    for i=1:1:belements
        [cli cdi]=Calc_Cl_Cd(Re10(i,k),re10_array,un_afdata,alpha(i,k),aoa);
        cl_cd(i,k,1)=cli;
        cl_cd(i,k,2)=cdi;
    end
end
cl=cl_cd(:,:,1);
cd=cl_cd(:,:,2);

pn50=SF*factor.*( (cl.*cos(phi_rep_deg))+(cd.*sin(phi_rep_deg)) );
pt50=SF*factor.*( (cl.*sin(phi_rep_deg))-(cd.*cos(phi_rep_deg)) );

%M - Root Flap Bending Moment
nu_rep=repmat(nu,1,lenphi);
rtflapbend=trapz((nu*R),(nu_rep.*pn50),1);
rtedgebend=trapz((nu*R),(nu_rep.*pt50),1);

if (round(abs(min(rtflapbend))*10^6)/10^6)<=(round((max(rtflapbend))*10^6)/10^6)
    %[o1_max j1_max]=max(rtflapbend);
   [ o1, j1 ]=max(rtflapbend);
else
    %[o1_min j1_min]=min(rtflapbend);
   [ o1, j1 ]=min(rtflapbend);
end

if (round(abs(min(rtedgebend))*10^6)/10^6)<=(round((max(rtedgebend))*10^6)/10^6)
    %[o2_max j2_max]=max(rtedgebend);
   [ o2, j2 ]=max(rtedgebend);
else
    %[o2_min j2_min]=min(rtedgebend);
   [ o2, j2 ]=min(rtedgebend);
end

%[o1 j1]=max(abs(rtflapbend));
%[o2 j2]=max(abs(rtedgebend));
%[o1 j1]=max(rtflapbend);
%[o2 j2]=max(rtedgebend);

pn_worst=pn50(:,j1);
pt_worst=pt50(:,j2);
end

function [pn_worst,pt_worst] = overspeed(pitch,tsr,twists,chords,afdata,tempdata,nu,SF,rpm_max)

%Run the BEM code for given pitch & tsr ranges
dcp_dct=bem(pitch,tsr,twists,chords,tempdata,afdata,nu);

%Load Results
[cp,ct,rootflapbend,pn,pt]=loads(dcp_dct,tempdata,nu,rpm_max);

pn_worst=SF*pn;
pt_worst=SF*pt;
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
for i=2:length(X)
    if xi<=X(i)
        I=i;
        break
    end
end
yi=Y(I-1)+( (Y(I)-Y(I-1))*((xi-X(I-1))/(X(I)-X(I-1))) );
end