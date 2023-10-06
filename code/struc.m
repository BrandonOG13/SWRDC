function [mass_cons, dmass, fq,layouts_dF]=struc(chords,twists,pn,pt,data,nu,pitch,pvalstruct,rpm_max,add_structure)
% credit for this source code goes to:
% Matias Sessarego
belements=length(nu);
areavio=trapz((nu*data.R),chords);
%% 1st Flapwise Lowest Eigenfrequency constraint limits
p=rpm_max*(pi/30); p4=data.L*p; p5=data.U*p;

%% Profile Data load
limit=data.limit;
yu=data.yu;
yl=data.yl;
xu=data.xu;
xl=data.xl;

xy_profs=cat(2,xu,yu,yl);

%Profiles for entire blade
xy_con=repmat(xy_profs,[1,1,belements]);

%% Slope Calculation
%Slopes
dydx=zeros(limit,2);
%Forward Difference
dydx(1,1)=(yu(2)-yu(1))./(xu(2)-xu(1));
dydx(1,2)=(yl(2)-yl(1))./(xl(2)-xl(1));
%Backward Difference
dydx(end,1)=(yu(end)-yu(limit-1))./(xu(end)-xu(limit-1));
dydx(end,2)=(yl(end)-yl(limit-1))./(xl(end)-xl(limit-1));
%Centered Difference
for i=2:1:(limit-1)
    dydx(i,1)=(yu(i+1)-yu(i-1))./(xu(i+1)-xu(i-1));
    dydx(i,2)=(yl(i+1)-yl(i-1))./(xl(i+1)-xl(i-1));
end

%Inner Coordinates Construction
tdx_sign=sign(dydx);
tdx_sign(:,2)=-tdx_sign(:,2);
tdy_sign=ones(limit,2);
tdy_sign(:,1)=-tdy_sign(:,1);
angle=atan(abs(dydx.^(-1)));

%% Structural Layup
%Determine Max Spar Cap Thickness
%initial guess
ydiff=(xy_profs(:,2)-xy_profs(:,3))/2;
[ymax1,loc]=max(ydiff);
angdiff=angle(:,1);

ymax=ymax1(end)*sin(angdiff(loc));

maxsparTh=0.95*ymax;

%% Check Strain at Max Spar Cap Thickness
%max strain
ep_max_tens=data.ep_max_tens/1000000;

numlayers=2;
chords3D=repmat(chords',(3*limit*numlayers),1);
chords_mat=reshape(chords3D,[limit 3 (belements*numlayers)]);

    layoutsSpar1=trim(maxsparTh,xy_profs(:,1),xy_profs(:,1),xy_profs(:,2),xy_profs(:,3),angle,tdx_sign,tdy_sign,limit);
    %plot(xy_profs(:,1),xy_profs(:,2),xy_profs(:,1),xy_profs(:,3),layoutsSpar(:,1,i),layoutsSpar(:,2,i),layoutsSpar(:,1,i),layoutsSpar(:,3,i))
    layoutsSpar=repmat(layoutsSpar1,[1 1 belements]);
    
layouts=cat(2,xy_con,layoutsSpar);
layouts_d=reshape(layouts,[limit 3 (numlayers*belements)]).*chords_mat;
stiffs_mass_cg=inert(layouts_d,data,numlayers);
[b_moms_def,fq]=beam(twists,pitch,pn,pt,stiffs_mass_cg,data,nu);
[in_val,maxStrains]=strains(xy_con,stiffs_mass_cg,b_moms_def,chords,0.9999*ep_max_tens);
strainvio=max(maxStrains);
    
if (strainvio>ep_max_tens)||(add_structure==0)
    if (strainvio>ep_max_tens)
    disp('Max Strain at Max Thickness Exceeded')
    end
    strainviof=(strainvio-ep_max_tens)/ep_max_tens;
    mass=trapz((nu*data.R),stiffs_mass_cg(:,4));
    freqvio1=(p4-fq)/p4; freqvio2=(fq-p5)/p5;
    mass_cons=[mass strainviof freqvio1 freqvio2];
    dmass=stiffs_mass_cg(:,4);
    layouts_dF=layouts_d;
    if pvalstruct==1
        [~]=structuralplot(nu,maxsparTh,chords,stiffs_mass_cg,b_moms_def);
    end
    return
end

%% Strain Constraint / Spar Cap Thickness Distribution

[mass,strainvio,sparTh,stiffs_mass_cg,b_moms_def,fq,layouts_dF] = refinemass(maxsparTh,belements,limit,xy_profs,tdx_sign,tdy_sign,...
    angle,data,xy_con,chords_mat,chords,twists,pitch,pn,pt,nu,numlayers);

strainviof=(strainvio-ep_max_tens)/ep_max_tens;
freqvio1=(p4-fq)/p4; freqvio2=(fq-p5)/p5;

%% Plot/Save Results

%dlmwrite('out_struc.txt',stiffs_mass_cg, 'delimiter', '\t','precision', 6) 

if pvalstruct==1
    [~]=structuralplot(nu,sparTh,chords,stiffs_mass_cg,b_moms_def);
end

%% Output Results
mass_cons=[mass strainviof freqvio1 freqvio2];
dmass=stiffs_mass_cg(:,4);
end

function [h1]=structuralplot(nu,sparTh,chords,stiffs_mass_cg,b_moms_def)
h1=figure(2);
subplot(3,2,1);
plot(nu,(sparTh.*chords),'-k.',nu,sparTh,'-g.')
xlabel('r/R(dimensionless)');   ylabel('Shell Thickness');
%legendoptions1=legend('SCT (m)','SCT (t/c)');
%set(legendoptions1,'Location','NorthWest');
title('Shell Thickness vs. r/R')

subplot(3,2,2);
plot(nu,stiffs_mass_cg(:,1)/(10^9),'-r.',nu,stiffs_mass_cg(:,2)/(10^9),'-b.')
xlabel('r/R(dimensionless)');   ylabel('Stiffness EI (GPa m^2)');
legendoptions1=legend('EI_1','EI_2');
set(legendoptions1,'Location','NorthEast');
title('Stiffness EI vs. r/R')

subplot(3,2,3);
plot(nu,-b_moms_def(:,1)/(10^6),'-r.',nu,-b_moms_def(:,2)/(10^6),'-b.')
xlabel('r/R(dimensionless)');   ylabel('Bending Moment (MPa m)');
legendoptions1=legend('M1','M2');
set(legendoptions1,'Location','NorthEast');
title('Bending Moment vs. r/R')

subplot(3,2,4);
plot(nu,b_moms_def(:,3),'-r.',nu,b_moms_def(:,4),'-b.')
xlabel('r/R(dimensionless)');   ylabel('Deflection (m)');
legendoptions1=legend('u_n','u_t');
set(legendoptions1,'Location','NorthWest');
title('Deflection vs. r/R')

subplot(3,2,5);
plot(nu,b_moms_def(:,6),'-r.',nu,b_moms_def(:,7),'-b.')
xlabel('r/R(dimensionless)');   ylabel('1st Flapwise Eigenmode (-)');
legendoptions1=legend('U_{y,flap1}','U_{z,flap1}');
set(legendoptions1,'Location','NorthWest');
title('1st Flapwise Eigenmode vs. r/R')

set(h1, 'units', 'centimeters', 'pos', [0 2 18 18])
end