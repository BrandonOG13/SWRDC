function [in_val,maxStrains] = strains(xy_con,stiffs_mass_cg,b_moms_def,chords,ep_max_tens)
% credit for this source code goes to:
% Matias Sessarego
np=size(xy_con,1);

%Strain from bending moments
M1=b_moms_def(:,1);
M2=b_moms_def(:,2);

EI1=stiffs_mass_cg(:,1);
EI2=stiffs_mass_cg(:,2);

%Strain from normal force
N=b_moms_def(:,8);
EA=stiffs_mass_cg(:,11);
N_EA=N./EA;

%Determine coordinates
cos_theta=repmat(cos(stiffs_mass_cg(:,3))',np,1);
sin_theta=repmat(sin(stiffs_mass_cg(:,3))',np,1);

XEvec=repmat((stiffs_mass_cg(:,7)./chords)',np,1);
YEvec=repmat((stiffs_mass_cg(:,8)./chords)',np,1);

dx =squeeze(xy_con(:,1,:))-XEvec;
dyu=squeeze(xy_con(:,2,:))-YEvec;
dyl=squeeze(xy_con(:,3,:))-YEvec;

dxcos=(dx.*cos_theta);
dxsin=(dx.*sin_theta);

xu=(dyu.*sin_theta)+dxcos; %Upper X
yu=(dyu.*cos_theta)-dxsin; %Upper Y
xl=(dyl.*sin_theta)+dxcos; %Lower X
yl=(dyl.*cos_theta)-dxsin; %Lower Y 

x=cat(1,xu,xl);
y=cat(1,yu,yl);

% Use max y & x to determine strain (1) or combined moments on x & y (2)
method=2;
if method==1
    max_y=abs(max(y,[],1)').*chords;
    min_y=abs(min(y,[],1)').*chords;
    
    max_x=abs(max(x,[],1)').*chords;
    min_x=abs(min(x,[],1)').*chords;
    
    ep_f=cat(2,abs(M1.*max_y)./EI1,abs(M1.*min_y)./EI1);
    ep_e=cat(2,abs(M2.*max_x)./EI2,abs(M2.*min_x)./EI2);
    
    strains_mat=cat(2,ep_f,ep_e);
    M_EI=max(strains_mat,[],2);
    
    %Total strain
    maxStrains=M_EI+N_EA;
elseif method==2
    N_EA_mat=repmat(N_EA',size(y,1),1);
    
    M1_mat=repmat(M1',size(y,1),1);
    M2_mat=repmat(M2',size(y,1),1);
    
    EI1_mat=repmat(EI1',size(y,1),1);
    EI2_mat=repmat(EI2',size(y,1),1);
    
    chords_mat=repmat(chords',size(y,1),1);
    
    strainsB=((M1_mat.*y.*chords_mat)./EI1_mat)-((M2_mat.*x.*chords_mat)./EI2_mat)+N_EA_mat;
    
    %Total strain
    maxStrains = max(strainsB,[],1)';
end

in_val=(maxStrains>=ep_max_tens);