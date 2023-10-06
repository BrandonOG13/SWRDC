function [mass,strainvio,sparThNew,stiffs_mass_cg,b_moms_def,fq,layouts_d] = refinemass(maxsparTh,belements,limit,xy_profs,tdx_sign,tdy_sign,...
    angle,data,xy_con,chords_mat,chords,twists,pitch,pn,pt,nu,numlayers)
% credit for this source code goes to: Matias Sessarego
%% CALCULATE INITIAL GUESS FOR CAP THICKNESS DISTRIBUTION
ep_max_tens=data.ep_max_tens/1000000;

%load outer profile coordinates
xy1=xy_profs(:,1);  xy2=xy_profs(:,2);  xy3=xy_profs(:,3);

%inital point & boundary conditions
lb=data.MinsparTh;
ub=maxsparTh;
x0=zeros(belements,1)+lb;
scacc=ub/1000;
for j=1:50
    [x0,in_val,strainvio,scacc,flag,stiffs_mass_cg,b_moms_def,fq,layouts_d] = calc_spar(x0,chords,chords_mat,twists,pitch,pn,pt,nu,data,numlayers,xy_con,scacc,xy1,xy2,xy3,angle,tdx_sign,tdy_sign,limit,belements,lb,ub,ep_max_tens);
    if flag==1
        break
    end
    scacc=scacc/10;
end
sparThNew=x0;
mass=trapz((nu*data.R),stiffs_mass_cg(:,4));

end

%% Algorithm for calculating spar cap thickness distribution
function [x0,in_val,strainvio,scacc,flag,stiffs_mass_cg,b_moms_def,fq,layouts_d] = calc_spar(x0in,chords,chords_mat,twists,pitch,pn,pt,nu,data,numlayers,xy_con,scacc,xy1,xy2,xy3,angle,tdx_sign,tdy_sign,limit,belements,lb,ub,ep_max_tens)
maxiter=500;
flag=0;

x0=x0in;
layoutsSpar=zeros(limit,3,belements);
in_val=true(belements,1);
for w1=1:1:maxiter
    count=find(in_val);
    lencount=length(count);
    for k=1:lencount
        ki=count(k);
        layoutsSpar(:,:,ki)=trim(x0(ki),xy1,xy1,xy2,xy3,angle,tdx_sign,tdy_sign,limit);
    end
    layouts=cat(2,xy_con,layoutsSpar);
    layouts_d=reshape(layouts,[limit 3 (numlayers*belements)]).*chords_mat;
    stiffs_mass_cg=inert(layouts_d,data,numlayers);
    
    %     figure(32);plot(layouts_d(:,1,33),layouts_d(:,2,33),'-k',layouts_d(:,1,33),layouts_d(:,3,33),'-k',...
    %         layouts_d(:,1,34),layouts_d(:,2,34),'-k',layouts_d(:,1,34),layouts_d(:,3,34),'-k',...
    %         stiffs_mass_cg(17,7),stiffs_mass_cg(17,8),'s'); axis equal
    %     axis([-0.001 max(layouts_d(:,1,33))+0.001 0 1])
    %     axis 'auto y'
    %     axis off
    
    [b_moms_def,fq]=beam(twists,pitch,pn,pt,stiffs_mass_cg,data,nu);
    [in_val,maxStrains]=strains(xy_con,stiffs_mass_cg,b_moms_def,chords,ep_max_tens);
    
    %       disp(in_val)
    
    check=any(in_val);
    if check==0
        break
    end
    %in_val_e=in_val;
    x0(in_val)=x0(in_val)+scacc;
    if any(x0>ub)
        x0(x0>ub)=ub;
        break
    end
end

if scacc<1e-9
    flag=1;
    strainvio=max(maxStrains);
    return
end

strainvio=max(maxStrains);
%x0(in_val_e)=x0(in_val_e)-scacc;
end