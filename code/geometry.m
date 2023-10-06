function [chords,twists,thicks] = geometry(x,data,nu,pval)
%Extract variables from data base
thick=data.thick;
N=data.N;
hubR_rotR=data.hubR/data.R;

numcps_ch_tw=((data.N+1)*2)-1;
%apply tip twist constraint
cps_tw=cat(1,x(1:N),data.tw_tip);
cps_ch=x((N+1):numcps_ch_tw);

cp_aero=zeros((N+1),2);
cp_aero(:,1)=cps_tw*(pi/180);
cp_aero(:,2)=cps_ch;
cp_aero=sort(cp_aero,1,'descend');

%% Blade Geometry
CP_r=([0 1]*(1-hubR_rotR))+hubR_rotR;

%Generate Chord Distribution
%aerodynamic blade portion
cp_ch_r_aero = (data.r_space*(1-CP_r(1)))+CP_r(1);
cp_ch_aero = cp_aero(:,2);

cp_ch=cp_ch_aero;

%Generate Twist distribution
cp_tw_aero = cp_aero(:,1);

chords=bezier(cp_ch,N,nu,cp_ch_r_aero);
twists=bezier(cp_tw_aero,N,nu,cp_ch_r_aero);
thicks=chords.*thick;

if pval==1
    figure(1)
    subplot(2,2,1);
    plot(nu,(twists*(180/pi)),'-b.',cp_ch_r_aero,cp_tw_aero*(180/pi),'ks')
    xlabel('r/R(dimensionless)');   ylabel('Twist (degrees)');
    legendoptions1=legend('Twist','CP');
    set(legendoptions1,'Location','NorthEast');
    title('Twist vs. r/R')
    
    subplot(2,2,2);
    plot(nu,chords,'-r.',cp_ch_r_aero,cp_ch,'ks')
    xlabel('r/R(dimensionless)');   ylabel('Chord (m)');
    legendoptions1=legend('Chord','CP');
    set(legendoptions1,'Location','NorthEast');
    title('Chord vs. r/R')
    
    subplot(2,2,3);
    plot(nu,(thick.*ones(size(nu))),'-g.',CP_r,(thick.*ones(size(CP_r))),'ks')
    xlabel('r/R(dimensionless)');   ylabel('Thickness/Chord (%)');
    legendoptions1=legend('Thickness','CP');
    set(legendoptions1,'Location','NorthEast');
    title('% Thickness vs. r/R')
    
    subplot(2,2,4);
    plot(nu,chords.*thick,'-k.')
    xlabel('r/R(dimensionless)');   ylabel('Thickness (m)');
    legendoptions1=legend('Thickness');
    set(legendoptions1,'Location','NorthEast');
    title('Dimensional Thickness vs. r/R')
end

end

