% ============================ START PROGRAM ======================= %

                          % ================== %
                          %    Main Program    %
                          % ================== %
             
% << All Variables >> %

% alpha         -->  Angle of attack
% a_ground      -->  The radius of the tower at the root
% a_top         -->  The radius of the tower at the top
% A_weighting   -->  Standard frequency weighting filter corrected with the
%                    sensitivity of the human auditory system
% c             -->  Segment chord length
% co            -->  Sound speed in air
% Chord         -->  Chord length
% Density       -->  The air density
% dis           -->  Distance between observer and rotor center
% dspl          -->  Intermediate parameter
% dsplw         -->  Intermediate parameter
% espl          -->  Intermediate parameter
% esplw         -->  Intermediate parameter
% Freq          -->  1/3 octave frequencies
% Gamma         -->  Paremeter denote the amount of shear,range between
%                    0.1-0.25
% h             -->  Trailing edge bluntness
% ho            -->  Altitude of the observer referenced to the tower base 
% Htower        -->  The height of the tower
% L             -->  Blade segment span
% logi_trip     -->  Condition for boundary layer tripping:
%                    tripped=0, untripped=1, partially tripped=2
% Nb            -->  Total number of blades
% Nfre          -->  Total number of frequencies
% Nwing         -->  Total number of azimuthal angles of blade
% Npsi          -->  Total number of observer view angles
% Nseg          -->  Total number of blade segments
% omega         -->  The angular rotating speed of the blade 
% P             -->  Total sound pressure
% P1            -->  Sound pressure due to TBLTE at pressure side
% P2            -->  Sound pressure due to TBLTE at suction side
% P3            -->  Sound pressure due to TBLTE at nonzero angle
% P4            -->  Total sound pressure due to TBLTE 
% P5            -->  Sound pressure due to LBLVS 
% P6            -->  Sound pressure due to TIP
% P7            -->  Sound pressure due to TEBVS 
% P8            -->  Sound pressure due to TBINF
% psi           -->  The observer angle reference to the downstream
%                    direction
% PSI           -->  The observer angle reference to the downstream
%                    direction
% r             -->  The blade radius at each segment
% R             -->  All the blade radius at each segment
% Ratio         -->  The coffecient to refine alpha_tip
% ro            -->  Distance between the observer and the tower center
% Round         -->  The shape of the blade tip, Round = 0 correspond to
%                    the rounded tip, else the tip is rectangular
% shaftlength   -->  The shaft length of the turbine
% SPL           -->  Sound pressure level
% SPLa          -->  Sound pressure level due to TBLTE at nonzero angle 
% SPLINF        -->  Sound pressure level due to TBINF
% SPLLBL        -->  Sound pressure level due to LBLVS
% SPLp          -->  Sound pressure level due to TBLTE at pressure side
% SPLs          -->  Sound pressure level due to TBLTE at suction side
% SPLTBL        -->  Total sound pressure level due to TBLTTE
% SPLTEB        -->  Sound pressure level due to TEBVS
% SPLTIP        -->  Sound pressure level due to TIP
% SPLtot        -->  Total equivilent sound pressure level  without A-weighting
% SPLtotw       -->  Total equivilent sound pressure level with A-weighting
% SPLw_         -->  Sound power level
% SPLwA_        -->  A-weighted sound power level
% Sum_dspl      -->  Intermediate parameter
% Sum_dsplw     -->  Intermediate parameter
% Sum_espl      -->  Intermediate parameter
% Sum_esplw     -->  Intermediate parameter
% Sum_P         -->  Intermediate parameter
% Sum_P1        -->  Intermediate parameter
% Sum_P2        -->  Intermediate parameter
% Sum_P3        -->  Intermediate parameter
% Sum_P4        -->  Intermediate parameter
% Sum_P5        -->  Intermediate parameter
% Sum_P6        -->  Intermediate parameter
% Sum_P7        -->  Intermediate parameter
% Sum_P8        -->  Intermediate parameter                
% TE_angle      -->  Trailing edge angle
% theta_cone    -->  Cone angle
% theta_tilt    -->  Tilt angle
% theta_twist   -->  Twist angle
% theta_wing    -->  Wing angle
% theta_yaw     -->  Yaw angle
% Tinten        -->  Turbulence intensity
% Tscale        -->  Turbulence length scale
% vis           -->  Kinematic viscosity
% Vo            -->  The mean wind speed at the hub hight
% Vrel          -->  The relative wind velovity
% wing          -->  Different wing operating conditions
% WING          -->  Wing angle for one of the blade
% wSPLtot       -->  Total equivilent sound power level  without A-weighting
% wSPLtot       -->  Total equivilent sound power level  with A-weighting
% zo            -->  Surface roughness

% <<All Mechanisums>> %
% TBLTE         -->  Turbulent boundary layer trailing edge noise
% LBLVS         -->  Laminar boundary layer vortex shedding noise
% TIP           -->  Tip vortex formation noise
% TEBVS         -->  Trailing edge bluntness vortex shedding noise
% TBINF         -->  Turbulent inflow noise

function [SPLw,SPLw1,SPLw2,SPLw3,SPLw4,SPLw5,SPLw6,SPLw7,SPLw8,...
          SPLwA,SPLwA1,SPLwA2,SPLwA3,SPLwA4,SPLwA5,SPLwA6,SPLwA7,SPLwA8,...
          SPL,SPL1,SPL2,SPL3,SPL4,SPL5,SPL6,SPL7,SPL8,...
        scale,inten,ro,ho,PSI,Freq,A_weighting,SPLtot,SPLtotw,wSPLtot,wSPLtotw]=MAIN(Airfoil,Vo,Density,vis,...
    Tscale,Tinten,ro,ho,PSI,co,Ratio,theta_yaw,theta_tilt,theta_cone,wing,...
    a_top,a_ground,Htower,shaftlength,omega,Gamma,Round,logi_trip,R,Chord,Twist,Hblunt,TEANGLE,zo,...
    noisedatabase,Rref,Aref)



Center_frequency =...
  [  20                             -50.5
     25                             -44.7
     31.5                           -39.4
     40                             -34.6
     50                             -30.2
     63                             -26.2
     80                             -22.5
     100                            -19.1
     125                            -16.1
     160                            -13.4
     200                            -10.9
     250                            -8.6
     315                            -6.6
     400                            -4.8
     500                            -3.2
     630                            -1.9
     800                            -0.8
     1000                            0
     1250                            0.6
     1600                            1.0
     2000                            1.2
     2500                            1.3
     3150                            1.2
     4000                            1.0
     5000                            0.5
     6300                           -0.1
     8000                           -1.1
     10000                          -2.5];
    %12500                          -4.3
    %16000                          -6.6
    %20000                          -9.3 ];
%      |                               |     
%  Center frequency [Hz]          A-weighting [dB] 


Freq = Center_frequency(:,1);
A_weighting = Center_frequency(:,2);


Nwing = length(wing);
Nb = 3;
Nseg = length(R);
Nfre = length(Freq);
Npsi = length(PSI);
dis = sqrt(ro^2+(Htower-ho)^2);
% ro
% dis
for i = 1:Npsi    % L 10    % To compute the cases for different observer angles.
    psi = PSI(i);

Sum_P1 = zeros(1,Nfre);
Sum_P2 = zeros(1,Nfre);
Sum_P3 = zeros(1,Nfre);
Sum_P4 = zeros(1,Nfre);
Sum_P5 = zeros(1,Nfre);
Sum_P6 = zeros(1,Nfre);
Sum_P7 = zeros(1,Nfre);
Sum_P8 = zeros(1,Nfre);
Sum_P = zeros(1,Nfre);
Sum_dspl = 0;
Sum_dsplw = 0;
Sum_espl = 0;
Sum_esplw = 0;


    for m = 1:Nwing   % L 20  % To compute the cases for different blade angles.
P = zeros(1,Nfre);
P1 = zeros(1,Nfre);
P2 = zeros(1,Nfre);
P3 = zeros(1,Nfre);
P4 = zeros(1,Nfre);
P5 = zeros(1,Nfre);
P6 = zeros(1,Nfre);
P7 = zeros(1,Nfre);
P8 = zeros(1,Nfre);
   
time0 = cputime;

    for n = 1:Nb   % L 30 
        tn0 = cputime;
         WING=wing(m);
        theta_wing = WING + (360/Nb)*(n-1);
        for k = 1:Nseg-1  % L 40
            tk0 = cputime;
            L=R(k+1)-R(k);
            r=R(k+1)-0.5*L+1;
            c=(Chord(k+1)+Chord(k))*0.5;
            theta_twist=(Twist(k+1)+Twist(k))*0.5;
            h = (Hblunt(k+1)+Hblunt(k))*0.5;
            TE_angle = (TEANGLE(k+1)+TEANGLE(k))*0.5;
            
            % without effect from induced velocity
             [Vrel,alpha]=V_relative(theta_yaw,theta_tilt,theta_cone,theta_wing,theta_twist,a_top,a_ground,...
                          shaftlength,Gamma,Htower,Vo,omega,r,zo);
                  % with effect from induced velocity, run BEM code first     
%            [Vrel,alpha]=V_relative_WK(theta_yaw,theta_tilt,theta_cone,theta_wing,theta_twist,a_top,a_ground,...
%                          shaftlength,Gamma,Htower,Vo,omega,r,zo,R);
           %attack(k)=alpha
           % Compute BL thickness and BL displacement thickness from XFOIL(or from experiments):
           t00 = cputime;
           [delt_p,delt_star_p,delt_s,delt_star_s]=Thickness_xfoil(Airfoil,Vrel,R,r,c,vis,co,logi_trip,alpha,...
               noisedatabase,Rref,Aref); 
           t_thick(n,k) = cputime -t00;
           % [delt_p,delt_star_p,delt_s,delt_star_s]=Thickness(Vrel,c,vis,co,logi_trip,alpha); % Experiment
           % delta(k)=delt_p;  % plot the thickness vs. blade length...
           % Compute directivity:
                      t00 = cputime;

            [THETA,PHI,Distance]=Trans_direc(ro,ho,psi,theta_yaw,theta_tilt,theta_cone,theta_wing,theta_twist,...
             shaftlength,Htower,r);
            [Di_hi_fr,Di_lo_fr]=Directivity(Vrel,co,THETA,PHI);
            t_direct(n,k) = cputime -t00;
% Distance
           % Compute TBLTE noise:
           t00 = cputime;
            [SPLTBL,SPLp,SPLs,SPLa]=TBLTE(delt_p,delt_star_p,delt_s,delt_star_s,Di_hi_fr,Di_lo_fr,Distance,...
                Freq,Vrel,c,vis,co,logi_trip,alpha,ro,ho,psi,theta_yaw,theta_tilt,...
                 theta_cone,theta_wing,theta_twist,shaftlength,Htower,r,L);
             t_tblte(n,k)=cputime -t00;
             
                 
                 t00 = cputime;
            % Compute LBLVS noise:
            [SPLLBL]=LBLVS(delt_p,delt_star_p,delt_s,delt_star_s,Di_hi_fr,Di_lo_fr,Distance,...
                Freq,Vrel,c,vis,co,logi_trip,alpha,ro,ho,psi,theta_yaw,theta_tilt,...
                 theta_cone,theta_wing,theta_twist,shaftlength,Htower,r,L);
            t_lblvs(n,k)=cputime -t00;
                 
                 t00 = cputime;
            % Compute TIP noise only for the last segment:
            if k==Nseg-1
                [SPLTIP]=TIP(Di_hi_fr,Di_lo_fr,Distance,...
                    Freq,Vrel,co,c,alpha,ro,ho,psi,theta_yaw,theta_tilt,theta_cone,theta_wing,theta_twist,...
                shaftlength,Htower,r,Round,Ratio);
            end
            t_tip(n,k)=cputime-t00;
            
            
            t00 = cputime;
             % Compute TEBVS noise:
              [SPLTEB]=TEBVS(delt_p,delt_star_p,delt_s,delt_star_s,Di_hi_fr,Di_lo_fr,Distance,...
                  Freq,h,TE_angle,Vrel,c,vis,co,logi_trip,alpha,ro,ho,psi,theta_yaw,theta_tilt,...
           theta_cone,theta_wing,theta_twist,shaftlength,Htower,r,L);
       t_tebvs(n,k)=cputime -t00;
           
           t00 = cputime;
             % Compute TBLINF noise:
           [SPLINF,Tbscale,Tbinten]=TBINF(Di_hi_fr,Di_lo_fr,Distance,...
               Freq,Vo,Vrel,c,L,vis,co,ro,ho,psi,theta_yaw,theta_tilt,theta_cone,theta_wing,theta_twist,...
            shaftlength,Htower,r,Density,Tscale,Tinten,zo,Gamma);
          scale(i,m,n,k) = Tbscale;
          inten(i,m,n,k) = Tbinten;
          t_tblinf(n,k)=cputime -t00;
          
          t00 = cputime;
            for f = 1:Nfre
                % sum up the sound pressure from all blade segments 
                P1(f) = P1(f) + 10^(SPLp(f)/10);
                P2(f) = P2(f) + 10^(SPLs(f)/10);
                P3(f) = P3(f) + 10^(SPLa(f)/10);
                P4(f) = P1(f) + P2(f) + P3(f);
                
%                 P5(f) = P5(f) + 10^(SPLLBL(f)/10);
% % % set to zero if this effect is not included
                  fac = 0.001;
                  P5(f) = P5(f) + 10^(SPLLBL(f)/10)*fac;
                
                if k==Nseg-1
                P6(f) = P6(f) + 10^(SPLTIP(f)/10);
                end
                
                
                P7(f) = P7(f) + 10^(SPLTEB(f)/10);
                
                P8(f) = P8(f) + 10^(SPLINF(f)/10);
                
                P(f) = P4(f) + P5(f) + P6(f) + P7(f) + P8(f);  
            end  
            t_sum(n,k)=cputime -t00;
             tk(k) = cputime - tk0;
          end   % L 40
          tn(n)=cputime-tn0;
       end   %  L 30
%        cputime-time0
%        tk
%        sum(tk)
%        tn
%        sum(tn)
%        t_thick
%        t_direct
%        t_tblte
%        t_lblvs
%        t_tip
%        t_tebvs
%        t_tblinf
%        t_sum
%        sum(tn)
%        sum(tn)/0.0156
    % Store the pressure values at different rotor-azimuth angles
    Sum_P1 = Sum_P1 + P1;
    Sum_P2 = Sum_P2 + P2;
    Sum_P3 = Sum_P3 + P3;
    Sum_P4 = Sum_P4 + P4;
    Sum_P5 = Sum_P5 + P5;
    Sum_P6 = Sum_P6 + P6;
    Sum_P7 = Sum_P7 + P7;
    Sum_P8 = Sum_P8 + P8;
    Sum_P = Sum_P + P;
     end  % L 20 




% Compute the average pressure values at different rotor-azimuth angles
    P1 = Sum_P1/Nwing;
    P2 = Sum_P2/Nwing;
    P3 = Sum_P3/Nwing;
    P4 = Sum_P4/Nwing;
    P5 = Sum_P5/Nwing;
    P6 = Sum_P6/Nwing;
    P7 = Sum_P7/Nwing;
    P8 = Sum_P8/Nwing;
    P = Sum_P/Nwing;

% Compute the sound pressure level for all mechanisms
  for f = 1:Nfre  %/
    if P(f)~=0
    SPL(f) = 10*log10(P(f));
    SPLw(f) = SPL(f)+10*log10(4*pi*dis^2);
    SPLwA(f) = SPL(f)+10*log10(4*pi*dis^2)+A_weighting(f);
    end
    if P1(f)~=0
    SPL1(f) = 10*log10(P1(f));
    SPLw1(f) = SPL1(f)+10*log10(4*pi*dis^2);
    SPLwA1(f) = SPL1(f)+10*log10(4*pi*dis^2)+A_weighting(f);
    end
    if P2(f)~=0
    SPL2(f) = 10*log10(P2(f));
    SPLw2(f) = SPL2(f)+10*log10(4*pi*dis^2);
    SPLwA2(f) = SPL2(f)+10*log10(4*pi*dis^2)+A_weighting(f);
    end
    if P3(f)~=0
    SPL3(f) = 10*log10(P3(f));
    SPLw3(f) = SPL3(f)+10*log10(4*pi*dis^2);
    SPLwA3(f) = SPL3(f)+10*log10(4*pi*dis^2)+A_weighting(f);
    end
    if P4(f)~=0
    SPL4(f) = 10*log10(P4(f));
    SPLw4(f) = SPL4(f)+10*log10(4*pi*dis^2);
    SPLwA4(f) = SPL4(f)+10*log10(4*pi*dis^2)+A_weighting(f);
    end 
    if P5(f)~=0
    SPL5(f) = 10*log10(P5(f));
    SPLw5(f) = SPL5(f)+10*log10(4*pi*dis^2);
    SPLwA5(f) = SPL5(f)+10*log10(4*pi*dis^2)+A_weighting(f);
    end
    if P6(f)~=0
    SPL6(f) = 10*log10(P6(f));
    SPLw6(f) = SPL6(f)+10*log10(4*pi*dis^2);
    SPLwA6(f) = SPL6(f)+10*log10(4*pi*dis^2)+A_weighting(f);
    end
    if P7(f)~=0
    SPL7(f) = 10*log10(P7(f));
    SPLw7(f) = SPL7(f)+10*log10(4*pi*dis^2);
    SPLwA7(f) = SPL7(f)+10*log10(4*pi*dis^2)+A_weighting(f);
    end
    if P8(f)~=0
    SPL8(f) = 10*log10(P8(f));
    SPLw8(f) = SPL8(f)+10*log10(4*pi*dis^2);
    SPLwA8(f) = SPL8(f)+10*log10(4*pi*dis^2)+A_weighting(f);
    end
  end  %/

% Compute the total sound pressure level(include the human response 
% of standard A-weighted filter in one-third octive bands)
   for f = 1:Nfre
       dspl(f) = SPL(f)/10;
       Sum_dspl = Sum_dspl + 10^dspl(f); % Sum up the sound pressure levels
       dsplw(f) = (SPL(f) + A_weighting(f))/10;
       Sum_dsplw = Sum_dsplw + 10^dsplw(f); % Sum up the A-weighted sound pressure levels
       espl(f) = SPLw(f)/10;
       Sum_espl = Sum_espl + 10^espl(f); % Sum up the sound power levels
       esplw(f) = (SPLw(f) + A_weighting(f))/10;
       Sum_esplw = Sum_esplw + 10^esplw(f); % Sum up the A-weighted sound power levels
   end
% --------------------------------------------------
SPLtotw(i)=10*log10(Sum_dsplw);  % A-weighted sound pressure level 
SPLtot(i)=10*log10(Sum_dspl);    % Sound pressure level without A-weighted 
wSPLtotw(i)=10*log10(Sum_esplw);  % A-weighted sound power level 
wSPLtot(i)=10*log10(Sum_espl);    % Sound power level without A-weighted 


end  % L 10


%plot(SPL,'r*-')    % 'd' , 'h' ,'v', ,'p', 's' ,'o' ,'x','+' , '*'
%grid on
%hold on
%plot(SPL8,'k+-')
%hold on
%plot(SPL7,'gv-')
%hold on
%plot(SPL6,'cd-')
%hold on
%plot(SPL5,'mh-')
%hold on
%plot(SPL4,'bp-')
%hold on

% polar(PSI*pi/180,SPLtot,'*-')
%x=SPLtot(:).*cos(PSI(:)*pi/180);
%y=SPLtot(:).*sin(PSI(:)*pi/180);
%plot(x,y,'*-'),axis equal
%grid on
%hold on
%plot(0,0,'p')
%hold on


% n=[1:Nfre];
% plot(n,SPL,'r*-',n,SPL1,'bo-',n,SPL2,'bs-',n,SPL3,'bx-',n,SPL4,'bp-',n,SPL5,'mh-',n,SPL6,'cd-',n,SPL7,'gv-',n,SPL8,'k+-')
% legend('SPL-TOTAL','SPL-TBLTE(Pressure side)','SPL-TBLTE(Suction side)','SPL-SEPARATION','SPL-TBLTE(Total)',...
%     'SPL-LBLVS','SPL-TIP','SPL-TEBLUNT','SPL-TBINFLOW')
% grid on,axis([1 Nfre 0 max(SPL)+10])
% xlabel('Frequency (Hz)'),ylabel('Sound Pressure level (dB)')
% set(gca,'Xtick',[1:2:Nfre])
% set(gca,'XtickLabel',[Freq(1:2:Nfre)])       

%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

              % ==================================================== %
              %                    < TBL-TE >                        %
              %     Turbulent-Boundary-Layer-Trailing-Edge Noise     %
              %              and Separated Flow Noise                %
              % ==================================================== %
                                 
                         % << All Variables >> %
% a             -->  Ratio of Strouhal number
% A             -->  Shape function
% alpha         -->  Angle of attack
% b             -->  Ratio of Strouhal number
% B             -->  Shape function
% c             -->  Segment chord length
% co            -->  Sound speed in air
% coeff         -->  Intemediate parameter
% delt_p        -->  Boundary layer thickness on pressure side
% delt_star_p   -->  Boundary layer displacement thickness on pressure side
% delt_s        -->  Boundary layer thickness on suction side
% delt_star_s   -->  Boundary layer displacement thickness on suction side
% Freq          -->  1/3 octave frequencies
% h             -->  Trailing edge bluntness
% ho            -->  Altitude of the observer referenced to the tower base 
% Htower        -->  The hight of the tower
% K1,K2,dK1     -->  Shape functions
% L             -->  Blade segment span
% logi_trip     -->  Condition for boundary layer tripping:
%                    tripped=0, untripped=1, partially tripped=2
% Ma            -->  Mach number
% P1            -->  Sound pressure due to TBLTE at pressure side
% P2            -->  Sound pressure due to TBLTE at suction side
% P3            -->  Sound pressure due to TBLTE at nonzero angle
% psi           -->  The observer angle reference to the downstream
%                    direction
% PSI           -->  The observer angle reference to the downstream
%                    direction
% r             -->  The blade radius at each segment
% Rc            -->  Reynolds number based on chord
% Rp            -->  Reynolds number based on the pressure  side
%                    displacement thickness
% ro            -->  Distance between the observer and the tower center
% shaftlength   -->  The shaft length of the turbine
% SPL           -->  The total sound pressure level vs.frequencies
% SPLa          -->  Sound pressure level due to TBLTE at nonzero angle 
% SPLp          -->  Sound pressure level due to TBLTE at pressure side
% SPLs          -->  Sound pressure level due to TBLTE at suction side
% SPLTBL        -->  Sound pressure level due to TBLTE 
% St1           -->  Peak Strouhal number
% St2           -->  Strouhal number
% Stp           -->  Strouhal number based on pressure side displacement
%                    thickness
% Sts           -->  Strouhal number based on suction side displacement
%                    thickness

% theta_cone    -->  Cone angle
% theta_tilt    -->  Tilt angle
% theta_twist   -->  Twist angle
% theta_wing    -->  Wing angle
% theta_yaw     -->  Yaw angle
% vis           -->  Kinematic viscosity
% Vo            -->  The mean wind speed at the hub center
% Vrel          -->  The relative wind velovity


function[SPLTBL,SPLp,SPLs,SPLa]=TBLTE(delt_p,delt_star_p,delt_s,delt_star_s,Di_hi_fr,Di_lo_fr,Distance,...
    Freq,Vrel,c,vis,co,logi_trip,alpha,ro,ho,psi,theta_yaw,theta_tilt,...
                 theta_cone,theta_wing,theta_twist,shaftlength,Htower,r,L)
             

% Compute Mach number and Reynolds number based on chord length 
Rc=Vrel*c/vis;
Ma=Vrel/co;

   

% Compute boundary layer thickness
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%[delt_p,delt_star_p,delt_s,delt_star_s]=Thickness(Vrel,c,vis,co,logi_trip,alpha);

% Compute the Reynolds number based on the pressure  side
% displacement thickness.
Rp = delt_star_p * Vrel / vis;

% Compute sound directivity
% [THETA,PHI,Distance]=Trans_direc(ro,ho,psi,theta_yaw,theta_tilt,theta_cone,theta_wing,theta_twist,...
%     shaftlength,Htower,r);
% [Di_hi_fr,Di_lo_fr]=Directivity(Vrel,co,THETA,PHI);

% Compute peak Strouhal number
St1 = 0.02 * Ma^(-0.6);
if alpha<1.33
    coeff = 1;
elseif alpha>=1.33 & alpha<=12.5
    coeff = 10^(0.0054*(alpha-1.33)^2);
else
    coeff = 4.72;
end
St2 = coeff * St1;
%St1_bar = (St1 + St2)/2;

% Compute K1,K2,dK1
[K1]=fun_K1(Rc);
[K2,gamma0,gamma]=fun_K2(Ma,Rc,alpha);
[dK1]=fun_dK1(alpha,Rp);

% Compute the pressure levels for each center frequency

for i = 1:length(Freq)
    
Stp(i)=Freq(i)*delt_star_p/Vrel;
Sts(i)=Freq(i)*delt_star_s/Vrel;


% Compute SPL for pressure side
a = abs(log10(Stp(i)/St1));
[A]=fun_A(a,Rc);
SPLp(i) = A + (K1 - 3) + dK1 + 10 * log10(delt_star_p*Ma^5*L*Di_hi_fr/Distance^2);

% Compute SPL for suction side
a = abs(log10(Sts(i)/St1));
[A]=fun_A(a,Rc);
SPLs(i) = A + (K1 - 3) + 10 * log10(delt_star_s*Ma^5*L*Di_hi_fr/Distance^2);

% Compute SPL for nonzero angle of attack
b = abs(log10(Sts(i)/St2));
[B]=fun_B(b,Rc);
SPLa(i) = B + K2 + 10 * log10(delt_star_s*Ma^5*L*Di_hi_fr/Distance^2);

if alpha>=gamma0 | alpha>12.5
   SPLp(i) = 10 * log10(delt_star_s*Ma^5*L*Di_lo_fr/Distance^2);  % !!!!!!!!!!!!!
   SPLs(i) = 10 * log10(delt_star_s*Ma^5*L*Di_lo_fr/Distance^2);
   
   [B]=fun_A(b,3*Rc);
   SPLa(i) = B + K2 + 10 * log10(delt_star_s*Ma^5*L*Di_lo_fr/Distance^2);
end

if SPLp(i)<-100
    SPLp(i)=-100;
end
if SPLs(i)<-100
    SPLs(i)=-100;
end
if SPLa(i)<-100
    SPLa(i)=-100;
end    

P1(i) = 10^(SPLp(i)/10);
P2(i) = 10^(SPLs(i)/10);
P3(i) = 10^(SPLa(i)/10);

SPLTBL(i) = 10 * log10(P1(i) + P2(i) + P3(i));

end

%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

               % ==================================================== %
               %                    < LBL-VS >                        %
               %     Laminar-Boundary-Layer-Vortex-Shedding Noise     %
               % ==================================================== %
               
               
               % << All Variables >> %
% alpha         -->  Angle of attack
% c             -->  Chord length
% co            -->  Sound speed in air
% d             -->  Ratio of Reynolds number
% delt_p        -->  Pressure side boundary layer thickness
% delt_s        -->  Suction side boundary layer thickness
% delt_star_p   -->  Pressure side boundary layer displacement thickness
% delt_star_s   -->  Suction side boundary layer displacment thickness
% Distance      -->  The distance from observer to the blade segment
% Di_hi_fr      -->  High frequency directivity
% Di_lo_fr      -->  Low frequency directivity
% e             -->  Ratio of Strouhal number
% Freq          -->  1/3 octave frequencies
% G1            -->  Shape function
% G2            -->  Shape function
% G3            -->  Shape function
% ho            -->  Altitude of the observer referenced to the tower base 
% Htower        -->  The tower hight
% L             -->  The span of each blade segment
% logi_trip     -->  Condition for boundary layer tripping:
%                    tripped=0, untripped=1, partially tripped=2
% Ma            -->  Mach number
% PHI           -->  The sound directivity angle
% psi           -->  The angle of the observer referenced to the
%                    downstream direction
% r             -->  Blade element radius
% ro            -->  Distance between the observer and the tower center
% Rc            -->  Reynolds number based on chord
% Rco           -->  Reference Reynolds number
% shaftlength   -->  The shaft length of the wind turbine
% SPLLBL        -->  Sound pressure level caused by laminar boundary layer vortex shedding
% St1_prim      -->  Strouhal number
% Stpeak_prim   -->  Strouhal number
% THETA         -->  The sound directivity angle
% theta_yaw     -->  The yaw angle of the turbine
% theta_tilt    -->  The tilt angle of the turbine
% theta_cone    -->  The cone angle of the turbine
% theta_wing    -->  The wing angle of the turbine
% theta_twist   -->  The twist angle of the turbine
% Vrel          -->  The relative velocity seen by the blade segment

function[SPLLBL]=LBLVS(delt_p,delt_star_p,delt_s,delt_star_s,Di_hi_fr,Di_lo_fr,Distance,...
    Freq,Vrel,c,vis,co,logi_trip,alpha,ro,ho,psi,theta_yaw,theta_tilt,...
                 theta_cone,theta_wing,theta_twist,shaftlength,Htower,r,L)

% Compute Mach number and Reynolds number based on chord length 
Rc=Vrel*c/vis;
Ma=Vrel/co;

   

% Compute boundary layer thickness
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%[delt_p,delt_star_p,delt_s,delt_star_s]=Thickness(Vrel,c,vis,co,logi_trip,alpha);

% Compute sound directivity
% [THETA,PHI,Distance]=Trans_direc(ro,ho,psi,theta_yaw,theta_tilt,theta_cone,theta_wing,theta_twist,...
%     shaftlength,Htower,r);
% [Di_hi_fr,Di_lo_fr]=Directivity(Vrel,co,THETA,PHI);

% Compute Strouhal number
if Rc<=1.3e05
    St1_prim = 0.18;
elseif Rc>1.3e05 & Rc<=4.0e05
    St1_prim = 0.001756*Rc^0.3931;
else
    St1_prim = 0.28;
end

Stpeak_prim = St1_prim * 10^(-0.04*alpha);


% Compute the reference Reynolds number
if alpha <=3.0
    Rco = 10^(0.215*alpha+4.978);
else
    Rco = 10^(0.120*alpha+5.263);
end
% Compute the shape function 'G2'
d = Rc/Rco;
G2 = fun_G2(d);

% Compute the angle correction function 'G3'
G3 = 171.04 - 3.03 *alpha;

for i = 1:length(Freq)
    St_prim(i) = Freq(i) * delt_p/Vrel;
    e = St_prim(i)/Stpeak_prim;
    G1 = fun_G1(e);
    
    SPLLBL(i) = 10*log10(delt_p*Ma^5*L*Di_hi_fr/Distance^2)  + G1 + G2 + G3;
    
end
    
%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

               % =================================== %
               %               < TIP >               %
               %      Tip-Vortex Formation Noise     %
               % =================================== %
               
               
               % << All Variables >> %
               
% alpha         -->  Angle of attack
% alpha_        -->  Corrected angle of attack
% c             -->  Segment chord length
% co            -->  Sound speed in air
% Freq          -->  1/3 octave frequencies
% h             -->  Trailing edge bluntness
% ho            -->  Altitude of the observer referenced to the tower base 
% Htower        -->  The hight of the tower
% L             -->  Blade segment span
% Ma            -->  Mach number
% Mmax          -->  Maximum Mach number
% psi           -->  The observer angle reference to the downstream
%                    direction
% PSI           -->  The observer angle reference to the downstream
%                    direction
% r             -->  The blade radius at each segment
% Ratio         -->  The coffecient to refine alpha_tip
% Round         -->  Logical parameter. Rounded tip=0, squared tip=1
% Rc            -->  Reynolds number based on chord
% ro            -->  Distance between the observer and the tower center
% shaftlength   -->  The shaft length of the turbine
% SPLTIP        -->  Sound pressure level due to TIP
% St_prim2      -->  Strouhal number
% theta_cone    -->  Cone angle
% theta_tilt    -->  Tilt angle
% theta_twist   -->  Twist angle
% theta_wing    -->  Wing angle
% theta_yaw     -->  Yaw angle
% vis           -->  Kinematic viscosity
% Vmax          -->  Maximum velocity
% Vo            -->  The mean wind speed at the hub center
% Vrel          -->  The relative wind velovity

function[SPLTIP]=TIP(Di_hi_fr,Di_lo_fr,Distance,...
    Freq,Vrel,co,c,alpha,ro,ho,psi,theta_yaw,theta_tilt,theta_cone,theta_wing,theta_twist,...
    shaftlength,Htower,r,Round,Ratio)
Ma=Vrel/co;
alpha_ = Ratio * alpha;

% Compute sound directivity
% [THETA,PHI,Distance]=Trans_direc(ro,ho,psi,theta_yaw,theta_tilt,theta_cone,theta_wing,theta_twist,...
%     shaftlength,Htower,r);
% [Di_hi_fr,Di_lo_fr]=Directivity(Vrel,co,THETA,PHI);

if Round == 0    % If the tip is rounded.
    l = 0.008 * alpha_ *c;
else
    if alpha_>=0 & alpha_<=2
        l = 0.0230 + 0.0169 * alpha_;
    else
        l = 0.0378 + 0.0095* alpha_;
    end
end
Mmax = Ma * (1 + 0.036 * alpha);
Vmax = co * Mmax;
for i = 1:length(Freq)
    Stprim2(i) = Freq(i) * l/Vmax;
    SPLTIP(i) = 10*log10(Ma^2*Mmax^3*l^2*Di_hi_fr/Distance^2) - ...
    30.5 * (log10(Stprim2(i))+0.3)^2 +126;
end


%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

               % ==================================================== %
               %                    < TEB-VS >                        %
               %     Trailing-Edge-Bluntness-Vortex-Shedding Noise     %
               % ==================================================== %
               
               
               % << All Variables >> %
% alpha         -->  Angle of attack
% c             -->  Segment chord length
% co            -->  Sound speed in air
% delt_p        -->  Boundary layer thickness on pressure side
% delt_star_p   -->  Boundary layer displacement thickness on pressure side
% delt_s        -->  Boundary layer thickness on suction side
% delt_star_s   -->  Boundary layer displacement thickness on suction side
% delt_star_avg -->  Average boundary layer displacement thickness
% Freq          -->  1/3 octave frequencies
% G4            -->  Specific shape function
% G5            -->  Specific shape function
% G5_0          -->  Specific shape function at TE-angle of zero degree
% G5_14         -->  Specific shape function at TE-angle of 14 degree
% G5min         -->  Specific shape function of minimum h_delt 
% G5min_0       -->  Specific shape function of minimum h_delt at TE-angle
%                    of zero degree
% G5min_14      -->  Specific shape function of minimum h_delt at TE-angle
%                    of 14 degree
% h             -->  Trailing edge bluntness
% h_delt        -->  Bluntness over average displacment thickness
% h_delt_prim   -->  Modified h_delt
% ho            -->  Altitude of the observer referenced to the tower base 
% Htower        -->  The height of the tower
% L             -->  Blade segment span
% logi_trip     -->  Condition for boundary layer tripping:
%                    tripped=0, untripped=1, partially tripped=2
% Ma            -->  Mach number
% psi           -->  The observer angle reference to the downstream
%                    direction
% PSI           -->  The observer angle reference to the downstream
%                    direction
% r             -->  The blade radius at each segment
% rat_St        -->  Ratio of Strouhal number
% Rc            -->  Reynolds number based on chord
% ro            -->  Distance between the observer and the tower center
% shaftlength   -->  The shaft length of the turbine
% SPLTEB        -->  Sound pressure level due to TEBVS
% Stpeak_prim3  -->  Peak Strouhal number
% theta_cone    -->  Cone angle
% theta_tilt    -->  Tilt angle
% theta_twist   -->  Twist angle
% theta_wing    -->  Wing angle
% theta_yaw     -->  Yaw angle
% vis           -->  Kinematic viscosity
% Vo            -->  The mean wind speed at the hub center
% Vrel          -->  The relative wind velovity


function[SPLTEB]=TEBVS(delt_p,delt_star_p,delt_s,delt_star_s,Di_hi_fr,Di_lo_fr,Distance,...
    Freq,h,TE_angle,Vrel,c,vis,co,logi_trip,alpha,ro,ho,psi,theta_yaw,theta_tilt,...
    theta_cone,theta_wing,theta_twist,shaftlength,Htower,r,L)

% Compute Mach number and Reynolds number based on chord length 
Rc=Vrel*c/vis;
Ma=Vrel/co;

% Compute boundary layer thickness
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
%[delt_p,delt_star_p,delt_s,delt_star_s]=Thickness(Vrel,c,vis,co,logi_trip,alpha);
delt_star_avg = (delt_star_p + delt_star_s)/2;

h_delt = h / delt_star_avg;
h_delt_prim = 6.724  * h_delt^2 -4.019 * h_delt + 1.107;
% Compute sound directivity
% [THETA,PHI,Distance]=Trans_direc(ro,ho,psi,theta_yaw,theta_tilt,theta_cone,theta_wing,theta_twist,...
%     shaftlength,Htower,r);
% [Di_hi_fr,Di_lo_fr]=Directivity(Vrel,co,THETA,PHI);

% Compute peak Strouhal number
if h_delt>=0.2
    Stpeak_prim3 = (0.212-0.0045*TE_angle)/(1+0.235*(h_delt)^(-1)-0.0132*h_delt^(-2));
else
    Stpeak_prim3 = 0.1 * h_delt + 0.095 - 0.00243 * TE_angle;
end

G4 = fun_G4(h_delt,TE_angle);

for i = 1:length(Freq)
    
    St_prim3(i) = Freq(i) * h / Vrel;
    rat_St = St_prim3(i)/Stpeak_prim3;
    
    G5_14 = fun_G5(h_delt,rat_St);
    
    G5_0 = fun_G5(h_delt_prim,rat_St);
    
    G5 = G5_0 + 0.0714*TE_angle*(G5_14 - G5_0);
    
    G5min_14 = fun_G5(0.25,rat_St);
    
    G5min_0 = fun_G5(0.25,rat_St);
    
    G5min = G5min_0 + 0.0714*TE_angle*(G5min_14 - G5min_0);
    
    if G5>0
        G5 = 0;
    end
    
    if G5>G5min
        G5 = G5min;
    end
    
    SPLTEB(i) = 10*log10(h*Ma^5.5*L*Di_hi_fr/Distance^2) + G4 + G5 ;
    
end

%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

               % ================================= %
               %              < TB-INF >           %
               %       Turbulent Inflow  Noise     %
               % ================================= %
               
               
               % << All Variables >> %
% alpha         -->  Angle of attack
% Bsq           -->  Intemediate parameter
% c             -->  Chord length
% co            -->  Sound speed in air
% Density       -->  Air density
% Distance      -->  The distance between observer and the blade segment
% Di_hi_fr      -->  Directivity at high frequency
% Di_lo_fr      -->  Directivity at low frequency
% Freq          -->  1/3 octave frequencies
% ho            -->  Altitude of the observer referenced to the tower base
% Htower        -->  The height of the tower
% K             -->  Wave number for low frequencies
% Kh            -->  Wave number for high frequencies
% LFC           -->  Intemediate parameter
% Ma            -->  Mach number
% n             -->  Power law factor(or shear factor)
% Nfreq         -->  The total number of frequencies
% PHI           -->  Directivity angle
% psi           -->  Observer angle referenced to the downstream direction
% r             -->  The blade radius at each segment
% Rc            -->  Reynolds number based on chord
% ro            -->  Distance between the observer and the tower center
% shaftlength   -->  Shaft length of the turbine
% SPLINF        -->  Sound pressure level due turbulent inflow
% SPLINFh       -->  Sound pressure level at high frequencies due turbulent inflow
% SPLINFl       -->  Sound pressure level at low frequencies due turbulent inflow
% SPLINFmid     -->  Sound pressure level with a smooth curve
% Ssq           -->  Intemediate parameter
% THETA         -->  Directivity angle
% theta_cone    -->  Cone angle
% theta_tilt    -->  Tilt angle
% theta_twist   -->  Twist angle
% theta_wing    -->  Wing angle
% theta_yaw     -->  Yaw angle
% Tinten        -->  Turbulence intensity
% Tscale        -->  Turbulence length scale
% vis           -->  Kinematic viscosity
% Vo            -->  The mean wind speed at the hub center
% Vrel          -->  The relative wind velovity
% zo            -->  Surface roughness
% Z             -->  Elevation above ground

function[SPLINF,Tbscale,Tbinten]=TBINF(Di_hi_fr,Di_lo_fr,Distance,...
    Freq,Vo,Vrel,c,L,vis,co,ro,ho,psi,theta_yaw,theta_tilt,theta_cone,theta_wing,theta_twist,...
    shaftlength,Htower,r,Density,Tscale,Tinten,zo,Gamma)
% Compute Mach number and Reynolds number based on chord length 
Rc=Vrel*c/vis;
Ma=Vrel/co;

if Tscale==0    % If length scale is not given in the input window, do the calculation as following!
Z=Htower + r*cos(theta_wing*pi/180);
%Tscale=300*(Z/300)^(0.46+0.074*log(zo));  % Gives lager length scale
Tscale=25*Z^(0.35)*zo^(-0.063);             % Gives smaller length scale
end
Tbscale=Tscale;

if Tinten==0   % If intensity is not given in the input window, do the calculation as following!
Z=Htower + r*cos(theta_wing*pi/180); 
%Gamma=0.24+0.096*log10(zo)+0.016*(log10(zo))^2; % Can be computed or just give it a value.
Tinten=Gamma*log(30/zo)/log(Z/zo);
end
Tbinten=Tinten;

% Compute sound directivity
% [THETA,PHI,Distance]=Trans_direc(ro,ho,psi,theta_yaw,theta_tilt,theta_cone,theta_wing,theta_twist,...
%     shaftlength,Htower,r);
% [Di_hi_fr,Di_lo_fr]=Directivity(Vrel,co,THETA,PHI);


Nfreq = length(Freq);
% SPLINFh = zeros(1,Nfreq);
% SPLINFmid = zeros(1,Nfreq);
% SPLINFl = zeros(1,Nfreq);

for i = 1:Nfreq
     Kh(i) = 8*pi*Freq(i)*Tscale/(3*Vrel);
     K(i) = pi*Freq(i)*c/Vrel;
     Bsq(i) = 1 - Ma^2;
     Ssq(i) = (2*pi*K(i)/Bsq(i) + (1+2.4*K(i)/Bsq(i))^(-1))^(-1);
     LFC(i) = 10 * Ssq(i)*Ma*K(i)^2*Bsq(i)^(-1);

     SPLINF(i) = 10*log10(Di_hi_fr*Density^2*co^2*Tscale*L/2*Ma^3*Tinten^2*K(i)^3*(1+K(i)^2)^(-7/3)/Distance^2) + 58.4...
             + 10*log10(LFC(i)/(1+LFC(i)));   % % this is another choice!
%          SPLINF(i) = 10*log10(Di_hi_fr*Tscale*0.5*L*Ma^5*Tinten^2*Kh(i)^3*(1+Kh(i)^2)^(-7/3)/Distance^2) +...
%              10*log10(10^18.13) + 10*log10(Di_lo_fr*LFC(i)/(1+LFC(i)));
%      SPLINF(i) = SPLINFmid(i);
   
end

%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

              % ========================================================= %
              %    Coordinate system tranformation between local blade    %
              %      coordinate system and global coordinate system       %
              % ========================================================= %
                                 
                         % << All Variables >> %

                         
% a12           -->  Tranformation matrix of system 1 and system 2
% a23           -->  Tranformation matrix of system 2 and system 3
% a34           -->  Tranformation matrix of system 3 and system 4
% Htower        -->  The height of the tower
% r             -->  Blade element radius
% r1            -->  Vectors given in coordinate system 1
% r2            -->  Vectors given in coordinate system 2.
% r3            -->  Vectors given in coordinate system 3
% r4            -->  Vectors given in coordinate system 4
% rs            -->  Vectors given in coordinate system 2, denote the shaft
% rt            -->  Vectors given in coordinate system 1, denote the tower
% shaftlength   -->  The shaftlength
% theta_cone    -->  Cone angle
% theta_tilt    -->  Tile angle
% theta_wing    -->  Wing angle
% theta_yaw     -->  Yaw angle

% Hint: Used to transform any vectors from coordinate system 4 to coordinate
% system 1.(From local blade coordinate system to global coordinate system.) 

% << Coordinate systems definition >>
% System 1: --> Placed at the root of the tower
% System 2: --> Placed in the nacelle(non-rotation)
% System 3: --> Fixed to the rotating shaft
% System 4: --> Aligned with one of the blades

function[rpoint]=Sys4_sys1(theta_yaw,theta_tilt,theta_cone,theta_wing,shaftlength,Htower,r);

 a1=[1 0 0;0 cos(theta_yaw*(pi/180)) sin(theta_yaw*(pi/180));0 -sin(theta_yaw*(pi/180)) cos(theta_yaw*(pi/180))];
 a2=[cos(theta_tilt*(pi/180)) 0 sin(theta_tilt*(pi/180));0 1 0;-sin(theta_tilt*(pi/180)) 0 cos(theta_tilt*(pi/180))];
 a3=[1 0 0;0 1 0;0 0 1];
 a12=a3*a2*a1;
 a23=[cos(theta_wing*(pi/180)) sin(theta_wing*(pi/180)) 0;-sin(theta_wing*(pi/180)) cos(theta_wing*(pi/180)) 0;0 0 1];
 a34=[cos(theta_cone*(pi/180)) 0 sin(theta_cone*(pi/180));0 1 0;-sin(theta_cone*(pi/180)) 0 cos(theta_cone*(pi/180))];
    
 r2=[0;0;-shaftlength];
 r1=a12\r2;
 rs=r1;
 r4=[r;0;0];
 r3=a34\r4;  
 r2=a23\r3;  
 r1=a12\r2;
 rt=[Htower;0;0];
 rpoint=(rt+rs+r1)';
     
%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

              % ========================================= %
              %        Transformation of directivity      %
              % ========================================= %
             
% << All Variables >> %

% a12           -->  Tranformation matrix of system 1 and system 2
% a23           -->  Tranformation matrix of system 2 and system 3
% a34           -->  Tranformation matrix of system 3 and system 4
% a45           -->  Tranformation matrix of system 4 and system 5
% Htower        -->  The height of the tower
% n1,n2         -->  Unit vectors
% r             -->  Blade element radius
% r1,r2,r3,r4   -->  Vectors given in coordinate system 5
% ro            -->  Distance between the observer and the tower center
% ho            -->  Altitude of the observer referenced to the tower base 
% Obs_1         -->  Vector of observer position in system 1
% Obs_5         -->  Vector of observer position in system 5
% Obs_xz_5       -->  Project Obs_5 into xz plane
% PHI           -->  Directivity angle relative to x axis in system 5
% psi           -->  The angle of observer relative to the downstream
%                    direction
% shaftlength   -->  The shaftlength 
% THETA         -->  Directivity angle relative to y axis in system 5
% theta_cone    -->  Cone angle
% theta_tilt    -->  Tile angle
% theta_twist   -->  Twist angle
% theta_wing    -->  Wing angle
% theta_yaw     -->  Yaw angle

% Hint: The observer position is transformed from the globle coordinate
% system 1 into the local coordinate system 5. Then the directivity angles
% of THETA and PHI are calculated in the system 5. The details are discussed
% in the final report.

% << Coordinate systems definition >>
% System 1: --> Placed at the root of the tower
% System 2: --> Placed in the nacelle(non-rotation)
% System 3: --> Fixed to the rotating shaft
% System 4: --> Aligned with one of the blades
% System 5: --> Aligned with one of the segments of the blade

function[THETA,PHI,Distance]=Trans_direc(ro,ho,psi,theta_yaw,theta_tilt,theta_cone,theta_wing,theta_twist,...
    shaftlength,Htower,r)
psi = psi*pi/180;
x = ho;
y = ro*sin(psi);
z = ro*cos(psi);
Obs_1 = [x;y;z];

a1=[1 0 0;0 cos(theta_yaw*(pi/180)) sin(theta_yaw*(pi/180));0 -sin(theta_yaw*(pi/180)) cos(theta_yaw*(pi/180))];
a2=[cos(theta_tilt*(pi/180)) 0 sin(theta_tilt*(pi/180));0 1 0;-sin(theta_tilt*(pi/180)) 0 cos(theta_tilt*(pi/180))];
a3=[1 0 0;0 1 0;0 0 1];
a12=a3*a2*a1;
a23=[cos(theta_wing*(pi/180)) sin(theta_wing*(pi/180)) 0;-sin(theta_wing*(pi/180)) cos(theta_wing*(pi/180)) 0;0 0 1];
a34=[cos(theta_cone*(pi/180)) 0 sin(theta_cone*(pi/180));0 1 0;-sin(theta_cone*(pi/180)) 0 cos(theta_cone*(pi/180))];
a45=[1 0 0;0 cos(theta_twist*(pi/180)) sin(theta_twist*(pi/180));0 -sin(theta_twist*(pi/180)) cos(theta_twist*(pi/180))];

r1 = a45*a34*a23*a12*Obs_1;
r2 = a45*a34*a23*[-Htower*cos(theta_tilt*pi/180);0;Htower*sin(theta_tilt*pi/180)];
r3 = a45*a34*[shaftlength*sin(theta_cone*pi/180);0;shaftlength*cos(theta_cone*pi/180)];
r4 = a45*[-r;0;0];

Obs_5 = r1+r2+r3+r4;
Distance = norm(Obs_5);
n1 = [0,1,0];
THETA = acos(abs(n1*Obs_5)/(norm(Obs_5)))*180/pi;
n2 = [1,0,0];
Obs_xz_5 = [Obs_5(1);0;Obs_5(3)];
PHI = acos(abs(n2*Obs_xz_5)/(norm(Obs_xz_5)))*180/pi;

%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

              % ======================================= %
              %        Sound directivity parameter      %
              % ======================================= %
                                 
                         % << All Variables >> %
                         
% co            -->  Sound speed in the air under normal condition
% Di_hi_fr      -->  Sound directivity at high frequency
% Di_lo_fr      -->  Sound directivity at low  frequency
% ho            -->  Altitude of the observer referenced to the tower base
% Htower        -->  The height of the tower
% Ma            -->  Mach number 
% Mc            -->  Mach number for the flow past the plate trailing edge
% psi           -->  The angle of observer relative to the downstream
%                    direction
% PHI           -->  Directivity angle relative to x axis in system 5
% ro            -->  Distance between the observer and the tower base
% shaftlength   -->  The shaftlength 
% THETA         -->  Directivity angle relative to y axis in system 5
% theta_yaw     -->  Yaw angle
% theta_tilt    -->  Tile angle
% theta_cone    -->  Cone angle
% theta_wing    -->  Wing angle
% theta_twist   -->  Twist angle
% Vrel          -->  Relative velocity

function[Di_hi_fr,Di_lo_fr]=Directivity(Vrel,co,THETA,PHI)
Ma = Vrel/co;
Mc = 0.8*Ma;

% Compute directivity angles THETA and PHI
THETA = THETA * pi/180;
PHI = PHI * pi/180;

Di_hi_fr = 2*(sin(0.5*THETA))^2*(sin(PHI))^2/((1+Ma*cos(THETA))*(1+(Ma-Mc)*cos(THETA))^2);

Di_lo_fr = (sin(THETA))^2*(sin(PHI))^2/(1+Ma*cos(THETA))^4;

if Di_hi_fr==0  
    Di_hi_fr = 1e-5;
end
if Di_lo_fr==0
    Di_lo_fr = 1e-5;
end
    
%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

             % ============================================================= %
             %    Boundary layer thickness parameter at the trailing edge    %
             % ============================================================= %
             
                          % << All Variables >> %
% alpha         -->  Angle of attack             
% c             -->  Chord length
% co            -->  Sound speed
% delt_0        -->  Boundary layer thickness of zero angle of attack
% delt_star_0   -->  Boundary layer displacement thickness of zero angle of attack
% delt_p        -->  Boundary layer thickness on pressure side of nonzero angle of attack
% delt_s        -->  Boundary layer thickness on suction side of nonzero angle of attack
% delt_star_p   -->  Boundary layer displacement thickness on pressure side of nonzero angle of attack
% delt_star_s   -->  Boundary layer displacement thickness on suction side of nonzero angle of attack
% logi_trip     -->  Condition for boundary layer tripping:
%                    tripped=0, untripped=1, partially tripped=2
% Ma            -->  Mach number
% Rc            -->  Reynolds number based on chord
% vis           -->  kinematic viscosity
% Vrel          -->  Relative velocity
function[delt_p,delt_star_p,delt_s,delt_star_s]=Thickness(Vrel,c,vis,co,logi_trip,alpha)
Rc=Vrel*c/vis;
Ma=Vrel/co;

%% Compute boundary layer thickness of zero angle of attack: 
% (1) For the heavily tripped boundary layer
delt_0 = c * 10^(1.892-0.9045*log10(Rc)+0.0596*(log10(Rc))^2);
% (2) For the untripped boundary layer
if logi_trip==1
    delt_0 = c * 10^(1.6569-0.9045*log10(Rc)+0.0596*(log10(Rc))^2);
% (3) For the partially tripped boundary layer
elseif logi_trip==2
    delt_0 = 0.6 * c * 10^(1.6569-0.9045*log10(Rc)+0.0596*(log10(Rc))^2);
end

%% Compute boundary layer displacement thickness of zero angle of attack:
if logi_trip==0 
    if Rc<=0.3e06
        delt_star_0 = c * 0.0601 * Rc^(-0.114);
    elseif Rc>=0.3e06 
        delt_star_0 = c * 10^(3.411-1.5397*log10(Rc)+0.1059*(log10(Rc))^2);
    end
elseif logi_trip==2
   if Rc<=0.3e06
        delt_star_0 = c * 0.0601 * Rc^(-0.114);
    elseif Rc>=0.3e06 
        delt_star_0 = c * 10^(3.411-1.5397*log10(Rc)+0.1059*(log10(Rc))^2);
    end
    delt_star_0 = 0.6 * delt_star_0;
elseif logi_trip==1
    delt_star_0 = c * 10^(3.0187-1.5397*log10(Rc)+0.1059*(log10(Rc))^2);
end

%% Compute boundary layer (displacement)thickness of nonzero angle of attack:
% Compute pressure side boundary layer (displacment)thickness:
delt_p = delt_0 * 10^(-0.04175*alpha+0.00106*alpha^2);
delt_star_p = delt_star_0 * 10^(-0.0432*alpha+0.00113*alpha^2);
% Compute suction side boundary layer (displacment)thickness:
if logi_trip==0  % tripped
    if alpha>=0 & alpha<=5
        delt_s = delt_0 * 10^(0.0311*alpha);
        delt_star_s = delt_star_0 * 10^(0.0679*alpha);
    elseif alpha>5 & alpha<=12.5
        delt_s = delt_0 * 0.3468 *10^(0.1231*alpha);
        delt_star_s = delt_star_0 * 0.381 * 10^(0.1516*alpha);
    else%if alpha>12.5 & alpha<=25
        delt_s = delt_0 * 5.718 *10^(0.0258*alpha);
        delt_star_s = delt_star_0 * 14.296 * 10^(0.0258*alpha);
    end
else
    if alpha>=0 & alpha<=7.5
        delt_s = delt_0 * 10^(0.03114*alpha);
        delt_star_s = delt_star_0 * 10^(0.0679*alpha);
    elseif alpha>7.5 & alpha<=12.5
        delt_s = delt_0 * 0.0303 * 10^(0.2336*alpha);
        delt_star_s = delt_star_0 * 0.0162 * 10^(0.3066*alpha);
    else%if alpha>12.5 & alpha<=25
        delt_s = delt_0 * 12 *10^(0.0258*alpha);
        delt_star_s = delt_star_0 * 52.42 * 10^(0.0258*alpha);
    end
end

%   [delt_p,delt_star_p,delt_s,delt_star_s]=thickness(10,1,1.82*1e-5,331,0,5)

%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
              % ======================================= %
              %      The spectral shape function A      %
              % ======================================= %
               
              % << All Variables >> %
% a             -->  Absolute value of the logarithm of Strouhal number
% A             -->  'A' is the spectral shape function
% Amax          -->  Defined for the interpolation
% Amin          -->  Defined for the interpolation
% alpha         -->  Angle of attack
% ao            -->  Defined for the interpolation
% Ar            -->  Interpolation factor
% Rc            -->  Reynolds number based on chord

% Hint: Function 'A' is an important parameter related with the trubulent
% boundary layer trailing edge noise and separated flow noise.

function [A]=fun_A(a,Rc)

if Rc<9.52e04
    ao = 0.57;
elseif Rc>=9.52e04 & Rc<=8.57e05
    ao = (-9.57e-13)*(Rc-8.57e05)^2 + 1.13;
else
    ao = 1.13;
end

if ao<0.204
    Amin = sqrt(67.552-886.788*ao^2)-8.219;
elseif ao>=0.204 & ao<=0.244
    Amin = -32.665*ao + 3.981;
else
    Amin = -142.795*ao^3 + 103.656*ao^2 - 57.757*ao + 6.006;
end

if ao<0.13
    Amax = sqrt(67.552-886.788*ao^2)-8.219;
elseif ao>=0.13 & ao<=0.321
    Amax = -15.901*ao + 1.098;
else
    Amax = -4.669*ao^3 + 3.491*ao^2 - 16.699*ao + 1.149;
end
Ar = (-20-Amin)/(Amax-Amin);

%----------------------------------------
if a<0.204
    Amin = sqrt(67.552-886.788*a^2)-8.219;
elseif a>=0.204 & a<=0.244
    Amin = -32.665*a + 3.981;
else
    Amin = -142.795*a^3 + 103.656*a^2 - 57.757*a + 6.006;
end

if a<0.13
    Amax = sqrt(67.552-886.788*a^2)-8.219;
elseif a>=0.13 & a<=0.321
    Amax = -15.901*a + 1.098;
else
    Amax = -4.669*a^3 + 3.491*a^2 - 16.699*a + 1.149;
end

A = Amin + Ar * (Amax-Amin);

%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

              % ======================================= %
              %      The spectral shape function B      %
              % ======================================= %
               
              % << All Variables >> %
% b             -->  Absolute value of the logarithm of Strouhal number
% B             -->  'B' is the spectral shape function
% Bmax          -->  Defined for the interpolation
% Bmin          -->  Defined for the interpolation
% alpha         -->  Angle of attack
% bo            -->  Defined for the interpolation
% Br            -->  Interpolation factor
% Rc            -->  Reynolds number based on chord

% Hint: Function 'B' is an important parameter related with the trubulent
% boundary layer trailing edge noise and separated flow noise.

function [B]=fun_B(b,Rc)

if Rc<9.52e04
    bo = 0.30;
elseif Rc>=9.52e04 & Rc<=8.57e05
    bo = (-4.48e-13)*(Rc-8.57e05)^2 + 0.56;
else
    bo = 0.56;
end

if bo<0.13
    Bmin = sqrt(16.888-886.788*bo^2)-4.109;
elseif bo>=0.13 & bo<=0.145
    Bmin = -83.607*bo + 8.138;
else
    Bmin = -817.810*bo^3 + 355.210*bo^2 - 135.024*bo + 10.619;
end

if bo<0.10
    Bmax = sqrt(16.888-886.788*bo^2)-4.109;
elseif bo>=0.10 & bo<=0.187
    Bmax = -31.330*bo + 1.854;
else
    Bmax = -80.541*bo^3 + 44.174*bo^2 - 39.381*bo + 2.344;
end
Br = (-20-Bmin)/(Bmax-Bmin);

%----------------------------------------
if b<0.13
    Bmin = sqrt(16.888-886.788*b^2)-4.109;
elseif b>=0.13 & b<=0.145
    Bmin = -83.607*b + 8.138;
else
    Bmin = -817.810*b^3 + 355.210*b^2 - 135.024*b + 10.619;
end

if b<0.10
    Bmax = sqrt(16.888-886.788*b^2)-4.109;
elseif b>=0.10 & b<=0.187
    Bmax = -31.330*b + 1.854;
else
    Bmax = -80.541*b^3 + 44.174*b^2 - 39.381*b + 2.344;
end

B = Bmin + Br * (Bmax-Bmin);

%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

              % ==================================== %
              %      The adjustment function dK1     %
              % ==================================== %
               
              % << All Variables >> %
% alpha         -->  Angle of attack               
% dK1           -->  Adjustment function
% Rp            -->  Reynolds number based on the pressure-side
%                    displacement thickness

% Hint: dK1 represents the adjustment for the pressure side contribution for
% nonzero angles of attack in case of trubulent
% boundary layer trailing edge noise and separated flow noise.

function [dK1]=fun_dK1(alpha,Rp)

if Rp<=5000
    dK1 = alpha * (1.43 * log10(Rp) - 5.29);
else
    dK1 = 0;
end

%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

              % ======================================= %
              %      The spectral shape function G1     %
              % ======================================= %
              
             % << All Variables >> %

% e             -->  Strouhal number ratio
% G1            -->  Shape function

% Hint: Function 'G1' is the function related with the laminar
% boundary layer vortex shedding noise.

function[G1]=fun_G1(e)

if e<=0.5974
    G1 = 39.8*log10(e) - 11.12;
elseif e>0.5974 & e<=0.8545
    G1 = 98.409*log10(e) + 2.0;
elseif e>0.8545 & e<=1.17
    G1 = -5.076 + sqrt(2.484 - 506.25*(log10(e))^2);
elseif e>1.17 & e<=1.674
    G1 = -98.409*log10(e) + 2.0;
else
    G1 = -39.8*log10(e) - 11.12;
end

%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

              % ======================================= %
              %      The spectral shape function G2     %
              % ======================================= %
              
             % << All Variables >> %

% d             -->  Reynolds number ratio
% G2            -->  Shape function


% Hint: Function 'G2' is the function related with the laminar
% boundary layer vortex shedding noise.

function[G2]=fun_G2(d)

if d<=0.3237
    G2 = 77.852*log10(d) + 15.328;
elseif d>0.3237 & d<=0.5689
    G2 = 65.188*log10(d) + 9.125;
elseif d>0.5689 & d<=1.7579
    G2 = -114.052*(log10(d))^2;
elseif d>1.7579 & d<=3.0889
    G2 = -65.188*log10(d) + 9.125;
else
    G2 = -77.852*log10(d) + 15.328;
end

%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

              % ======================================= %
              %      The spectral shape function G4     %
              % ======================================= %
              
             % << All Variables >> %

% G4            -->  Shape function
% h_delt        -->  Bluntness over average displacement thickness
% TE_angle      -->  Trailing edge angle

% Hint: Function 'G4' is the function related with the trailing
% edge bluntness vortex shedding noise.

function[G4]=fun_G4(h_delt,TE_angle)

if h_delt<=5
    G4 = 17.5*log10(h_delt) + 157.5 - 1.114*TE_angle;
else
    G4 = 169.7 - 1.114*TE_angle;
end

%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

              % ======================================= %
              %      The spectral shape function G5     %
              % ======================================= %
              
             % << All Variables >> %

% eta           -->  Strouhal number ratio
% eta0          -->  Reference Strouhal number ratio
% G5            -->  Shape funtion
% h_delt        -->  Bluntness over average displacement thickness
% m,k,mu        -->  Intermediate parameter
% rat_St        -->  Ratio of Strouhal number

% Hint: Function 'G5' is the function related with the trailing
% edge bluntness vortex shedding noise.

function[G5x]=fun_G5(h_delt,rat_St)

eta = log10(rat_St);

% Compute 'mu'
if h_delt<0.25
    mu = 0.1221;
elseif h_delt>=0.25 & h_delt<0.62
    mu = -0.2175 * h_delt + 0.1755;
elseif h_delt>=0.62 & h_delt<1.15
    mu = -0.0308 * h_delt + 0.0596;
else
    mu = 0.0242;
end

% Compute 'm'
if h_delt<=0.02
    m = 0;
elseif h_delt>0.02 & h_delt<=0.5
    m = 68.724 * h_delt - 1.35;
elseif h_delt>0.5 & h_delt<=0.62
    m = 308.475 * h_delt - 121.23;
elseif h_delt>0.62 & h_delt<=1.15
    m = 224.811 * h_delt - 69.35;
elseif h_delt>1.15 & h_delt<=1.2
    m = 1583.28 * h_delt - 1631.59;
else
    m = 268.344;
end
      
eta0 = -sqrt(m^2*mu^4/(6.25+m^2*mu^2));
k = 2.5*sqrt(1 - (eta0/mu)^2) - 2.5 - m*eta0;

if eta<eta0
    G5x = m*eta + k;
elseif eta>=eta0 & eta<0
    G5x = 2.5 * sqrt(1 - (eta/mu)^2) - 2.5;
elseif eta>=0 & eta<0.03616
    G5x = sqrt(1.5625 - 1194.99*eta^2) - 1.25;
else
    G5x = -155.543*eta + 4.375;
end

%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

              % =================================== %
              %      The amplitude function K1      %
              % ==================================== %
               
              % << All Variables >> %
% K1            -->  Amplitude function
% Rc            -->  Reynolds number based on chord

% Hint: K1 is the amplitude function related with the trubulent
% boundary layer trailing edge noise and separated flow noise.

function[K1]=fun_K1(Rc)

if Rc<2.47e05
    K1 = -4.31*log10(Rc) + 156.3;
elseif Rc>=2.47e05 & Rc<=8.0e05
    K1 = -9.0*log10(Rc) + 181.6;
else
    K1 = 128.5;
end

%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

              % =================================== %
              %      The amplitude function K2      %
              % ==================================== %
               
              % << All Variables >> %
% alpha         -->  Angle of attack
% beta          -->  Intermediate parameter
% beta0         -->  Intermediate parameter
% gamma         -->  Intermediate parameter
% gamma0        -->  Intermediate parameter
% K1            -->  Amplitude function
% K2            -->  Amplitude function
% Ma            -->  Mach number
% Rc            -->  Reynolds number based on chord

% Hint: K2 is the amplitude function related with the trubulent
% boundary layer trailing edge noise and separated flow noise.

function[K2,gamma0,gamma]=fun_K2(Ma,Rc,alpha)

gamma = 27.094 * Ma + 3.31;
gamma0 = 23.43 * Ma + 4.651;
beta = 72.65 * Ma + 10.74;
beta0 = -34.19 * Ma - 13.82;

[K1]=fun_K1(Rc);

if alpha<gamma0-gamma
    K2 = K1 - 1000;
elseif alpha>=gamma0-gamma & alpha<=gamma0+gamma
    K2 = K1 + sqrt(beta^2-(beta/gamma)^2*(alpha-gamma0)^2) + beta0;
else
    K2 = K1 - 12;
end

%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

                          % ================================ %
                          %    The relative velocity ver1    %
                          % ================================ %
                          % no wake
% << All Variables >> %

% alpha         -->  Angle of attack
% a1,...        -->  The transformation matrix
% Blade_geometry-->  The input data file gives the geometry informations of
%                    the blade 
% Htower        -->  The hight of the tower
% Gamma         -->  Paremeter denote the amount of shear,range between
%                    0.1-0.25
% M             -->  The total number of blade segments
% N             -->  The total number of azimuth positions which will be tested
% omega         -->  Angular totating speed of the blade
% phi           -->  The flow angle
% r             -->  Blade element radius
% R             -->  The array denotes the radius of all the blade segments
% radius        -->  The section radius of the tower at a certain hight
% r_distance    -->  The distance to the tower center line.
% r1,...        -->  The vectors in different coordinate systems
% rpoint        -->  The point in coordinate system 4
% shaftlength   -->  The shaftlength 
% theta_cone    -->  Cone angle
% theta_tilt    -->  Tile angle
% theta_wing    -->  Wing angle
% theta_yaw     -->  Yaw angle
% Vo            -->  Mean velocity at the hub hight
% Vr            -->  The wind velocity in radical direction.
% Vrel          -->  The length of relative velocity vector 
% V_rel         -->  The relative velocity vector 
% Vrely         -->  The y-component of relative velocity
% Vrelz         -->  The z-component of relative velocity
% Vrote         -->  The rotating speed seen by the blade
% Vshear        -->  The vector form of the oncoming wind velocity
% Vshear_z      -->  Atmosperic boundary layer wind speed in z-axis
%                    direction (in coordinate system 5). The other two
%                    components of wind speed are zero simply because we
%                    have defined the z-axis aligned with the wind
%                    direction.
% Vtheta        -->  The wind velocity in tangential direction.
% Vtower,Vwind  -->  The wind velocity vecotor effected by the tower
% Vx            -->  The wind velocity in y-axis direction
% Vy            -->  The wind velocity in y-axis direction
% Vz            -->  The wind velocity in z-axis direction
% wing          -->  The array denotes the wing angles
% x,y,z         -->  The position of vector 'rpoint' in three axis

% Hint: In this part, the atmospheric boundary layer wind speed is
% considered. The wind is also influenced by the tower, a model for the
% inluence of the tower is to assume potential flow. Further more, the
% tower shape model is needed as one of the input parameters.

function[Vrel,alpha]=V_relative(theta_yaw,theta_tilt,theta_cone,theta_wing,theta_twist,a_top,a_ground,...
                          shaftlength,Gamma,Htower,Vo,omega,r,zo);

   top_height = 0.05*Htower;
%    if Gamma==0  % If the shear factor is not given in the input window, do following computation!
%    Gamma=0.24+0.096*log10(zo)+0.016*(log10(zo))^2; % This is one way to  calculate Gamma.
%    end
%load Blade_geometry.dat;
%R = Blade_geometry(:,2);
%M = length(R);
%N = length(wing);

%for n=1:N   
%     theta_wing=wing(n);
%   for m=1:M        
%        r=R(m);
        rpoint=Sys4_sys1(theta_yaw,theta_tilt,theta_cone,theta_wing,shaftlength,Htower,r);
        x=rpoint(1);
        y=rpoint(2);
        z=rpoint(3);
        Vshear_z=Vo*(x/Htower)^Gamma;
        Vshear=[0;0;Vshear_z];
        
        r_distance=sqrt(z^2+y^2);
        if x<Htower    % the effect of the tower    
           radius=a_ground*(Htower+40-x)/(Htower+40);  
           
            Vr=Vshear_z*(1-(radius/r_distance)^2)*(z/r_distance);
            Vtheta=-Vshear_z*(1+(radius/r_distance)^2)*(-y/r_distance);
            Vz=Vr*(z/r_distance)-Vtheta*(-y/r_distance);
            Vy=-Vr*(-y/r_distance)-Vtheta*(z/r_distance);
            Vtower=[0;Vy;Vz];
            Vwind = Vtower; 
                  
           elseif x<Htower + top_height & x>=Htower  
            radius=a_top*(top_height-(x-Htower))/top_height;
             
            Vr=Vshear_z*(1-(radius/r_distance)^2)*(z/r_distance);
            Vtheta=-Vshear_z*(1+(radius/r_distance)^2)*(-y/r_distance);
            Vz=Vr*(z/r_distance)-Vtheta*(-y/r_distance);
            Vy=-Vr*(-y/r_distance)-Vtheta*(z/r_distance);         
            Vtower=[0;Vy;Vz];
            Vwind = Vtower;
                       
            else  
             
             Vwind = Vshear;
             
            end
             Vwind_x=Vwind(1);
             Vwind_y=Vwind(2);
             Vwind_z=Vwind(3);
             
a1=[1 0 0;0 cos(theta_yaw*(pi/180)) sin(theta_yaw*(pi/180));0 -sin(theta_yaw*(pi/180)) cos(theta_yaw*(pi/180))];
a2=[cos(theta_tilt*(pi/180)) 0 sin(theta_tilt*(pi/180));0 1 0;-sin(theta_tilt*(pi/180)) 0 cos(theta_tilt*(pi/180))];
a3=[1 0 0;0 1 0;0 0 1];
a12=a3*a2*a1;
a23=[cos(theta_wing*(pi/180)) sin(theta_wing*(pi/180)) 0;-sin(theta_wing*(pi/180)) cos(theta_wing*(pi/180)) 0;0 0 1];
a34=[cos(theta_cone*(pi/180)) 0 sin(theta_cone*(pi/180));0 1 0;-sin(theta_cone*(pi/180)) 0 cos(theta_cone*(pi/180))];
  
r1=[Vwind_x;Vwind_y;Vwind_z];
r2=a12*r1;
r3=a23*r2;
r4=a34*r3;
rp=r4;

Vw=rp;
Vx=rp(1);
Vy=rp(2);
Vz=rp(3);   
Vrote=[0;-r*omega*cos(theta_cone);0];
V_rel = Vw + Vrote;     
Vrely = V_rel(2);
Vrelz = V_rel(3);
phi = atan(-Vrelz/Vrely)*180/pi;
alpha = phi - theta_twist;
Vrel=norm(V_rel);
%   end
%end


% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

                          % ================================ %
                          %    The relative velocity ver2    %
                          % ================================ %
                          % consider Wn, Wt
% << All Variables >> %

% alpha         -->  Angle of attack
% a1,...        -->  The transformation matrix
% Blade_geometry-->  The input data file gives the geometry informations of
%                    the blade 
% Htower        -->  The hight of the tower
% Gamma         -->  Paremeter denote the amount of shear,range between
%                    0.1-0.25
% M             -->  The total number of blade segments
% N             -->  The total number of azimuth positions which will be tested
% omega         -->  Angular totating speed of the blade
% phi           -->  The flow angle
% r             -->  Blade element radius
% R             -->  The array denotes the radius of all the blade segments
% radius        -->  The section radius of the tower at a certain hight
% r_distance    -->  The distance to the tower center line.
% r1,...        -->  The vectors in different coordinate systems
% rpoint        -->  The point in coordinate system 4
% shaftlength   -->  The shaftlength 
% theta_cone    -->  Cone angle
% theta_tilt    -->  Tile angle
% theta_wing    -->  Wing angle
% theta_yaw     -->  Yaw angle
% Vo            -->  Mean velocity at the hub hight
% Vr            -->  The wind velocity in radical direction.
% Vrel          -->  The length of relative velocity vector 
% V_rel         -->  The relative velocity vector 
% Vrely         -->  The y-component of relative velocity
% Vrelz         -->  The z-component of relative velocity
% Vrote         -->  The rotating speed seen by the blade
% Vshear        -->  The vector form of the oncoming wind velocity
% Vshear_z      -->  Atmosperic boundary layer wind speed in z-axis
%                    direction (in coordinate system 5). The other two
%                    components of wind speed are zero simply because we
%                    have defined the z-axis aligned with the wind
%                    direction.
% Vtheta        -->  The wind velocity in tangential direction.
% Vtower,Vwind  -->  The wind velocity vecotor effected by the tower
% Vx            -->  The wind velocity in y-axis direction
% Vy            -->  The wind velocity in y-axis direction
% Vz            -->  The wind velocity in z-axis direction
% wing          -->  The array denotes the wing angles
% x,y,z         -->  The position of vector 'rpoint' in three axis

% Hint: In this part, the atmospheric boundary layer wind speed is
% considered. The wind is also influenced by the tower, a model for the
% inluence of the tower is to assume potential flow. Further more, the
% tower shape model is needed as one of the input parameters.

function[Vrel,alpha]=V_relative_WK(theta_yaw,theta_tilt,theta_cone,theta_wing,theta_twist,a_top,a_ground,...
                          shaftlength,Gamma,Htower,Vo,omega,r,zo,R);

   top_height = 0.05*Htower;
%    if Gamma==0  % If the shear factor is not given in the input window, do following computation!
%    Gamma=0.24+0.096*log10(zo)+0.016*(log10(zo))^2; % This is one way to  calculate Gamma.
%    end
%load Blade_geometry.dat;
%R = Blade_geometry(:,2);
%M = length(R);
%N = length(wing);

%for n=1:N   
%     theta_wing=wing(n);
%   for m=1:M        
%        r=R(m);
        rpoint=Sys4_sys1(theta_yaw,theta_tilt,theta_cone,theta_wing,shaftlength,Htower,r);
        x=rpoint(1);
        y=rpoint(2);
        z=rpoint(3);
        Vshear_z=Vo*(x/Htower)^Gamma;
        Vshear=[0;0;Vshear_z];
        
        r_distance=sqrt(z^2+y^2);
        if x<Htower    % the effect of the tower    
           radius=a_ground*(Htower+40-x)/(Htower+40);  
           
            Vr=Vshear_z*(1-(radius/r_distance)^2)*(z/r_distance);
            Vtheta=-Vshear_z*(1+(radius/r_distance)^2)*(-y/r_distance);
            Vz=Vr*(z/r_distance)-Vtheta*(-y/r_distance);
            Vy=-Vr*(-y/r_distance)-Vtheta*(z/r_distance);
            Vtower=[0;Vy;Vz];
            Vwind = Vtower; 
                  
           elseif x<Htower + top_height & x>=Htower  
            radius=a_top*(top_height-(x-Htower))/top_height;
             
            Vr=Vshear_z*(1-(radius/r_distance)^2)*(z/r_distance);
            Vtheta=-Vshear_z*(1+(radius/r_distance)^2)*(-y/r_distance);
            Vz=Vr*(z/r_distance)-Vtheta*(-y/r_distance);
            Vy=-Vr*(-y/r_distance)-Vtheta*(z/r_distance);         
            Vtower=[0;Vy;Vz];
            Vwind = Vtower;
                       
            else  
             
             Vwind = Vshear;
             
            end
             Vwind_x=Vwind(1);
             Vwind_y=Vwind(2);
             Vwind_z=Vwind(3);
             
a1=[1 0 0;0 cos(theta_yaw*(pi/180)) sin(theta_yaw*(pi/180));0 -sin(theta_yaw*(pi/180)) cos(theta_yaw*(pi/180))];
a2=[cos(theta_tilt*(pi/180)) 0 sin(theta_tilt*(pi/180));0 1 0;-sin(theta_tilt*(pi/180)) 0 cos(theta_tilt*(pi/180))];
a3=[1 0 0;0 1 0;0 0 1];
a12=a3*a2*a1;
a23=[cos(theta_wing*(pi/180)) sin(theta_wing*(pi/180)) 0;-sin(theta_wing*(pi/180)) cos(theta_wing*(pi/180)) 0;0 0 1];
a34=[cos(theta_cone*(pi/180)) 0 sin(theta_cone*(pi/180));0 1 0;-sin(theta_cone*(pi/180)) 0 cos(theta_cone*(pi/180))];
  
r1=[Vwind_x;Vwind_y;Vwind_z];
r2=a12*r1;
r3=a23*r2;
r4=a34*r3;
rp=r4;

%%%%% read in the dynamic wake data in z and y direction %%%%%
fid=fopen('wakedata','r');
res=fscanf(fid,'%g \n',[3,inf]);
fclose(fid);
rr = res(1,:);
Wn = res(2,:);
Wt = res(3,:);
%plot(Wn),hold on,plot(Wt)
% perform linear interpolation at r position
%Wn(1:round(0.33*length(Wn))) = 0;
%%%% read done %%%%
if r<=max(rr) & r>=min(rr)
    Wz = interp1(rr,Wn,r);
    Wy = interp1(rr,Wt,r);
elseif r>max(rr)
    Wz = Wn(length(Wn));; 
    Wy = Wt(length(Wt));;
elseif r<min(rr)
    Wz = Wn(1); 
    Wy = Wt(1);
end
Vw=rp;
Vx=rp(1);
Vy=rp(2);
Vz=rp(3);   
Vrote=[0;-r*omega*cos(theta_cone);0];
V_rel = Vw + Vrote;     
Vrely = V_rel(2)-Wy;
%Vrelz = V_rel(3)-Wz;
Vrelz = V_rel(3)*(1-Wz/Vo);
phi = atan(-Vrelz/Vrely)*180/pi;
alpha = phi - theta_twist;
Vrel=norm(V_rel);
%   end
%end


% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
             % ============================================================= %
             %    Boundary layer thickness parameter at the trailing edge    %
             % ============================================================= %
             
                          % << All Variables >> %
% alpha         -->  Angle of attack             
% c             -->  Chord length
% co            -->  Sound speed
% delt_0        -->  Boundary layer thickness of zero angle of attack
% delt_star_0   -->  Boundary layer displacement thickness of zero angle of attack
% delt_p        -->  Boundary layer thickness on pressure side of nonzero angle of attack
% delt_s        -->  Boundary layer thickness on suction side of nonzero angle of attack
% delt_star_p   -->  Boundary layer displacement thickness on pressure side of nonzero angle of attack
% delt_star_s   -->  Boundary layer displacement thickness on suction side of nonzero angle of attack
% logi_trip     -->  Condition for boundary layer tripping:
%                    tripped=0, untripped=1, partially tripped=2
% Rc            -->  Reynolds number based on chord
% vis           -->  kinematic viscosity
% Vrel          -->  Relative velocity
function[delt_p,delt_star_p,delt_s,delt_star_s]=Thickness_xfoil(Airfoil,Vrel,R,r,c,vis,co,logi_trip,alpha,...
    noisedatabase,Rref,Aref)
% clc
% clear

%format long
% logi_trip=1;
%Airfoil='NACA6320';
% Vrel=100;
% c=1;
% vis=1.8e-5;
% co=343;
% alpha=7;

Rc=Vrel*c/vis;
[rp,rs]=Thickness_ratio(logi_trip,Rc,alpha);

for k=1:3
if Airfoil(k)==1         
    if logi_trip==1
      [P,S]=P1S1NACA0012;
%     load P1NACA0012.mat;
%     load S1NACA0012.mat;
    elseif logi_trip==0
      [P,S]=P0S0NACA0012;
%     load P0NACA0012.mat;
%     load S0NACA0012.mat;   
    end
elseif Airfoil(k)==2          
    if logi_trip==1
       [P,S]=P1S1NACA6320;
%     load P1NACA6320.mat;
%     load S1NACA6320.mat;
    elseif logi_trip==0
       [P,S]=P0S0NACA6320;
%     load P0NACA6320.mat;
%     load S0NACA6320.mat;
    end
elseif Airfoil(k)==3          
    if logi_trip==1
       [P,S]=P1S1NACA63212;
%     load P1NACA63212.mat;
%     load S1NACA63212.mat;
    elseif logi_trip==0
       [P,S]=P0S0NACA63212;
%     load P0NACA63212.mat;
%     load S0NACA63212.mat;
    end
elseif Airfoil(k)==4
    if logi_trip==1
        [P,S]=P1S1NACA63215;
%     load P1NACA63215.mat;
%     load S1NACA63215.mat;
    elseif logi_trip==0
        [P,S]=P0S0NACA63215;
%     load P0NACA63215.mat;
%     load S0NACA63215.mat;
    end  
elseif Airfoil(k)==5
    if logi_trip==1
        [P,S]=P1S1NACA63218;
%     load P1NACA63218.mat;
%     load S1NACA63218.mat;
    elseif logi_trip==0
        [P,S]=P0S0NACA63218;
%     load P0NACA63218.mat;
%     load S0NACA63218.mat;
    end  
elseif Airfoil(k)==6
    if logi_trip==1
        [P,S]=P1S1NACA63221;
%     load P1NACA63221.mat;
%     load S1NACA63221.mat;
    elseif logi_trip==0
        [P,S]=P0S0NACA63221;
%     load P0NACA63221.mat;
%     load S0NACA63221.mat;
    end  
elseif Airfoil(k)==7
    if logi_trip==1
        [P,S]=P1S1NACA63412;
%     load P1NACA63412.mat;
%     load S1NACA63412.mat;
    elseif logi_trip==0
        [P,S]=P0S0NACA63412;
%     load P0NACA63412.mat;
%     load S0NACA63412.mat;
    end
elseif Airfoil(k)==8
    if logi_trip==1
        [P,S]=P1S1NACA63415;
%     load P1NACA63415.mat;
%     load S1NACA63415.mat;
    elseif logi_trip==0
        [P,S]=P0S0NACA63415;
%     load P0NACA63415.mat;
%     load S0NACA63415.mat;
    end
elseif Airfoil(k)==9
    if logi_trip==1
        [P,S]=P1S1NACA63418;
%     load P1NACA63418.mat;
%     load S1NACA63418.mat;
    elseif logi_trip==0
        [P,S]=P0S0NACA63418;
%     load P0NACA63418.mat;
%     load S0NACA63418.mat;
    end
elseif Airfoil(k)==10
    if logi_trip==1
        [P,S]=P1S1NACA63421;
%     load P1NACA63421.mat;
%     load S1NACA63421.mat;
    elseif logi_trip==0
        [P,S]=P0S0NACA63421;
%     load P0NACA63421.mat;
%     load S0NACA63421.mat;
    end
elseif Airfoil(k)==11
    if logi_trip==1
        [P,S]=P1S1FFA_W3_211;
%     load P1FFA-W3-211.mat;
%     load S1FFA-W3-211.mat;
    elseif logi_trip==0
        [P,S]=P0S0FFA_W3_211;
%     load P0FFA-W3-211.mat;
%     load S0FFA-W3-211.mat;
    end
elseif Airfoil(k)==12
    
    noisedata_dat=importdata(noisedatabase);
    no_RE=size(noisedata_dat,1)/2;
    P=noisedata_dat(1:no_RE,:);
    S=noisedata_dat((no_RE+1):end,:);
    
end
                   
for j=1:length(Aref)
   for i=1:length(Rref)-1
       if Rc<=Rref(i+1) & Rc>Rref(i)
          deltp(j)=linear(Rref(i),Rref(i+1),P(i,j),P(i+1,j),Rc);
          delts(j)=linear(Rref(i),Rref(i+1),S(i,j),S(i+1,j),Rc);
      end
   end
   
      if Rc<=min(Rref)
          deltp(j)=P(1,j);
          delts(j)=S(1,j);
        %  'Small Rc !!!'
      end
      
       if Rc>=max(Rref)
          deltp(j)=P(length(Rref),j);
          delts(j)=S(length(Rref),j);
         % 'Large Rc !!!!'
       end
       
end

  for j=1:length(Aref)-1
       if alpha<=Aref(j+1) & alpha>Aref(j)
          delt_star_p=linear(Aref(j),Aref(j+1),deltp(j),deltp(j+1),alpha);
          delt_star_s=linear(Aref(j),Aref(j+1),delts(j),delts(j+1),alpha);
       end
   end

    if alpha<=min(Aref)
          delt_star_p=deltp(1);
          delt_star_s=delts(1);
         % 'Small alpha !!!'
      end
      
       if alpha>=max(Aref)
          delt_star_p=deltp(length(Aref));
          delt_star_s=delts(length(Aref));
         % 'Large alpha !!!!'
       end
 D_star_p(k)=c*delt_star_p;
 D_star_s(k)=c*delt_star_s;
 D_p(k) = rp*D_star_p(k);
 D_s(k) = rs*D_star_s(k);
end
% 

% Perform the linear interpolation along the blade:
if r<=0.5*max(R)
    delt_star_p=linear(0,0.5*max(R),D_star_p(1),D_star_p(2),r);
    delt_star_s=linear(0,0.5*max(R),D_star_s(1),D_star_s(2),r);
    delt_p=linear(0,0.5*max(R),D_p(1),D_p(2),r);
    delt_s=linear(0,0.5*max(R),D_s(1),D_s(2),r);
elseif r>0.5*max(R)
    delt_star_p=linear(0.5*max(R),max(R),D_star_p(2),D_star_p(3),r);
    delt_star_s=linear(0.5*max(R),max(R),D_star_s(2),D_star_s(3),r);
    delt_p=linear(0.5*max(R),max(R),D_p(2),D_p(3),r);
    delt_s=linear(0.5*max(R),max(R),D_s(2),D_s(3),r);
end

function[y]=linear(a1,a2,b1,b2,x)
y = (b1-b2)/(a1-a2)*x + (a1*b2-a2*b1)/(a1-a2); 
              


             % ========================================================= %
             %    Ratio of BL thickness and BL displacement thickness    %
             % ========================================================= %
             
                          % << All Variables >> %
% alpha         -->  Angle of attack             
% c             -->  Chord length
% co            -->  Sound speed
% delt_0        -->  Boundary layer thickness of zero angle of attack
% delt_star_0   -->  Boundary layer displacement thickness of zero angle of attack
% delt_p        -->  Boundary layer thickness on pressure side of nonzero angle of attack
% delt_s        -->  Boundary layer thickness on suction side of nonzero angle of attack
% delt_star_p   -->  Boundary layer displacement thickness on pressure side of nonzero angle of attack
% delt_star_s   -->  Boundary layer displacement thickness on suction side of nonzero angle of attack
% logi_trip     -->  Condition for boundary layer tripping:
%                    tripped=0, untripped=1, partially tripped=2
% Ma            -->  Mach number
% Rc            -->  Reynolds number based on chord
% vis           -->  kinematic viscosity
% Vrel          -->  Relative velocity
%function[delt_p,delt_star_p,delt_s,delt_star_s]=Thickness(Vrel,c,vis,co,logi_trip,alpha)
%Rc=Vrel*c/vis;
%Ma=Vrel/co;

function[rp,rs]=Thickness_ratio(logi_trip,Rc,alpha)
c=1;
%logi_trip=0;
% kk(1)=0.1e5;
% Rc(1)=0.1e6;
% for ii=1:30
%     kk(ii+1)=1.2*kk(ii);
%     Rc(ii+1)=Rc(ii)+kk(ii+1);
% end
%alpha=[0:1:80];

for i=1:length(Rc)
    for j=1:length(alpha)
%% Compute boundary layer thickness of zero angle of attack: 
% (1) For the heavily tripped boundary layer
delt_0(i,j) = c * 10^(1.892-0.9045*log10(Rc(i))+0.0596*(log10(Rc(i)))^2);
% (2) For the untripped boundary layer
if logi_trip==1
    delt_0(i,j) = c * 10^(1.6569-0.9045*log10(Rc(i))+0.0596*(log10(Rc(i)))^2);
% (3) For the partially tripped boundary layer
elseif logi_trip==2
    delt_0(i,j) = 0.6 * c * 10^(1.892-0.9045*log10(Rc(i))+0.0596*(log10(Rc(i)))^2);
end

%% Compute boundary layer displacement thickness of zero angle of attack:
if logi_trip==0 
    if Rc(i)<=0.3e06
        delt_star_0(i,j) = c * 0.0601 * Rc(i)^(-0.114);
    elseif Rc(i)>0.3e06 
        delt_star_0(i,j) = c * 10^(3.411-1.5397*log10(Rc(i))+0.1059*(log10(Rc(i)))^2);
    end
elseif logi_trip==2
   if Rc(i)<=0.3e06
        delt_star_0(i,j) = c * 0.0601 * Rc(i)^(-0.114);
    elseif Rc(i)>0.3e06 
        delt_star_0(i,j) = c * 10^(3.411-1.5397*log10(Rc(i))+0.1059*(log10(Rc(i)))^2);
    end
    delt_star_0(i,j) = 0.6 * delt_star_0(i,j);
elseif logi_trip==1
    delt_star_0(i,j) = c * 10^(3.0187-1.5397*log10(Rc(i))+0.1059*(log10(Rc(i)))^2);
end
k0(i)=delt_0(i,j)/delt_star_0(i,j);

%% Compute boundary layer (displacement)thickness of nonzero angle of attack:
% Compute pressure side boundary layer (displacment)thickness:
delt_p(i,j) = delt_0(i,j) * 10^(-0.04175*alpha(j)+0.00106*alpha(j)^2);
delt_star_p(i,j) = delt_star_0(i,j) * 10^(-0.0432*alpha(j)+0.00113*alpha(j)^2);
% % % Compute suction side boundary layer (displacment)thickness:
if logi_trip==0  % tripped
    if alpha(j)>=0 & alpha(j)<=5
        delt_s(i,j) = delt_0(i,j) * 10^(0.0311*alpha(j));
        delt_star_s(i,j) = delt_star_0(i,j) * 10^(0.0679*alpha(j));
    elseif alpha(j)>5 & alpha(j)<=12.5
        delt_s(i,j) = delt_0(i,j) * 0.3468 *10^(0.1231*alpha(j));
        delt_star_s(i,j) = delt_star_0(i,j) * 0.381 * 10^(0.1516*alpha(j));
    else%if alpha>12.5 & alpha<=25
        delt_s(i,j) = delt_0(i,j) * 5.718 *10^(0.0258*alpha(j));
        delt_star_s(i,j) = delt_star_0(i,j) * 14.296 * 10^(0.0258*alpha(j));
    end
else
    if alpha(j)>=0 & alpha(j)<=7.5
        delt_s(i,j) = delt_0(i,j) * 10^(0.03114*alpha(j));
        delt_star_s(i,j) = delt_star_0(i,j) * 10^(0.0679*alpha(j));
    elseif alpha(j)>7.5 & alpha(j)<=12.5
        delt_s(i,j) = delt_0(i,j) * 0.0303 * 10^(0.2336*alpha(j));
        delt_star_s(i,j) = delt_star_0(i,j) * 0.0162 * 10^(0.3066*alpha(j));
    else%if alpha>12.5 & alpha<=25
        delt_s(i,j) = delt_0(i,j) * 12 *10^(0.0258*alpha(j));
        delt_star_s(i,j) = delt_star_0(i,j) * 52.42 * 10^(0.0258*alpha(j));
    end
end
rp(i,j)=delt_p(i,j)/delt_star_p(i,j);
rs(i,j)=delt_s(i,j)/delt_star_s(i,j);
end
end
%surf(rp)
%surf(rs)
%save('')


%========== END PROGRAM ============ %




