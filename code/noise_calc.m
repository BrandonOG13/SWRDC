function [noise] = noise_calc(chords,twists,nu,data,rpm_design)
%Noise Calculation

%Load noise database file
noisedatabase=data.noisedatabase;

%Single variable
Airfoil = [12,12,12];
Vo = data.vo_design;
Density = data.rho;
vis = data.vis;
Tscale = data.Tscale;
Tinten = data.Tinten;
ro = data.r0;
ho = data.h0;
PSI = data.PSI;
co = 340;
Ratio = data.Ratio;
theta_yaw  = 0;
theta_tilt  = 0;
theta_cone = 0;
wing = 0;
a_top = data.a_top;
a_ground = data.a_ground;
Htower = data.Htower;
shaftlength = data.shaftlength;
omega =rpm_design*(pi/30);
Gamma = data.Gamma;
Round = data.round;%Rounded tip=0, squared tip=1
h=data.hblunt;
TEangle=data.TEangle;
logi_trip = data.logi_trip; %tripped=0 (Ncrit=4), untripped=1 (Ncrit=9), partially tripped=2
R = nu'*data.R;

Rref=data.Rref;  % Reference Reynolds number   
Aref=data.Aref;  % Reference alfa number

%Vectors
Chord =  chords';
Twist =  twists'*(180/pi);
Hblunt =  ones(size(Chord))*h;
TEANGLE = ones(size(Chord))*TEangle;
zo = data.z0;

%Main Noise Function
[SPLw,SPLw1,SPLw2,SPLw3,SPLw4,SPLw5,SPLw6,SPLw7,SPLw8,...
    SPLwA,SPLwA1,SPLwA2,SPLwA3,SPLwA4,SPLwA5,SPLwA6,SPLwA7,SPLwA8,...
    SPL,SPL1,SPL2,SPL3,SPL4,SPL5,SPL6,SPL7,SPL8,...
    scale,inten,ro,ho,PSI,Freq,A_weighting,SPLtot,SPLtotw,wSPLtot,wSPLtotw]=MAIN(Airfoil,Vo,Density,vis,...
    Tscale,Tinten,ro,ho,PSI,co,Ratio,theta_yaw,theta_tilt,theta_cone,wing,...
    a_top,a_ground,Htower,shaftlength,omega,Gamma,Round,logi_trip,R,Chord,Twist,Hblunt,TEANGLE,zo,...
    noisedatabase,Rref,Aref);

%Total sound pressure level including the human response of standard
%A-weighted filter in one-third octive bands
noise = SPLtotw;
end