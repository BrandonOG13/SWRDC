function [data] = setup(data,varargin)
% credit for this source code goes to:
% Matias Sessarego, M.Sc., msessare@ucalgary.ca, msessare@gmail.com
if isempty(varargin)
%% BEM Code Setup
%Max Iterations
data.iter=250;

%Number of Blade Elements
data.belements=30;

%Relaxation Factor
data.relax=0.25;

%Blade Element Spacing
data.spacingBE=1;

%% Plotting
%blade geometry
data.pval=0;
%structural performance plot
data.pvalstruct=0;
%optimization progress (3 objectives only)
data.optprogress=0;
%final 3D blade plot
data.valBplot=0;
%print final blade coordinates
data.printcoords=0;

%Display Cp-TSR curve
data.tsrmax=10;
data.tsrmin=1;
data.cptsr=1;

%% Wind Data
%Air Density (kg/m^3) - Sea Level=1.225
data.rho=1.225;
%Air Viscosity
data.vis=1.78e-5;

%% AEP data
data.optAEP=0;
data.winddata='wind_data\Marrowstone_Island_C5.dat';
data.mws=8; %(m/s)
data.probtype=1;
data.vcutin=1; %(m/s)
data.vcutout=20; %(m/s)
data.k_weibull=1.9;
data.A_weibull=6.8; %(m/s)

%% Turbine Data
%Generator Rated Power (kW) - Rotor Rated Power = data.rat_pow/data.gagbe
data.maxpower=754;

%Generator and Gearbox Efficiency (/1.0)
%data.gagbe=0.944;

%Rated Rotor Speed
data.rat_rpm=550;

%Number of Blades
data.numBlades=3;

%Rotor and Hub Radii (m)
data.R=(3/2);
data.hubR=13.33/100;

%% Starting Data
%Wind speed (m/s) at starting
data.u_start=5;

%Tip speed ratio to complete starting
data.lambda_start=1;

%Moment of inertia from generator & drive train (kg.m^2)
data.J_gen=0.00;

% Resistive torque of generator and drive train (N.m)
data.tor_resis=0.0;

%% Variable Speed & Fixed Pitch Inputs
%Design wind speed and tip speed ratio at power extraction
data.vo_design=10.5; %10.5 implies IEC Class III
data.lambda_design=8;

%Pitch (degrees) & TSR Ranges
data.pitch=0;

%% Aerodynamic Data
data.airfoildatastring='airfoildatabase\sg6043_liftdrag_AR14_5.txt';

%Reynolds number range
data.re_array = [100000 150000 200000 300000 500000];

%% Structural Data
data.profiledatastring='airfoildatabase\sg6043_profile.txt';

%airfoil relative thickness (t/c) - c is chord length
data.thick=0.10;

%Material Elastic Moduli (Pa) & Densities (kg/m^3)
data.E=31.0e9;    %Uniaxial [0]
data.rh=550;

%Minimum Cap Thickness (relative to chord)
data.MinsparTh=0.0005;

%Safety Factors
% GL Safety Factors: 
% 1.35 (general material factor), 1.50 (influence of aging)
% 1.10 (temperature effect), 1.20 (hand lay-up laminate)
% 1.10 (non-post-cured laminate)
% TOTAL=2.94
data.materialSF=1.35*1.50*1.10*1.20*1.10;
%Safety Factor on Loads
data.loadSF=1.35;

%Load Case (Parked Condition: 2, Rotor Overspeed: 3)
data.loadcase=2;
data.rotoroverspeed=1.12; %in percentage (%)

%Load Case
%Designed for Class III, 
%parked/stationary conditions with 50-year extreme wind speed,
data.Ve50=1.4*(37.5); %m/s

%Maximum Allowable Strain (micro-strain) - axial (or tensile) direction 
%caused by normal and tangential bending moments
%7586 from Table 8 WindPACT study (7586/2.94=2580)
data.ep_failure=7586;
data.ep_max_tens=data.ep_failure/data.materialSF;

%Flapwise Lowest Eigenfrequency (1st)
%Lower (L*p) & Upper (U*p) bound where p 
%is rotor rotational frequency
data.L=3.5;
data.U=9999;

%% Noise Data
% Constants for noise model
data.zo = 10.00;        %  data.z0 = surface roughness
data.Ratio = 1.0;       %  data.Ratio = 1.0;
data.a_top = 0.2;       %  data.a_top = tower width at top of tower;
data.a_ground = 0.3;    %  data.a_ground = tower width at bottom of tower;
data.Htower = 15.0;     %  data.Htower = 15.0; Tower height
data.shaftlength = 1.0; %  data.shaftlength = 1.0; Turbine shaft length
data.Gamma = 0.15;      %  data.Gamma = 0.0; 
data.Tscale = 0.0;      %  data.Tscale = Turbulent length scale, set to 0 to calculate in program
data.Tinten = 0.0;      %  data.Tinten = Turbulent intensity scale, set to 0 to calculate in program
data.r0 = data.Htower;  %  data.r0 = horizontal distance from observer to turbine
data.round = 1;         %  data.round = 0; flag, 0 = rounded tip, 1 = square
data.hblunt = 0.001;
data.TEangle=14;
data.h0=0;
data.PSI=0;
data.add_noise=0;
data.estimate=0;

data.Rref=[0.1,0.4,0.8,1.2,2,2.5,3,3.5,4,5,6,7,8,10,15,20]*10^6;  % Reference Reynolds number
data.Aref=[-5,-4,-3,-2,-1,0,1,2,3,4,6,8,10,12,14,16,20,25,30];    % Reference alfa number

%% XFOIL Data
data.xfoilexecute=1;
data.ITER=200;
data.Ncrit=1;

%% Blade Optimization Data
%Degree of Bezier Curve for Post-Blade Root Region
data.N=5;

%Radial Location of Control Points 
%half-cosine spacing (towards hub)
%data.r_space = sort(cos((pi/2):(pi/(data.N*2)):pi)+1);
%half-cosine spacing (towards tip)
%data.r_space = sort(cos(0:(pi/(data.N*2)):(pi/2)));
%full-cosine spacing
%data.r_space = sort((cos(0:(pi/(data.N)):pi)+1)/2);
%uniform spacing
%data.r_space = (0:(1/data.N):1);
data.spacing=1;

%Gene Bounds
%twist upper and lower bounds:
data.uptwist=30;
data.dotwist=0;
%chord upper and lower bounds:
data.upchord=0.3;
data.dochord=0.045;

% data.upbnd=cat(2,[10 10 10 10 10],[0.35 0.35 0.35 0.35 0.35 0.35]);
% data.lobnd=cat(2,[0 0 0 0 0],[0.15 0.15 0.05 0.05 0.05 0.05]);

%% Baseline Blade
%Specify if a baseline wind turbine blade is provided
data.runBaseline=1;
data.turbdata='turbinedatabase\anderson.txt';

%% Determine Weights
data.weights=cat(2,0.25,0.25,0.25,0.25);

%% Genetic Algorithm Inputs
%Population Size (MUST BE EVEN)
data.pop=50;
%Number of Generations
data.mingen=100;

%Crossover Type
%SBX (val=1) | Uniform (val=2)
data.val=1;
%Simulated binary crossover (SBX) operator index (POSITIVE INTEGER ONLY)
data.nc=15;
%Crossover probability
data.pc=0.9;

%Parameter-based mutation operator index (POSITIVE INTEGER ONLY)
data.nm=20;
%Mutation probability
data.pm=1/18;

%Load initial population
data.flag3=0;

%% Execute Code using Parallel Computing Toolbox?
data.parallel=0;

%% Default Excel Results Output filename and directory
data.saveTextbox=[pwd '\output\results.xls'];

end

%% CALCULATE ADDITIONAL INPUTS
%----------------------------------------------------------------------
%----------------------------------------------------------------------
%----------------------------------------------------------------------
%----------------------------------------------------------------------

% Remove objectives and weights that user has deselected
zind=false(1,4);
w=data.weights;
if data.add_noise==0;
    zind(4)=1;
    w(4)=NaN;
end
if data.perform_structure==0;
    zind(3)=1;
    w(3)=NaN;
end
if data.add_starting==0;
    zind(2)=1;
    w(2)=NaN;
end
w(isnan(w))=[];
data.weights=w;
data.zind=zind;

%Maximum strain
data.ep_max_tens=data.ep_failure/data.materialSF;

%Control point data
data.upbnd=cat(2,ones(1,data.N)*data.uptwist,ones(1,data.N+1)*data.upchord);
data.lobnd=cat(2,ones(1,data.N)*data.dotwist,ones(1,data.N+1)*data.dochord);

%radial location of control points 
spacing=data.spacing;
switch spacing
    case 2        %uniform spacing
        data.r_space = (0:(1/data.N):1);
    case 3        %half-cosine spacing (finer towards hub)
        data.r_space = sort(cos((pi/2):(pi/(data.N*2)):pi)+1);
    case 4        %half-cosine spacing (finer towards tip)
        data.r_space = sort(cos(0:(pi/(data.N*2)):(pi/2)));
    case 5        %full-cosine spacing
        data.r_space = sort((cos(0:(pi/(data.N)):pi)+1)/2);
    otherwise     %uniform spacing
        data.r_space = (0:(1/data.N):1);
end

%Blade Element r/R Distribution (dimensionless)
spacingBE=data.spacingBE;
hubR_rotR=data.hubR/data.R;
switch spacingBE
    case 3        %half-cosine spacing (finer towards hub)
        space=sort(cos((pi/2):(pi/((data.belements-1)*2)):pi)+1)';
    case 4        %full-cosine spacing
        space=sort((cos(0:(pi/(data.belements-1)):pi)+1)/2)';
    case 5        %uniform spacing
        space =(0:(1/(data.belements-1)):1)';
    otherwise     %half-cosine spacing (finer towards tip)
        space=sort(cos(0:(pi/((data.belements-1)*2)):(pi/2)))';
end
data.nu=(space*(1-hubR_rotR))+hubR_rotR;

%Blade tip twist constraint (degrees)
%(input for missing control point at N+1 for twist)
data.tw_tip = 0.0;

%Airfoil Data
[data] = airfoildata(data);
%--------------------------------------------

%Pitch
data.pitch=(data.pitch*(pi/180));

%Starting
data.t_max=200;

%Noise Data
%surface roughness
zo=data.zo;
switch zo
    case 1;        z0=0.01;     %Very smooth, ice or mud
    case 2;        z0=0.20;     %Calm open sea
    case 3;        z0=0.50;     %Blown sea
    case 4;        z0=3.00;     %Snow surface
    case 5;        z0=8.00;     %Lawn grass
    case 6;        z0=10.00;    %Rough pasture
    case 7;        z0=30.00;    %Fallow field
    case 8;        z0=50.00;    %Crops
    case 9;        z0=100.00;   %Few trees
    case 10;       z0=250.00;   %Many trees, hedges
    case 11;       z0=500.00;   %Forest and woodlands
    case 12;       z0=1500.00;  %Suburbs
    case 13;       z0=3000.00;  %Centers of cities with tall buildings
end
data.z0=z0/1000; %convert mm to m
%power law factor
if data.estimate==1
    zotemp=data.z0;
    data.Gamma=0.24+0.096*log10(zotemp/1000)+0.016*(log10(zotemp/1000))^2;
end

%XFOIL Data
xfoilexecute=data.xfoilexecute;
Ncrit=9;
data.logi_trip=1;
%logi_trip --> Condition for boundary layer tripping: tripped=0, untripped=1, partially tripped=2
if data.Ncrit==3
    Ncrit=4;
    data.logi_trip=0;
end
itersXFOIL=data.ITER;
Rref=data.Rref;  % Reference Reynolds number   
Aref=data.Aref;  % Reference alfa number
if xfoilexecute==1
    noise_FileName='TEMPnoisedatabase.dat';
    saveNoiseDataName=[pwd '\noisedatabase\' noise_FileName];
    data.noisedatabase=saveNoiseDataName;
    [~]=xfoil_execution(Rref,Aref,itersXFOIL,Ncrit,saveNoiseDataName);
else
    delete('airfoilXFOIL.txt')
end

%Utopia Point
data.res=[(16/27) 0 0];

%% Prepare Excel file and output
%Copy the template for the Excel output file to the new output file
try
    copyfile([pwd '\output\TEMPLATE_Do_Not_Delete.xls'],data.saveTextbox);
catch
    disp('ABORT: Excel is operating in the background. Please close it through the Task Manager.')
    disp('Or an incorrect path was specified for the results file. Returning to GUI...')
    data=NaN;
    return
end
%% Display Disclaimer
fprintf('\n')
fprintf('Small Wind-turbine Rotor Design Code (SWRDC) \nby M. Sessarego \n')
fprintf('Department of Mechanical and Manufacturing Engineering, University of Calgary, \n')
fprintf('2500 University Dr. NW, Calgary, Alberta, Canada, T2N 1N4. \n')
fprintf('E-mail: msessare@ucalgary.ca \n')
fprintf('Copyright (c) 2013, Matias Sessarego \nAll rights reserved. \n')
fprintf('BY INSTALLING OR USING THIS SOFTWARE, YOU AGREE TO HAVE READ THE LICENSE TERMS AS DESCRIBED IN THE README_FIRST_LICENSE.txt file \n')

%% Baseline Turbine
if data.runBaseline==1
    [~] = runbaseline(data);
end

end