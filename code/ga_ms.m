function [P,I2]=ga_ms(data)
% credit for this source code goes to:
% Matias Sessarego, M.Sc., msessare@ucalgary.ca, msessare@gmail.com
%Genetic Algorithm
format compact
format shortG
%% INPUTS
%Load Wind Turbine Optimization Data
wtdata=data;

%Population Size (MUST BE EVEN)
pop=wtdata.pop;
%Number of Generations
mingen=wtdata.mingen;
%Generation Start
t0=1;

%Load Gene Bounds
upbnd=wtdata.upbnd;
lobnd=wtdata.lobnd;

%Crossover Type
%SBX (val=1) | Uniform (val=2)
val=wtdata.val;
%Simulated binary crossover (SBX) operator index (POSITIVE INTEGER ONLY)
nc=wtdata.nc;
%Crossover probability
pc=wtdata.pc;

%Parameter-based mutation operator index (POSITIVE INTEGER ONLY)
nm=wtdata.nm;
%Mutation probability
pm=wtdata.pm;

%Load initial population (1), load chros (2), set random (0)
flag3=data.flag3;

%% ALGORITHM
bnds=cat(1,upbnd,lobnd);
numgenes=size(bnds,2);

%Initial Population
[P,chros]=inipop(pop,bnds,wtdata,flag3);

if size(chros,1)~=pop
    %Size of loaded initial population does not match
    %population size specified in GUI.
    %Optimization aborted.'
    I2=NaN;
    return
end
        
f=[P.o]';

%Construct a wait bar
h = waitbar(0,'Generation = 0','Name','Optimizing...','CreateCancelBtn','setappdata(gcbf,''canceling'',1)');
setappdata(h,'canceling',0)

%Main Loop
Ptot=[];
disp(' ')
disp(cat(2,'Generation = ',num2str(t0-1)))
disp(' ')
for t=t0:1:mingen
    %Check if cancel button is pressed
    if getappdata(h,'canceling')
        disp('OPTIMIZATION CANCELLED')
        disp(' ')
        break
    end
    
    disp('Cp / AEP  Starting  Mass (kg)  C.1: Strain  C.2: Freq.L  C.3: Freq.U  Noise (dB)  1stFlapEigenFreq (rad/s)')
    
    Ptot=cat(2,Ptot,f);
    [R]=makenewpop(P,pop,bnds,nc,nm,pc,pm,numgenes,wtdata,val);
    P=R;
    f=[P.o]';
    
    save('endPop.mat','P')
    save('PopTot.mat','Ptot')
    
    %Check if population converged to optimum solution
    if all(isnan(f))
        I2=1;
        disp(' ')
        disp('ABORT: All individuals in population are identical')
        disp(' ')
        break
    end
    
    %Determine best individual
    I0=[P.o]~=Inf;
    I1=([P.cons]==min([P(I0).cons]));
    I2=find([P.o]==min([P(I1).o]));
    I2=I2(end);
    bestObjValues=P(I2).objvals;
    
    %Display result on command window
    disp(' ')
    disp('Best Value  Best Constraint')
    fprintf('%5.6f   % 2.6f\n',cat(2,P(I2).o,P(I2).cons))
    disp(' ')
    disp(cat(2,'Generation = ',num2str(t)))
    disp(' ')
    
    %Plot results
    plotpop(P,t,numgenes,wtdata,bestObjValues)
    
    %Update on optimization progress
    waitbar(t/mingen,h,sprintf('Generation %i of %i',t,mingen))
end
waitbar(t/mingen,h,'Please wait...')

%% EXPORT RESULTS
disp('Cp / AEP  Starting  Mass (kg)  C.1: Strain  C.2: Freq.L  C.3: Freq.U  Noise (dB)  1stFlapEigenFreq (rad/s)')
[~,~,~,~] = finalresults(P(I2).chr,data);
disp(' ')

save('endPop.mat','P')
Ptot=cat(2,Ptot,f);
save('PopuGen.mat','Ptot')

Best=P(I2);
save('finalblade.mat','Best')
delete(h)
disp('FINISHED')
disp(' ')
