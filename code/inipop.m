function [P,chros] = inipop(pop,bnds,wtdata,flag3)
% credit for this source code goes to:
% Matias Sessarego, M.Sc., msessare@ucalgary.ca, msessare@gmail.com
% initialize population
numgenes=size(bnds,2);
upbnd=bnds(1,:);
lobnd=bnds(2,:);
dibnd=upbnd-lobnd;

P(1:pop)=struct('chr',[],'o',[],'cons',[],'objvals',[]);

disp(' ')
fprintf('Loading and calculating initial population...\n')
disp(' ')
disp('Cp / AEP  Starting  Mass (kg)  C.1: Strain  C.2: Freq.L  C.3: Freq.U  Noise (dB)  1stFlapEigenFreq (rad/s)')

switch flag3
    case 0
        %Uniform random initial population
        chros=repmat(lobnd,pop,1)+(repmat(dibnd,pop,1).*rand(pop,numgenes));
    case 1
        %Load initial population
        load endPop.mat
        chros=reshape([P.chr],numgenes,[])';
        %return
        if size(chros,1)~=pop
            disp('Size of loaded initial population does not match')
            disp('population size specified in GUI.')
            disp('Optimization aborted.')
            return
        end
    case 2
        %load chros
        load chros.mat
        %chros=H;
        if size(chros,1)~=pop
            disp('Size of loaded initial population does not match')
            disp('population size specified in the GUI.')
            disp('Optimization aborted.')
            return
        end
end

numobjcons=5;
objcons=zeros(pop,numobjcons);

clear e0
if wtdata.parallel==1;
    parfor e0=1:pop
        objcons(e0,:)=objFcalc(chros(e0,:),wtdata);
    end
else
    for e0=1:pop
        objcons(e0,:)=objFcalc(chros(e0,:),wtdata);
    end
end

objcons(isinf(objcons))=NaN;
zmax=max(objcons(:,1:(end-1)));
zmin=min(objcons(:,1:(end-1)));
objcons(isnan(objcons))=Inf;

o_c=zeros(pop,2);
for e1=1:pop
    o_c(e1,:)=objF(wtdata,objcons(e1,:),zmax,zmin);
end

for e2=1:pop
    P(e2).chr=chros(e2,:);
    P(e2).o=o_c(e2,1:(end-1));
    P(e2).cons=o_c(e2,end);
    P(e2).objvals=objcons(e2,1:(end-1));
end
disp(' ')
disp('Initial population calculation complete.')