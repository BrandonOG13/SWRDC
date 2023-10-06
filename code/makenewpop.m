function [R] = makenewpop(P,pop,bnds,nc,nm,pc,pm,numgenes,wtdata,val)
%{
Credit goes to: Multi-objective NSGA-II code in C. by Kanpur Genetic Algorithms Laboratory
Revision 1.1.6 (08 July 2011) (for Linux only- 64-bit bug for binary coding fixed): NSGA-II in C with gnuplot (Real + Binary + Constraint Handling)
http://www.iitk.ac.in/kangal/codes.shtml
%}
%% binary tournament selection
matpool=pop/2;
ParS=bitournament(P,pop,matpool);

Par=reshape([ParS.chr],numgenes,[])';

%% crossover
EPS=1e-14;
childgenes=zeros(matpool,numgenes);
for i=1:2:(matpool-1)
	parAgenes=Par(i,:);
	parBgenes=Par((i+1),:);
    if (rand(1)<=pc)
        for j=1:numgenes
            switch val
                case 1
                    [c1,c2] = sbx(nc,bnds,parAgenes,parBgenes,j,EPS);
                case 2
                    [c1,c2] = unc(parAgenes,parBgenes,j,EPS);
            end
            childgenes(i,j)=c1;
            childgenes((i+1),j)=c2;
        end
    else
        childgenes(i,:)=parAgenes;
        childgenes((i+1),:)=parBgenes;
    end
end
%% parameter-based mutation
for i=1:matpool
    for j=1:numgenes
        if rand(1)<=pm
            x=childgenes(i,j);
            xl=bnds(2,j);
            xu=bnds(1,j);
            u=rand(1);
            if u<=0.5
                delta1=(1-((x-xl)/(xu-xl)))^(nm+1);
                deltabar=(((2*u)+((1-(2*u))*delta1))^(1/(nm+1)))-1;
            else
                delta2=(1-((xu-x)/(xu-xl)))^(nm+1);
                deltabar=1-(((2*(1-u))+ (2*(u-0.5)*delta2))^(1/(nm+1)));
            end
        	childgenes(i,j)=x+(deltabar*(xu-xl));
            childgenes(i,j)=max(childgenes(i,j),xl);
            childgenes(i,j)=min(childgenes(i,j),xu);
        end
    end
end
%% combine parent & offspring population
Q(1:matpool)=struct('chr',[],'o',[],'cons',[],'objvals',[]);
numobjcons=length(ParS(1).objvals)+1;

objcons=zeros(matpool,numobjcons);

clear e0
if wtdata.parallel==1;
    parfor e0=1:matpool
        objcons(e0,:)=objFcalc(childgenes(e0,:),wtdata);
    end
else
    for e0=1:matpool
        objcons(e0,:)=objFcalc(childgenes(e0,:),wtdata);
    end
end

for e1=1:matpool
    Q(e1).chr=childgenes(e1,:);
    Q(e1).cons=objcons(e1,end);
    Q(e1).objvals=objcons(e1,1:(end-1));
end
R=cat(2,ParS,Q);

objcons2=cat(2,reshape([R.objvals],numobjcons-1,[])',reshape([R.cons],1,[])');
objcons2(isinf(objcons2))=NaN;
zmax=max(objcons2(:,1:(end-1)));
zmin=min(objcons2(:,1:(end-1)));
objcons2(isnan(objcons2))=Inf;
o_c=zeros(pop,2);
for e2=1:pop
    o_c(e2,:)=objF(wtdata,objcons2(e2,:),zmax,zmin);
    R(e2).o=o_c(e2,1:(end-1));
end

end

function [c1,c2] = sbx(nc,bnds,parAgenes,parBgenes,j,EPS)
%simulated binary crossover (SBX)
if (rand(1)<=0.5) && (abs(parAgenes(j) - parBgenes(j)) > EPS)
    u=rand(1);
    if parAgenes(j)<=parBgenes(j)
        x1=parAgenes(j);
        x2=parBgenes(j);
    else
        x1=parBgenes(j);
        x2=parAgenes(j);
    end
    xl=bnds(2,j);
    xu=bnds(1,j);
    B=1+((2*(x1-xl))/(x2-x1));
    alpha=2-(B^(-(nc+1)));
    if u<=(1/alpha)
        Bbar=(alpha*u)^(1/(nc+1));
    else
        Bbar=(1/(2-(alpha*u)))^(1/(nc+1));
    end
    c1=0.5*((x1+x2)-(Bbar*abs(x2-x1)));
    B=1+((2*(xu-x2))/(x2-x1));
    alpha=2-(B^(-(nc+1)));
    if u<=(1/alpha)
        Bbar=(alpha*u)^(1/(nc+1));
    else
        Bbar=(1/(2-(alpha*u)))^(1/(nc+1));
    end
    c2=0.5*((x1+x2)+(Bbar*abs(x2-x1)));
    c1=max(c1,xl);
    c1=min(c1,xu);
    c2=max(c2,xl);
    c2=min(c2,xu);
else
    c1=parAgenes(j);
    c2=parBgenes(j);
end
end

function [c1,c2] = unc(parAgenes,parBgenes,j,EPS)
%uniform crossover
if (rand(1)<=0.5) && (abs(parAgenes(j) - parBgenes(j)) > EPS)
    c1=parBgenes(j);
    c2=parAgenes(j);
else
    c1=parAgenes(j);
    c2=parBgenes(j);
end
end