function [Pmatpool]=bitournament(P,pop,matpool)
% credit for this source code goes to:
% Matias Sessarego, M.Sc., msessare@ucalgary.ca, msessare@gmail.com
%% binary tournament selection
Pmatpool(1:matpool)=struct('chr',[],'o',[],'cons',[],'objvals',[]);
for i=1:matpool
    Pt=P;
    popt=pop;
    indv1=randi(popt);
    best=Pt(indv1);
    
    Pt(indv1)=[];
    popt=popt-1;
    indv2=randi(popt);
    next=Pt(indv2);
    
    if next.cons==0 && best.cons>0
        best=next;
    elseif next.cons==0 && best.cons==0
        if next.o<best.o
            best=next;
        elseif (next.o==best.o) && (rand()<=0.5)
            best=next;
        end
    end
    Pmatpool(i)=best;
end