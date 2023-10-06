function plotpop(P,t,numgenes,data,fx_consBest)
% credit for this source code goes to:
% Matias Sessarego, M.Sc., msessare@ucalgary.ca, msessare@gmail.com
Pchromo=reshape([P.chr],numgenes,[])';
numobj_p_cons=length(data.weights);
fx_consBest(data.zind)=[];

G=zeros(size(Pchromo,1),numobj_p_cons);
for i=1:size(G,1)
    G(i,:)=P(i).objvals(~data.zind);
end

if data.optprogress==1
    if numobj_p_cons==3
        G1=1./G(:,1);
        G2=G(:,2);
        G3=G(:,3);
        
        h1=figure(30);
        plot3(G1,G2,G3,'r.','MarkerSize',5)
        xlabel('Obj.1','FontSize',12)
        ylabel('Obj.2','FontSize',12)
        zlabel('Obj.3','FontSize',12)
        view(120,30)
        %set(h1, 'units', 'centimeters', 'pos', [0 14 15 11])
        hold on
        
        if data.optAEP==0
            plot3(data.res(:,1),data.res(:,2),data.res(:,3),'b.','MarkerSize',5)
            text(data.res(:,1),data.res(:,2),data.res(:,3),'Utopia')
        end
        
        plot3((1/fx_consBest(1)),fx_consBest(2),fx_consBest(3),'pm')
        text((1/fx_consBest(1)),fx_consBest(2),fx_consBest(3),'Best')
        hold off
        
    elseif numobj_p_cons==2
        G1=1./G(:,1);
        G2=G(:,2);
        
        h1=figure(30);
        plot(G1,G2,'r.','MarkerSize',5)
        xlabel('Obj. 1','FontSize',12)
        ylabel('Obj. 2','FontSize',12)
        hold on
        
        if data.optAEP==0
            plot(data.res(:,1),data.res(:,2),'b.','MarkerSize',5)
            text(data.res(:,1),data.res(:,2),'Utopia')
        end
        
        plot((1/fx_consBest(1)),fx_consBest(2),'pm')
        text((1/fx_consBest(1)),fx_consBest(2),'Best')
        %set(h1, 'units', 'centimeters', 'pos', [0 14 15 11])
        hold off
        %axis([0 5 0 10])
    end
end

end