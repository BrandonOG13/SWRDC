function [inner_coords] = trim(th,xu,xl,yu,yl,angle,tdx_sign,tdy_sign,limit)
% credit for this source code goes to:
% Matias Sessarego
thdx=(th*cos(angle)).*tdx_sign;
thdy=(th*sin(angle)).*tdy_sign;
in_x1=cat(2,xu,xl)+thdx;
in_y1=cat(2,yu,yl)+thdy;

in_x2=cat(1,in_x1(2:end,:),in_x1(1,:));
in_y2=cat(1,in_y1(2:end,:),in_y1(1,:));
%[B,index]=sortrows(A,...)
x1=in_x1(:,1);
y1=in_y1(:,1);
x3=in_x1(:,2);
y3=in_y1(:,2);

x2=in_x2(:,1);
y2=in_y2(:,1);
x4=in_x2(:,2);
y4=in_y2(:,2);

x2_x1=x2-x1;
x4_x3=x4-x3;

denom=(((y4-y3).*x2_x1)-(x4_x3.*(y2-y1)));

num_a=((x4_x3.*(y1-y3))-((y4-y3).*(x1-x3)));
num_b=((x2_x1.*(y1-y3))-((y2-y1).*(x1-x3)));

ua=num_a./denom;
ub=num_b./denom;

xia=x1+(ua.*x2_x1);
xib=x3+(ub.*x4_x3);

gfa=(xia.*((xia>=th)&(xia<=(1-th))&(abs(ua)>=0)&(abs(ua)<=1)));
gfa(isnan(gfa))=0;
gfb=(xib.*((xib>=th)&(xib<=(1-th))&(abs(ub)>=0)&(abs(ub)<=1)));
gfb(isnan(gfb))=0;
a1=find(gfa,1,'last');
a2=find(gfa,1,'first');
b1=find(gfb,1,'last');
b2=find(gfb,1,'first');

if isempty(a1);    a1=[];    end
if isempty(a2);    a2=[];    end
if isempty(b1);    b1=[];    end
if isempty(b2);    b2=[];    end

% xa=max(cat(2,a2,b2));
% xb=min(cat(2,a1,b1));
xa=min(cat(2,a2,b2));
xb=max(cat(2,a1,b1));

xif=cat(1,x1(xa)+(ua(xa).*x2_x1(xa)),x3(xb)+(ub(xb).*x4_x3(xb)));

xis=min(xif);
xie=max(xif);
    
spacing=sort((cos(0:(pi/(limit-1)):pi)+1)/2)*(xie-xis);
xcos=(spacing+xis)';

in_yu_fin=interp1(in_x1(:,1),in_y1(:,1),xcos);
in_yl_fin=interp1(in_x1(:,2),in_y1(:,2),xcos);

inner_coords=cat(2,xcos,in_yu_fin,in_yl_fin);