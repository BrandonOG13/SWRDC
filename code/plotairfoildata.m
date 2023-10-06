function [h1] = plotairfoildata(data)
%This function will collect the airfoil data and display 
%it in an organized figure
[dataTEMP] = airfoildata(data);
yu=dataTEMP.yu;
yl=dataTEMP.yl;
xu=dataTEMP.xu;
xl=dataTEMP.xl;
re_array=data.re_array;
delete('airfoilXFOIL.txt')

h1=figure(50);
subplot(2,2,1);
plot(xu,yu,'.-b',xl,yl,'.-r')
title('Profile Coordinates')
xlabel('\it{x/c}');   ylabel('\it{y/c}');
set(gca,'XTick',0:0.25:1)
axis equal

afdata=dataTEMP.afdata;
d=size(afdata,3);

% Define your own Color Order Matrix
% For this example I used the default.
% It should be a 3-column matrix with entries in [0,1]
ColOrd = get(gca,'ColorOrder');

% Determine the number of colors in
% the matrix
[m,n] = size(ColOrd);

subplot(2,2,3);
hold on
for i=1:1:d
    % Determine which row to use in the
    % Color Order Matrix
    ColRow = rem(i,m);
    if ColRow == 0
        ColRow = m;
    end
    % Get the color
    Col = ColOrd(ColRow,:);
    
    plot(afdata(:,1,i),afdata(:,2,i),'Color',Col)
end
axis([-180 180 -2 2])
axis 'auto y'
title('Lift Data')
ylabel('\it{C_l}');   xlabel('Angle of Attack ({\circ})');
set(gca,'XTick',-180:90:180)

hold off

subplot(2,2,4);
hold on
for i=1:1:d
    % Determine which row to use in the
    % Color Order Matrix
    ColRow = rem(i,m);
    if ColRow == 0
        ColRow = m;
    end
    % Get the color
    Col = ColOrd(ColRow,:);
    
    plot(afdata(:,1,i),afdata(:,3,i),'Color',Col)
end
axis([-180 180 -2 2])
axis 'auto y'
title('Drag Data')
ylabel('\it{C_d}');   xlabel('Angle of Attack ({\circ})');
set(gca,'XTick',-180:90:180)

hold off

subplot(2,2,2);
legendtext =cell(d,1);
hold on
for i=1:1:d
    % Determine which row to use in the
    % Color Order Matrix
    ColRow = rem(i,m);
    if ColRow == 0
        ColRow = m;
    end
    % Get the color
    Col = ColOrd(ColRow,:);
    legendtext{i} = num2str(re_array(i));
    
    plot([-0.5 -0.5],[-0.5 -0.5],'Color',Col)
end
axis([-1 1 -1 1])
set(gca,'visible','off')
legendoptions1=legend(legendtext);
v = get(legendoptions1,'title');
set(v,'string','Reynolds Number');
set(legendoptions1,'Location','SouthWest');

hold off

set(h1, 'units', 'centimeters', 'pos', [2 2 20 16])

end