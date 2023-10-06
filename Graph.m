x = (3:7);
y = [1071.145337, 1087.2581838916, 1085.707473056, 1087.64921784555, 1095.69405643246];

plot(x, y, "-square", "LineWidth",2, "Color",[0.4660 0.6740 0.1880]);
xlabel("Cantidad de aspas");
ylabel("Producción anual de energía (kW/h)");
xlim([3 7]);
grid on;
box off;


legend("TSR=4", 'Location', 'eastoutside');



