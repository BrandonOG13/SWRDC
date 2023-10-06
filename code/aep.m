%% AEP Calculation
function [AEP]=aep(cp,data,plotting)
%Load AEP data
probtype=data.probtype;
vcutin=data.vcutin;%m/s
vcutout=data.vcutout;%m/s
vo_rng=linspace(vcutin,vcutout,200);

% Calculate the Power Curve
if cp~=0
    pow_rng=0.5*data.rho*(vo_rng.^3)*pi*(data.R^2)*cp;
    
    %If power exceeds maximum
    if any(pow_rng>data.maxpower)
        vo_maxpower=interp1(pow_rng,vo_rng,data.maxpower);
        vo_rng=cat(2,vo_rng(vo_rng<vo_maxpower),vo_maxpower,vo_rng(vo_rng>vo_maxpower));
        pow_rng=0.5*data.rho*(vo_rng.^3)*pi*(data.R^2)*cp;
        pow_rng(pow_rng>data.maxpower)=0;
    end
end

% Calculate Probability Distribution
switch probtype
    case 3
        %Weibull Probability Distribution
        k=data.k_weibull;
        A=data.A_weibull;
        prob=(k/A)*((vo_rng/A).^(k-1)).*exp(-((vo_rng/A).^k));
    case 4
        %User-Defined Probability Distribution
        winddata=importdata(data.winddata);
        prob=interp1(winddata.data(:,1),winddata.data(:,2),vo_rng);
        if any(isnan(prob))
            disp('User-defined probability distribution does not contain cut-in and/or cut-out windspeed(s)')
            prob(isnan(prob))=0;
        end
    otherwise
        %Rayleigh Probability Distribution
        Vbar=data.mws;%m/s
        prob=(pi/2)*(vo_rng/(Vbar^2)).*exp(-(pi/4)*((vo_rng/Vbar).^2));
end

%Plot probabiltiy distribution
if plotting==1
    h18=figure(50);
    plot(vo_rng,prob,'-b.')
    title('Probability distribution - h(V_o)')
    xlabel('V_o: Windspeed (m/s)');   ylabel('h(V_o)');
end

%% Calculate the AEP
if cp~=0
    AEP=(trapz(vo_rng,pow_rng.*prob)*8760)/1000; %kWh/year
else
    AEP=0;
end

end