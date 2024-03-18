tic
clear;close all;clc

%% Computational Paramters

t0 = 1;
years = 150;
tf = years*365*24*60*60;
tspan = linspace(t0,tf,years);
n = 40;

%% Initial Conditions

wl0 = 5000;
hm0 = (0.7167*1.4-0.483)/2;
wb0 = 10000;
cl0 = 1.2;
m0 = 0;
mflux0 = 0;

%% Model Parameters

Co = linspace(0,150/1000,n);
RSLR = linspace(0,15*(10^-3)/(3600*24*365),n);
rhos = 1000;
period = 12.5*3600*1;
ws = 0.5*10^-3;
tcr = 0.1;
wind = 6;
ka = 2;
ke = 0.16/(365*24*3600);
amp = 1.4/2;
rhom = 1000;
lamda = 0.0001;
beta = 10^100;
k = 0.012/(24*60*60);
wtidal = 20*(wb0-wl0);

%% Solving System of Ordinary Differential Equations

options = odeset('AbsTol',10^-6,'RelTol',10^-6,'Events', @behavior); 
options_2 = optimset('Algorithm','Levenberg-Marquardt','TolFun',10^-28,'TolX',10^-28,'MaxFunEvals',10000);
outcome = zeros(n,n);
tdata = NaN(years,7,4);
endpoint = NaN(7,n,n);

for i=1:n

    for j=1:n

        p = [Co(j) RSLR(i) rhos period ws tcr wind ka ke amp rhom lamda beta k wtidal];
        hl0 = fsolve(@(y) lagoon_depth(y,p,wl0),2,options_2);

        if hl0 < 0.5
            
            hl0 = 0.5;

        end
        
        y0 = [wl0 hl0 hm0 wb0 cl0 m0 mflux0];
        [t,y,ie] = ode23s(@(t,y) system_equations(t,y,p), tspan, y0, options);
        endpoint(:,i,j) = y(end,:);

        if i == 1 && j == 35

            tdata(1:length(t),:,1) = y; % fill

        elseif i == 15 && j == 35

            tdata(1:length(t),:,2) = y; % stable

        elseif i == 20 && j ==15

            tdata(1:length(t),:,3) = y; % retreat

        elseif i == 35 && j == 10

            tdata(1:length(t),:,4) = y; % drown

        end

        if size(ie) == size(outcome(i,j))

            outcome(i,j) = ie;

        end
        
    end

end

%% Organizing Data

tdata(:,6:7,:) = carbon_to_methane(tdata(:,6:7,:));
mflux = NaN(years,4);
mflux(1:end-1,1) = diff(tdata(:,7,1))./(365*24*60*60); % fill
mflux(1:end-1,2) = diff(tdata(:,7,2))./(365*24*60*60); % stable
mflux(1:end-1,3) = diff(tdata(:,7,3))./(365*24*60*60); % retreat
mflux(1:end-1,4) = diff(tdata(:,7,4))./(365*24*60*60); % drown
wl = NaN(n,n);
hl = NaN(n,n);
wm = NaN(n,n);
hm = NaN(n,n);
ml = NaN(n,n);
wl(:,:) = endpoint(1,:,:);
hl(:,:) = endpoint(2,:,:);
wm(:,:) = endpoint(4,:,:)-endpoint(1,:,:);
hm(:,:) = endpoint(3,:,:);
m(:,:) = endpoint(6,:,:);
limits = NaN(years,2);
limits(:,1) = (17.7+9.7)/(3600*10^6);
limits(:,2) = (17.7-9.7)/(3600*10^6);
optdepth = NaN(years,1);
optdepth(:,1) = hm0;
tspanyr = linspace(t0,years,years);

%% Plotting Results

figure(1)

tiledlayout(2,4);

nexttile([2,2]);
contourf(Co*1000,RSLR.*(3600*24*365)./10^-3,m)
title("Methane Emissions After 150 Years");
xlabel("Reference Sediment Concentration (mg/mL)");
ylabel("Rate of Sea Level Rise (mm/yr)");
colormap("turbo")
cb = colorbar();
cb.Label.String = "Kilograms";

nexttile([1,1]);
contourf(Co*1000,RSLR.*(3600*24*365)./10^-3,wl)
title("Width of the Lagoon After 150 Years");
colormap("turbo")
cb = colorbar();
cb.Label.String = "Meters";

nexttile([1,1]);
contourf(Co*1000,RSLR.*(3600*24*365)./10^-3,hl)
title("Depth of the Lagoon After 150 Years");
colormap("turbo")
cb = colorbar();
cb.Label.String = "Meters";

nexttile([1,1]);
contourf(Co*1000,RSLR.*(3600*24*365)./10^-3,wm)
title("Width of the Marsh Platform After 150 Years");
colormap("turbo")
cb = colorbar();
cb.Label.String = "Meters";

nexttile([1,1]);
contourf(Co*1000,RSLR.*(3600*24*365)./10^-3,hm)
title("Depth of the Marsh Platform Below MHW After 100");
colormap("turbo")
cb = colorbar();
cb.Label.String = "Meters";

figure(2)

tiledlayout(4,3);

nexttile([2,3]);
hold on
plot(tspanyr,tdata(:,6,1))
plot(tspanyr,tdata(:,6,2))
plot(tspanyr,tdata(:,6,3)) 
plot(tspanyr,tdata(:,6,4)) 
title("Methane Emissions");
ylabel("Kilograms of Methane");
xlabel("Years")
legend("Marsh Expansion to Full","Marsh Expansion to Equilibrium","Marsh Retreat (Lateral Erosion)","Marsh Vertical Drowning",'Location','northwest')
hold off

nexttile([2,1]);
hold on
plot(tspanyr,tdata(:,4,1)-tdata(:,1,1))
plot(tspanyr,tdata(:,4,2)-tdata(:,1,2))
plot(tspanyr,tdata(:,4,3)-tdata(:,1,3))
plot(tspanyr,tdata(:,4,4)-tdata(:,1,4))
title("Width of the Marsh");
ylabel("Meters");
hold off

nexttile([2,1]);
hold on
plot(tspanyr,tdata(:,3,1))
plot(tspanyr,tdata(:,3,2))
plot(tspanyr,tdata(:,3,3)) 
plot(tspanyr,tdata(:,3,4))
l = plot(tspanyr,optdepth,"--g");
legend(l,"Optimal Depth for Biomass Production","Location","southoutside")
title("Depth of the Marsh Below MHW");
ylabel("Meters");
hold off

nexttile([2,1]);
hold on
plot(tspanyr,mflux(:,1))
plot(tspanyr,mflux(:,2))
plot(tspanyr,mflux(:,3))
plot(tspanyr,mflux(:,4))
ll = plot(tspanyr,limits(:,1),"--k");
plot(tspanyr,limits(:,2),"--k")
legend(ll,"Comer-Warner et al. 2022","Location","southoutside")
title("Methane Flux");
ylabel("Kilograms of Methane per Meter^2 per Second");
hold off

toc