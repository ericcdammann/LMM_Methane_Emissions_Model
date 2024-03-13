tic 
clc; clear variables;

%% Computational Parameters

years = 100;
n_iter = 10000;
scaled_n_iter = n_iter/250;
t_init = 0;
t_fin = years*365*24*60*60;
h = (t_fin-t_init)/n_iter;
t_sec = linspace(t_init,t_fin,n_iter);
t = t_sec/(60*60*24*365);

%% Initial Conditions

wm_init = 5000;
wl_init = 5000;
hm_init = 0.25;
hl_init = 3.25;
wb_init = wl_init+wm_init;
ctot_init = 0;
clab_init = 0.45;
cref_init = 0.45;
cdec_init = 0;

%% Model Parameters

rhos = 1000;
P = 12.5*3600*1;
ws = 0.5*10^-3;
tcr = 0.1;
Co = 0.075;
wind = linspace(0,20,scaled_n_iter);
ka = 2;
ke = 0.16/(365*24*3600);
amp = 1.4/2;
rslr = linspace(0,15,scaled_n_iter);
RSLR = rslr*(10^-3)/(3600*24*365);
rhom = 1000;
lamda = 0.0001;
beta = 10^10;
kk = 0.012/(24*60*60);
wtidal = wm_init*7;

%% Preallocationg Arrays

WL = zeros(n_iter,scaled_n_iter,scaled_n_iter);
HM = zeros(n_iter,scaled_n_iter,scaled_n_iter);
HL = zeros(n_iter,scaled_n_iter,scaled_n_iter);
WB = zeros(n_iter,scaled_n_iter,scaled_n_iter);
CTOT = zeros(n_iter,scaled_n_iter,scaled_n_iter);
CLAB = zeros(n_iter,scaled_n_iter,scaled_n_iter);
CREF = zeros(n_iter,scaled_n_iter,scaled_n_iter);
MFLUX = zeros(n_iter,scaled_n_iter,scaled_n_iter);
MTOT = zeros(n_iter,scaled_n_iter,scaled_n_iter);
WL(1,:,:) = wl_init;
HM(1,:,:) = hm_init;
HL(1,:,:) = hl_init;
WB(1,:,:) = wb_init;
CTOT(1,:,:) = ctot_init;
CLAB(1,:,:) = clab_init;
CREF(1,:,:) = cref_init;

%% Behavior

behavior = ones(scaled_n_iter,scaled_n_iter);
stable = 1;
fill = 2;
drown = 3;
erode = 4;

%% Solving Equations

for j=1:scaled_n_iter % wind

    for k=1:scaled_n_iter  % RSLR

        for l=1:n_iter-1

            fetch = WL(l,k,j);
            hl = HL(l,k,j);
            hm = HM(l,k,j);
            wb = WB(l,k,j);
            wm = wb-fetch;
            ctot = CTOT(l,k,j);
            clab = CLAB(l,k,j);
            cref = CREF(l,k,j);
            mtot = MTOT(l,k,j);

            % Organic Sediment Production

            BMax = 0.900;
            Dmin = 0;
            Dmax = 0.7167*2*amp-0.483;
            AA = 0.25*(-Dmin-Dmax)*(Dmax-3*Dmin);
            Bpeak = BMax*(hm-Dmax)*(hm-Dmin)/AA;

            if Bpeak <= 1e-3

                Bpeak = 0;

            end

            Bfrac = (Bpeak/BMax);
            nuGp = 0.0138;
            AMC = (365/2)*Bpeak*(nuGp)/(365*24*60*60);
            por = 1000/2650;
            chiref = 0.15;
            Rref = AMC*chiref;
            po = 1000;
            O = (1/por)*(Rref/po);
       
            % Average Depths

            fac = min(1,hl/(2*amp));
            fac2 = min(1,hm/(2*amp));
            Hf = (hl+(hl-fac*2*amp))/2;
            Hm = (hm+(hm-fac2*2*amp))/2;

            % Sediment Concentration in the Lagoon

            tau_1 =  lagoon_bed_shear_stress(fetch,wind(j),Hf);
            tau_lagoon = max((tau_1-tcr)/tcr,0)*lamda;
            Cr = rhos*tau_lagoon/(1+tau_lagoon);

            % Sediment Concentration in the Marsh

            if Hm > 1e-4
    
                tau_2 =  marsh_bed_shear_stress(fetch,wind(j),Hm,Bfrac);
    
            else 
        
                tau_2 = 0;
    
            end
    
            tau_marsh = max((tau_2-tcr)/tcr,0)*lamda;
            Cm = rhos*tau_marsh/(1+tau_marsh);

            % Sediment Flux Between the Lagoon and Marsh

            Fm = (Cr-Cm)*min(2*amp,hm)/P/rhom;

            % Sediment Flux Between the Open Ocean and Lagoon

            Fc = (Cr-Co)*(fac*2*amp)/P/rhom;

            % Marsh Platform Erosion rate

            dist = 10;
            hs = hm+(hl-hm)*(1-exp(-dist*0.1/hl));
            W = wave_power(amp,wind(j),fetch,hs);
            Be = ke*W/(hs-hm);

            % Marsh Platform Progradation Rate

            Ba = ka*Cr*ws/rhom;

            % Marsh Platform Upland Migration Rate

            Umig = (Fm+O)/beta;

            % Labile Carbon Deposition Rate and Biogeochemicial Model Equations
        
            chilab = 1-chiref;
            Rlab  = AMC*chilab;
            phi = wm/(2*wtidal);

            % Carbon Dynamics Equations

            dCTOTdt = AMC-(kk*clab);
            dCLABdt = Rlab-(kk*clab);
            dCREFdt = Rref;
            mflux = phi*kk*clab;

            % Morphodynamic Equations

            dWLdt = Be-Ba;
            dHMdt = -Fm-O+RSLR(k);
            dHLdt = -(Be-Ba)*(hl-hm)/fetch+Fm*wm/fetch+Fc+RSLR(k);
            dWBdt = Umig;

            % Foward Euler Method

            if hm+h*dHMdt >= hl+h*dHLdt || hm+h*dHMdt >= 2*amp
                
                behavior(j,k) = drown;

                break

            elseif (wb+h*dWBdt)-(fetch+h*dWLdt) <=1 && hm+h*dHMdt <= 2*amp

                behavior(j,k) = erode;

                break

            elseif (wb+h*dWBdt) <= (wb+h*dWBdt)-(fetch+h*dWLdt)

                behavior(j,k) = fill;

                break

            else

                WL(l+1,k,j) = fetch+h*dWLdt;
                HM(l+1,k,j) = hm+h*dHMdt;
                HL(l+1,k,j) = hl+h*dHLdt;
                WB(l+1,k,j) = wb+h*dWBdt;
                CTOT(l+1,k,j) = ctot+h*dCTOTdt;
                CLAB(l+1,k,j) = clab+h*dCLABdt;
                CREF(l+1,k,j) = cref+h*dCREFdt*wm;
                MTOT(l+1,k,j) = mtot+h*mflux*wm;
                MFLUX(l+1,k,j) = mflux;

            end

        end

    end

end

WM = WB - WL;
WM( WM <= 0 ) = NaN;
WL( WL <= 0 ) = NaN;
HM( HM <= 0 ) = NaN;
HL( HL <= 0 ) = NaN;
CTOT( CTOT <= 0 ) = NaN;
CREF( CREF <= 0 ) = NaN;
CLAB( CLAB <= 0 ) = NaN;
MTOT( MTOT <= 0 ) = NaN;
MFLUX( MFLUX <= 0 ) = NaN;

%% Plotting 

a_data = zeros(scaled_n_iter,scaled_n_iter);
b_data = zeros(scaled_n_iter,scaled_n_iter);
c_data = zeros(scaled_n_iter,scaled_n_iter);
d_data = zeros(scaled_n_iter,scaled_n_iter);
e_data = zeros(scaled_n_iter,scaled_n_iter);
f_data = zeros(scaled_n_iter,scaled_n_iter);

for m=1:scaled_n_iter

    a_data(:,m) = WM(n_iter,m,:);
    b_data(:,m) = WL(n_iter,m,:);
    c_data(:,m) = HM(n_iter,m,:);
    d_data(:,m) = HL(n_iter,m,:);
    e_data(:,m) = CREF(n_iter,m,:);
    f_data(:,m) = MTOT(n_iter,m,:);

end

figure(1)

contourf(rslr,wind,a_data);
title("Width of the Marsh Platform After 100 Years Under Different Environmental Parameters");
xlabel("Rate of Sea Level Rise (mm/yr)");
ylabel("Average Wind Speed (m/s)");
colormap("turbo")
cb = colorbar();
cb.Label.String = "Meters";

figure(2)

contourf(rslr,wind,b_data);
title("Width of the Lagoon After 100 Years Under Different Environmental Parameters");
xlabel("Rate of Sea Level Rise (mm/yr)");
ylabel("Average Wind Speed (m/s)");
colormap("turbo")
cb = colorbar();
cb.Label.String = "Meters";

figure(3)

contourf(rslr,wind,c_data);
title("Depth of the Marsh Below MHW After 100 Years Under Different Environmental Parameters");
xlabel("Rate of Sea Level Rise (mm/yr)");
ylabel("Average Wind Speed (m/s)");
colormap("turbo")
cb = colorbar();
cb.Label.String = "Meters";

figure(4)

contourf(rslr,wind,d_data);
title("Depth of the Lagoon After 100 Years Under Different Environmental Parameters");
xlabel("Rate of Sea Level Rise (mm/yr)");
ylabel("Average Wind Speed (m/s)");
colormap("turbo")
cb = colorbar();
cb.Label.String = "Meters";

figure(5)

contourf(rslr,wind,e_data);
title("Carbon Stored per m^2 after 100 Years Under Different Environmental Parameters");
xlabel("Rate of Sea Level Rise (mm/yr)");
ylabel("Average Wind Speed (m/s)");
colormap("turbo")
cb = colorbar();
cb.Label.String = "Kilograms of Carbon per Meter^2";

figure(6)

contourf(rslr,wind,f_data);
title("Methane Emissions after 100 Years Under Different Environmental Parameters");
xlabel("Rate of Sea Level Rise (mm/yr)");
ylabel("Average Wind Speed (m/s)");
colormap("turbo")
cb = colorbar();
cb.Label.String = "Kilograms of Carbon";

figure(7)

tl = tiledlayout(4,3);
title(tl,"Time Series Analysis of Critical Model State Variables Under Different Wind Conditions")

nexttile([2,3]);
hold on
plot(t,carbon_to_methane(MTOT(:,18,3))) % Drown
plot(t,carbon_to_methane(MTOT(:,12,10))) % Erode
plot(t,carbon_to_methane(MTOT(:,5,5))) % Stable
plot(t,carbon_to_methane(MTOT(:,1,5))) % Fill
title("Methane Emissions");
ylabel("Kilograms of Methane");
xlabel("Years")
legend("Marsh Vertical Drowning","Marsh Retreat (Lateral Erosion)","Marsh Expansion to Equilibrium","Marsh Expansion to Full",'Location','northwest')
hold off

nexttile([2,1]);
hold on
plot(t,WM(:,18,3)) % Drown
plot(t,WM(:,12,10)) % Erode
plot(t,WM(:,5,5)) % Stable
plot(t,WM(:,1,5)) % Fill
title("Width of the Marsh");
ylabel("Meters");
hold off

nexttile([2,1]);
hold on
plot(t,HM(:,18,3)) % Drown
plot(t,HM(:,12,10)) % Erode
plot(t,HM(:,5,5)) % Stable
plot(t,HM(:,1,5)) % Fill
title("Depth of the Marsh Below MHW");
ylabel("Meters");
hold off

nexttile([2,1]);

upperlimit = zeros(1,n_iter) + ((17.7+9.7)/(3600*10^6));
lowerlimit = zeros(1,n_iter) + ((17.7-9.7)/(3600*10^6));

hold on
plot(t,carbon_to_methane(MFLUX(:,18,3))) % Drown
plot(t,carbon_to_methane(MFLUX(:,12,10))) % Erode
plot(t,carbon_to_methane(MFLUX(:,5,5))) % Stable
plot(t,carbon_to_methane(MFLUX(:,1,5))) % Fill
lit = plot(t,upperlimit,"--k",t,lowerlimit,"--k");
title("Methane Flux");
ylabel("Kilograms of Methane per Meter^2 per Second");
hold off

legend(lit,"Comer-Warner et al. 2022", "Location","southoutside")

toc