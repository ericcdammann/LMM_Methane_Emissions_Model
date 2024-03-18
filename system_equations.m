function dydt = system_equations(~,y,p)

%% State Variables

fetch = y(1);
hl = y(2);
hm = y(3);
wb = y(4);
cl = y(5);
wm = wb-fetch;

%% Model Parameters

Co = p(1);
RSLR = p(2);
rhos = p(3);
period = p(4);
ws = p(5);
tcr = p(6);
wind = p(7);
ka = p(8);
ke = p(9);
amp = p(10);
rhom = p(11);
lamda = p(12);
beta = p(13);
k = p(14);
wtidal = p(15);

%% Organic Sediment Production

BMax = 2.500;
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
       
%% Average Depths

fac = min(1,hl/(2*amp));
fac2 = min(1,hm/(2*amp));
Hf = (hl+(hl-fac*2*amp))/2;
Hm = (hm+(hm-fac2*2*amp))/2;

%% Sediment Concentration in the Lagoon

tau_1 =  lagoon_bed_shear_stress(fetch,wind,Hf);
tau_lagoon = max((tau_1-tcr)/tcr,0)*lamda;
Cr = rhos*tau_lagoon/(1+tau_lagoon);

%% Sediment Concentration in the Marsh

if Hm > 1e-4
    
   tau_2 =  marsh_bed_shear_stress(fetch,wind,Hm,Bfrac);
    
else 
        
   tau_2 = 0;
    
end
    
tau_marsh = max((tau_2-tcr)/tcr,0)*lamda;
Cm = rhos*tau_marsh/(1+tau_marsh);

%% Sediment Flux Between the Lagoon and Marsh

Fm = (Cr-Cm)*min(2*amp,hm)/period/rhom;

%% Sediment Flux Between the Open Ocean and Lagoon

Fc = (Cr-Co)*(fac*2*amp)/period/rhom;

%% Marsh Platform Erosion rate

dist = 10;
hs = hm+(hl-hm)*(1-exp(-dist*0.1/hl));
W = wave_power(amp,wind,fetch,hs);
Be = ke*W/(hs-hm);

%% Marsh Platform Progradation Rate

Ba = ka*Cr*ws/rhom;

%% Marsh Platform Upland Migration Rate

Umig = (Fm+O)/beta;

%% Biogeochemical Submodule

chilab = 1-chiref;
Rlab  = AMC*chilab;
phi = wm/(2*wtidal);
mflux = phi*k*cl;

%% Morphodynamic Equations

dydt = zeros(7,1);
dydt(1) = Be-Ba;
dydt(2) = -(Be-Ba)*(hl-hm)/fetch+Fm*wm/fetch+Fc+RSLR;
dydt(3) = -Fm-O+RSLR;
dydt(4) = Umig;
dydt(5) = (Rlab-(k*cl));
dydt(6) = mflux*wm;
dydt(7) = mflux;

end