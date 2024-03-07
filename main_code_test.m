clc; clear variables;
tic 

%% Computational Parameters

years = 100;
n_iter = 1000;
scaled_n_iter = n_iter/25;
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
Co = 0.05;
wind = 4.58974;
ka = 2;
ke = 0.16/(365*24*3600);
amp = 1.4/2;
RSLR = 7.82051*(10^-3)/(3600*24*365);
rhom = 1000;
lamda = 0.0001;
beta = 10^10;
kk = 0.012/(24*60*60);
wtidal = 5*wm_init;

%% Preallocationg Arrays

WL = zeros(n_iter,1);
HM = zeros(n_iter,1);
HL = zeros(n_iter,1);
WB = zeros(n_iter,1);
CTOT = zeros(n_iter,1);
CLAB = zeros(n_iter,1);
CREF = zeros(n_iter,1);
MFLUX = zeros(n_iter,1);
MTOT = zeros(n_iter,1);
WL(1) = wl_init;
HM(1) = hm_init;
HL(1) = hl_init;
WB(1) = wb_init;
CTOT(1) = ctot_init;
CLAB(1) = clab_init;
CREF(1) = cref_init;

%% Solving Equations

for j=1:n_iter-1

    fetch = WL(j);
    hl = HL(j);
    hm = HM(j);
    wb = WB(j);
    wm = wb-fetch;
    ctot = CTOT(j);
    clab = CLAB(j);
    cref = CREF(j);
    mtot = MTOT(j);

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

    tau_1 =  lagoon_bed_shear_stress(fetch,wind,Hf);
    tau_lagoon = max((tau_1-tcr)/tcr,0)*lamda;
    Cr = rhos*tau_lagoon/(1+tau_lagoon);

    % Sediment Concentration in the Marsh

    if Hm > 1e-4
    
       tau_2 =  marsh_bed_shear_stress(fetch,wind,Hm,Bfrac);
    
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
    W = wave_power(amp,wind,fetch,hs);
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
    dHMdt = -Fm-O+RSLR;
    dHLdt = -(Be-Ba)*(hl-hm)/fetch+Fm*wm/fetch+Fc+RSLR;
    dWBdt = Umig;

    % Foward Euler Method

    if (wb+h*dWBdt)-(fetch+h*dWLdt) <=1 || hm+h*dHMdt>=2*amp || hm+h*dHMdt >= hl+h*dHLdt

       break

    else

       WL(j+1) = fetch+h*dWLdt;
       HM(j+1) = hm+h*dHMdt;
       HL(j+1) = hl+h*dHLdt;
       WB(j+1) = wb+h*dWBdt;
       CTOT(j+1) = ctot+h*dCTOTdt;
       CLAB(j+1) = clab+h*dCLABdt;
       CREF(j+1) = cref+h*dCREFdt;
       MTOT(j+1) = mtot+h*mflux*wm; 

    end

end