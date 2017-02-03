function finalArray = iterateOneSimStep(rhoValue,eValue,sValue,MeqValue,fileID,path)

s       = sValue;
e       = eValue;       % homing rate
rho     = rhoValue;     % NHEJ rate
M_eq    = MeqValue;     % equilibrium adult population size
[rho,e,s,M_eq];
% Mosquito life cycle parameters:

beta     = 32;           % eggs laid per day by female mosquito
beta_HH  = 0;            % eggs laid per day by HH female mosquito
beta_Hh  = beta;         % eggs laid per day by Hh female mosquito
beta_HR  = beta;         % eggs laid per day by HR female mosquito
beta_hh  = beta;         % eggs laid per day by hh female mosquito
beta_hR  = beta;         % eggs laid per day by hR female mosquito
beta_RR  = beta;         % eggs laid per day by RR female mosquito
mu_m     = 0.123;        % adult mosquito daily mortality
mu_m_HH  = mu_m;
mu_m_Hh  = mu_m;
mu_m_HR  = mu_m;
mu_m_hh  = mu_m;
mu_m_hR  = mu_m;
mu_m_RR  = mu_m;
rm       = 1.096;        % daily population growth rate
T0       = 1;
TL       = 14;
TP       = 1;

% Derived mosquito life cycle parameters:

g       = T0 + TL + TP + 1/mu_m;  % average generation time
Rm      = rm^g;
mu_l    = 1 - (Rm*mu_m/((beta/2)*(1-mu_m)))^(1/(T0+TL+TP));        % larval daily mortality
theta0  = (1 - mu_l)^T0;          % probability of surviving egg stage
thetaL  = (1 - mu_l)^TL;          % probability of surviving larval stage
thetaP  = (1 - mu_l)^TP;          % probability of surviving pupal stage
%M_eq    = 10000;          % equilibrium adult population size
migrate = 0.01/g;
alpha   = (beta*theta0*(M_eq/2)/(Rm-1))*((1-thetaL/Rm)/(1 - (thetaL/Rm)^(1/TL)));

% Homing parameters:

%e       = 0.98;     % homing rate
%rho     = 0.01;  % NHEJ rate

% Equilibrium population sizes:
L_eq    = round(alpha*(Rm-1));

% Initial conditions (population 1):

LHH1(1:50)    = 0;
LHh1(1:50)    = 0;
LHR1(1:50)    = 0;
Lhh1(1:50)    = L_eq;
LhR1(1:50)    = 0;
LRR1(1:50)    = 0;
L1(1:50)      = L_eq;

MHH1(1:50)    = 0;
MHh1(1:50)    = 0;
MHR1(1:50)    = 0;
Mhh1(1:50)    = M_eq/2;
MhR1(1:50)    = 0;
MRR1(1:50)    = 0;

MHH_HH1(1:50)    = 0;
MHH_Hh1(1:50)    = 0;
MHH_HR1(1:50)    = 0;
MHH_hh1(1:50)    = 0;
MHH_hR1(1:50)    = 0;
MHH_RR1(1:50)    = 0;

MHh_HH1(1:50)    = 0;
MHh_Hh1(1:50)    = 0;
MHh_HR1(1:50)    = 0;
MHh_hh1(1:50)    = 0;
MHh_hR1(1:50)    = 0;
MHh_RR1(1:50)    = 0;

MHR_HH1(1:50)    = 0;
MHR_Hh1(1:50)    = 0;
MHR_HR1(1:50)    = 0;
MHR_hh1(1:50)    = 0;
MHR_hR1(1:50)    = 0;
MHR_RR1(1:50)    = 0;

Mhh_HH1(1:50)    = 0;
Mhh_Hh1(1:50)    = 0;
Mhh_HR1(1:50)    = 0;
Mhh_hh1(1:50)    = M_eq/2;
Mhh_hR1(1:50)    = 0;
Mhh_RR1(1:50)    = 0;

MhR_HH1(1:50)    = 0;
MhR_Hh1(1:50)    = 0;
MhR_HR1(1:50)    = 0;
MhR_hh1(1:50)    = 0;
MhR_hR1(1:50)    = 0;
MhR_RR1(1:50)    = 0;

MRR_HH1(1:50)    = 0;
MRR_Hh1(1:50)    = 0;
MRR_HR1(1:50)    = 0;
MRR_hh1(1:50)    = 0;
MRR_hR1(1:50)    = 0;
MRR_RR1(1:50)    = 0;

% Initial conditions (population 2):

LHH2(1:50)    = 0;
LHh2(1:50)    = 0;
LHR2(1:50)    = 0;
Lhh2(1:50)    = L_eq;
LhR2(1:50)    = 0;
LRR2(1:50)    = 0;
L2(1:50)      = L_eq;

MHH2(1:50)    = 0;
MHh2(1:50)    = 0;
MHR2(1:50)    = 0;
Mhh2(1:50)    = M_eq/2;
MhR2(1:50)    = 0;
MRR2(1:50)    = 0;

MHH_HH2(1:50)    = 0;
MHH_Hh2(1:50)    = 0;
MHH_HR2(1:50)    = 0;
MHH_hh2(1:50)    = 0;
MHH_hR2(1:50)    = 0;
MHH_RR2(1:50)    = 0;

MHh_HH2(1:50)    = 0;
MHh_Hh2(1:50)    = 0;
MHh_HR2(1:50)    = 0;
MHh_hh2(1:50)    = 0;
MHh_hR2(1:50)    = 0;
MHh_RR2(1:50)    = 0;

MHR_HH2(1:50)    = 0;
MHR_Hh2(1:50)    = 0;
MHR_HR2(1:50)    = 0;
MHR_hh2(1:50)    = 0;
MHR_hR2(1:50)    = 0;
MHR_RR2(1:50)    = 0;

Mhh_HH2(1:50)    = 0;
Mhh_Hh2(1:50)    = 0;
Mhh_HR2(1:50)    = 0;
Mhh_hh2(1:50)    = M_eq/2;
Mhh_hR2(1:50)    = 0;
Mhh_RR2(1:50)    = 0;

MhR_HH2(1:50)    = 0;
MhR_Hh2(1:50)    = 0;
MhR_HR2(1:50)    = 0;
MhR_hh2(1:50)    = 0;
MhR_hR2(1:50)    = 0;
MhR_RR2(1:50)    = 0;

MRR_HH2(1:50)    = 0;
MRR_Hh2(1:50)    = 0;
MRR_HR2(1:50)    = 0;
MRR_hh2(1:50)    = 0;
MRR_hR2(1:50)    = 0;
MRR_RR2(1:50)    = 0;

% Initial conditions (population 3):

LHH3(1:50)    = 0;
LHh3(1:50)    = 0;
LHR3(1:50)    = 0;
Lhh3(1:50)    = L_eq;
LhR3(1:50)    = 0;
LRR3(1:50)    = 0;
L3(1:50)      = L_eq;

MHH3(1:50)    = 0;
MHh3(1:50)    = 0;
MHR3(1:50)    = 0;
Mhh3(1:50)    = M_eq/2;
MhR3(1:50)    = 0;
MRR3(1:50)    = 0;

MHH_HH3(1:50)    = 0;
MHH_Hh3(1:50)    = 0;
MHH_HR3(1:50)    = 0;
MHH_hh3(1:50)    = 0;
MHH_hR3(1:50)    = 0;
MHH_RR3(1:50)    = 0;

MHh_HH3(1:50)    = 0;
MHh_Hh3(1:50)    = 0;
MHh_HR3(1:50)    = 0;
MHh_hh3(1:50)    = 0;
MHh_hR3(1:50)    = 0;
MHh_RR3(1:50)    = 0;

MHR_HH3(1:50)    = 0;
MHR_Hh3(1:50)    = 0;
MHR_HR3(1:50)    = 0;
MHR_hh3(1:50)    = 0;
MHR_hR3(1:50)    = 0;
MHR_RR3(1:50)    = 0;

Mhh_HH3(1:50)    = 0;
Mhh_Hh3(1:50)    = 0;
Mhh_HR3(1:50)    = 0;
Mhh_hh3(1:50)    = M_eq/2;
Mhh_hR3(1:50)    = 0;
Mhh_RR3(1:50)    = 0;

MhR_HH3(1:50)    = 0;
MhR_Hh3(1:50)    = 0;
MhR_HR3(1:50)    = 0;
MhR_hh3(1:50)    = 0;
MhR_hR3(1:50)    = 0;
MhR_RR3(1:50)    = 0;

MRR_HH3(1:50)    = 0;
MRR_Hh3(1:50)    = 0;
MRR_HR3(1:50)    = 0;
MRR_hh3(1:50)    = 0;
MRR_hR3(1:50)    = 0;
MRR_RR3(1:50)    = 0;

% Initial conditions (population 4):

LHH4(1:50)    = 0;
LHh4(1:50)    = 0;
LHR4(1:50)    = 0;
Lhh4(1:50)    = L_eq;
LhR4(1:50)    = 0;
LRR4(1:50)    = 0;
L4(1:50)      = L_eq;

MHH4(1:50)    = 0;
MHh4(1:50)    = 0;
MHR4(1:50)    = 0;
Mhh4(1:50)    = M_eq/2;
MhR4(1:50)    = 0;
MRR4(1:50)    = 0;

MHH_HH4(1:50)    = 0;
MHH_Hh4(1:50)    = 0;
MHH_HR4(1:50)    = 0;
MHH_hh4(1:50)    = 0;
MHH_hR4(1:50)    = 0;
MHH_RR4(1:50)    = 0;

MHh_HH4(1:50)    = 0;
MHh_Hh4(1:50)    = 0;
MHh_HR4(1:50)    = 0;
MHh_hh4(1:50)    = 0;
MHh_hR4(1:50)    = 0;
MHh_RR4(1:50)    = 0;

MHR_HH4(1:50)    = 0;
MHR_Hh4(1:50)    = 0;
MHR_HR4(1:50)    = 0;
MHR_hh4(1:50)    = 0;
MHR_hR4(1:50)    = 0;
MHR_RR4(1:50)    = 0;

Mhh_HH4(1:50)    = 0;
Mhh_Hh4(1:50)    = 0;
Mhh_HR4(1:50)    = 0;
Mhh_hh4(1:50)    = M_eq/2;
Mhh_hR4(1:50)    = 0;
Mhh_RR4(1:50)    = 0;

MhR_HH4(1:50)    = 0;
MhR_Hh4(1:50)    = 0;
MhR_HR4(1:50)    = 0;
MhR_hh4(1:50)    = 0;
MhR_hR4(1:50)    = 0;
MhR_RR4(1:50)    = 0;

MRR_HH4(1:50)    = 0;
MRR_Hh4(1:50)    = 0;
MRR_HR4(1:50)    = 0;
MRR_hh4(1:50)    = 0;
MRR_hR4(1:50)    = 0;
MRR_RR4(1:50)    = 0;

% Daily difference equations:

max_t = 1000;

for t = 51:max_t

    % Population 1:
    thetaLA = thetaL;
    for i=1:TL
        thetaLA = thetaLA*((alpha/(alpha+L1(t-i)))^(1/TL));
    end
    
    L1(t) = max(0,binornd(L1(t-1),(1-mu_l)*((alpha/(alpha+L1(t-1)))^(1/TL))) + binornd(poissrnd( beta_HH*MHH_HH1(t-T0) + beta_HH*MHH_Hh1(t-T0) + beta_HH*MHH_HR1(t-T0) + beta_HH*MHH_hh1(t-T0) + beta_HH*MHH_hR1(t-T0) + beta_HH*MHH_RR1(t-T0)  +  beta_Hh*MHh_HH1(t-T0) + beta_Hh*MHh_Hh1(t-T0) + beta_Hh*MHh_HR1(t-T0) + beta_Hh*MHh_hh1(t-T0) + beta_Hh*MHh_hR1(t-T0) + beta_Hh*MHh_RR1(t-T0)  +  beta_HR*MHR_HH1(t-T0) + beta_HR*MHR_Hh1(t-T0) + beta_HR*MHR_HR1(t-T0) + beta_HR*MHR_hh1(t-T0) + beta_HR*MHR_hR1(t-T0) + beta_HR*MHR_RR1(t-T0)  +  beta_hh*Mhh_HH1(t-T0) + beta_hh*Mhh_Hh1(t-T0) + beta_hh*Mhh_HR1(t-T0) + beta_hh*Mhh_hh1(t-T0) + beta_hh*Mhh_hR1(t-T0) + beta_hh*Mhh_RR1(t-T0)  +  beta_hR*MhR_HH1(t-T0) + beta_hR*MhR_Hh1(t-T0) + beta_hR*MhR_HR1(t-T0) + beta_hR*MhR_hh1(t-T0) + beta_hR*MhR_hR1(t-T0) + beta_hR*MhR_RR1(t-T0)  +  beta_RR*MRR_HH1(t-T0) + beta_RR*MRR_Hh1(t-T0) + beta_RR*MRR_HR1(t-T0) + beta_RR*MRR_hh1(t-T0) + beta_RR*MRR_hR1(t-T0) + beta_RR*MRR_RR1(t-T0)  ),theta0) - binornd(poissrnd(beta_HH*MHH_HH1(t-T0-TL) + beta_HH*MHH_Hh1(t-T0-TL) + beta_HH*MHH_HR1(t-T0-TL) + beta_HH*MHH_hh1(t-T0-TL) + beta_HH*MHH_hR1(t-T0-TL) + beta_HH*MHH_RR1(t-T0-TL)  +  beta_Hh*MHh_HH1(t-T0-TL) + beta_Hh*MHh_Hh1(t-T0-TL) + beta_Hh*MHh_HR1(t-T0-TL) + beta_Hh*MHh_hh1(t-T0-TL) + beta_Hh*MHh_hR1(t-T0-TL) + beta_Hh*MHh_RR1(t-T0-TL)  +  beta_HR*MHR_HH1(t-T0-TL) + beta_HR*MHR_Hh1(t-T0-TL) + beta_HR*MHR_HR1(t-T0-TL) + beta_HR*MHR_hh1(t-T0-TL) + beta_HR*MHR_hR1(t-T0-TL) + beta_HR*MHR_RR1(t-T0-TL)  +  beta_hh*Mhh_HH1(t-T0-TL) + beta_hh*Mhh_Hh1(t-T0-TL) + beta_hh*Mhh_HR1(t-T0-TL) + beta_hh*Mhh_hh1(t-T0-TL) + beta_hh*Mhh_hR1(t-T0-TL) + beta_hh*Mhh_RR1(t-T0-TL)  +  beta_hR*MhR_HH1(t-T0-TL) + beta_hR*MhR_Hh1(t-T0-TL) + beta_hR*MhR_HR1(t-T0-TL) + beta_hR*MhR_hh1(t-T0-TL) + beta_hR*MhR_hR1(t-T0-TL) + beta_hR*MhR_RR1(t-T0-TL)  +  beta_RR*MRR_HH1(t-T0-TL) + beta_RR*MRR_Hh1(t-T0-TL) + beta_RR*MRR_HR1(t-T0-TL) + beta_RR*MRR_hh1(t-T0-TL) + beta_RR*MRR_hR1(t-T0-TL) + beta_RR*MRR_RR1(t-T0-TL) ),theta0*thetaLA));
    
    thetaLA = thetaL;
    for i=1:TL
        thetaLA = thetaLA*((alpha/(alpha+L1(t-i-TP)))^(1/TL));
    end
    
    % Eggs produced by HH females:
    
    E_HHHH1 = poissrnd(beta_HH*MHH_HH1(t-T0-TL-TP));
    E_HHHH_HH1 = E_HHHH1;

    E_HHHh1 = mnrnd(poissrnd(beta_HH*MHH_Hh1(t-T0-TL-TP)),[((1+e)/2),((1-e-rho)/2),(rho/2)]);
    E_HHHh_HH1 = E_HHHh1(1); E_HHHh_Hh1 = E_HHHh1(2); E_HHHh_HR1 = E_HHHh1(3);

    E_HHHR1 = mnrnd(poissrnd(beta_HH*MHH_HR1(t-T0-TL-TP)),[(1/2),(1/2)]);
    E_HHHR_HH1 = E_HHHR1(1); E_HHHR_HR1 = E_HHHR1(2);

    E_HHhh1 = poissrnd(beta_HH*MHH_hh1(t-T0-TL-TP));
    E_HHhh_Hh1 = E_HHhh1;

    E_HHhR1 = mnrnd(poissrnd(beta_HH*MHH_hR1(t-T0-TL-TP)),[(1/2),(1/2)]);
    E_HHhR_Hh1 = E_HHhR1(1); E_HHhR_HR1 = E_HHhR1(2);

    E_HHRR1 = poissrnd(beta_HH*MHH_RR1(t-T0-TL-TP));
    E_HHRR_HR1 = E_HHRR1;

    % Eggs produced by Hh females:
    
    E_HhHH1 = mnrnd(poissrnd(beta_Hh*MHh_HH1(t-T0-TL-TP)),[((1+e)/2),((1-e-rho)/2),(rho/2)]);
    E_HhHH_HH1 = E_HhHH1(1); E_HhHH_Hh1 = E_HhHH1(2); E_HhHH_HR1 = E_HhHH1(3);

    E_HhHh1 = mnrnd(poissrnd(beta_Hh*MHh_Hh1(t-T0-TL-TP)),[(((1+e)^2)/4),((1+e)*(1-e-rho)/2),((1+e)*rho/2),(((1-e-rho)^2)/4),((1-e-rho)*rho/2),((rho^2)/4)]);
    E_HhHh_HH1 = E_HhHh1(1); E_HhHh_Hh1 = E_HhHh1(2); E_HhHh_HR1 = E_HhHh1(3); E_HhHh_hh1 = E_HhHh1(4); E_HhHh_hR1 = E_HhHh1(5); E_HhHh_RR1 = E_HhHh1(6);

    E_HhHR1 = mnrnd(poissrnd(beta_Hh*MHh_HR1(t-T0-TL-TP)),[((1+e)/4),((1-e-rho)/4),((1+e+rho)/4),((1-e-rho)/4),(rho/4)]);
    E_HhHR_HH1 = E_HhHR1(1); E_HhHR_Hh1 = E_HhHR1(2); E_HhHR_HR1 = E_HhHR1(3); E_HhHR_hR1 = E_HhHR1(4); E_HhHR_RR1 = E_HhHR1(5);

    E_Hhhh1 = mnrnd(poissrnd(beta_Hh*MHh_hh1(t-T0-TL-TP)),[((1+e)/2),((1-e-rho)/2),(rho/2)]);
    E_Hhhh_Hh1 = E_Hhhh1(1); E_Hhhh_hh1 = E_Hhhh1(2); E_Hhhh_hR1 = E_Hhhh1(3);

    E_HhhR1 = mnrnd(poissrnd(beta_Hh*MHh_hR1(t-T0-TL-TP)),[((1+e)/4),((1+e)/4),((1-e-rho)/4),((1-e)/4),(rho/4)]);
    E_HhhR_Hh1 = E_HhhR1(1); E_HhhR_HR1 = E_HhhR1(2); E_HhhR_hh1 = E_HhhR1(3); E_HhhR_hR1 = E_HhhR1(4); E_HhhR_RR1 = E_HhhR1(5);
    
    E_HhRR1 = mnrnd(poissrnd(beta_Hh*MHh_RR1(t-T0-TL-TP)),[((1+e)/2),((1-e-rho)/2),(rho/2)]);
    E_HhRR_HR1 = E_HhRR1(1); E_HhRR_hR1 = E_HhRR1(2); E_HhRR_RR1 = E_HhRR1(3);
    
    % Eggs produced by HR females:
    
    E_HRHH1 = mnrnd(poissrnd(beta_HR*MHR_HH1(t-T0-TL-TP)),[(1/2),(1/2)]);
    E_HRHH_HH1 = E_HRHH1(1); E_HRHH_HR1 = E_HRHH1(2);

    E_HRHh1 = mnrnd(poissrnd(beta_HR*MHR_Hh1(t-T0-TL-TP)),[((1+e)/4),((1-e-rho)/4),((1+e+rho)/4),((1-e-rho)/4),(rho/4)]);
    E_HRHh_HH1 = E_HRHh1(1); E_HRHh_Hh1 = E_HRHh1(2); E_HRHh_HR1 = E_HRHh1(3); E_HRHh_hR1 = E_HRHh1(4); E_HRHh_RR1 = E_HRHh1(5);

    E_HRHR1 = mnrnd(poissrnd(beta_HR*MHR_HR1(t-T0-TL-TP)),[(1/4),(1/2),(1/4)]);
    E_HRHR_HH1 = E_HRHR1(1); E_HRHR_HR1 = E_HRHR1(2); E_HRHR_RR1 = E_HRHR1(3);

    E_HRhh1 = mnrnd(poissrnd(beta_HR*MHR_hh1(t-T0-TL-TP)),[(1/2),(1/2)]);
    E_HRhh_Hh1 = E_HRhh1(1); E_HRhh_hR1 = E_HRhh1(2);

    E_HRhR1 = mnrnd(poissrnd(beta_HR*MHR_hR1(t-T0-TL-TP)),[(1/4),(1/4),(1/4),(1/4)]);
    E_HRhR_Hh1 = E_HRhR1(1); E_HRhR_HR1 = E_HRhR1(2); E_HRhR_hR1 = E_HRhR1(3); E_HRhR_RR1 = E_HRhR1(4);
    
    E_HRRR1 = mnrnd(poissrnd(beta_HR*MHR_RR1(t-T0-TL-TP)),[(1/2),(1/2)]);
    E_HRRR_HR1 = E_HRRR1(1); E_HRRR_RR1 = E_HRRR1(2);
    
    % Eggs produced by hh females:
    
    E_hhHH1 = poissrnd(beta_hh*Mhh_HH1(t-T0-TL-TP));
    E_hhHH_Hh1 = E_hhHH1;

    E_hhHh1 = mnrnd(poissrnd(beta_hh*Mhh_Hh1(t-T0-TL-TP)),[((1+e)/2),((1-e-rho)/2),(rho/2)]);
    E_hhHh_Hh1 = E_hhHh1(1); E_hhHh_hh1 = E_hhHh1(2); E_hhHh_hR1 = E_hhHh1(3);

    E_hhHR1 = mnrnd(poissrnd(beta_hh*Mhh_HR1(t-T0-TL-TP)),[(1/2),(1/2)]);
    E_hhHR_Hh1 = E_hhHR1(1); E_hhHR_hR1 = E_hhHR1(2);

    E_hhhh1 = poissrnd(beta_hh*Mhh_hh1(t-T0-TL-TP));
    E_hhhh_hh1 = E_hhhh1;

    E_hhhR1 = mnrnd(poissrnd(beta_hh*Mhh_hR1(t-T0-TL-TP)),[(1/2),(1/2)]);
    E_hhhR_hh1 = E_hhhR1(1); E_hhhR_hR1 = E_hhhR1(2);

    E_hhRR1 = poissrnd(beta_hh*Mhh_RR1(t-T0-TL-TP));
    E_hhRR_hR1 = E_hhRR1;
    
    % Eggs produced by hR females:
    
    E_hRHH1 = mnrnd(poissrnd(beta_hR*MhR_HH1(t-T0-TL-TP)),[(1/2),(1/2)]);
    E_hRHH_Hh1 = E_hRHH1(1); E_hRHH_HR1 = E_hRHH1(2);

    E_hRHh1 = mnrnd(poissrnd(beta_hR*MhR_Hh1(t-T0-TL-TP)),[((1+e)/4),((1+e)/4),((1-e-rho)/4),((1-e)/4),(rho/4)]);
    E_hRHh_Hh1 = E_hRHh1(1); E_hRHh_HR1 = E_hRHh1(2); E_hRHh_hh1 = E_hRHh1(3); E_hRHh_hR1 = E_hRHh1(4); E_hRHh_RR1 = E_hRHh1(5);

    E_hRHR1 = mnrnd(poissrnd(beta_hR*MhR_HR1(t-T0-TL-TP)),[(1/4),(1/4),(1/4),(1/4)]);
    E_hRHR_Hh1 = E_hRHR1(1); E_hRHR_HR1 = E_hRHR1(2); E_hRHR_hR1 = E_hRHR1(3); E_hRHR_RR1 = E_hRHR1(4);

    E_hRhh1 = mnrnd(poissrnd(beta_hR*MHR_hh1(t-T0-TL-TP)),[(1/2),(1/2)]);
    E_hRhh_hh1 = E_hRhh1(1); E_hRhh_hR1 = E_hRhh1(2);

    E_hRhR1 = mnrnd(poissrnd(beta_hR*MhR_hR1(t-T0-TL-TP)),[(1/4),(1/2),(1/4)]);
    E_hRhR_hh1 = E_hRhR1(1); E_hRhR_hR1 = E_hRhR1(2); E_hRhR_RR1 = E_hRhR1(3);
    
    E_hRRR1 = mnrnd(poissrnd(beta_hR*MhR_RR1(t-T0-TL-TP)),[(1/2),(1/2)]);
    E_hRRR_hR1 = E_hRRR1(1); E_hRRR_RR1 = E_hRRR1(2);
    
    % Eggs produced by RR females:
    
    E_RRHH1 = poissrnd(beta_RR*MRR_HH1(t-T0-TL-TP));
    E_RRHH_HR1 = E_RRHH1;

    E_RRHh1 = mnrnd(poissrnd(beta_RR*MRR_Hh1(t-T0-TL-TP)),[((1+e)/2),((1-e-rho)/2),(rho/2)]);
    E_RRHh_HR1 = E_RRHh1(1); E_RRHh_hR1 = E_RRHh1(2); E_RRHh_RR1 = E_RRHh1(3);

    E_RRHR1 = mnrnd(poissrnd(beta_RR*MRR_HR1(t-T0-TL-TP)),[(1/2),(1/2)]);
    E_RRHR_HR1 = E_RRHR1(1); E_RRHR_RR1 = E_RRHR1(2);

    E_RRhh1 = poissrnd(beta_RR*MRR_hh1(t-T0-TL-TP));
    E_RRhh_hR1 = E_RRhh1;

    E_RRhR1 = mnrnd(poissrnd(beta_RR*MRR_hR1(t-T0-TL-TP)),[(1/2),(1/2)]);
    E_RRhR_hR1 = E_RRhR1(1); E_RRhR_RR1 = E_RRhR1(2);

    E_RRRR1 = poissrnd(beta_RR*MRR_RR1(t-T0-TL-TP));
    E_RRRR_RR1 = E_RRRR1;
    
    % Adult male & female newborns:
    
    MHH_Births1 = mnrnd(binornd((E_HHHH_HH1 + E_HHHh_HH1 + E_HHHR_HH1 + E_HhHH_HH1 + E_HhHh_HH1 + E_HhHR_HH1 + E_HRHH_HH1 + E_HRHh_HH1 + E_HRHR_HH1),theta0*thetaLA*thetaP*(1-mu_m_HH)),[(1/2),(1/2)]);
    MHh_Births1 = mnrnd(binornd((E_HHHh_Hh1 + E_HHhh_Hh1 + E_HHhR_Hh1 + E_HhHH_Hh1 + E_HhHh_Hh1 + E_HhHR_Hh1 + E_Hhhh_Hh1 + E_HhhR_Hh1 + E_HRHh_Hh1 + E_HRhh_Hh1 + E_HRhR_Hh1 + E_hhHH_Hh1 + E_hhHh_Hh1 + E_hhHR_Hh1 + E_hRHH_Hh1 + E_hRHh_Hh1 + E_hRHR_Hh1),theta0*thetaLA*thetaP*(1-mu_m_Hh)),[(1/2),(1/2)]);
    MHR_Births1 = mnrnd(binornd((E_HHHh_HR1 + E_HHHR_HR1 + E_HHhR_HR1 + E_HHRR_HR1 + E_HhHH_HR1 + E_HhHh_HR1 + E_HhHR_HR1 + E_HhhR_HR1 + E_HhRR_HR1 + E_HRHH_HR1 + E_HRHh_HR1 + E_HRHR_HR1 + E_HRhR_HR1 + E_HRRR_HR1 + E_hRHH_HR1 + E_hRHh_HR1 + E_hRHR_HR1 + E_RRHH_HR1 + E_RRHh_HR1 + E_RRHR_HR1),theta0*thetaLA*thetaP*(1-mu_m_HR)),[(1/2),(1/2)]);
    Mhh_Births1 = mnrnd(binornd((E_HhHh_hh1 + E_Hhhh_hh1 + E_HhhR_hh1 + E_hhHh_hh1 + E_hhhh_hh1 + E_hhhR_hh1 + E_hRHh_hh1 + E_hRhh_hh1 + E_hRhR_hh1),theta0*thetaLA*thetaP*(1-mu_m_hh)),[(1/2),(1/2)]);
    MhR_Births1 = mnrnd(binornd((E_HhHh_hR1 + E_HhHR_hR1 + E_Hhhh_hR1 + E_HhhR_hR1 + E_HhRR_hR1 + E_HRHh_hR1 + E_HRhh_hR1 + E_HRhR_hR1 + E_hhHh_hR1 + E_hhHR_hR1 + E_hhhR_hR1 + E_hhRR_hR1 + E_hRHh_hR1 + E_hRHR_hR1 + E_hRhh_hR1 + E_hRhR_hR1 + E_hRRR_hR1 + E_RRHh_hR1 + E_RRhh_hR1 + E_RRhR_hR1),theta0*thetaLA*thetaP*(1-mu_m_hR)),[(1/2),(1/2)]);
    MRR_Births1 = mnrnd(binornd((E_HhHh_RR1 + E_HhHR_RR1 + E_HhhR_RR1 + E_HhRR_RR1 + E_HRHh_RR1 + E_HRHR_RR1 + E_HRhR_RR1 + E_HRRR_RR1 + E_hRHh_RR1 + E_hRHR_RR1 + E_hRhR_RR1 + E_hRRR_RR1 + E_RRHh_RR1 + E_RRHR_RR1 + E_RRhR_RR1 + E_RRRR_RR1),theta0*thetaLA*thetaP*(1-mu_m_RR)),[(1/2),(1/2)]);

    % Adult male genotypes:
    
    MHH1(t) = max(0,binornd(MHH1(t-1),(1-mu_m_HH)) + MHH_Births1(1));
    MHh1(t) = max(0,binornd(MHh1(t-1),(1-mu_m_Hh)) + MHh_Births1(1));
    MHR1(t) = max(0,binornd(MHR1(t-1),(1-mu_m_HR)) + MHR_Births1(1));
    Mhh1(t) = max(0,binornd(Mhh1(t-1),(1-mu_m_hh)) + Mhh_Births1(1));
    MhR1(t) = max(0,binornd(MhR1(t-1),(1-mu_m_hR)) + MhR_Births1(1));
    MRR1(t) = max(0,binornd(MRR1(t-1),(1-mu_m_RR)) + MRR_Births1(1));
    
    % Adult female genotypes:
    
    X_HH1 = mnrnd(MHH_Births1(2) , [ (MHH1(t-1)/(MHH1(t-1)+MHh1(t-1)+MHR1(t-1)+Mhh1(t-1)+MhR1(t-1)+MRR1(t-1))) , (MHh1(t-1)/(MHH1(t-1)+MHh1(t-1)+MHR1(t-1)+Mhh1(t-1)+MhR1(t-1)+MRR1(t-1))) , (MHR1(t-1)/(MHH1(t-1)+MHh1(t-1)+MHR1(t-1)+Mhh1(t-1)+MhR1(t-1)+MRR1(t-1))) , (Mhh1(t-1)/(MHH1(t-1)+MHh1(t-1)+MHR1(t-1)+Mhh1(t-1)+MhR1(t-1)+MRR1(t-1))) , (MhR1(t-1)/(MHH1(t-1)+MHh1(t-1)+MHR1(t-1)+Mhh1(t-1)+MhR1(t-1)+MRR1(t-1))) , (MRR1(t-1)/(MHH1(t-1)+MHh1(t-1)+MHR1(t-1)+Mhh1(t-1)+MhR1(t-1)+MRR1(t-1))) ]);
    X_Hh1 = mnrnd(MHh_Births1(2) , [ (MHH1(t-1)/(MHH1(t-1)+MHh1(t-1)+MHR1(t-1)+Mhh1(t-1)+MhR1(t-1)+MRR1(t-1))) , (MHh1(t-1)/(MHH1(t-1)+MHh1(t-1)+MHR1(t-1)+Mhh1(t-1)+MhR1(t-1)+MRR1(t-1))) , (MHR1(t-1)/(MHH1(t-1)+MHh1(t-1)+MHR1(t-1)+Mhh1(t-1)+MhR1(t-1)+MRR1(t-1))) , (Mhh1(t-1)/(MHH1(t-1)+MHh1(t-1)+MHR1(t-1)+Mhh1(t-1)+MhR1(t-1)+MRR1(t-1))) , (MhR1(t-1)/(MHH1(t-1)+MHh1(t-1)+MHR1(t-1)+Mhh1(t-1)+MhR1(t-1)+MRR1(t-1))) , (MRR1(t-1)/(MHH1(t-1)+MHh1(t-1)+MHR1(t-1)+Mhh1(t-1)+MhR1(t-1)+MRR1(t-1))) ]);
    X_HR1 = mnrnd(MHR_Births1(2) , [ (MHH1(t-1)/(MHH1(t-1)+MHh1(t-1)+MHR1(t-1)+Mhh1(t-1)+MhR1(t-1)+MRR1(t-1))) , (MHh1(t-1)/(MHH1(t-1)+MHh1(t-1)+MHR1(t-1)+Mhh1(t-1)+MhR1(t-1)+MRR1(t-1))) , (MHR1(t-1)/(MHH1(t-1)+MHh1(t-1)+MHR1(t-1)+Mhh1(t-1)+MhR1(t-1)+MRR1(t-1))) , (Mhh1(t-1)/(MHH1(t-1)+MHh1(t-1)+MHR1(t-1)+Mhh1(t-1)+MhR1(t-1)+MRR1(t-1))) , (MhR1(t-1)/(MHH1(t-1)+MHh1(t-1)+MHR1(t-1)+Mhh1(t-1)+MhR1(t-1)+MRR1(t-1))) , (MRR1(t-1)/(MHH1(t-1)+MHh1(t-1)+MHR1(t-1)+Mhh1(t-1)+MhR1(t-1)+MRR1(t-1))) ]);
    X_hh1 = mnrnd(Mhh_Births1(2) , [ (MHH1(t-1)/(MHH1(t-1)+MHh1(t-1)+MHR1(t-1)+Mhh1(t-1)+MhR1(t-1)+MRR1(t-1))) , (MHh1(t-1)/(MHH1(t-1)+MHh1(t-1)+MHR1(t-1)+Mhh1(t-1)+MhR1(t-1)+MRR1(t-1))) , (MHR1(t-1)/(MHH1(t-1)+MHh1(t-1)+MHR1(t-1)+Mhh1(t-1)+MhR1(t-1)+MRR1(t-1))) , (Mhh1(t-1)/(MHH1(t-1)+MHh1(t-1)+MHR1(t-1)+Mhh1(t-1)+MhR1(t-1)+MRR1(t-1))) , (MhR1(t-1)/(MHH1(t-1)+MHh1(t-1)+MHR1(t-1)+Mhh1(t-1)+MhR1(t-1)+MRR1(t-1))) , (MRR1(t-1)/(MHH1(t-1)+MHh1(t-1)+MHR1(t-1)+Mhh1(t-1)+MhR1(t-1)+MRR1(t-1))) ]);
    X_hR1 = mnrnd(MhR_Births1(2) , [ (MHH1(t-1)/(MHH1(t-1)+MHh1(t-1)+MHR1(t-1)+Mhh1(t-1)+MhR1(t-1)+MRR1(t-1))) , (MHh1(t-1)/(MHH1(t-1)+MHh1(t-1)+MHR1(t-1)+Mhh1(t-1)+MhR1(t-1)+MRR1(t-1))) , (MHR1(t-1)/(MHH1(t-1)+MHh1(t-1)+MHR1(t-1)+Mhh1(t-1)+MhR1(t-1)+MRR1(t-1))) , (Mhh1(t-1)/(MHH1(t-1)+MHh1(t-1)+MHR1(t-1)+Mhh1(t-1)+MhR1(t-1)+MRR1(t-1))) , (MhR1(t-1)/(MHH1(t-1)+MHh1(t-1)+MHR1(t-1)+Mhh1(t-1)+MhR1(t-1)+MRR1(t-1))) , (MRR1(t-1)/(MHH1(t-1)+MHh1(t-1)+MHR1(t-1)+Mhh1(t-1)+MhR1(t-1)+MRR1(t-1))) ]);    
    X_RR1 = mnrnd(MRR_Births1(2) , [ (MHH1(t-1)/(MHH1(t-1)+MHh1(t-1)+MHR1(t-1)+Mhh1(t-1)+MhR1(t-1)+MRR1(t-1))) , (MHh1(t-1)/(MHH1(t-1)+MHh1(t-1)+MHR1(t-1)+Mhh1(t-1)+MhR1(t-1)+MRR1(t-1))) , (MHR1(t-1)/(MHH1(t-1)+MHh1(t-1)+MHR1(t-1)+Mhh1(t-1)+MhR1(t-1)+MRR1(t-1))) , (Mhh1(t-1)/(MHH1(t-1)+MHh1(t-1)+MHR1(t-1)+Mhh1(t-1)+MhR1(t-1)+MRR1(t-1))) , (MhR1(t-1)/(MHH1(t-1)+MHh1(t-1)+MHR1(t-1)+Mhh1(t-1)+MhR1(t-1)+MRR1(t-1))) , (MRR1(t-1)/(MHH1(t-1)+MHh1(t-1)+MHR1(t-1)+Mhh1(t-1)+MhR1(t-1)+MRR1(t-1))) ]);

    MHH1_Temp = [binornd(MHH_HH1(t-1),(1-mu_m_HH)) , binornd(MHH_Hh1(t-1),(1-mu_m_HH)) , binornd(MHH_HR1(t-1),(1-mu_m_HH)) , binornd(MHH_hh1(t-1),(1-mu_m_HH)) , binornd(MHH_hR1(t-1),(1-mu_m_HH)) , binornd(MHH_RR1(t-1),(1-mu_m_HH))] + X_HH1;
    MHH_HH1(t) = MHH1_Temp(1); MHH_Hh1(t) = MHH1_Temp(2); MHH_HR1(t) = MHH1_Temp(3); MHH_hh1(t) = MHH1_Temp(4); MHH_hR1(t) = MHH1_Temp(5); MHH_RR1(t) = MHH1_Temp(6);

    MHh1_Temp = [binornd(MHh_HH1(t-1),(1-mu_m_Hh)) , binornd(MHh_Hh1(t-1),(1-mu_m_Hh)) , binornd(MHh_HR1(t-1),(1-mu_m_Hh)) , binornd(MHh_hh1(t-1),(1-mu_m_Hh)) , binornd(MHh_hR1(t-1),(1-mu_m_Hh)) , binornd(MHh_RR1(t-1),(1-mu_m_Hh))] + X_Hh1;
    MHh_HH1(t) = MHh1_Temp(1); MHh_Hh1(t) = MHh1_Temp(2); MHh_HR1(t) = MHh1_Temp(3); MHh_hh1(t) = MHh1_Temp(4); MHh_hR1(t) = MHh1_Temp(5); MHh_RR1(t) = MHh1_Temp(6);

    MHR1_Temp = [binornd(MHR_HH1(t-1),(1-mu_m_HR)) , binornd(MHR_Hh1(t-1),(1-mu_m_HR)) , binornd(MHR_HR1(t-1),(1-mu_m_HR)) , binornd(MHR_hh1(t-1),(1-mu_m_HR)) , binornd(MHR_hR1(t-1),(1-mu_m_HR)) , binornd(MHR_RR1(t-1),(1-mu_m_HR))] + X_HR1;
    MHR_HH1(t) = MHR1_Temp(1); MHR_Hh1(t) = MHR1_Temp(2); MHR_HR1(t) = MHR1_Temp(3); MHR_hh1(t) = MHR1_Temp(4); MHR_hR1(t) = MHR1_Temp(5); MHR_RR1(t) = MHR1_Temp(6);

    Mhh1_Temp = [binornd(Mhh_HH1(t-1),(1-mu_m_hh)) , binornd(Mhh_Hh1(t-1),(1-mu_m_hh)) , binornd(Mhh_HR1(t-1),(1-mu_m_hh)) , binornd(Mhh_hh1(t-1),(1-mu_m_hh)) , binornd(Mhh_hR1(t-1),(1-mu_m_hh)) , binornd(Mhh_RR1(t-1),(1-mu_m_hh))] + X_hh1;
    Mhh_HH1(t) = Mhh1_Temp(1); Mhh_Hh1(t) = Mhh1_Temp(2); Mhh_HR1(t) = Mhh1_Temp(3); Mhh_hh1(t) = Mhh1_Temp(4); Mhh_hR1(t) = Mhh1_Temp(5); Mhh_RR1(t) = Mhh1_Temp(6);

    MhR1_Temp = [binornd(MhR_HH1(t-1),(1-mu_m_hR)) , binornd(MhR_Hh1(t-1),(1-mu_m_hR)) , binornd(MhR_HR1(t-1),(1-mu_m_hR)) , binornd(MhR_hh1(t-1),(1-mu_m_hR)) , binornd(MhR_hR1(t-1),(1-mu_m_hR)) , binornd(MhR_RR1(t-1),(1-mu_m_hR))] + X_hR1;
    MhR_HH1(t) = MhR1_Temp(1); MhR_Hh1(t) = MhR1_Temp(2); MhR_HR1(t) = MhR1_Temp(3); MhR_hh1(t) = MhR1_Temp(4); MhR_hR1(t) = MhR1_Temp(5); MhR_RR1(t) = MhR1_Temp(6);

    MRR1_Temp = [binornd(MRR_HH1(t-1),(1-mu_m_RR)) , binornd(MRR_Hh1(t-1),(1-mu_m_RR)) , binornd(MRR_HR1(t-1),(1-mu_m_RR)) , binornd(MRR_hh1(t-1),(1-mu_m_RR)) , binornd(MRR_hR1(t-1),(1-mu_m_RR)) , binornd(MRR_RR1(t-1),(1-mu_m_RR))] + X_RR1;
    MRR_HH1(t) = MRR1_Temp(1); MRR_Hh1(t) = MRR1_Temp(2); MRR_HR1(t) = MRR1_Temp(3); MRR_hh1(t) = MRR1_Temp(4); MRR_hR1(t) = MRR1_Temp(5); MRR_RR1(t) = MRR1_Temp(6);
    
    MHH_HH1(t) = max(0,MHH_HH1(t)); MHH_Hh1(t) = max(0,MHH_Hh1(t)); MHH_HR1(t) = max(0,MHH_HR1(t)); MHH_hh1(t) = max(0,MHH_hh1(t)); MHH_hR1(t) = max(0,MHH_hR1(t)); MHH_RR1(t) = max(0,MHH_RR1(t));
    MHh_HH1(t) = max(0,MHh_HH1(t)); MHh_Hh1(t) = max(0,MHh_Hh1(t)); MHh_HR1(t) = max(0,MHh_HR1(t)); MHh_hh1(t) = max(0,MHh_hh1(t)); MHh_hR1(t) = max(0,MHh_hR1(t)); MHh_RR1(t) = max(0,MHh_RR1(t));
    MHR_HH1(t) = max(0,MHR_HH1(t)); MHR_Hh1(t) = max(0,MHR_Hh1(t)); MHR_HR1(t) = max(0,MHR_HR1(t)); MHR_hh1(t) = max(0,MHR_hh1(t)); MHR_hR1(t) = max(0,MHR_hR1(t)); MHR_RR1(t) = max(0,MHR_RR1(t));
    Mhh_HH1(t) = max(0,Mhh_HH1(t)); Mhh_Hh1(t) = max(0,Mhh_Hh1(t)); Mhh_HR1(t) = max(0,Mhh_HR1(t)); Mhh_hh1(t) = max(0,Mhh_hh1(t)); Mhh_hR1(t) = max(0,Mhh_hR1(t)); Mhh_RR1(t) = max(0,Mhh_RR1(t));
    MhR_HH1(t) = max(0,MhR_HH1(t)); MhR_Hh1(t) = max(0,MhR_Hh1(t)); MhR_HR1(t) = max(0,MhR_HR1(t)); MhR_hh1(t) = max(0,MhR_hh1(t)); MhR_hR1(t) = max(0,MhR_hR1(t)); MhR_RR1(t) = max(0,MhR_RR1(t));
    MRR_HH1(t) = max(0,MRR_HH1(t)); MRR_Hh1(t) = max(0,MRR_Hh1(t)); MRR_HR1(t) = max(0,MRR_HR1(t)); MRR_hh1(t) = max(0,MRR_hh1(t)); MRR_hR1(t) = max(0,MRR_hR1(t)); MRR_RR1(t) = max(0,MRR_RR1(t));
     
    if (t==51)
        MHH1(t) = MHH1(t) + M_eq;
        % MHH_HH1(t) = MHH_HH1(t) + M_eq/2;
    end
    
    % Population 2:
    thetaLA = thetaL;
    for i=1:TL
        thetaLA = thetaLA*((alpha/(alpha+L2(t-i)))^(1/TL));
    end
    
    L2(t) = max(0,binornd(L2(t-1),(1-mu_l)*((alpha/(alpha+L2(t-1)))^(1/TL))) + binornd(poissrnd( beta_HH*MHH_HH2(t-T0) + beta_HH*MHH_Hh2(t-T0) + beta_HH*MHH_HR2(t-T0) + beta_HH*MHH_hh2(t-T0) + beta_HH*MHH_hR2(t-T0) + beta_HH*MHH_RR2(t-T0)  +  beta_Hh*MHh_HH2(t-T0) + beta_Hh*MHh_Hh2(t-T0) + beta_Hh*MHh_HR2(t-T0) + beta_Hh*MHh_hh2(t-T0) + beta_Hh*MHh_hR2(t-T0) + beta_Hh*MHh_RR2(t-T0)  +  beta_HR*MHR_HH2(t-T0) + beta_HR*MHR_Hh2(t-T0) + beta_HR*MHR_HR2(t-T0) + beta_HR*MHR_hh2(t-T0) + beta_HR*MHR_hR2(t-T0) + beta_HR*MHR_RR2(t-T0)  +  beta_hh*Mhh_HH2(t-T0) + beta_hh*Mhh_Hh2(t-T0) + beta_hh*Mhh_HR2(t-T0) + beta_hh*Mhh_hh2(t-T0) + beta_hh*Mhh_hR2(t-T0) + beta_hh*Mhh_RR2(t-T0)  +  beta_hR*MhR_HH2(t-T0) + beta_hR*MhR_Hh2(t-T0) + beta_hR*MhR_HR2(t-T0) + beta_hR*MhR_hh2(t-T0) + beta_hR*MhR_hR2(t-T0) + beta_hR*MhR_RR2(t-T0)  +  beta_RR*MRR_HH2(t-T0) + beta_RR*MRR_Hh2(t-T0) + beta_RR*MRR_HR2(t-T0) + beta_RR*MRR_hh2(t-T0) + beta_RR*MRR_hR2(t-T0) + beta_RR*MRR_RR2(t-T0)  ),theta0) - binornd(poissrnd(beta_HH*MHH_HH2(t-T0-TL) + beta_HH*MHH_Hh2(t-T0-TL) + beta_HH*MHH_HR2(t-T0-TL) + beta_HH*MHH_hh2(t-T0-TL) + beta_HH*MHH_hR2(t-T0-TL) + beta_HH*MHH_RR2(t-T0-TL)  +  beta_Hh*MHh_HH2(t-T0-TL) + beta_Hh*MHh_Hh2(t-T0-TL) + beta_Hh*MHh_HR2(t-T0-TL) + beta_Hh*MHh_hh2(t-T0-TL) + beta_Hh*MHh_hR2(t-T0-TL) + beta_Hh*MHh_RR2(t-T0-TL)  +  beta_HR*MHR_HH2(t-T0-TL) + beta_HR*MHR_Hh2(t-T0-TL) + beta_HR*MHR_HR2(t-T0-TL) + beta_HR*MHR_hh2(t-T0-TL) + beta_HR*MHR_hR2(t-T0-TL) + beta_HR*MHR_RR2(t-T0-TL)  +  beta_hh*Mhh_HH2(t-T0-TL) + beta_hh*Mhh_Hh2(t-T0-TL) + beta_hh*Mhh_HR2(t-T0-TL) + beta_hh*Mhh_hh2(t-T0-TL) + beta_hh*Mhh_hR2(t-T0-TL) + beta_hh*Mhh_RR2(t-T0-TL)  +  beta_hR*MhR_HH2(t-T0-TL) + beta_hR*MhR_Hh2(t-T0-TL) + beta_hR*MhR_HR2(t-T0-TL) + beta_hR*MhR_hh2(t-T0-TL) + beta_hR*MhR_hR2(t-T0-TL) + beta_hR*MhR_RR2(t-T0-TL)  +  beta_RR*MRR_HH2(t-T0-TL) + beta_RR*MRR_Hh2(t-T0-TL) + beta_RR*MRR_HR2(t-T0-TL) + beta_RR*MRR_hh2(t-T0-TL) + beta_RR*MRR_hR2(t-T0-TL) + beta_RR*MRR_RR2(t-T0-TL) ),theta0*thetaLA));
    
    thetaLA = thetaL;
    for i=1:TL
        thetaLA = thetaLA*((alpha/(alpha+L2(t-i-TP)))^(1/TL));
    end
    
    % Eggs produced by HH females:
    
    E_HHHH2 = poissrnd(beta_HH*MHH_HH2(t-T0-TL-TP));
    E_HHHH_HH2 = E_HHHH2;

    E_HHHh2 = mnrnd(poissrnd(beta_HH*MHH_Hh2(t-T0-TL-TP)),[((1+e)/2),((1-e-rho)/2),(rho/2)]);
    E_HHHh_HH2 = E_HHHh2(1); E_HHHh_Hh2 = E_HHHh2(2); E_HHHh_HR2 = E_HHHh2(3);

    E_HHHR2 = mnrnd(poissrnd(beta_HH*MHH_HR2(t-T0-TL-TP)),[(1/2),(1/2)]);
    E_HHHR_HH2 = E_HHHR2(1); E_HHHR_HR2 = E_HHHR2(2);

    E_HHhh2 = poissrnd(beta_HH*MHH_hh2(t-T0-TL-TP));
    E_HHhh_Hh2 = E_HHhh2;

    E_HHhR2 = mnrnd(poissrnd(beta_HH*MHH_hR2(t-T0-TL-TP)),[(1/2),(1/2)]);
    E_HHhR_Hh2 = E_HHhR2(1); E_HHhR_HR2 = E_HHhR2(2);

    E_HHRR2 = poissrnd(beta_HH*MHH_RR2(t-T0-TL-TP));
    E_HHRR_HR2 = E_HHRR2;

    % Eggs produced by Hh females:
    
    E_HhHH2 = mnrnd(poissrnd(beta_Hh*MHh_HH2(t-T0-TL-TP)),[((1+e)/2),((1-e-rho)/2),(rho/2)]);
    E_HhHH_HH2 = E_HhHH2(1); E_HhHH_Hh2 = E_HhHH2(2); E_HhHH_HR2 = E_HhHH2(3);

    E_HhHh2 = mnrnd(poissrnd(beta_Hh*MHh_Hh2(t-T0-TL-TP)),[(((1+e)^2)/4),((1+e)*(1-e-rho)/2),((1+e)*rho/2),(((1-e-rho)^2)/4),((1-e-rho)*rho/2),((rho^2)/4)]);
    E_HhHh_HH2 = E_HhHh2(1); E_HhHh_Hh2 = E_HhHh2(2); E_HhHh_HR2 = E_HhHh2(3); E_HhHh_hh2 = E_HhHh2(4); E_HhHh_hR2 = E_HhHh2(5); E_HhHh_RR2 = E_HhHh2(6);

    E_HhHR2 = mnrnd(poissrnd(beta_Hh*MHh_HR2(t-T0-TL-TP)),[((1+e)/4),((1-e-rho)/4),((1+e+rho)/4),((1-e-rho)/4),(rho/4)]);
    E_HhHR_HH2 = E_HhHR2(1); E_HhHR_Hh2 = E_HhHR2(2); E_HhHR_HR2 = E_HhHR2(3); E_HhHR_hR2 = E_HhHR2(4); E_HhHR_RR2 = E_HhHR2(5);

    E_Hhhh2 = mnrnd(poissrnd(beta_Hh*MHh_hh2(t-T0-TL-TP)),[((1+e)/2),((1-e-rho)/2),(rho/2)]);
    E_Hhhh_Hh2 = E_Hhhh2(1); E_Hhhh_hh2 = E_Hhhh2(2); E_Hhhh_hR2 = E_Hhhh2(3);

    E_HhhR2 = mnrnd(poissrnd(beta_Hh*MHh_hR2(t-T0-TL-TP)),[((1+e)/4),((1+e)/4),((1-e-rho)/4),((1-e)/4),(rho/4)]);
    E_HhhR_Hh2 = E_HhhR2(1); E_HhhR_HR2 = E_HhhR2(2); E_HhhR_hh2 = E_HhhR2(3); E_HhhR_hR2 = E_HhhR2(4); E_HhhR_RR2 = E_HhhR2(5);
    
    E_HhRR2 = mnrnd(poissrnd(beta_Hh*MHh_RR2(t-T0-TL-TP)),[((1+e)/2),((1-e-rho)/2),(rho/2)]);
    E_HhRR_HR2 = E_HhRR2(1); E_HhRR_hR2 = E_HhRR2(2); E_HhRR_RR2 = E_HhRR2(3);
    
    % Eggs produced by HR females:
    
    E_HRHH2 = mnrnd(poissrnd(beta_HR*MHR_HH2(t-T0-TL-TP)),[(1/2),(1/2)]);
    E_HRHH_HH2 = E_HRHH2(1); E_HRHH_HR2 = E_HRHH2(2);

    E_HRHh2 = mnrnd(poissrnd(beta_HR*MHR_Hh2(t-T0-TL-TP)),[((1+e)/4),((1-e-rho)/4),((1+e+rho)/4),((1-e-rho)/4),(rho/4)]);
    E_HRHh_HH2 = E_HRHh2(1); E_HRHh_Hh2 = E_HRHh2(2); E_HRHh_HR2 = E_HRHh2(3); E_HRHh_hR2 = E_HRHh2(4); E_HRHh_RR2 = E_HRHh2(5);

    E_HRHR2 = mnrnd(poissrnd(beta_HR*MHR_HR2(t-T0-TL-TP)),[(1/4),(1/2),(1/4)]);
    E_HRHR_HH2 = E_HRHR2(1); E_HRHR_HR2 = E_HRHR2(2); E_HRHR_RR2 = E_HRHR2(3);

    E_HRhh2 = mnrnd(poissrnd(beta_HR*MHR_hh2(t-T0-TL-TP)),[(1/2),(1/2)]);
    E_HRhh_Hh2 = E_HRhh2(1); E_HRhh_hR2 = E_HRhh2(2);

    E_HRhR2 = mnrnd(poissrnd(beta_HR*MHR_hR2(t-T0-TL-TP)),[(1/4),(1/4),(1/4),(1/4)]);
    E_HRhR_Hh2 = E_HRhR2(1); E_HRhR_HR2 = E_HRhR2(2); E_HRhR_hR2 = E_HRhR2(3); E_HRhR_RR2 = E_HRhR2(4);
    
    E_HRRR2 = mnrnd(poissrnd(beta_HR*MHR_RR2(t-T0-TL-TP)),[(1/2),(1/2)]);
    E_HRRR_HR2 = E_HRRR2(1); E_HRRR_RR2 = E_HRRR2(2);
    
    % Eggs produced by hh females:
    
    E_hhHH2 = poissrnd(beta_hh*Mhh_HH2(t-T0-TL-TP));
    E_hhHH_Hh2 = E_hhHH2;

    E_hhHh2 = mnrnd(poissrnd(beta_hh*Mhh_Hh2(t-T0-TL-TP)),[((1+e)/2),((1-e-rho)/2),(rho/2)]);
    E_hhHh_Hh2 = E_hhHh2(1); E_hhHh_hh2 = E_hhHh2(2); E_hhHh_hR2 = E_hhHh2(3);

    E_hhHR2 = mnrnd(poissrnd(beta_hh*Mhh_HR2(t-T0-TL-TP)),[(1/2),(1/2)]);
    E_hhHR_Hh2 = E_hhHR2(1); E_hhHR_hR2 = E_hhHR2(2);

    E_hhhh2 = poissrnd(beta_hh*Mhh_hh2(t-T0-TL-TP));
    E_hhhh_hh2 = E_hhhh2;

    E_hhhR2 = mnrnd(poissrnd(beta_hh*Mhh_hR2(t-T0-TL-TP)),[(1/2),(1/2)]);
    E_hhhR_hh2 = E_hhhR2(1); E_hhhR_hR2 = E_hhhR2(2);

    E_hhRR2 = poissrnd(beta_hh*Mhh_RR2(t-T0-TL-TP));
    E_hhRR_hR2 = E_hhRR2;
    
    % Eggs produced by hR females:
    
    E_hRHH2 = mnrnd(poissrnd(beta_hR*MhR_HH2(t-T0-TL-TP)),[(1/2),(1/2)]);
    E_hRHH_Hh2 = E_hRHH2(1); E_hRHH_HR2 = E_hRHH2(2);

    E_hRHh2 = mnrnd(poissrnd(beta_hR*MhR_Hh2(t-T0-TL-TP)),[((1+e)/4),((1+e)/4),((1-e-rho)/4),((1-e)/4),(rho/4)]);
    E_hRHh_Hh2 = E_hRHh2(1); E_hRHh_HR2 = E_hRHh2(2); E_hRHh_hh2 = E_hRHh2(3); E_hRHh_hR2 = E_hRHh2(4); E_hRHh_RR2 = E_hRHh2(5);

    E_hRHR2 = mnrnd(poissrnd(beta_hR*MhR_HR2(t-T0-TL-TP)),[(1/4),(1/4),(1/4),(1/4)]);
    E_hRHR_Hh2 = E_hRHR2(1); E_hRHR_HR2 = E_hRHR2(2); E_hRHR_hR2 = E_hRHR2(3); E_hRHR_RR2 = E_hRHR2(4);

    E_hRhh2 = mnrnd(poissrnd(beta_hR*MHR_hh2(t-T0-TL-TP)),[(1/2),(1/2)]);
    E_hRhh_hh2 = E_hRhh2(1); E_hRhh_hR2 = E_hRhh2(2);

    E_hRhR2 = mnrnd(poissrnd(beta_hR*MhR_hR2(t-T0-TL-TP)),[(1/4),(1/2),(1/4)]);
    E_hRhR_hh2 = E_hRhR2(1); E_hRhR_hR2 = E_hRhR2(2); E_hRhR_RR2 = E_hRhR2(3);
    
    E_hRRR2 = mnrnd(poissrnd(beta_hR*MhR_RR2(t-T0-TL-TP)),[(1/2),(1/2)]);
    E_hRRR_hR2 = E_hRRR2(1); E_hRRR_RR2 = E_hRRR2(2);
    
    % Eggs produced by RR females:
    
    E_RRHH2 = poissrnd(beta_RR*MRR_HH2(t-T0-TL-TP));
    E_RRHH_HR2 = E_RRHH2;

    E_RRHh2 = mnrnd(poissrnd(beta_RR*MRR_Hh2(t-T0-TL-TP)),[((1+e)/2),((1-e-rho)/2),(rho/2)]);
    E_RRHh_HR2 = E_RRHh2(1); E_RRHh_hR2 = E_RRHh2(2); E_RRHh_RR2 = E_RRHh2(3);

    E_RRHR2 = mnrnd(poissrnd(beta_RR*MRR_HR2(t-T0-TL-TP)),[(1/2),(1/2)]);
    E_RRHR_HR2 = E_RRHR2(1); E_RRHR_RR2 = E_RRHR2(2);

    E_RRhh2 = poissrnd(beta_RR*MRR_hh2(t-T0-TL-TP));
    E_RRhh_hR2 = E_RRhh2;

    E_RRhR2 = mnrnd(poissrnd(beta_RR*MRR_hR2(t-T0-TL-TP)),[(1/2),(1/2)]);
    E_RRhR_hR2 = E_RRhR2(1); E_RRhR_RR2 = E_RRhR2(2);

    E_RRRR2 = poissrnd(beta_RR*MRR_RR2(t-T0-TL-TP));
    E_RRRR_RR2 = E_RRRR2;
    
    % Adult male & female newborns:
    
    MHH_Births2 = mnrnd(binornd((E_HHHH_HH2 + E_HHHh_HH2 + E_HHHR_HH2 + E_HhHH_HH2 + E_HhHh_HH2 + E_HhHR_HH2 + E_HRHH_HH2 + E_HRHh_HH2 + E_HRHR_HH2),theta0*thetaLA*thetaP*(1-mu_m_HH)),[(1/2),(1/2)]);
    MHh_Births2 = mnrnd(binornd((E_HHHh_Hh2 + E_HHhh_Hh2 + E_HHhR_Hh2 + E_HhHH_Hh2 + E_HhHh_Hh2 + E_HhHR_Hh2 + E_Hhhh_Hh2 + E_HhhR_Hh2 + E_HRHh_Hh2 + E_HRhh_Hh2 + E_HRhR_Hh2 + E_hhHH_Hh2 + E_hhHh_Hh2 + E_hhHR_Hh2 + E_hRHH_Hh2 + E_hRHh_Hh2 + E_hRHR_Hh2),theta0*thetaLA*thetaP*(1-mu_m_Hh)),[(1/2),(1/2)]);
    MHR_Births2 = mnrnd(binornd((E_HHHh_HR2 + E_HHHR_HR2 + E_HHhR_HR2 + E_HHRR_HR2 + E_HhHH_HR2 + E_HhHh_HR2 + E_HhHR_HR2 + E_HhhR_HR2 + E_HhRR_HR2 + E_HRHH_HR2 + E_HRHh_HR2 + E_HRHR_HR2 + E_HRhR_HR2 + E_HRRR_HR2 + E_hRHH_HR2 + E_hRHh_HR2 + E_hRHR_HR2 + E_RRHH_HR2 + E_RRHh_HR2 + E_RRHR_HR2),theta0*thetaLA*thetaP*(1-mu_m_HR)),[(1/2),(1/2)]);
    Mhh_Births2 = mnrnd(binornd((E_HhHh_hh2 + E_Hhhh_hh2 + E_HhhR_hh2 + E_hhHh_hh2 + E_hhhh_hh2 + E_hhhR_hh2 + E_hRHh_hh2 + E_hRhh_hh2 + E_hRhR_hh2),theta0*thetaLA*thetaP*(1-mu_m_hh)),[(1/2),(1/2)]);
    MhR_Births2 = mnrnd(binornd((E_HhHh_hR2 + E_HhHR_hR2 + E_Hhhh_hR2 + E_HhhR_hR2 + E_HhRR_hR2 + E_HRHh_hR2 + E_HRhh_hR2 + E_HRhR_hR2 + E_hhHh_hR2 + E_hhHR_hR2 + E_hhhR_hR2 + E_hhRR_hR2 + E_hRHh_hR2 + E_hRHR_hR2 + E_hRhh_hR2 + E_hRhR_hR2 + E_hRRR_hR2 + E_RRHh_hR2 + E_RRhh_hR2 + E_RRhR_hR2),theta0*thetaLA*thetaP*(1-mu_m_hR)),[(1/2),(1/2)]);
    MRR_Births2 = mnrnd(binornd((E_HhHh_RR2 + E_HhHR_RR2 + E_HhhR_RR2 + E_HhRR_RR2 + E_HRHh_RR2 + E_HRHR_RR2 + E_HRhR_RR2 + E_HRRR_RR2 + E_hRHh_RR2 + E_hRHR_RR2 + E_hRhR_RR2 + E_hRRR_RR2 + E_RRHh_RR2 + E_RRHR_RR2 + E_RRhR_RR2 + E_RRRR_RR2),theta0*thetaLA*thetaP*(1-mu_m_RR)),[(1/2),(1/2)]);

    % Adult male genotypes:
    
    MHH2(t) = max(0,binornd(MHH2(t-1),(1-mu_m_HH)) + MHH_Births2(1));
    MHh2(t) = max(0,binornd(MHh2(t-1),(1-mu_m_Hh)) + MHh_Births2(1));
    MHR2(t) = max(0,binornd(MHR2(t-1),(1-mu_m_HR)) + MHR_Births2(1));
    Mhh2(t) = max(0,binornd(Mhh2(t-1),(1-mu_m_hh)) + Mhh_Births2(1));
    MhR2(t) = max(0,binornd(MhR2(t-1),(1-mu_m_hR)) + MhR_Births2(1));
    MRR2(t) = max(0,binornd(MRR2(t-1),(1-mu_m_RR)) + MRR_Births2(1));
    
    % Adult female genotypes:
    
    X_HH2 = mnrnd(MHH_Births2(2) , [ (MHH2(t-1)/(MHH2(t-1)+MHh2(t-1)+MHR2(t-1)+Mhh2(t-1)+MhR2(t-1)+MRR2(t-1))) , (MHh2(t-1)/(MHH2(t-1)+MHh2(t-1)+MHR2(t-1)+Mhh2(t-1)+MhR2(t-1)+MRR2(t-1))) , (MHR2(t-1)/(MHH2(t-1)+MHh2(t-1)+MHR2(t-1)+Mhh2(t-1)+MhR2(t-1)+MRR2(t-1))) , (Mhh2(t-1)/(MHH2(t-1)+MHh2(t-1)+MHR2(t-1)+Mhh2(t-1)+MhR2(t-1)+MRR2(t-1))) , (MhR2(t-1)/(MHH2(t-1)+MHh2(t-1)+MHR2(t-1)+Mhh2(t-1)+MhR2(t-1)+MRR2(t-1))) , (MRR2(t-1)/(MHH2(t-1)+MHh2(t-1)+MHR2(t-1)+Mhh2(t-1)+MhR2(t-1)+MRR2(t-1))) ]);
    X_Hh2 = mnrnd(MHh_Births2(2) , [ (MHH2(t-1)/(MHH2(t-1)+MHh2(t-1)+MHR2(t-1)+Mhh2(t-1)+MhR2(t-1)+MRR2(t-1))) , (MHh2(t-1)/(MHH2(t-1)+MHh2(t-1)+MHR2(t-1)+Mhh2(t-1)+MhR2(t-1)+MRR2(t-1))) , (MHR2(t-1)/(MHH2(t-1)+MHh2(t-1)+MHR2(t-1)+Mhh2(t-1)+MhR2(t-1)+MRR2(t-1))) , (Mhh2(t-1)/(MHH2(t-1)+MHh2(t-1)+MHR2(t-1)+Mhh2(t-1)+MhR2(t-1)+MRR2(t-1))) , (MhR2(t-1)/(MHH2(t-1)+MHh2(t-1)+MHR2(t-1)+Mhh2(t-1)+MhR2(t-1)+MRR2(t-1))) , (MRR2(t-1)/(MHH2(t-1)+MHh2(t-1)+MHR2(t-1)+Mhh2(t-1)+MhR2(t-1)+MRR2(t-1))) ]);
    X_HR2 = mnrnd(MHR_Births2(2) , [ (MHH2(t-1)/(MHH2(t-1)+MHh2(t-1)+MHR2(t-1)+Mhh2(t-1)+MhR2(t-1)+MRR2(t-1))) , (MHh2(t-1)/(MHH2(t-1)+MHh2(t-1)+MHR2(t-1)+Mhh2(t-1)+MhR2(t-1)+MRR2(t-1))) , (MHR2(t-1)/(MHH2(t-1)+MHh2(t-1)+MHR2(t-1)+Mhh2(t-1)+MhR2(t-1)+MRR2(t-1))) , (Mhh2(t-1)/(MHH2(t-1)+MHh2(t-1)+MHR2(t-1)+Mhh2(t-1)+MhR2(t-1)+MRR2(t-1))) , (MhR2(t-1)/(MHH2(t-1)+MHh2(t-1)+MHR2(t-1)+Mhh2(t-1)+MhR2(t-1)+MRR2(t-1))) , (MRR2(t-1)/(MHH2(t-1)+MHh2(t-1)+MHR2(t-1)+Mhh2(t-1)+MhR2(t-1)+MRR2(t-1))) ]);
    X_hh2 = mnrnd(Mhh_Births2(2) , [ (MHH2(t-1)/(MHH2(t-1)+MHh2(t-1)+MHR2(t-1)+Mhh2(t-1)+MhR2(t-1)+MRR2(t-1))) , (MHh2(t-1)/(MHH2(t-1)+MHh2(t-1)+MHR2(t-1)+Mhh2(t-1)+MhR2(t-1)+MRR2(t-1))) , (MHR2(t-1)/(MHH2(t-1)+MHh2(t-1)+MHR2(t-1)+Mhh2(t-1)+MhR2(t-1)+MRR2(t-1))) , (Mhh2(t-1)/(MHH2(t-1)+MHh2(t-1)+MHR2(t-1)+Mhh2(t-1)+MhR2(t-1)+MRR2(t-1))) , (MhR2(t-1)/(MHH2(t-1)+MHh2(t-1)+MHR2(t-1)+Mhh2(t-1)+MhR2(t-1)+MRR2(t-1))) , (MRR2(t-1)/(MHH2(t-1)+MHh2(t-1)+MHR2(t-1)+Mhh2(t-1)+MhR2(t-1)+MRR2(t-1))) ]);
    X_hR2 = mnrnd(MhR_Births2(2) , [ (MHH2(t-1)/(MHH2(t-1)+MHh2(t-1)+MHR2(t-1)+Mhh2(t-1)+MhR2(t-1)+MRR2(t-1))) , (MHh2(t-1)/(MHH2(t-1)+MHh2(t-1)+MHR2(t-1)+Mhh2(t-1)+MhR2(t-1)+MRR2(t-1))) , (MHR2(t-1)/(MHH2(t-1)+MHh2(t-1)+MHR2(t-1)+Mhh2(t-1)+MhR2(t-1)+MRR2(t-1))) , (Mhh2(t-1)/(MHH2(t-1)+MHh2(t-1)+MHR2(t-1)+Mhh2(t-1)+MhR2(t-1)+MRR2(t-1))) , (MhR2(t-1)/(MHH2(t-1)+MHh2(t-1)+MHR2(t-1)+Mhh2(t-1)+MhR2(t-1)+MRR2(t-1))) , (MRR2(t-1)/(MHH2(t-1)+MHh2(t-1)+MHR2(t-1)+Mhh2(t-1)+MhR2(t-1)+MRR2(t-1))) ]);    
    X_RR2 = mnrnd(MRR_Births2(2) , [ (MHH2(t-1)/(MHH2(t-1)+MHh2(t-1)+MHR2(t-1)+Mhh2(t-1)+MhR2(t-1)+MRR2(t-1))) , (MHh2(t-1)/(MHH2(t-1)+MHh2(t-1)+MHR2(t-1)+Mhh2(t-1)+MhR2(t-1)+MRR2(t-1))) , (MHR2(t-1)/(MHH2(t-1)+MHh2(t-1)+MHR2(t-1)+Mhh2(t-1)+MhR2(t-1)+MRR2(t-1))) , (Mhh2(t-1)/(MHH2(t-1)+MHh2(t-1)+MHR2(t-1)+Mhh2(t-1)+MhR2(t-1)+MRR2(t-1))) , (MhR2(t-1)/(MHH2(t-1)+MHh2(t-1)+MHR2(t-1)+Mhh2(t-1)+MhR2(t-1)+MRR2(t-1))) , (MRR2(t-1)/(MHH2(t-1)+MHh2(t-1)+MHR2(t-1)+Mhh2(t-1)+MhR2(t-1)+MRR2(t-1))) ]);

    MHH2_Temp = [binornd(MHH_HH2(t-1),(1-mu_m_HH)) , binornd(MHH_Hh2(t-1),(1-mu_m_HH)) , binornd(MHH_HR2(t-1),(1-mu_m_HH)) , binornd(MHH_hh2(t-1),(1-mu_m_HH)) , binornd(MHH_hR2(t-1),(1-mu_m_HH)) , binornd(MHH_RR2(t-1),(1-mu_m_HH))] + X_HH2;
    MHH_HH2(t) = MHH2_Temp(1); MHH_Hh2(t) = MHH2_Temp(2); MHH_HR2(t) = MHH2_Temp(3); MHH_hh2(t) = MHH2_Temp(4); MHH_hR2(t) = MHH2_Temp(5); MHH_RR2(t) = MHH2_Temp(6);

    MHh2_Temp = [binornd(MHh_HH2(t-1),(1-mu_m_Hh)) , binornd(MHh_Hh2(t-1),(1-mu_m_Hh)) , binornd(MHh_HR2(t-1),(1-mu_m_Hh)) , binornd(MHh_hh2(t-1),(1-mu_m_Hh)) , binornd(MHh_hR2(t-1),(1-mu_m_Hh)) , binornd(MHh_RR2(t-1),(1-mu_m_Hh))] + X_Hh2;
    MHh_HH2(t) = MHh2_Temp(1); MHh_Hh2(t) = MHh2_Temp(2); MHh_HR2(t) = MHh2_Temp(3); MHh_hh2(t) = MHh2_Temp(4); MHh_hR2(t) = MHh2_Temp(5); MHh_RR2(t) = MHh2_Temp(6);

    MHR2_Temp = [binornd(MHR_HH2(t-1),(1-mu_m_HR)) , binornd(MHR_Hh2(t-1),(1-mu_m_HR)) , binornd(MHR_HR2(t-1),(1-mu_m_HR)) , binornd(MHR_hh2(t-1),(1-mu_m_HR)) , binornd(MHR_hR2(t-1),(1-mu_m_HR)) , binornd(MHR_RR2(t-1),(1-mu_m_HR))] + X_HR2;
    MHR_HH2(t) = MHR2_Temp(1); MHR_Hh2(t) = MHR2_Temp(2); MHR_HR2(t) = MHR2_Temp(3); MHR_hh2(t) = MHR2_Temp(4); MHR_hR2(t) = MHR2_Temp(5); MHR_RR2(t) = MHR2_Temp(6);

    Mhh2_Temp = [binornd(Mhh_HH2(t-1),(1-mu_m_hh)) , binornd(Mhh_Hh2(t-1),(1-mu_m_hh)) , binornd(Mhh_HR2(t-1),(1-mu_m_hh)) , binornd(Mhh_hh2(t-1),(1-mu_m_hh)) , binornd(Mhh_hR2(t-1),(1-mu_m_hh)) , binornd(Mhh_RR2(t-1),(1-mu_m_hh))] + X_hh2;
    Mhh_HH2(t) = Mhh2_Temp(1); Mhh_Hh2(t) = Mhh2_Temp(2); Mhh_HR2(t) = Mhh2_Temp(3); Mhh_hh2(t) = Mhh2_Temp(4); Mhh_hR2(t) = Mhh2_Temp(5); Mhh_RR2(t) = Mhh2_Temp(6);

    MhR2_Temp = [binornd(MhR_HH2(t-1),(1-mu_m_hR)) , binornd(MhR_Hh2(t-1),(1-mu_m_hR)) , binornd(MhR_HR2(t-1),(1-mu_m_hR)) , binornd(MhR_hh2(t-1),(1-mu_m_hR)) , binornd(MhR_hR2(t-1),(1-mu_m_hR)) , binornd(MhR_RR2(t-1),(1-mu_m_hR))] + X_hR2;
    MhR_HH2(t) = MhR2_Temp(1); MhR_Hh2(t) = MhR2_Temp(2); MhR_HR2(t) = MhR2_Temp(3); MhR_hh2(t) = MhR2_Temp(4); MhR_hR2(t) = MhR2_Temp(5); MhR_RR2(t) = MhR2_Temp(6);

    MRR2_Temp = [binornd(MRR_HH2(t-1),(1-mu_m_RR)) , binornd(MRR_Hh2(t-1),(1-mu_m_RR)) , binornd(MRR_HR2(t-1),(1-mu_m_RR)) , binornd(MRR_hh2(t-1),(1-mu_m_RR)) , binornd(MRR_hR2(t-1),(1-mu_m_RR)) , binornd(MRR_RR2(t-1),(1-mu_m_RR))] + X_RR2;
    MRR_HH2(t) = MRR2_Temp(1); MRR_Hh2(t) = MRR2_Temp(2); MRR_HR2(t) = MRR2_Temp(3); MRR_hh2(t) = MRR2_Temp(4); MRR_hR2(t) = MRR2_Temp(5); MRR_RR2(t) = MRR2_Temp(6);
    
    MHH_HH2(t) = max(0,MHH_HH2(t)); MHH_Hh2(t) = max(0,MHH_Hh2(t)); MHH_HR2(t) = max(0,MHH_HR2(t)); MHH_hh2(t) = max(0,MHH_hh2(t)); MHH_hR2(t) = max(0,MHH_hR2(t)); MHH_RR2(t) = max(0,MHH_RR2(t));
    MHh_HH2(t) = max(0,MHh_HH2(t)); MHh_Hh2(t) = max(0,MHh_Hh2(t)); MHh_HR2(t) = max(0,MHh_HR2(t)); MHh_hh2(t) = max(0,MHh_hh2(t)); MHh_hR2(t) = max(0,MHh_hR2(t)); MHh_RR2(t) = max(0,MHh_RR2(t));
    MHR_HH2(t) = max(0,MHR_HH2(t)); MHR_Hh2(t) = max(0,MHR_Hh2(t)); MHR_HR2(t) = max(0,MHR_HR2(t)); MHR_hh2(t) = max(0,MHR_hh2(t)); MHR_hR2(t) = max(0,MHR_hR2(t)); MHR_RR2(t) = max(0,MHR_RR2(t));
    Mhh_HH2(t) = max(0,Mhh_HH2(t)); Mhh_Hh2(t) = max(0,Mhh_Hh2(t)); Mhh_HR2(t) = max(0,Mhh_HR2(t)); Mhh_hh2(t) = max(0,Mhh_hh2(t)); Mhh_hR2(t) = max(0,Mhh_hR2(t)); Mhh_RR2(t) = max(0,Mhh_RR2(t));
    MhR_HH2(t) = max(0,MhR_HH2(t)); MhR_Hh2(t) = max(0,MhR_Hh2(t)); MhR_HR2(t) = max(0,MhR_HR2(t)); MhR_hh2(t) = max(0,MhR_hh2(t)); MhR_hR2(t) = max(0,MhR_hR2(t)); MhR_RR2(t) = max(0,MhR_RR2(t));
    MRR_HH2(t) = max(0,MRR_HH2(t)); MRR_Hh2(t) = max(0,MRR_Hh2(t)); MRR_HR2(t) = max(0,MRR_HR2(t)); MRR_hh2(t) = max(0,MRR_hh2(t)); MRR_hR2(t) = max(0,MRR_hR2(t)); MRR_RR2(t) = max(0,MRR_RR2(t));
    
    if (t==51)
        MHH2(t) = MHH2(t) + M_eq;
        % MHH_HH2(t) = MHH_HH(t) + M_eq/2;
    end
        
    % Population 3:
    thetaLA = thetaL;
    for i=1:TL
        thetaLA = thetaLA*((alpha/(alpha+L3(t-i)))^(1/TL));
    end
    
    L3(t) = max(0,binornd(L3(t-1),(1-mu_l)*((alpha/(alpha+L3(t-1)))^(1/TL))) + binornd(poissrnd( beta_HH*MHH_HH3(t-T0) + beta_HH*MHH_Hh3(t-T0) + beta_HH*MHH_HR3(t-T0) + beta_HH*MHH_hh3(t-T0) + beta_HH*MHH_hR3(t-T0) + beta_HH*MHH_RR3(t-T0)  +  beta_Hh*MHh_HH3(t-T0) + beta_Hh*MHh_Hh3(t-T0) + beta_Hh*MHh_HR3(t-T0) + beta_Hh*MHh_hh3(t-T0) + beta_Hh*MHh_hR3(t-T0) + beta_Hh*MHh_RR3(t-T0)  +  beta_HR*MHR_HH3(t-T0) + beta_HR*MHR_Hh3(t-T0) + beta_HR*MHR_HR3(t-T0) + beta_HR*MHR_hh3(t-T0) + beta_HR*MHR_hR3(t-T0) + beta_HR*MHR_RR3(t-T0)  +  beta_hh*Mhh_HH3(t-T0) + beta_hh*Mhh_Hh3(t-T0) + beta_hh*Mhh_HR3(t-T0) + beta_hh*Mhh_hh3(t-T0) + beta_hh*Mhh_hR3(t-T0) + beta_hh*Mhh_RR3(t-T0)  +  beta_hR*MhR_HH3(t-T0) + beta_hR*MhR_Hh3(t-T0) + beta_hR*MhR_HR3(t-T0) + beta_hR*MhR_hh3(t-T0) + beta_hR*MhR_hR3(t-T0) + beta_hR*MhR_RR3(t-T0)  +  beta_RR*MRR_HH3(t-T0) + beta_RR*MRR_Hh3(t-T0) + beta_RR*MRR_HR3(t-T0) + beta_RR*MRR_hh3(t-T0) + beta_RR*MRR_hR3(t-T0) + beta_RR*MRR_RR3(t-T0)  ),theta0) - binornd(poissrnd(beta_HH*MHH_HH3(t-T0-TL) + beta_HH*MHH_Hh3(t-T0-TL) + beta_HH*MHH_HR3(t-T0-TL) + beta_HH*MHH_hh3(t-T0-TL) + beta_HH*MHH_hR3(t-T0-TL) + beta_HH*MHH_RR3(t-T0-TL)  +  beta_Hh*MHh_HH3(t-T0-TL) + beta_Hh*MHh_Hh3(t-T0-TL) + beta_Hh*MHh_HR3(t-T0-TL) + beta_Hh*MHh_hh3(t-T0-TL) + beta_Hh*MHh_hR3(t-T0-TL) + beta_Hh*MHh_RR3(t-T0-TL)  +  beta_HR*MHR_HH3(t-T0-TL) + beta_HR*MHR_Hh3(t-T0-TL) + beta_HR*MHR_HR3(t-T0-TL) + beta_HR*MHR_hh3(t-T0-TL) + beta_HR*MHR_hR3(t-T0-TL) + beta_HR*MHR_RR3(t-T0-TL)  +  beta_hh*Mhh_HH3(t-T0-TL) + beta_hh*Mhh_Hh3(t-T0-TL) + beta_hh*Mhh_HR3(t-T0-TL) + beta_hh*Mhh_hh3(t-T0-TL) + beta_hh*Mhh_hR3(t-T0-TL) + beta_hh*Mhh_RR3(t-T0-TL)  +  beta_hR*MhR_HH3(t-T0-TL) + beta_hR*MhR_Hh3(t-T0-TL) + beta_hR*MhR_HR3(t-T0-TL) + beta_hR*MhR_hh3(t-T0-TL) + beta_hR*MhR_hR3(t-T0-TL) + beta_hR*MhR_RR3(t-T0-TL)  +  beta_RR*MRR_HH3(t-T0-TL) + beta_RR*MRR_Hh3(t-T0-TL) + beta_RR*MRR_HR3(t-T0-TL) + beta_RR*MRR_hh3(t-T0-TL) + beta_RR*MRR_hR3(t-T0-TL) + beta_RR*MRR_RR3(t-T0-TL) ),theta0*thetaLA));
    
    thetaLA = thetaL;
    for i=1:TL
        thetaLA = thetaLA*((alpha/(alpha+L3(t-i-TP)))^(1/TL));
    end
    
    % Eggs produced by HH females:
    
    E_HHHH3 = poissrnd(beta_HH*MHH_HH3(t-T0-TL-TP));
    E_HHHH_HH3 = E_HHHH3;

    E_HHHh3 = mnrnd(poissrnd(beta_HH*MHH_Hh3(t-T0-TL-TP)),[((1+e)/2),((1-e-rho)/2),(rho/2)]);
    E_HHHh_HH3 = E_HHHh3(1); E_HHHh_Hh3 = E_HHHh3(2); E_HHHh_HR3 = E_HHHh3(3);

    E_HHHR3 = mnrnd(poissrnd(beta_HH*MHH_HR3(t-T0-TL-TP)),[(1/2),(1/2)]);
    E_HHHR_HH3 = E_HHHR3(1); E_HHHR_HR3 = E_HHHR3(2);

    E_HHhh3 = poissrnd(beta_HH*MHH_hh3(t-T0-TL-TP));
    E_HHhh_Hh3 = E_HHhh3;

    E_HHhR3 = mnrnd(poissrnd(beta_HH*MHH_hR3(t-T0-TL-TP)),[(1/2),(1/2)]);
    E_HHhR_Hh3 = E_HHhR3(1); E_HHhR_HR3 = E_HHhR3(2);

    E_HHRR3 = poissrnd(beta_HH*MHH_RR3(t-T0-TL-TP));
    E_HHRR_HR3 = E_HHRR3;

    % Eggs produced by Hh females:
    
    E_HhHH3 = mnrnd(poissrnd(beta_Hh*MHh_HH3(t-T0-TL-TP)),[((1+e)/2),((1-e-rho)/2),(rho/2)]);
    E_HhHH_HH3 = E_HhHH3(1); E_HhHH_Hh3 = E_HhHH3(2); E_HhHH_HR3 = E_HhHH3(3);

    E_HhHh3 = mnrnd(poissrnd(beta_Hh*MHh_Hh3(t-T0-TL-TP)),[(((1+e)^2)/4),((1+e)*(1-e-rho)/2),((1+e)*rho/2),(((1-e-rho)^2)/4),((1-e-rho)*rho/2),((rho^2)/4)]);
    E_HhHh_HH3 = E_HhHh3(1); E_HhHh_Hh3 = E_HhHh3(2); E_HhHh_HR3 = E_HhHh3(3); E_HhHh_hh3 = E_HhHh3(4); E_HhHh_hR3 = E_HhHh3(5); E_HhHh_RR3 = E_HhHh3(6);

    E_HhHR3 = mnrnd(poissrnd(beta_Hh*MHh_HR3(t-T0-TL-TP)),[((1+e)/4),((1-e-rho)/4),((1+e+rho)/4),((1-e-rho)/4),(rho/4)]);
    E_HhHR_HH3 = E_HhHR3(1); E_HhHR_Hh3 = E_HhHR3(2); E_HhHR_HR3 = E_HhHR3(3); E_HhHR_hR3 = E_HhHR3(4); E_HhHR_RR3 = E_HhHR3(5);

    E_Hhhh3 = mnrnd(poissrnd(beta_Hh*MHh_hh3(t-T0-TL-TP)),[((1+e)/2),((1-e-rho)/2),(rho/2)]);
    E_Hhhh_Hh3 = E_Hhhh3(1); E_Hhhh_hh3 = E_Hhhh3(2); E_Hhhh_hR3 = E_Hhhh3(3);

    E_HhhR3 = mnrnd(poissrnd(beta_Hh*MHh_hR3(t-T0-TL-TP)),[((1+e)/4),((1+e)/4),((1-e-rho)/4),((1-e)/4),(rho/4)]);
    E_HhhR_Hh3 = E_HhhR3(1); E_HhhR_HR3 = E_HhhR3(2); E_HhhR_hh3 = E_HhhR3(3); E_HhhR_hR3 = E_HhhR3(4); E_HhhR_RR3 = E_HhhR3(5);
    
    E_HhRR3 = mnrnd(poissrnd(beta_Hh*MHh_RR3(t-T0-TL-TP)),[((1+e)/2),((1-e-rho)/2),(rho/2)]);
    E_HhRR_HR3 = E_HhRR3(1); E_HhRR_hR3 = E_HhRR3(2); E_HhRR_RR3 = E_HhRR3(3);
    
    % Eggs produced by HR females:
    
    E_HRHH3 = mnrnd(poissrnd(beta_HR*MHR_HH3(t-T0-TL-TP)),[(1/2),(1/2)]);
    E_HRHH_HH3 = E_HRHH3(1); E_HRHH_HR3 = E_HRHH3(2);

    E_HRHh3 = mnrnd(poissrnd(beta_HR*MHR_Hh3(t-T0-TL-TP)),[((1+e)/4),((1-e-rho)/4),((1+e+rho)/4),((1-e-rho)/4),(rho/4)]);
    E_HRHh_HH3 = E_HRHh3(1); E_HRHh_Hh3 = E_HRHh3(2); E_HRHh_HR3 = E_HRHh3(3); E_HRHh_hR3 = E_HRHh3(4); E_HRHh_RR3 = E_HRHh3(5);

    E_HRHR3 = mnrnd(poissrnd(beta_HR*MHR_HR3(t-T0-TL-TP)),[(1/4),(1/2),(1/4)]);
    E_HRHR_HH3 = E_HRHR3(1); E_HRHR_HR3 = E_HRHR3(2); E_HRHR_RR3 = E_HRHR3(3);

    E_HRhh3 = mnrnd(poissrnd(beta_HR*MHR_hh3(t-T0-TL-TP)),[(1/2),(1/2)]);
    E_HRhh_Hh3 = E_HRhh3(1); E_HRhh_hR3 = E_HRhh3(2);

    E_HRhR3 = mnrnd(poissrnd(beta_HR*MHR_hR3(t-T0-TL-TP)),[(1/4),(1/4),(1/4),(1/4)]);
    E_HRhR_Hh3 = E_HRhR3(1); E_HRhR_HR3 = E_HRhR3(2); E_HRhR_hR3 = E_HRhR3(3); E_HRhR_RR3 = E_HRhR3(4);
    
    E_HRRR3 = mnrnd(poissrnd(beta_HR*MHR_RR3(t-T0-TL-TP)),[(1/2),(1/2)]);
    E_HRRR_HR3 = E_HRRR3(1); E_HRRR_RR3 = E_HRRR3(2);
    
    % Eggs produced by hh females:
    
    E_hhHH3 = poissrnd(beta_hh*Mhh_HH3(t-T0-TL-TP));
    E_hhHH_Hh3 = E_hhHH3;

    E_hhHh3 = mnrnd(poissrnd(beta_hh*Mhh_Hh3(t-T0-TL-TP)),[((1+e)/2),((1-e-rho)/2),(rho/2)]);
    E_hhHh_Hh3 = E_hhHh3(1); E_hhHh_hh3 = E_hhHh3(2); E_hhHh_hR3 = E_hhHh3(3);

    E_hhHR3 = mnrnd(poissrnd(beta_hh*Mhh_HR3(t-T0-TL-TP)),[(1/2),(1/2)]);
    E_hhHR_Hh3 = E_hhHR3(1); E_hhHR_hR3 = E_hhHR3(2);

    E_hhhh3 = poissrnd(beta_hh*Mhh_hh3(t-T0-TL-TP));
    E_hhhh_hh3 = E_hhhh3;

    E_hhhR3 = mnrnd(poissrnd(beta_hh*Mhh_hR3(t-T0-TL-TP)),[(1/2),(1/2)]);
    E_hhhR_hh3 = E_hhhR3(1); E_hhhR_hR3 = E_hhhR3(2);

    E_hhRR3 = poissrnd(beta_hh*Mhh_RR3(t-T0-TL-TP));
    E_hhRR_hR3 = E_hhRR3;
    
    % Eggs produced by hR females:
    
    E_hRHH3 = mnrnd(poissrnd(beta_hR*MhR_HH3(t-T0-TL-TP)),[(1/2),(1/2)]);
    E_hRHH_Hh3 = E_hRHH3(1); E_hRHH_HR3 = E_hRHH3(2);

    E_hRHh3 = mnrnd(poissrnd(beta_hR*MhR_Hh3(t-T0-TL-TP)),[((1+e)/4),((1+e)/4),((1-e-rho)/4),((1-e)/4),(rho/4)]);
    E_hRHh_Hh3 = E_hRHh3(1); E_hRHh_HR3 = E_hRHh3(2); E_hRHh_hh3 = E_hRHh3(3); E_hRHh_hR3 = E_hRHh3(4); E_hRHh_RR3 = E_hRHh3(5);

    E_hRHR3 = mnrnd(poissrnd(beta_hR*MhR_HR3(t-T0-TL-TP)),[(1/4),(1/4),(1/4),(1/4)]);
    E_hRHR_Hh3 = E_hRHR3(1); E_hRHR_HR3 = E_hRHR3(2); E_hRHR_hR3 = E_hRHR3(3); E_hRHR_RR3 = E_hRHR3(4);

    E_hRhh3 = mnrnd(poissrnd(beta_hR*MHR_hh3(t-T0-TL-TP)),[(1/2),(1/2)]);
    E_hRhh_hh3 = E_hRhh3(1); E_hRhh_hR3 = E_hRhh3(2);

    E_hRhR3 = mnrnd(poissrnd(beta_hR*MhR_hR3(t-T0-TL-TP)),[(1/4),(1/2),(1/4)]);
    E_hRhR_hh3 = E_hRhR3(1); E_hRhR_hR3 = E_hRhR3(2); E_hRhR_RR3 = E_hRhR3(3);
    
    E_hRRR3 = mnrnd(poissrnd(beta_hR*MhR_RR3(t-T0-TL-TP)),[(1/2),(1/2)]);
    E_hRRR_hR3 = E_hRRR3(1); E_hRRR_RR3 = E_hRRR3(2);
    
    % Eggs produced by RR females:
    
    E_RRHH3 = poissrnd(beta_RR*MRR_HH3(t-T0-TL-TP));
    E_RRHH_HR3 = E_RRHH3;

    E_RRHh3 = mnrnd(poissrnd(beta_RR*MRR_Hh3(t-T0-TL-TP)),[((1+e)/2),((1-e-rho)/2),(rho/2)]);
    E_RRHh_HR3 = E_RRHh3(1); E_RRHh_hR3 = E_RRHh3(2); E_RRHh_RR3 = E_RRHh3(3);

    E_RRHR3 = mnrnd(poissrnd(beta_RR*MRR_HR3(t-T0-TL-TP)),[(1/2),(1/2)]);
    E_RRHR_HR3 = E_RRHR3(1); E_RRHR_RR3 = E_RRHR3(2);

    E_RRhh3 = poissrnd(beta_RR*MRR_hh3(t-T0-TL-TP));
    E_RRhh_hR3 = E_RRhh3;

    E_RRhR3 = mnrnd(poissrnd(beta_RR*MRR_hR3(t-T0-TL-TP)),[(1/2),(1/2)]);
    E_RRhR_hR3 = E_RRhR3(1); E_RRhR_RR3 = E_RRhR3(2);

    E_RRRR3 = poissrnd(beta_RR*MRR_RR3(t-T0-TL-TP));
    E_RRRR_RR3 = E_RRRR3;
    
    % Adult male & female newborns:
    
    MHH_Births3 = mnrnd(binornd((E_HHHH_HH3 + E_HHHh_HH3 + E_HHHR_HH3 + E_HhHH_HH3 + E_HhHh_HH3 + E_HhHR_HH3 + E_HRHH_HH3 + E_HRHh_HH3 + E_HRHR_HH3),theta0*thetaLA*thetaP*(1-mu_m_HH)),[(1/2),(1/2)]);
    MHh_Births3 = mnrnd(binornd((E_HHHh_Hh3 + E_HHhh_Hh3 + E_HHhR_Hh3 + E_HhHH_Hh3 + E_HhHh_Hh3 + E_HhHR_Hh3 + E_Hhhh_Hh3 + E_HhhR_Hh3 + E_HRHh_Hh3 + E_HRhh_Hh3 + E_HRhR_Hh3 + E_hhHH_Hh3 + E_hhHh_Hh3 + E_hhHR_Hh3 + E_hRHH_Hh3 + E_hRHh_Hh3 + E_hRHR_Hh3),theta0*thetaLA*thetaP*(1-mu_m_Hh)),[(1/2),(1/2)]);
    MHR_Births3 = mnrnd(binornd((E_HHHh_HR3 + E_HHHR_HR3 + E_HHhR_HR3 + E_HHRR_HR3 + E_HhHH_HR3 + E_HhHh_HR3 + E_HhHR_HR3 + E_HhhR_HR3 + E_HhRR_HR3 + E_HRHH_HR3 + E_HRHh_HR3 + E_HRHR_HR3 + E_HRhR_HR3 + E_HRRR_HR3 + E_hRHH_HR3 + E_hRHh_HR3 + E_hRHR_HR3 + E_RRHH_HR3 + E_RRHh_HR3 + E_RRHR_HR3),theta0*thetaLA*thetaP*(1-mu_m_HR)),[(1/2),(1/2)]);
    Mhh_Births3 = mnrnd(binornd((E_HhHh_hh3 + E_Hhhh_hh3 + E_HhhR_hh3 + E_hhHh_hh3 + E_hhhh_hh3 + E_hhhR_hh3 + E_hRHh_hh3 + E_hRhh_hh3 + E_hRhR_hh3),theta0*thetaLA*thetaP*(1-mu_m_hh)),[(1/2),(1/2)]);
    MhR_Births3 = mnrnd(binornd((E_HhHh_hR3 + E_HhHR_hR3 + E_Hhhh_hR3 + E_HhhR_hR3 + E_HhRR_hR3 + E_HRHh_hR3 + E_HRhh_hR3 + E_HRhR_hR3 + E_hhHh_hR3 + E_hhHR_hR3 + E_hhhR_hR3 + E_hhRR_hR3 + E_hRHh_hR3 + E_hRHR_hR3 + E_hRhh_hR3 + E_hRhR_hR3 + E_hRRR_hR3 + E_RRHh_hR3 + E_RRhh_hR3 + E_RRhR_hR3),theta0*thetaLA*thetaP*(1-mu_m_hR)),[(1/2),(1/2)]);
    MRR_Births3 = mnrnd(binornd((E_HhHh_RR3 + E_HhHR_RR3 + E_HhhR_RR3 + E_HhRR_RR3 + E_HRHh_RR3 + E_HRHR_RR3 + E_HRhR_RR3 + E_HRRR_RR3 + E_hRHh_RR3 + E_hRHR_RR3 + E_hRhR_RR3 + E_hRRR_RR3 + E_RRHh_RR3 + E_RRHR_RR3 + E_RRhR_RR3 + E_RRRR_RR3),theta0*thetaLA*thetaP*(1-mu_m_RR)),[(1/2),(1/2)]);

    % Adult male genotypes:
    
    MHH3(t) = max(0,binornd(MHH3(t-1),(1-mu_m_HH)) + MHH_Births3(1));
    MHh3(t) = max(0,binornd(MHh3(t-1),(1-mu_m_Hh)) + MHh_Births3(1));
    MHR3(t) = max(0,binornd(MHR3(t-1),(1-mu_m_HR)) + MHR_Births3(1));
    Mhh3(t) = max(0,binornd(Mhh3(t-1),(1-mu_m_hh)) + Mhh_Births3(1));
    MhR3(t) = max(0,binornd(MhR3(t-1),(1-mu_m_hR)) + MhR_Births3(1));
    MRR3(t) = max(0,binornd(MRR3(t-1),(1-mu_m_RR)) + MRR_Births3(1));
    
    % Adult female genotypes:
    
    X_HH3 = mnrnd(MHH_Births3(2) , [ (MHH3(t-1)/(MHH3(t-1)+MHh3(t-1)+MHR3(t-1)+Mhh3(t-1)+MhR3(t-1)+MRR3(t-1))) , (MHh3(t-1)/(MHH3(t-1)+MHh3(t-1)+MHR3(t-1)+Mhh3(t-1)+MhR3(t-1)+MRR3(t-1))) , (MHR3(t-1)/(MHH3(t-1)+MHh3(t-1)+MHR3(t-1)+Mhh3(t-1)+MhR3(t-1)+MRR3(t-1))) , (Mhh3(t-1)/(MHH3(t-1)+MHh3(t-1)+MHR3(t-1)+Mhh3(t-1)+MhR3(t-1)+MRR3(t-1))) , (MhR3(t-1)/(MHH3(t-1)+MHh3(t-1)+MHR3(t-1)+Mhh3(t-1)+MhR3(t-1)+MRR3(t-1))) , (MRR3(t-1)/(MHH3(t-1)+MHh3(t-1)+MHR3(t-1)+Mhh3(t-1)+MhR3(t-1)+MRR3(t-1))) ]);
    X_Hh3 = mnrnd(MHh_Births3(2) , [ (MHH3(t-1)/(MHH3(t-1)+MHh3(t-1)+MHR3(t-1)+Mhh3(t-1)+MhR3(t-1)+MRR3(t-1))) , (MHh3(t-1)/(MHH3(t-1)+MHh3(t-1)+MHR3(t-1)+Mhh3(t-1)+MhR3(t-1)+MRR3(t-1))) , (MHR3(t-1)/(MHH3(t-1)+MHh3(t-1)+MHR3(t-1)+Mhh3(t-1)+MhR3(t-1)+MRR3(t-1))) , (Mhh3(t-1)/(MHH3(t-1)+MHh3(t-1)+MHR3(t-1)+Mhh3(t-1)+MhR3(t-1)+MRR3(t-1))) , (MhR3(t-1)/(MHH3(t-1)+MHh3(t-1)+MHR3(t-1)+Mhh3(t-1)+MhR3(t-1)+MRR3(t-1))) , (MRR3(t-1)/(MHH3(t-1)+MHh3(t-1)+MHR3(t-1)+Mhh3(t-1)+MhR3(t-1)+MRR3(t-1))) ]);
    X_HR3 = mnrnd(MHR_Births3(2) , [ (MHH3(t-1)/(MHH3(t-1)+MHh3(t-1)+MHR3(t-1)+Mhh3(t-1)+MhR3(t-1)+MRR3(t-1))) , (MHh3(t-1)/(MHH3(t-1)+MHh3(t-1)+MHR3(t-1)+Mhh3(t-1)+MhR3(t-1)+MRR3(t-1))) , (MHR3(t-1)/(MHH3(t-1)+MHh3(t-1)+MHR3(t-1)+Mhh3(t-1)+MhR3(t-1)+MRR3(t-1))) , (Mhh3(t-1)/(MHH3(t-1)+MHh3(t-1)+MHR3(t-1)+Mhh3(t-1)+MhR3(t-1)+MRR3(t-1))) , (MhR3(t-1)/(MHH3(t-1)+MHh3(t-1)+MHR3(t-1)+Mhh3(t-1)+MhR3(t-1)+MRR3(t-1))) , (MRR3(t-1)/(MHH3(t-1)+MHh3(t-1)+MHR3(t-1)+Mhh3(t-1)+MhR3(t-1)+MRR3(t-1))) ]);
    X_hh3 = mnrnd(Mhh_Births3(2) , [ (MHH3(t-1)/(MHH3(t-1)+MHh3(t-1)+MHR3(t-1)+Mhh3(t-1)+MhR3(t-1)+MRR3(t-1))) , (MHh3(t-1)/(MHH3(t-1)+MHh3(t-1)+MHR3(t-1)+Mhh3(t-1)+MhR3(t-1)+MRR3(t-1))) , (MHR3(t-1)/(MHH3(t-1)+MHh3(t-1)+MHR3(t-1)+Mhh3(t-1)+MhR3(t-1)+MRR3(t-1))) , (Mhh3(t-1)/(MHH3(t-1)+MHh3(t-1)+MHR3(t-1)+Mhh3(t-1)+MhR3(t-1)+MRR3(t-1))) , (MhR3(t-1)/(MHH3(t-1)+MHh3(t-1)+MHR3(t-1)+Mhh3(t-1)+MhR3(t-1)+MRR3(t-1))) , (MRR3(t-1)/(MHH3(t-1)+MHh3(t-1)+MHR3(t-1)+Mhh3(t-1)+MhR3(t-1)+MRR3(t-1))) ]);
    X_hR3 = mnrnd(MhR_Births3(2) , [ (MHH3(t-1)/(MHH3(t-1)+MHh3(t-1)+MHR3(t-1)+Mhh3(t-1)+MhR3(t-1)+MRR3(t-1))) , (MHh3(t-1)/(MHH3(t-1)+MHh3(t-1)+MHR3(t-1)+Mhh3(t-1)+MhR3(t-1)+MRR3(t-1))) , (MHR3(t-1)/(MHH3(t-1)+MHh3(t-1)+MHR3(t-1)+Mhh3(t-1)+MhR3(t-1)+MRR3(t-1))) , (Mhh3(t-1)/(MHH3(t-1)+MHh3(t-1)+MHR3(t-1)+Mhh3(t-1)+MhR3(t-1)+MRR3(t-1))) , (MhR3(t-1)/(MHH3(t-1)+MHh3(t-1)+MHR3(t-1)+Mhh3(t-1)+MhR3(t-1)+MRR3(t-1))) , (MRR3(t-1)/(MHH3(t-1)+MHh3(t-1)+MHR3(t-1)+Mhh3(t-1)+MhR3(t-1)+MRR3(t-1))) ]);    
    X_RR3 = mnrnd(MRR_Births3(2) , [ (MHH3(t-1)/(MHH3(t-1)+MHh3(t-1)+MHR3(t-1)+Mhh3(t-1)+MhR3(t-1)+MRR3(t-1))) , (MHh3(t-1)/(MHH3(t-1)+MHh3(t-1)+MHR3(t-1)+Mhh3(t-1)+MhR3(t-1)+MRR3(t-1))) , (MHR3(t-1)/(MHH3(t-1)+MHh3(t-1)+MHR3(t-1)+Mhh3(t-1)+MhR3(t-1)+MRR3(t-1))) , (Mhh3(t-1)/(MHH3(t-1)+MHh3(t-1)+MHR3(t-1)+Mhh3(t-1)+MhR3(t-1)+MRR3(t-1))) , (MhR3(t-1)/(MHH3(t-1)+MHh3(t-1)+MHR3(t-1)+Mhh3(t-1)+MhR3(t-1)+MRR3(t-1))) , (MRR3(t-1)/(MHH3(t-1)+MHh3(t-1)+MHR3(t-1)+Mhh3(t-1)+MhR3(t-1)+MRR3(t-1))) ]);

    MHH3_Temp = [binornd(MHH_HH3(t-1),(1-mu_m_HH)) , binornd(MHH_Hh3(t-1),(1-mu_m_HH)) , binornd(MHH_HR3(t-1),(1-mu_m_HH)) , binornd(MHH_hh3(t-1),(1-mu_m_HH)) , binornd(MHH_hR3(t-1),(1-mu_m_HH)) , binornd(MHH_RR3(t-1),(1-mu_m_HH))] + X_HH3;
    MHH_HH3(t) = MHH3_Temp(1); MHH_Hh3(t) = MHH3_Temp(2); MHH_HR3(t) = MHH3_Temp(3); MHH_hh3(t) = MHH3_Temp(4); MHH_hR3(t) = MHH3_Temp(5); MHH_RR3(t) = MHH3_Temp(6);

    MHh3_Temp = [binornd(MHh_HH3(t-1),(1-mu_m_Hh)) , binornd(MHh_Hh3(t-1),(1-mu_m_Hh)) , binornd(MHh_HR3(t-1),(1-mu_m_Hh)) , binornd(MHh_hh3(t-1),(1-mu_m_Hh)) , binornd(MHh_hR3(t-1),(1-mu_m_Hh)) , binornd(MHh_RR3(t-1),(1-mu_m_Hh))] + X_Hh3;
    MHh_HH3(t) = MHh3_Temp(1); MHh_Hh3(t) = MHh3_Temp(2); MHh_HR3(t) = MHh3_Temp(3); MHh_hh3(t) = MHh3_Temp(4); MHh_hR3(t) = MHh3_Temp(5); MHh_RR3(t) = MHh3_Temp(6);

    MHR3_Temp = [binornd(MHR_HH3(t-1),(1-mu_m_HR)) , binornd(MHR_Hh3(t-1),(1-mu_m_HR)) , binornd(MHR_HR3(t-1),(1-mu_m_HR)) , binornd(MHR_hh3(t-1),(1-mu_m_HR)) , binornd(MHR_hR3(t-1),(1-mu_m_HR)) , binornd(MHR_RR3(t-1),(1-mu_m_HR))] + X_HR3;
    MHR_HH3(t) = MHR3_Temp(1); MHR_Hh3(t) = MHR3_Temp(2); MHR_HR3(t) = MHR3_Temp(3); MHR_hh3(t) = MHR3_Temp(4); MHR_hR3(t) = MHR3_Temp(5); MHR_RR3(t) = MHR3_Temp(6);

    Mhh3_Temp = [binornd(Mhh_HH3(t-1),(1-mu_m_hh)) , binornd(Mhh_Hh3(t-1),(1-mu_m_hh)) , binornd(Mhh_HR3(t-1),(1-mu_m_hh)) , binornd(Mhh_hh3(t-1),(1-mu_m_hh)) , binornd(Mhh_hR3(t-1),(1-mu_m_hh)) , binornd(Mhh_RR3(t-1),(1-mu_m_hh))] + X_hh3;
    Mhh_HH3(t) = Mhh3_Temp(1); Mhh_Hh3(t) = Mhh3_Temp(2); Mhh_HR3(t) = Mhh3_Temp(3); Mhh_hh3(t) = Mhh3_Temp(4); Mhh_hR3(t) = Mhh3_Temp(5); Mhh_RR3(t) = Mhh3_Temp(6);

    MhR3_Temp = [binornd(MhR_HH3(t-1),(1-mu_m_hR)) , binornd(MhR_Hh3(t-1),(1-mu_m_hR)) , binornd(MhR_HR3(t-1),(1-mu_m_hR)) , binornd(MhR_hh3(t-1),(1-mu_m_hR)) , binornd(MhR_hR3(t-1),(1-mu_m_hR)) , binornd(MhR_RR3(t-1),(1-mu_m_hR))] + X_hR3;
    MhR_HH3(t) = MhR3_Temp(1); MhR_Hh3(t) = MhR3_Temp(2); MhR_HR3(t) = MhR3_Temp(3); MhR_hh3(t) = MhR3_Temp(4); MhR_hR3(t) = MhR3_Temp(5); MhR_RR3(t) = MhR3_Temp(6);

    MRR3_Temp = [binornd(MRR_HH3(t-1),(1-mu_m_RR)) , binornd(MRR_Hh3(t-1),(1-mu_m_RR)) , binornd(MRR_HR3(t-1),(1-mu_m_RR)) , binornd(MRR_hh3(t-1),(1-mu_m_RR)) , binornd(MRR_hR3(t-1),(1-mu_m_RR)) , binornd(MRR_RR3(t-1),(1-mu_m_RR))] + X_RR3;
    MRR_HH3(t) = MRR3_Temp(1); MRR_Hh3(t) = MRR3_Temp(2); MRR_HR3(t) = MRR3_Temp(3); MRR_hh3(t) = MRR3_Temp(4); MRR_hR3(t) = MRR3_Temp(5); MRR_RR3(t) = MRR3_Temp(6);
    
    MHH_HH3(t) = max(0,MHH_HH3(t)); MHH_Hh3(t) = max(0,MHH_Hh3(t)); MHH_HR3(t) = max(0,MHH_HR3(t)); MHH_hh3(t) = max(0,MHH_hh3(t)); MHH_hR3(t) = max(0,MHH_hR3(t)); MHH_RR3(t) = max(0,MHH_RR3(t));
    MHh_HH3(t) = max(0,MHh_HH3(t)); MHh_Hh3(t) = max(0,MHh_Hh3(t)); MHh_HR3(t) = max(0,MHh_HR3(t)); MHh_hh3(t) = max(0,MHh_hh3(t)); MHh_hR3(t) = max(0,MHh_hR3(t)); MHh_RR3(t) = max(0,MHh_RR3(t));
    MHR_HH3(t) = max(0,MHR_HH3(t)); MHR_Hh3(t) = max(0,MHR_Hh3(t)); MHR_HR3(t) = max(0,MHR_HR3(t)); MHR_hh3(t) = max(0,MHR_hh3(t)); MHR_hR3(t) = max(0,MHR_hR3(t)); MHR_RR3(t) = max(0,MHR_RR3(t));
    Mhh_HH3(t) = max(0,Mhh_HH3(t)); Mhh_Hh3(t) = max(0,Mhh_Hh3(t)); Mhh_HR3(t) = max(0,Mhh_HR3(t)); Mhh_hh3(t) = max(0,Mhh_hh3(t)); Mhh_hR3(t) = max(0,Mhh_hR3(t)); Mhh_RR3(t) = max(0,Mhh_RR3(t));
    MhR_HH3(t) = max(0,MhR_HH3(t)); MhR_Hh3(t) = max(0,MhR_Hh3(t)); MhR_HR3(t) = max(0,MhR_HR3(t)); MhR_hh3(t) = max(0,MhR_hh3(t)); MhR_hR3(t) = max(0,MhR_hR3(t)); MhR_RR3(t) = max(0,MhR_RR3(t));
    MRR_HH3(t) = max(0,MRR_HH3(t)); MRR_Hh3(t) = max(0,MRR_Hh3(t)); MRR_HR3(t) = max(0,MRR_HR3(t)); MRR_hh3(t) = max(0,MRR_hh3(t)); MRR_hR3(t) = max(0,MRR_hR3(t)); MRR_RR3(t) = max(0,MRR_RR3(t));

    if (t==51)
        MHH3(t) = MHH3(t) + M_eq;
        % MHH_HH3(t) = MHH_HH(t) + M_eq/2;
    end
    
    % Population 4:
    thetaLA = thetaL;
    for i=1:TL
        thetaLA = thetaLA*((alpha/(alpha+L4(t-i)))^(1/TL));
    end
    
    L4(t) = max(0,binornd(L4(t-1),(1-mu_l)*((alpha/(alpha+L4(t-1)))^(1/TL))) + binornd(poissrnd( beta_HH*MHH_HH4(t-T0) + beta_HH*MHH_Hh4(t-T0) + beta_HH*MHH_HR4(t-T0) + beta_HH*MHH_hh4(t-T0) + beta_HH*MHH_hR4(t-T0) + beta_HH*MHH_RR4(t-T0)  +  beta_Hh*MHh_HH4(t-T0) + beta_Hh*MHh_Hh4(t-T0) + beta_Hh*MHh_HR4(t-T0) + beta_Hh*MHh_hh4(t-T0) + beta_Hh*MHh_hR4(t-T0) + beta_Hh*MHh_RR4(t-T0)  +  beta_HR*MHR_HH4(t-T0) + beta_HR*MHR_Hh4(t-T0) + beta_HR*MHR_HR4(t-T0) + beta_HR*MHR_hh4(t-T0) + beta_HR*MHR_hR4(t-T0) + beta_HR*MHR_RR4(t-T0)  +  beta_hh*Mhh_HH4(t-T0) + beta_hh*Mhh_Hh4(t-T0) + beta_hh*Mhh_HR4(t-T0) + beta_hh*Mhh_hh4(t-T0) + beta_hh*Mhh_hR4(t-T0) + beta_hh*Mhh_RR4(t-T0)  +  beta_hR*MhR_HH4(t-T0) + beta_hR*MhR_Hh4(t-T0) + beta_hR*MhR_HR4(t-T0) + beta_hR*MhR_hh4(t-T0) + beta_hR*MhR_hR4(t-T0) + beta_hR*MhR_RR4(t-T0)  +  beta_RR*MRR_HH4(t-T0) + beta_RR*MRR_Hh4(t-T0) + beta_RR*MRR_HR4(t-T0) + beta_RR*MRR_hh4(t-T0) + beta_RR*MRR_hR4(t-T0) + beta_RR*MRR_RR4(t-T0)  ),theta0) - binornd(poissrnd(beta_HH*MHH_HH4(t-T0-TL) + beta_HH*MHH_Hh4(t-T0-TL) + beta_HH*MHH_HR4(t-T0-TL) + beta_HH*MHH_hh4(t-T0-TL) + beta_HH*MHH_hR4(t-T0-TL) + beta_HH*MHH_RR4(t-T0-TL)  +  beta_Hh*MHh_HH4(t-T0-TL) + beta_Hh*MHh_Hh4(t-T0-TL) + beta_Hh*MHh_HR4(t-T0-TL) + beta_Hh*MHh_hh4(t-T0-TL) + beta_Hh*MHh_hR4(t-T0-TL) + beta_Hh*MHh_RR4(t-T0-TL)  +  beta_HR*MHR_HH4(t-T0-TL) + beta_HR*MHR_Hh4(t-T0-TL) + beta_HR*MHR_HR4(t-T0-TL) + beta_HR*MHR_hh4(t-T0-TL) + beta_HR*MHR_hR4(t-T0-TL) + beta_HR*MHR_RR4(t-T0-TL)  +  beta_hh*Mhh_HH4(t-T0-TL) + beta_hh*Mhh_Hh4(t-T0-TL) + beta_hh*Mhh_HR4(t-T0-TL) + beta_hh*Mhh_hh4(t-T0-TL) + beta_hh*Mhh_hR4(t-T0-TL) + beta_hh*Mhh_RR4(t-T0-TL)  +  beta_hR*MhR_HH4(t-T0-TL) + beta_hR*MhR_Hh4(t-T0-TL) + beta_hR*MhR_HR4(t-T0-TL) + beta_hR*MhR_hh4(t-T0-TL) + beta_hR*MhR_hR4(t-T0-TL) + beta_hR*MhR_RR4(t-T0-TL)  +  beta_RR*MRR_HH4(t-T0-TL) + beta_RR*MRR_Hh4(t-T0-TL) + beta_RR*MRR_HR4(t-T0-TL) + beta_RR*MRR_hh4(t-T0-TL) + beta_RR*MRR_hR4(t-T0-TL) + beta_RR*MRR_RR4(t-T0-TL) ),theta0*thetaLA));
    
    thetaLA = thetaL;
    for i=1:TL
        thetaLA = thetaLA*((alpha/(alpha+L4(t-i-TP)))^(1/TL));
    end
    
    % Eggs produced by HH females:
    
    E_HHHH4 = poissrnd(beta_HH*MHH_HH4(t-T0-TL-TP));
    E_HHHH_HH4 = E_HHHH4;

    E_HHHh4 = mnrnd(poissrnd(beta_HH*MHH_Hh4(t-T0-TL-TP)),[((1+e)/2),((1-e-rho)/2),(rho/2)]);
    E_HHHh_HH4 = E_HHHh4(1); E_HHHh_Hh4 = E_HHHh4(2); E_HHHh_HR4 = E_HHHh4(3);

    E_HHHR4 = mnrnd(poissrnd(beta_HH*MHH_HR4(t-T0-TL-TP)),[(1/2),(1/2)]);
    E_HHHR_HH4 = E_HHHR4(1); E_HHHR_HR4 = E_HHHR4(2);

    E_HHhh4 = poissrnd(beta_HH*MHH_hh4(t-T0-TL-TP));
    E_HHhh_Hh4 = E_HHhh4;

    E_HHhR4 = mnrnd(poissrnd(beta_HH*MHH_hR4(t-T0-TL-TP)),[(1/2),(1/2)]);
    E_HHhR_Hh4 = E_HHhR4(1); E_HHhR_HR4 = E_HHhR4(2);

    E_HHRR4 = poissrnd(beta_HH*MHH_RR4(t-T0-TL-TP));
    E_HHRR_HR4 = E_HHRR4;

    % Eggs produced by Hh females:
    
    E_HhHH4 = mnrnd(poissrnd(beta_Hh*MHh_HH4(t-T0-TL-TP)),[((1+e)/2),((1-e-rho)/2),(rho/2)]);
    E_HhHH_HH4 = E_HhHH4(1); E_HhHH_Hh4 = E_HhHH4(2); E_HhHH_HR4 = E_HhHH4(3);

    E_HhHh4 = mnrnd(poissrnd(beta_Hh*MHh_Hh4(t-T0-TL-TP)),[(((1+e)^2)/4),((1+e)*(1-e-rho)/2),((1+e)*rho/2),(((1-e-rho)^2)/4),((1-e-rho)*rho/2),((rho^2)/4)]);
    E_HhHh_HH4 = E_HhHh4(1); E_HhHh_Hh4 = E_HhHh4(2); E_HhHh_HR4 = E_HhHh4(3); E_HhHh_hh4 = E_HhHh4(4); E_HhHh_hR4 = E_HhHh4(5); E_HhHh_RR4 = E_HhHh4(6);

    E_HhHR4 = mnrnd(poissrnd(beta_Hh*MHh_HR4(t-T0-TL-TP)),[((1+e)/4),((1-e-rho)/4),((1+e+rho)/4),((1-e-rho)/4),(rho/4)]);
    E_HhHR_HH4 = E_HhHR4(1); E_HhHR_Hh4 = E_HhHR4(2); E_HhHR_HR4 = E_HhHR4(3); E_HhHR_hR4 = E_HhHR4(4); E_HhHR_RR4 = E_HhHR4(5);

    E_Hhhh4 = mnrnd(poissrnd(beta_Hh*MHh_hh4(t-T0-TL-TP)),[((1+e)/2),((1-e-rho)/2),(rho/2)]);
    E_Hhhh_Hh4 = E_Hhhh4(1); E_Hhhh_hh4 = E_Hhhh4(2); E_Hhhh_hR4 = E_Hhhh4(3);

    E_HhhR4 = mnrnd(poissrnd(beta_Hh*MHh_hR4(t-T0-TL-TP)),[((1+e)/4),((1+e)/4),((1-e-rho)/4),((1-e)/4),(rho/4)]);
    E_HhhR_Hh4 = E_HhhR4(1); E_HhhR_HR4 = E_HhhR4(2); E_HhhR_hh4 = E_HhhR4(3); E_HhhR_hR4 = E_HhhR4(4); E_HhhR_RR4 = E_HhhR4(5);
    
    E_HhRR4 = mnrnd(poissrnd(beta_Hh*MHh_RR4(t-T0-TL-TP)),[((1+e)/2),((1-e-rho)/2),(rho/2)]);
    E_HhRR_HR4 = E_HhRR4(1); E_HhRR_hR4 = E_HhRR4(2); E_HhRR_RR4 = E_HhRR4(3);
    
    % Eggs produced by HR females:
    
    E_HRHH4 = mnrnd(poissrnd(beta_HR*MHR_HH4(t-T0-TL-TP)),[(1/2),(1/2)]);
    E_HRHH_HH4 = E_HRHH4(1); E_HRHH_HR4 = E_HRHH4(2);

    E_HRHh4 = mnrnd(poissrnd(beta_HR*MHR_Hh4(t-T0-TL-TP)),[((1+e)/4),((1-e-rho)/4),((1+e+rho)/4),((1-e-rho)/4),(rho/4)]);
    E_HRHh_HH4 = E_HRHh4(1); E_HRHh_Hh4 = E_HRHh4(2); E_HRHh_HR4 = E_HRHh4(3); E_HRHh_hR4 = E_HRHh4(4); E_HRHh_RR4 = E_HRHh4(5);

    E_HRHR4 = mnrnd(poissrnd(beta_HR*MHR_HR4(t-T0-TL-TP)),[(1/4),(1/2),(1/4)]);
    E_HRHR_HH4 = E_HRHR4(1); E_HRHR_HR4 = E_HRHR4(2); E_HRHR_RR4 = E_HRHR4(3);

    E_HRhh4 = mnrnd(poissrnd(beta_HR*MHR_hh4(t-T0-TL-TP)),[(1/2),(1/2)]);
    E_HRhh_Hh4 = E_HRhh4(1); E_HRhh_hR4 = E_HRhh4(2);

    E_HRhR4 = mnrnd(poissrnd(beta_HR*MHR_hR4(t-T0-TL-TP)),[(1/4),(1/4),(1/4),(1/4)]);
    E_HRhR_Hh4 = E_HRhR4(1); E_HRhR_HR4 = E_HRhR4(2); E_HRhR_hR4 = E_HRhR4(3); E_HRhR_RR4 = E_HRhR4(4);
    
    E_HRRR4 = mnrnd(poissrnd(beta_HR*MHR_RR4(t-T0-TL-TP)),[(1/2),(1/2)]);
    E_HRRR_HR4 = E_HRRR4(1); E_HRRR_RR4 = E_HRRR4(2);
    
    % Eggs produced by hh females:
    
    E_hhHH4 = poissrnd(beta_hh*Mhh_HH4(t-T0-TL-TP));
    E_hhHH_Hh4 = E_hhHH4;

    E_hhHh4 = mnrnd(poissrnd(beta_hh*Mhh_Hh4(t-T0-TL-TP)),[((1+e)/2),((1-e-rho)/2),(rho/2)]);
    E_hhHh_Hh4 = E_hhHh4(1); E_hhHh_hh4 = E_hhHh4(2); E_hhHh_hR4 = E_hhHh4(3);

    E_hhHR4 = mnrnd(poissrnd(beta_hh*Mhh_HR4(t-T0-TL-TP)),[(1/2),(1/2)]);
    E_hhHR_Hh4 = E_hhHR4(1); E_hhHR_hR4 = E_hhHR4(2);

    E_hhhh4 = poissrnd(beta_hh*Mhh_hh4(t-T0-TL-TP));
    E_hhhh_hh4 = E_hhhh4;

    E_hhhR4 = mnrnd(poissrnd(beta_hh*Mhh_hR4(t-T0-TL-TP)),[(1/2),(1/2)]);
    E_hhhR_hh4 = E_hhhR4(1); E_hhhR_hR4 = E_hhhR4(2);

    E_hhRR4 = poissrnd(beta_hh*Mhh_RR4(t-T0-TL-TP));
    E_hhRR_hR4 = E_hhRR4;
    
    % Eggs produced by hR females:
    
    E_hRHH4 = mnrnd(poissrnd(beta_hR*MhR_HH4(t-T0-TL-TP)),[(1/2),(1/2)]);
    E_hRHH_Hh4 = E_hRHH4(1); E_hRHH_HR4 = E_hRHH4(2);

    E_hRHh4 = mnrnd(poissrnd(beta_hR*MhR_Hh4(t-T0-TL-TP)),[((1+e)/4),((1+e)/4),((1-e-rho)/4),((1-e)/4),(rho/4)]);
    E_hRHh_Hh4 = E_hRHh4(1); E_hRHh_HR4 = E_hRHh4(2); E_hRHh_hh4 = E_hRHh4(3); E_hRHh_hR4 = E_hRHh4(4); E_hRHh_RR4 = E_hRHh4(5);

    E_hRHR4 = mnrnd(poissrnd(beta_hR*MhR_HR4(t-T0-TL-TP)),[(1/4),(1/4),(1/4),(1/4)]);
    E_hRHR_Hh4 = E_hRHR4(1); E_hRHR_HR4 = E_hRHR4(2); E_hRHR_hR4 = E_hRHR4(3); E_hRHR_RR4 = E_hRHR4(4);

    E_hRhh4 = mnrnd(poissrnd(beta_hR*MHR_hh4(t-T0-TL-TP)),[(1/2),(1/2)]);
    E_hRhh_hh4 = E_hRhh4(1); E_hRhh_hR4 = E_hRhh4(2);

    E_hRhR4 = mnrnd(poissrnd(beta_hR*MhR_hR4(t-T0-TL-TP)),[(1/4),(1/2),(1/4)]);
    E_hRhR_hh4 = E_hRhR4(1); E_hRhR_hR4 = E_hRhR4(2); E_hRhR_RR4 = E_hRhR4(3);
    
    E_hRRR4 = mnrnd(poissrnd(beta_hR*MhR_RR4(t-T0-TL-TP)),[(1/2),(1/2)]);
    E_hRRR_hR4 = E_hRRR4(1); E_hRRR_RR4 = E_hRRR4(2);
    
    % Eggs produced by RR females:
    
    E_RRHH4 = poissrnd(beta_RR*MRR_HH4(t-T0-TL-TP));
    E_RRHH_HR4 = E_RRHH4;

    E_RRHh4 = mnrnd(poissrnd(beta_RR*MRR_Hh4(t-T0-TL-TP)),[((1+e)/2),((1-e-rho)/2),(rho/2)]);
    E_RRHh_HR4 = E_RRHh4(1); E_RRHh_hR4 = E_RRHh4(2); E_RRHh_RR4 = E_RRHh4(3);

    E_RRHR4 = mnrnd(poissrnd(beta_RR*MRR_HR4(t-T0-TL-TP)),[(1/2),(1/2)]);
    E_RRHR_HR4 = E_RRHR4(1); E_RRHR_RR4 = E_RRHR4(2);

    E_RRhh4 = poissrnd(beta_RR*MRR_hh4(t-T0-TL-TP));
    E_RRhh_hR4 = E_RRhh4;

    E_RRhR4 = mnrnd(poissrnd(beta_RR*MRR_hR4(t-T0-TL-TP)),[(1/2),(1/2)]);
    E_RRhR_hR4 = E_RRhR4(1); E_RRhR_RR4 = E_RRhR4(2);

    E_RRRR4 = poissrnd(beta_RR*MRR_RR4(t-T0-TL-TP));
    E_RRRR_RR4 = E_RRRR4;
    
    % Adult male & female newborns:
    
    MHH_Births4 = mnrnd(binornd((E_HHHH_HH4 + E_HHHh_HH4 + E_HHHR_HH4 + E_HhHH_HH4 + E_HhHh_HH4 + E_HhHR_HH4 + E_HRHH_HH4 + E_HRHh_HH4 + E_HRHR_HH4),theta0*thetaLA*thetaP*(1-mu_m_HH)),[(1/2),(1/2)]);
    MHh_Births4 = mnrnd(binornd((E_HHHh_Hh4 + E_HHhh_Hh4 + E_HHhR_Hh4 + E_HhHH_Hh4 + E_HhHh_Hh4 + E_HhHR_Hh4 + E_Hhhh_Hh4 + E_HhhR_Hh4 + E_HRHh_Hh4 + E_HRhh_Hh4 + E_HRhR_Hh4 + E_hhHH_Hh4 + E_hhHh_Hh4 + E_hhHR_Hh4 + E_hRHH_Hh4 + E_hRHh_Hh4 + E_hRHR_Hh4),theta0*thetaLA*thetaP*(1-mu_m_Hh)),[(1/2),(1/2)]);
    MHR_Births4 = mnrnd(binornd((E_HHHh_HR4 + E_HHHR_HR4 + E_HHhR_HR4 + E_HHRR_HR4 + E_HhHH_HR4 + E_HhHh_HR4 + E_HhHR_HR4 + E_HhhR_HR4 + E_HhRR_HR4 + E_HRHH_HR4 + E_HRHh_HR4 + E_HRHR_HR4 + E_HRhR_HR4 + E_HRRR_HR4 + E_hRHH_HR4 + E_hRHh_HR4 + E_hRHR_HR4 + E_RRHH_HR4 + E_RRHh_HR4 + E_RRHR_HR4),theta0*thetaLA*thetaP*(1-mu_m_HR)),[(1/2),(1/2)]);
    Mhh_Births4 = mnrnd(binornd((E_HhHh_hh4 + E_Hhhh_hh4 + E_HhhR_hh4 + E_hhHh_hh4 + E_hhhh_hh4 + E_hhhR_hh4 + E_hRHh_hh4 + E_hRhh_hh4 + E_hRhR_hh4),theta0*thetaLA*thetaP*(1-mu_m_hh)),[(1/2),(1/2)]);
    MhR_Births4 = mnrnd(binornd((E_HhHh_hR4 + E_HhHR_hR4 + E_Hhhh_hR4 + E_HhhR_hR4 + E_HhRR_hR4 + E_HRHh_hR4 + E_HRhh_hR4 + E_HRhR_hR4 + E_hhHh_hR4 + E_hhHR_hR4 + E_hhhR_hR4 + E_hhRR_hR4 + E_hRHh_hR4 + E_hRHR_hR4 + E_hRhh_hR4 + E_hRhR_hR4 + E_hRRR_hR4 + E_RRHh_hR4 + E_RRhh_hR4 + E_RRhR_hR4),theta0*thetaLA*thetaP*(1-mu_m_hR)),[(1/2),(1/2)]);
    MRR_Births4 = mnrnd(binornd((E_HhHh_RR4 + E_HhHR_RR4 + E_HhhR_RR4 + E_HhRR_RR4 + E_HRHh_RR4 + E_HRHR_RR4 + E_HRhR_RR4 + E_HRRR_RR4 + E_hRHh_RR4 + E_hRHR_RR4 + E_hRhR_RR4 + E_hRRR_RR4 + E_RRHh_RR4 + E_RRHR_RR4 + E_RRhR_RR4 + E_RRRR_RR4),theta0*thetaLA*thetaP*(1-mu_m_RR)),[(1/2),(1/2)]);

    % Adult male genotypes:
    
    MHH4(t) = max(0,binornd(MHH4(t-1),(1-mu_m_HH)) + MHH_Births4(1));
    MHh4(t) = max(0,binornd(MHh4(t-1),(1-mu_m_Hh)) + MHh_Births4(1));
    MHR4(t) = max(0,binornd(MHR4(t-1),(1-mu_m_HR)) + MHR_Births4(1));
    Mhh4(t) = max(0,binornd(Mhh4(t-1),(1-mu_m_hh)) + Mhh_Births4(1));
    MhR4(t) = max(0,binornd(MhR4(t-1),(1-mu_m_hR)) + MhR_Births4(1));
    MRR4(t) = max(0,binornd(MRR4(t-1),(1-mu_m_RR)) + MRR_Births4(1));
    
    % Adult female genotypes:
    
    X_HH4 = mnrnd(MHH_Births4(2) , [ (MHH4(t-1)/(MHH4(t-1)+MHh4(t-1)+MHR4(t-1)+Mhh4(t-1)+MhR4(t-1)+MRR4(t-1))) , (MHh4(t-1)/(MHH4(t-1)+MHh4(t-1)+MHR4(t-1)+Mhh4(t-1)+MhR4(t-1)+MRR4(t-1))) , (MHR4(t-1)/(MHH4(t-1)+MHh4(t-1)+MHR4(t-1)+Mhh4(t-1)+MhR4(t-1)+MRR4(t-1))) , (Mhh4(t-1)/(MHH4(t-1)+MHh4(t-1)+MHR4(t-1)+Mhh4(t-1)+MhR4(t-1)+MRR4(t-1))) , (MhR4(t-1)/(MHH4(t-1)+MHh4(t-1)+MHR4(t-1)+Mhh4(t-1)+MhR4(t-1)+MRR4(t-1))) , (MRR4(t-1)/(MHH4(t-1)+MHh4(t-1)+MHR4(t-1)+Mhh4(t-1)+MhR4(t-1)+MRR4(t-1))) ]);
    X_Hh4 = mnrnd(MHh_Births4(2) , [ (MHH4(t-1)/(MHH4(t-1)+MHh4(t-1)+MHR4(t-1)+Mhh4(t-1)+MhR4(t-1)+MRR4(t-1))) , (MHh4(t-1)/(MHH4(t-1)+MHh4(t-1)+MHR4(t-1)+Mhh4(t-1)+MhR4(t-1)+MRR4(t-1))) , (MHR4(t-1)/(MHH4(t-1)+MHh4(t-1)+MHR4(t-1)+Mhh4(t-1)+MhR4(t-1)+MRR4(t-1))) , (Mhh4(t-1)/(MHH4(t-1)+MHh4(t-1)+MHR4(t-1)+Mhh4(t-1)+MhR4(t-1)+MRR4(t-1))) , (MhR4(t-1)/(MHH4(t-1)+MHh4(t-1)+MHR4(t-1)+Mhh4(t-1)+MhR4(t-1)+MRR4(t-1))) , (MRR4(t-1)/(MHH4(t-1)+MHh4(t-1)+MHR4(t-1)+Mhh4(t-1)+MhR4(t-1)+MRR4(t-1))) ]);
    X_HR4 = mnrnd(MHR_Births4(2) , [ (MHH4(t-1)/(MHH4(t-1)+MHh4(t-1)+MHR4(t-1)+Mhh4(t-1)+MhR4(t-1)+MRR4(t-1))) , (MHh4(t-1)/(MHH4(t-1)+MHh4(t-1)+MHR4(t-1)+Mhh4(t-1)+MhR4(t-1)+MRR4(t-1))) , (MHR4(t-1)/(MHH4(t-1)+MHh4(t-1)+MHR4(t-1)+Mhh4(t-1)+MhR4(t-1)+MRR4(t-1))) , (Mhh4(t-1)/(MHH4(t-1)+MHh4(t-1)+MHR4(t-1)+Mhh4(t-1)+MhR4(t-1)+MRR4(t-1))) , (MhR4(t-1)/(MHH4(t-1)+MHh4(t-1)+MHR4(t-1)+Mhh4(t-1)+MhR4(t-1)+MRR4(t-1))) , (MRR4(t-1)/(MHH4(t-1)+MHh4(t-1)+MHR4(t-1)+Mhh4(t-1)+MhR4(t-1)+MRR4(t-1))) ]);
    X_hh4 = mnrnd(Mhh_Births4(2) , [ (MHH4(t-1)/(MHH4(t-1)+MHh4(t-1)+MHR4(t-1)+Mhh4(t-1)+MhR4(t-1)+MRR4(t-1))) , (MHh4(t-1)/(MHH4(t-1)+MHh4(t-1)+MHR4(t-1)+Mhh4(t-1)+MhR4(t-1)+MRR4(t-1))) , (MHR4(t-1)/(MHH4(t-1)+MHh4(t-1)+MHR4(t-1)+Mhh4(t-1)+MhR4(t-1)+MRR4(t-1))) , (Mhh4(t-1)/(MHH4(t-1)+MHh4(t-1)+MHR4(t-1)+Mhh4(t-1)+MhR4(t-1)+MRR4(t-1))) , (MhR4(t-1)/(MHH4(t-1)+MHh4(t-1)+MHR4(t-1)+Mhh4(t-1)+MhR4(t-1)+MRR4(t-1))) , (MRR4(t-1)/(MHH4(t-1)+MHh4(t-1)+MHR4(t-1)+Mhh4(t-1)+MhR4(t-1)+MRR4(t-1))) ]);
    X_hR4 = mnrnd(MhR_Births4(2) , [ (MHH4(t-1)/(MHH4(t-1)+MHh4(t-1)+MHR4(t-1)+Mhh4(t-1)+MhR4(t-1)+MRR4(t-1))) , (MHh4(t-1)/(MHH4(t-1)+MHh4(t-1)+MHR4(t-1)+Mhh4(t-1)+MhR4(t-1)+MRR4(t-1))) , (MHR4(t-1)/(MHH4(t-1)+MHh4(t-1)+MHR4(t-1)+Mhh4(t-1)+MhR4(t-1)+MRR4(t-1))) , (Mhh4(t-1)/(MHH4(t-1)+MHh4(t-1)+MHR4(t-1)+Mhh4(t-1)+MhR4(t-1)+MRR4(t-1))) , (MhR4(t-1)/(MHH4(t-1)+MHh4(t-1)+MHR4(t-1)+Mhh4(t-1)+MhR4(t-1)+MRR4(t-1))) , (MRR4(t-1)/(MHH4(t-1)+MHh4(t-1)+MHR4(t-1)+Mhh4(t-1)+MhR4(t-1)+MRR4(t-1))) ]);    
    X_RR4 = mnrnd(MRR_Births4(2) , [ (MHH4(t-1)/(MHH4(t-1)+MHh4(t-1)+MHR4(t-1)+Mhh4(t-1)+MhR4(t-1)+MRR4(t-1))) , (MHh4(t-1)/(MHH4(t-1)+MHh4(t-1)+MHR4(t-1)+Mhh4(t-1)+MhR4(t-1)+MRR4(t-1))) , (MHR4(t-1)/(MHH4(t-1)+MHh4(t-1)+MHR4(t-1)+Mhh4(t-1)+MhR4(t-1)+MRR4(t-1))) , (Mhh4(t-1)/(MHH4(t-1)+MHh4(t-1)+MHR4(t-1)+Mhh4(t-1)+MhR4(t-1)+MRR4(t-1))) , (MhR4(t-1)/(MHH4(t-1)+MHh4(t-1)+MHR4(t-1)+Mhh4(t-1)+MhR4(t-1)+MRR4(t-1))) , (MRR4(t-1)/(MHH4(t-1)+MHh4(t-1)+MHR4(t-1)+Mhh4(t-1)+MhR4(t-1)+MRR4(t-1))) ]);

    MHH4_Temp = [binornd(MHH_HH4(t-1),(1-mu_m_HH)) , binornd(MHH_Hh4(t-1),(1-mu_m_HH)) , binornd(MHH_HR4(t-1),(1-mu_m_HH)) , binornd(MHH_hh4(t-1),(1-mu_m_HH)) , binornd(MHH_hR4(t-1),(1-mu_m_HH)) , binornd(MHH_RR4(t-1),(1-mu_m_HH))] + X_HH4;
    MHH_HH4(t) = MHH4_Temp(1); MHH_Hh4(t) = MHH4_Temp(2); MHH_HR4(t) = MHH4_Temp(4); MHH_hh4(t) = MHH4_Temp(4); MHH_hR4(t) = MHH4_Temp(5); MHH_RR4(t) = MHH4_Temp(6);

    MHh4_Temp = [binornd(MHh_HH4(t-1),(1-mu_m_Hh)) , binornd(MHh_Hh4(t-1),(1-mu_m_Hh)) , binornd(MHh_HR4(t-1),(1-mu_m_Hh)) , binornd(MHh_hh4(t-1),(1-mu_m_Hh)) , binornd(MHh_hR4(t-1),(1-mu_m_Hh)) , binornd(MHh_RR4(t-1),(1-mu_m_Hh))] + X_Hh4;
    MHh_HH4(t) = MHh4_Temp(1); MHh_Hh4(t) = MHh4_Temp(2); MHh_HR4(t) = MHh4_Temp(3); MHh_hh4(t) = MHh4_Temp(4); MHh_hR4(t) = MHh4_Temp(5); MHh_RR4(t) = MHh4_Temp(6);

    MHR4_Temp = [binornd(MHR_HH4(t-1),(1-mu_m_HR)) , binornd(MHR_Hh4(t-1),(1-mu_m_HR)) , binornd(MHR_HR4(t-1),(1-mu_m_HR)) , binornd(MHR_hh4(t-1),(1-mu_m_HR)) , binornd(MHR_hR4(t-1),(1-mu_m_HR)) , binornd(MHR_RR4(t-1),(1-mu_m_HR))] + X_HR4;
    MHR_HH4(t) = MHR4_Temp(1); MHR_Hh4(t) = MHR4_Temp(2); MHR_HR4(t) = MHR4_Temp(3); MHR_hh4(t) = MHR4_Temp(4); MHR_hR4(t) = MHR4_Temp(5); MHR_RR4(t) = MHR4_Temp(6);

    Mhh4_Temp = [binornd(Mhh_HH4(t-1),(1-mu_m_hh)) , binornd(Mhh_Hh4(t-1),(1-mu_m_hh)) , binornd(Mhh_HR4(t-1),(1-mu_m_hh)) , binornd(Mhh_hh4(t-1),(1-mu_m_hh)) , binornd(Mhh_hR4(t-1),(1-mu_m_hh)) , binornd(Mhh_RR4(t-1),(1-mu_m_hh))] + X_hh4;
    Mhh_HH4(t) = Mhh4_Temp(1); Mhh_Hh4(t) = Mhh4_Temp(2); Mhh_HR4(t) = Mhh4_Temp(3); Mhh_hh4(t) = Mhh4_Temp(4); Mhh_hR4(t) = Mhh4_Temp(5); Mhh_RR4(t) = Mhh4_Temp(6);

    MhR4_Temp = [binornd(MhR_HH4(t-1),(1-mu_m_hR)) , binornd(MhR_Hh4(t-1),(1-mu_m_hR)) , binornd(MhR_HR4(t-1),(1-mu_m_hR)) , binornd(MhR_hh4(t-1),(1-mu_m_hR)) , binornd(MhR_hR4(t-1),(1-mu_m_hR)) , binornd(MhR_RR4(t-1),(1-mu_m_hR))] + X_hR4;
    MhR_HH4(t) = MhR4_Temp(1); MhR_Hh4(t) = MhR4_Temp(2); MhR_HR4(t) = MhR4_Temp(3); MhR_hh4(t) = MhR4_Temp(4); MhR_hR4(t) = MhR4_Temp(5); MhR_RR4(t) = MhR4_Temp(6);

    MRR4_Temp = [binornd(MRR_HH4(t-1),(1-mu_m_RR)) , binornd(MRR_Hh4(t-1),(1-mu_m_RR)) , binornd(MRR_HR4(t-1),(1-mu_m_RR)) , binornd(MRR_hh4(t-1),(1-mu_m_RR)) , binornd(MRR_hR4(t-1),(1-mu_m_RR)) , binornd(MRR_RR4(t-1),(1-mu_m_RR))] + X_RR4;
    MRR_HH4(t) = MRR4_Temp(1); MRR_Hh4(t) = MRR4_Temp(2); MRR_HR4(t) = MRR4_Temp(3); MRR_hh4(t) = MRR4_Temp(4); MRR_hR4(t) = MRR4_Temp(5); MRR_RR4(t) = MRR4_Temp(6);
    
    MHH_HH4(t) = max(0,MHH_HH4(t)); MHH_Hh4(t) = max(0,MHH_Hh4(t)); MHH_HR4(t) = max(0,MHH_HR4(t)); MHH_hh4(t) = max(0,MHH_hh4(t)); MHH_hR4(t) = max(0,MHH_hR4(t)); MHH_RR4(t) = max(0,MHH_RR4(t));
    MHh_HH4(t) = max(0,MHh_HH4(t)); MHh_Hh4(t) = max(0,MHh_Hh4(t)); MHh_HR4(t) = max(0,MHh_HR4(t)); MHh_hh4(t) = max(0,MHh_hh4(t)); MHh_hR4(t) = max(0,MHh_hR4(t)); MHh_RR4(t) = max(0,MHh_RR4(t));
    MHR_HH4(t) = max(0,MHR_HH4(t)); MHR_Hh4(t) = max(0,MHR_Hh4(t)); MHR_HR4(t) = max(0,MHR_HR4(t)); MHR_hh4(t) = max(0,MHR_hh4(t)); MHR_hR4(t) = max(0,MHR_hR4(t)); MHR_RR4(t) = max(0,MHR_RR4(t));
    Mhh_HH4(t) = max(0,Mhh_HH4(t)); Mhh_Hh4(t) = max(0,Mhh_Hh4(t)); Mhh_HR4(t) = max(0,Mhh_HR4(t)); Mhh_hh4(t) = max(0,Mhh_hh4(t)); Mhh_hR4(t) = max(0,Mhh_hR4(t)); Mhh_RR4(t) = max(0,Mhh_RR4(t));
    MhR_HH4(t) = max(0,MhR_HH4(t)); MhR_Hh4(t) = max(0,MhR_Hh4(t)); MhR_HR4(t) = max(0,MhR_HR4(t)); MhR_hh4(t) = max(0,MhR_hh4(t)); MhR_hR4(t) = max(0,MhR_hR4(t)); MhR_RR4(t) = max(0,MhR_RR4(t));
    MRR_HH4(t) = max(0,MRR_HH4(t)); MRR_Hh4(t) = max(0,MRR_Hh4(t)); MRR_HR4(t) = max(0,MRR_HR4(t)); MRR_hh4(t) = max(0,MRR_hh4(t)); MRR_hR4(t) = max(0,MRR_hR4(t)); MRR_RR4(t) = max(0,MRR_RR4(t));

    if (t==51)
        MHH4(t) = MHH4(t) + M_eq;
        % MHH_HH4(t) = MHH_HH(t) + M_eq/2;
    end
    
    % Migrants pop 1 -> pop 2:
    MHH1to2(t) = poissrnd(MHH1(t)*migrate/3);
    MHh1to2(t) = poissrnd(MHh1(t)*migrate/3);
    MHR1to2(t) = poissrnd(MHR1(t)*migrate/3);
    Mhh1to2(t) = poissrnd(Mhh1(t)*migrate/3);
    MhR1to2(t) = poissrnd(MhR1(t)*migrate/3);
    MRR1to2(t) = poissrnd(MRR1(t)*migrate/3);

    MHH_HH1to2(t) = poissrnd(MHH_HH1(t)*migrate/3);
    MHH_Hh1to2(t) = poissrnd(MHH_Hh1(t)*migrate/3);
    MHH_HR1to2(t) = poissrnd(MHH_HR1(t)*migrate/3);
    MHH_hh1to2(t) = poissrnd(MHH_hh1(t)*migrate/3);
    MHH_hR1to2(t) = poissrnd(MHH_hR1(t)*migrate/3);
    MHH_RR1to2(t) = poissrnd(MHH_RR1(t)*migrate/3);

    MHh_HH1to2(t) = poissrnd(MHh_HH1(t)*migrate/3);
    MHh_Hh1to2(t) = poissrnd(MHh_Hh1(t)*migrate/3);
    MHh_HR1to2(t) = poissrnd(MHh_HR1(t)*migrate/3);
    MHh_hh1to2(t) = poissrnd(MHh_hh1(t)*migrate/3);
    MHh_hR1to2(t) = poissrnd(MHh_hR1(t)*migrate/3);
    MHh_RR1to2(t) = poissrnd(MHh_RR1(t)*migrate/3);

    MHR_HH1to2(t) = poissrnd(MHR_HH1(t)*migrate/3);
    MHR_Hh1to2(t) = poissrnd(MHR_Hh1(t)*migrate/3);
    MHR_HR1to2(t) = poissrnd(MHR_HR1(t)*migrate/3);
    MHR_hh1to2(t) = poissrnd(MHR_hh1(t)*migrate/3);
    MHR_hR1to2(t) = poissrnd(MHR_hR1(t)*migrate/3);
    MHR_RR1to2(t) = poissrnd(MHR_RR1(t)*migrate/3);

    Mhh_HH1to2(t) = poissrnd(Mhh_HH1(t)*migrate/3);
    Mhh_Hh1to2(t) = poissrnd(Mhh_Hh1(t)*migrate/3);
    Mhh_HR1to2(t) = poissrnd(Mhh_HR1(t)*migrate/3);
    Mhh_hh1to2(t) = poissrnd(Mhh_hh1(t)*migrate/3);
    Mhh_hR1to2(t) = poissrnd(Mhh_hR1(t)*migrate/3);
    Mhh_RR1to2(t) = poissrnd(Mhh_RR1(t)*migrate/3);

    MhR_HH1to2(t) = poissrnd(MhR_HH1(t)*migrate/3);
    MhR_Hh1to2(t) = poissrnd(MhR_Hh1(t)*migrate/3);
    MhR_HR1to2(t) = poissrnd(MhR_HR1(t)*migrate/3);
    MhR_hh1to2(t) = poissrnd(MhR_hh1(t)*migrate/3);
    MhR_hR1to2(t) = poissrnd(MhR_hR1(t)*migrate/3);
    MhR_RR1to2(t) = poissrnd(MhR_RR1(t)*migrate/3);

    MRR_HH1to2(t) = poissrnd(MRR_HH1(t)*migrate/3);
    MRR_Hh1to2(t) = poissrnd(MRR_Hh1(t)*migrate/3);
    MRR_HR1to2(t) = poissrnd(MRR_HR1(t)*migrate/3);
    MRR_hh1to2(t) = poissrnd(MRR_hh1(t)*migrate/3);
    MRR_hR1to2(t) = poissrnd(MRR_hR1(t)*migrate/3);
    MRR_RR1to2(t) = poissrnd(MRR_RR1(t)*migrate/3);

    % Migrants pop 1 -> pop 3:
    MHH1to3(t) = poissrnd(MHH1(t)*migrate/3);
    MHh1to3(t) = poissrnd(MHh1(t)*migrate/3);
    MHR1to3(t) = poissrnd(MHR1(t)*migrate/3);
    Mhh1to3(t) = poissrnd(Mhh1(t)*migrate/3);
    MhR1to3(t) = poissrnd(MhR1(t)*migrate/3);
    MRR1to3(t) = poissrnd(MRR1(t)*migrate/3);

    MHH_HH1to3(t) = poissrnd(MHH_HH1(t)*migrate/3);
    MHH_Hh1to3(t) = poissrnd(MHH_Hh1(t)*migrate/3);
    MHH_HR1to3(t) = poissrnd(MHH_HR1(t)*migrate/3);
    MHH_hh1to3(t) = poissrnd(MHH_hh1(t)*migrate/3);
    MHH_hR1to3(t) = poissrnd(MHH_hR1(t)*migrate/3);
    MHH_RR1to3(t) = poissrnd(MHH_RR1(t)*migrate/3);

    MHh_HH1to3(t) = poissrnd(MHh_HH1(t)*migrate/3);
    MHh_Hh1to3(t) = poissrnd(MHh_Hh1(t)*migrate/3);
    MHh_HR1to3(t) = poissrnd(MHh_HR1(t)*migrate/3);
    MHh_hh1to3(t) = poissrnd(MHh_hh1(t)*migrate/3);
    MHh_hR1to3(t) = poissrnd(MHh_hR1(t)*migrate/3);
    MHh_RR1to3(t) = poissrnd(MHh_RR1(t)*migrate/3);

    MHR_HH1to3(t) = poissrnd(MHR_HH1(t)*migrate/3);
    MHR_Hh1to3(t) = poissrnd(MHR_Hh1(t)*migrate/3);
    MHR_HR1to3(t) = poissrnd(MHR_HR1(t)*migrate/3);
    MHR_hh1to3(t) = poissrnd(MHR_hh1(t)*migrate/3);
    MHR_hR1to3(t) = poissrnd(MHR_hR1(t)*migrate/3);
    MHR_RR1to3(t) = poissrnd(MHR_RR1(t)*migrate/3);

    Mhh_HH1to3(t) = poissrnd(Mhh_HH1(t)*migrate/3);
    Mhh_Hh1to3(t) = poissrnd(Mhh_Hh1(t)*migrate/3);
    Mhh_HR1to3(t) = poissrnd(Mhh_HR1(t)*migrate/3);
    Mhh_hh1to3(t) = poissrnd(Mhh_hh1(t)*migrate/3);
    Mhh_hR1to3(t) = poissrnd(Mhh_hR1(t)*migrate/3);
    Mhh_RR1to3(t) = poissrnd(Mhh_RR1(t)*migrate/3);

    MhR_HH1to3(t) = poissrnd(MhR_HH1(t)*migrate/3);
    MhR_Hh1to3(t) = poissrnd(MhR_Hh1(t)*migrate/3);
    MhR_HR1to3(t) = poissrnd(MhR_HR1(t)*migrate/3);
    MhR_hh1to3(t) = poissrnd(MhR_hh1(t)*migrate/3);
    MhR_hR1to3(t) = poissrnd(MhR_hR1(t)*migrate/3);
    MhR_RR1to3(t) = poissrnd(MhR_RR1(t)*migrate/3);

    MRR_HH1to3(t) = poissrnd(MRR_HH1(t)*migrate/3);
    MRR_Hh1to3(t) = poissrnd(MRR_Hh1(t)*migrate/3);
    MRR_HR1to3(t) = poissrnd(MRR_HR1(t)*migrate/3);
    MRR_hh1to3(t) = poissrnd(MRR_hh1(t)*migrate/3);
    MRR_hR1to3(t) = poissrnd(MRR_hR1(t)*migrate/3);
    MRR_RR1to3(t) = poissrnd(MRR_RR1(t)*migrate/3);

    % Migrants pop 1 -> pop 4:
    MHH1to4(t) = poissrnd(MHH1(t)*migrate/3);
    MHh1to4(t) = poissrnd(MHh1(t)*migrate/3);
    MHR1to4(t) = poissrnd(MHR1(t)*migrate/3);
    Mhh1to4(t) = poissrnd(Mhh1(t)*migrate/3);
    MhR1to4(t) = poissrnd(MhR1(t)*migrate/3);
    MRR1to4(t) = poissrnd(MRR1(t)*migrate/3);

    MHH_HH1to4(t) = poissrnd(MHH_HH1(t)*migrate/3);
    MHH_Hh1to4(t) = poissrnd(MHH_Hh1(t)*migrate/3);
    MHH_HR1to4(t) = poissrnd(MHH_HR1(t)*migrate/3);
    MHH_hh1to4(t) = poissrnd(MHH_hh1(t)*migrate/3);
    MHH_hR1to4(t) = poissrnd(MHH_hR1(t)*migrate/3);
    MHH_RR1to4(t) = poissrnd(MHH_RR1(t)*migrate/3);

    MHh_HH1to4(t) = poissrnd(MHh_HH1(t)*migrate/3);
    MHh_Hh1to4(t) = poissrnd(MHh_Hh1(t)*migrate/3);
    MHh_HR1to4(t) = poissrnd(MHh_HR1(t)*migrate/3);
    MHh_hh1to4(t) = poissrnd(MHh_hh1(t)*migrate/3);
    MHh_hR1to4(t) = poissrnd(MHh_hR1(t)*migrate/3);
    MHh_RR1to4(t) = poissrnd(MHh_RR1(t)*migrate/3);

    MHR_HH1to4(t) = poissrnd(MHR_HH1(t)*migrate/3);
    MHR_Hh1to4(t) = poissrnd(MHR_Hh1(t)*migrate/3);
    MHR_HR1to4(t) = poissrnd(MHR_HR1(t)*migrate/3);
    MHR_hh1to4(t) = poissrnd(MHR_hh1(t)*migrate/3);
    MHR_hR1to4(t) = poissrnd(MHR_hR1(t)*migrate/3);
    MHR_RR1to4(t) = poissrnd(MHR_RR1(t)*migrate/3);

    Mhh_HH1to4(t) = poissrnd(Mhh_HH1(t)*migrate/3);
    Mhh_Hh1to4(t) = poissrnd(Mhh_Hh1(t)*migrate/3);
    Mhh_HR1to4(t) = poissrnd(Mhh_HR1(t)*migrate/3);
    Mhh_hh1to4(t) = poissrnd(Mhh_hh1(t)*migrate/3);
    Mhh_hR1to4(t) = poissrnd(Mhh_hR1(t)*migrate/3);
    Mhh_RR1to4(t) = poissrnd(Mhh_RR1(t)*migrate/3);

    MhR_HH1to4(t) = poissrnd(MhR_HH1(t)*migrate/3);
    MhR_Hh1to4(t) = poissrnd(MhR_Hh1(t)*migrate/3);
    MhR_HR1to4(t) = poissrnd(MhR_HR1(t)*migrate/3);
    MhR_hh1to4(t) = poissrnd(MhR_hh1(t)*migrate/3);
    MhR_hR1to4(t) = poissrnd(MhR_hR1(t)*migrate/3);
    MhR_RR1to4(t) = poissrnd(MhR_RR1(t)*migrate/3);

    MRR_HH1to4(t) = poissrnd(MRR_HH1(t)*migrate/3);
    MRR_Hh1to4(t) = poissrnd(MRR_Hh1(t)*migrate/3);
    MRR_HR1to4(t) = poissrnd(MRR_HR1(t)*migrate/3);
    MRR_hh1to4(t) = poissrnd(MRR_hh1(t)*migrate/3);
    MRR_hR1to4(t) = poissrnd(MRR_hR1(t)*migrate/3);
    MRR_RR1to4(t) = poissrnd(MRR_RR1(t)*migrate/3);

    % Migrants pop 2 -> pop 1:
    MHH2to1(t) = poissrnd(MHH2(t)*migrate/3);
    MHh2to1(t) = poissrnd(MHh2(t)*migrate/3);
    MHR2to1(t) = poissrnd(MHR2(t)*migrate/3);
    Mhh2to1(t) = poissrnd(Mhh2(t)*migrate/3);
    MhR2to1(t) = poissrnd(MhR2(t)*migrate/3);
    MRR2to1(t) = poissrnd(MRR2(t)*migrate/3);

    MHH_HH2to1(t) = poissrnd(MHH_HH2(t)*migrate/3);
    MHH_Hh2to1(t) = poissrnd(MHH_Hh2(t)*migrate/3);
    MHH_HR2to1(t) = poissrnd(MHH_HR2(t)*migrate/3);
    MHH_hh2to1(t) = poissrnd(MHH_hh2(t)*migrate/3);
    MHH_hR2to1(t) = poissrnd(MHH_hR2(t)*migrate/3);
    MHH_RR2to1(t) = poissrnd(MHH_RR2(t)*migrate/3);

    MHh_HH2to1(t) = poissrnd(MHh_HH2(t)*migrate/3);
    MHh_Hh2to1(t) = poissrnd(MHh_Hh2(t)*migrate/3);
    MHh_HR2to1(t) = poissrnd(MHh_HR2(t)*migrate/3);
    MHh_hh2to1(t) = poissrnd(MHh_hh2(t)*migrate/3);
    MHh_hR2to1(t) = poissrnd(MHh_hR2(t)*migrate/3);
    MHh_RR2to1(t) = poissrnd(MHh_RR2(t)*migrate/3);

    MHR_HH2to1(t) = poissrnd(MHR_HH2(t)*migrate/3);
    MHR_Hh2to1(t) = poissrnd(MHR_Hh2(t)*migrate/3);
    MHR_HR2to1(t) = poissrnd(MHR_HR2(t)*migrate/3);
    MHR_hh2to1(t) = poissrnd(MHR_hh2(t)*migrate/3);
    MHR_hR2to1(t) = poissrnd(MHR_hR2(t)*migrate/3);
    MHR_RR2to1(t) = poissrnd(MHR_RR2(t)*migrate/3);

    Mhh_HH2to1(t) = poissrnd(Mhh_HH2(t)*migrate/3);
    Mhh_Hh2to1(t) = poissrnd(Mhh_Hh2(t)*migrate/3);
    Mhh_HR2to1(t) = poissrnd(Mhh_HR2(t)*migrate/3);
    Mhh_hh2to1(t) = poissrnd(Mhh_hh2(t)*migrate/3);
    Mhh_hR2to1(t) = poissrnd(Mhh_hR2(t)*migrate/3);
    Mhh_RR2to1(t) = poissrnd(Mhh_RR2(t)*migrate/3);

    MhR_HH2to1(t) = poissrnd(MhR_HH2(t)*migrate/3);
    MhR_Hh2to1(t) = poissrnd(MhR_Hh2(t)*migrate/3);
    MhR_HR2to1(t) = poissrnd(MhR_HR2(t)*migrate/3);
    MhR_hh2to1(t) = poissrnd(MhR_hh2(t)*migrate/3);
    MhR_hR2to1(t) = poissrnd(MhR_hR2(t)*migrate/3);
    MhR_RR2to1(t) = poissrnd(MhR_RR2(t)*migrate/3);

    MRR_HH2to1(t) = poissrnd(MRR_HH2(t)*migrate/3);
    MRR_Hh2to1(t) = poissrnd(MRR_Hh2(t)*migrate/3);
    MRR_HR2to1(t) = poissrnd(MRR_HR2(t)*migrate/3);
    MRR_hh2to1(t) = poissrnd(MRR_hh2(t)*migrate/3);
    MRR_hR2to1(t) = poissrnd(MRR_hR2(t)*migrate/3);
    MRR_RR2to1(t) = poissrnd(MRR_RR2(t)*migrate/3);

    % Migrants pop 2 -> pop 3:
    MHH2to3(t) = poissrnd(MHH2(t)*migrate/3);
    MHh2to3(t) = poissrnd(MHh2(t)*migrate/3);
    MHR2to3(t) = poissrnd(MHR2(t)*migrate/3);
    Mhh2to3(t) = poissrnd(Mhh2(t)*migrate/3);
    MhR2to3(t) = poissrnd(MhR2(t)*migrate/3);
    MRR2to3(t) = poissrnd(MRR2(t)*migrate/3);

    MHH_HH2to3(t) = poissrnd(MHH_HH2(t)*migrate/3);
    MHH_Hh2to3(t) = poissrnd(MHH_Hh2(t)*migrate/3);
    MHH_HR2to3(t) = poissrnd(MHH_HR2(t)*migrate/3);
    MHH_hh2to3(t) = poissrnd(MHH_hh2(t)*migrate/3);
    MHH_hR2to3(t) = poissrnd(MHH_hR2(t)*migrate/3);
    MHH_RR2to3(t) = poissrnd(MHH_RR2(t)*migrate/3);

    MHh_HH2to3(t) = poissrnd(MHh_HH2(t)*migrate/3);
    MHh_Hh2to3(t) = poissrnd(MHh_Hh2(t)*migrate/3);
    MHh_HR2to3(t) = poissrnd(MHh_HR2(t)*migrate/3);
    MHh_hh2to3(t) = poissrnd(MHh_hh2(t)*migrate/3);
    MHh_hR2to3(t) = poissrnd(MHh_hR2(t)*migrate/3);
    MHh_RR2to3(t) = poissrnd(MHh_RR2(t)*migrate/3);

    MHR_HH2to3(t) = poissrnd(MHR_HH2(t)*migrate/3);
    MHR_Hh2to3(t) = poissrnd(MHR_Hh2(t)*migrate/3);
    MHR_HR2to3(t) = poissrnd(MHR_HR2(t)*migrate/3);
    MHR_hh2to3(t) = poissrnd(MHR_hh2(t)*migrate/3);
    MHR_hR2to3(t) = poissrnd(MHR_hR2(t)*migrate/3);
    MHR_RR2to3(t) = poissrnd(MHR_RR2(t)*migrate/3);

    Mhh_HH2to3(t) = poissrnd(Mhh_HH2(t)*migrate/3);
    Mhh_Hh2to3(t) = poissrnd(Mhh_Hh2(t)*migrate/3);
    Mhh_HR2to3(t) = poissrnd(Mhh_HR2(t)*migrate/3);
    Mhh_hh2to3(t) = poissrnd(Mhh_hh2(t)*migrate/3);
    Mhh_hR2to3(t) = poissrnd(Mhh_hR2(t)*migrate/3);
    Mhh_RR2to3(t) = poissrnd(Mhh_RR2(t)*migrate/3);

    MhR_HH2to3(t) = poissrnd(MhR_HH2(t)*migrate/3);
    MhR_Hh2to3(t) = poissrnd(MhR_Hh2(t)*migrate/3);
    MhR_HR2to3(t) = poissrnd(MhR_HR2(t)*migrate/3);
    MhR_hh2to3(t) = poissrnd(MhR_hh2(t)*migrate/3);
    MhR_hR2to3(t) = poissrnd(MhR_hR2(t)*migrate/3);
    MhR_RR2to3(t) = poissrnd(MhR_RR2(t)*migrate/3);

    MRR_HH2to3(t) = poissrnd(MRR_HH2(t)*migrate/3);
    MRR_Hh2to3(t) = poissrnd(MRR_Hh2(t)*migrate/3);
    MRR_HR2to3(t) = poissrnd(MRR_HR2(t)*migrate/3);
    MRR_hh2to3(t) = poissrnd(MRR_hh2(t)*migrate/3);
    MRR_hR2to3(t) = poissrnd(MRR_hR2(t)*migrate/3);
    MRR_RR2to3(t) = poissrnd(MRR_RR2(t)*migrate/3);

    % Migrants pop 2 -> pop 4:
    MHH2to4(t) = poissrnd(MHH2(t)*migrate/3);
    MHh2to4(t) = poissrnd(MHh2(t)*migrate/3);
    MHR2to4(t) = poissrnd(MHR2(t)*migrate/3);
    Mhh2to4(t) = poissrnd(Mhh2(t)*migrate/3);
    MhR2to4(t) = poissrnd(MhR2(t)*migrate/3);
    MRR2to4(t) = poissrnd(MRR2(t)*migrate/3);

    MHH_HH2to4(t) = poissrnd(MHH_HH2(t)*migrate/3);
    MHH_Hh2to4(t) = poissrnd(MHH_Hh2(t)*migrate/3);
    MHH_HR2to4(t) = poissrnd(MHH_HR2(t)*migrate/3);
    MHH_hh2to4(t) = poissrnd(MHH_hh2(t)*migrate/3);
    MHH_hR2to4(t) = poissrnd(MHH_hR2(t)*migrate/3);
    MHH_RR2to4(t) = poissrnd(MHH_RR2(t)*migrate/3);

    MHh_HH2to4(t) = poissrnd(MHh_HH2(t)*migrate/3);
    MHh_Hh2to4(t) = poissrnd(MHh_Hh2(t)*migrate/3);
    MHh_HR2to4(t) = poissrnd(MHh_HR2(t)*migrate/3);
    MHh_hh2to4(t) = poissrnd(MHh_hh2(t)*migrate/3);
    MHh_hR2to4(t) = poissrnd(MHh_hR2(t)*migrate/3);
    MHh_RR2to4(t) = poissrnd(MHh_RR2(t)*migrate/3);

    MHR_HH2to4(t) = poissrnd(MHR_HH2(t)*migrate/3);
    MHR_Hh2to4(t) = poissrnd(MHR_Hh2(t)*migrate/3);
    MHR_HR2to4(t) = poissrnd(MHR_HR2(t)*migrate/3);
    MHR_hh2to4(t) = poissrnd(MHR_hh2(t)*migrate/3);
    MHR_hR2to4(t) = poissrnd(MHR_hR2(t)*migrate/3);
    MHR_RR2to4(t) = poissrnd(MHR_RR2(t)*migrate/3);

    Mhh_HH2to4(t) = poissrnd(Mhh_HH2(t)*migrate/3);
    Mhh_Hh2to4(t) = poissrnd(Mhh_Hh2(t)*migrate/3);
    Mhh_HR2to4(t) = poissrnd(Mhh_HR2(t)*migrate/3);
    Mhh_hh2to4(t) = poissrnd(Mhh_hh2(t)*migrate/3);
    Mhh_hR2to4(t) = poissrnd(Mhh_hR2(t)*migrate/3);
    Mhh_RR2to4(t) = poissrnd(Mhh_RR2(t)*migrate/3);

    MhR_HH2to4(t) = poissrnd(MhR_HH2(t)*migrate/3);
    MhR_Hh2to4(t) = poissrnd(MhR_Hh2(t)*migrate/3);
    MhR_HR2to4(t) = poissrnd(MhR_HR2(t)*migrate/3);
    MhR_hh2to4(t) = poissrnd(MhR_hh2(t)*migrate/3);
    MhR_hR2to4(t) = poissrnd(MhR_hR2(t)*migrate/3);
    MhR_RR2to4(t) = poissrnd(MhR_RR2(t)*migrate/3);

    MRR_HH2to4(t) = poissrnd(MRR_HH2(t)*migrate/3);
    MRR_Hh2to4(t) = poissrnd(MRR_Hh2(t)*migrate/3);
    MRR_HR2to4(t) = poissrnd(MRR_HR2(t)*migrate/3);
    MRR_hh2to4(t) = poissrnd(MRR_hh2(t)*migrate/3);
    MRR_hR2to4(t) = poissrnd(MRR_hR2(t)*migrate/3);
    MRR_RR2to4(t) = poissrnd(MRR_RR2(t)*migrate/3);

    % Migrants pop 3 -> pop 1:
    MHH3to1(t) = poissrnd(MHH3(t)*migrate/3);
    MHh3to1(t) = poissrnd(MHh3(t)*migrate/3);
    MHR3to1(t) = poissrnd(MHR3(t)*migrate/3);
    Mhh3to1(t) = poissrnd(Mhh3(t)*migrate/3);
    MhR3to1(t) = poissrnd(MhR3(t)*migrate/3);
    MRR3to1(t) = poissrnd(MRR3(t)*migrate/3);

    MHH_HH3to1(t) = poissrnd(MHH_HH3(t)*migrate/3);
    MHH_Hh3to1(t) = poissrnd(MHH_Hh3(t)*migrate/3);
    MHH_HR3to1(t) = poissrnd(MHH_HR3(t)*migrate/3);
    MHH_hh3to1(t) = poissrnd(MHH_hh3(t)*migrate/3);
    MHH_hR3to1(t) = poissrnd(MHH_hR3(t)*migrate/3);
    MHH_RR3to1(t) = poissrnd(MHH_RR3(t)*migrate/3);

    MHh_HH3to1(t) = poissrnd(MHh_HH3(t)*migrate/3);
    MHh_Hh3to1(t) = poissrnd(MHh_Hh3(t)*migrate/3);
    MHh_HR3to1(t) = poissrnd(MHh_HR3(t)*migrate/3);
    MHh_hh3to1(t) = poissrnd(MHh_hh3(t)*migrate/3);
    MHh_hR3to1(t) = poissrnd(MHh_hR3(t)*migrate/3);
    MHh_RR3to1(t) = poissrnd(MHh_RR3(t)*migrate/3);

    MHR_HH3to1(t) = poissrnd(MHR_HH3(t)*migrate/3);
    MHR_Hh3to1(t) = poissrnd(MHR_Hh3(t)*migrate/3);
    MHR_HR3to1(t) = poissrnd(MHR_HR3(t)*migrate/3);
    MHR_hh3to1(t) = poissrnd(MHR_hh3(t)*migrate/3);
    MHR_hR3to1(t) = poissrnd(MHR_hR3(t)*migrate/3);
    MHR_RR3to1(t) = poissrnd(MHR_RR3(t)*migrate/3);

    Mhh_HH3to1(t) = poissrnd(Mhh_HH3(t)*migrate/3);
    Mhh_Hh3to1(t) = poissrnd(Mhh_Hh3(t)*migrate/3);
    Mhh_HR3to1(t) = poissrnd(Mhh_HR3(t)*migrate/3);
    Mhh_hh3to1(t) = poissrnd(Mhh_hh3(t)*migrate/3);
    Mhh_hR3to1(t) = poissrnd(Mhh_hR3(t)*migrate/3);
    Mhh_RR3to1(t) = poissrnd(Mhh_RR3(t)*migrate/3);

    MhR_HH3to1(t) = poissrnd(MhR_HH3(t)*migrate/3);
    MhR_Hh3to1(t) = poissrnd(MhR_Hh3(t)*migrate/3);
    MhR_HR3to1(t) = poissrnd(MhR_HR3(t)*migrate/3);
    MhR_hh3to1(t) = poissrnd(MhR_hh3(t)*migrate/3);
    MhR_hR3to1(t) = poissrnd(MhR_hR3(t)*migrate/3);
    MhR_RR3to1(t) = poissrnd(MhR_RR3(t)*migrate/3);

    MRR_HH3to1(t) = poissrnd(MRR_HH3(t)*migrate/3);
    MRR_Hh3to1(t) = poissrnd(MRR_Hh3(t)*migrate/3);
    MRR_HR3to1(t) = poissrnd(MRR_HR3(t)*migrate/3);
    MRR_hh3to1(t) = poissrnd(MRR_hh3(t)*migrate/3);
    MRR_hR3to1(t) = poissrnd(MRR_hR3(t)*migrate/3);
    MRR_RR3to1(t) = poissrnd(MRR_RR3(t)*migrate/3);

    % Migrants pop 3 -> pop 2:
    MHH3to2(t) = poissrnd(MHH3(t)*migrate/3);
    MHh3to2(t) = poissrnd(MHh3(t)*migrate/3);
    MHR3to2(t) = poissrnd(MHR3(t)*migrate/3);
    Mhh3to2(t) = poissrnd(Mhh3(t)*migrate/3);
    MhR3to2(t) = poissrnd(MhR3(t)*migrate/3);
    MRR3to2(t) = poissrnd(MRR3(t)*migrate/3);

    MHH_HH3to2(t) = poissrnd(MHH_HH3(t)*migrate/3);
    MHH_Hh3to2(t) = poissrnd(MHH_Hh3(t)*migrate/3);
    MHH_HR3to2(t) = poissrnd(MHH_HR3(t)*migrate/3);
    MHH_hh3to2(t) = poissrnd(MHH_hh3(t)*migrate/3);
    MHH_hR3to2(t) = poissrnd(MHH_hR3(t)*migrate/3);
    MHH_RR3to2(t) = poissrnd(MHH_RR3(t)*migrate/3);

    MHh_HH3to2(t) = poissrnd(MHh_HH3(t)*migrate/3);
    MHh_Hh3to2(t) = poissrnd(MHh_Hh3(t)*migrate/3);
    MHh_HR3to2(t) = poissrnd(MHh_HR3(t)*migrate/3);
    MHh_hh3to2(t) = poissrnd(MHh_hh3(t)*migrate/3);
    MHh_hR3to2(t) = poissrnd(MHh_hR3(t)*migrate/3);
    MHh_RR3to2(t) = poissrnd(MHh_RR3(t)*migrate/3);

    MHR_HH3to2(t) = poissrnd(MHR_HH3(t)*migrate/3);
    MHR_Hh3to2(t) = poissrnd(MHR_Hh3(t)*migrate/3);
    MHR_HR3to2(t) = poissrnd(MHR_HR3(t)*migrate/3);
    MHR_hh3to2(t) = poissrnd(MHR_hh3(t)*migrate/3);
    MHR_hR3to2(t) = poissrnd(MHR_hR3(t)*migrate/3);
    MHR_RR3to2(t) = poissrnd(MHR_RR3(t)*migrate/3);

    Mhh_HH3to2(t) = poissrnd(Mhh_HH3(t)*migrate/3);
    Mhh_Hh3to2(t) = poissrnd(Mhh_Hh3(t)*migrate/3);
    Mhh_HR3to2(t) = poissrnd(Mhh_HR3(t)*migrate/3);
    Mhh_hh3to2(t) = poissrnd(Mhh_hh3(t)*migrate/3);
    Mhh_hR3to2(t) = poissrnd(Mhh_hR3(t)*migrate/3);
    Mhh_RR3to2(t) = poissrnd(Mhh_RR3(t)*migrate/3);

    MhR_HH3to2(t) = poissrnd(MhR_HH3(t)*migrate/3);
    MhR_Hh3to2(t) = poissrnd(MhR_Hh3(t)*migrate/3);
    MhR_HR3to2(t) = poissrnd(MhR_HR3(t)*migrate/3);
    MhR_hh3to2(t) = poissrnd(MhR_hh3(t)*migrate/3);
    MhR_hR3to2(t) = poissrnd(MhR_hR3(t)*migrate/3);
    MhR_RR3to2(t) = poissrnd(MhR_RR3(t)*migrate/3);

    MRR_HH3to2(t) = poissrnd(MRR_HH3(t)*migrate/3);
    MRR_Hh3to2(t) = poissrnd(MRR_Hh3(t)*migrate/3);
    MRR_HR3to2(t) = poissrnd(MRR_HR3(t)*migrate/3);
    MRR_hh3to2(t) = poissrnd(MRR_hh3(t)*migrate/3);
    MRR_hR3to2(t) = poissrnd(MRR_hR3(t)*migrate/3);
    MRR_RR3to2(t) = poissrnd(MRR_RR3(t)*migrate/3);

    % Migrants pop 3 -> pop 4:
    MHH3to4(t) = poissrnd(MHH3(t)*migrate/3);
    MHh3to4(t) = poissrnd(MHh3(t)*migrate/3);
    MHR3to4(t) = poissrnd(MHR3(t)*migrate/3);
    Mhh3to4(t) = poissrnd(Mhh3(t)*migrate/3);
    MhR3to4(t) = poissrnd(MhR3(t)*migrate/3);
    MRR3to4(t) = poissrnd(MRR3(t)*migrate/3);

    MHH_HH3to4(t) = poissrnd(MHH_HH3(t)*migrate/3);
    MHH_Hh3to4(t) = poissrnd(MHH_Hh3(t)*migrate/3);
    MHH_HR3to4(t) = poissrnd(MHH_HR3(t)*migrate/3);
    MHH_hh3to4(t) = poissrnd(MHH_hh3(t)*migrate/3);
    MHH_hR3to4(t) = poissrnd(MHH_hR3(t)*migrate/3);
    MHH_RR3to4(t) = poissrnd(MHH_RR3(t)*migrate/3);

    MHh_HH3to4(t) = poissrnd(MHh_HH3(t)*migrate/3);
    MHh_Hh3to4(t) = poissrnd(MHh_Hh3(t)*migrate/3);
    MHh_HR3to4(t) = poissrnd(MHh_HR3(t)*migrate/3);
    MHh_hh3to4(t) = poissrnd(MHh_hh3(t)*migrate/3);
    MHh_hR3to4(t) = poissrnd(MHh_hR3(t)*migrate/3);
    MHh_RR3to4(t) = poissrnd(MHh_RR3(t)*migrate/3);

    MHR_HH3to4(t) = poissrnd(MHR_HH3(t)*migrate/3);
    MHR_Hh3to4(t) = poissrnd(MHR_Hh3(t)*migrate/3);
    MHR_HR3to4(t) = poissrnd(MHR_HR3(t)*migrate/3);
    MHR_hh3to4(t) = poissrnd(MHR_hh3(t)*migrate/3);
    MHR_hR3to4(t) = poissrnd(MHR_hR3(t)*migrate/3);
    MHR_RR3to4(t) = poissrnd(MHR_RR3(t)*migrate/3);

    Mhh_HH3to4(t) = poissrnd(Mhh_HH3(t)*migrate/3);
    Mhh_Hh3to4(t) = poissrnd(Mhh_Hh3(t)*migrate/3);
    Mhh_HR3to4(t) = poissrnd(Mhh_HR3(t)*migrate/3);
    Mhh_hh3to4(t) = poissrnd(Mhh_hh3(t)*migrate/3);
    Mhh_hR3to4(t) = poissrnd(Mhh_hR3(t)*migrate/3);
    Mhh_RR3to4(t) = poissrnd(Mhh_RR3(t)*migrate/3);

    MhR_HH3to4(t) = poissrnd(MhR_HH3(t)*migrate/3);
    MhR_Hh3to4(t) = poissrnd(MhR_Hh3(t)*migrate/3);
    MhR_HR3to4(t) = poissrnd(MhR_HR3(t)*migrate/3);
    MhR_hh3to4(t) = poissrnd(MhR_hh3(t)*migrate/3);
    MhR_hR3to4(t) = poissrnd(MhR_hR3(t)*migrate/3);
    MhR_RR3to4(t) = poissrnd(MhR_RR3(t)*migrate/3);

    MRR_HH3to4(t) = poissrnd(MRR_HH3(t)*migrate/3);
    MRR_Hh3to4(t) = poissrnd(MRR_Hh3(t)*migrate/3);
    MRR_HR3to4(t) = poissrnd(MRR_HR3(t)*migrate/3);
    MRR_hh3to4(t) = poissrnd(MRR_hh3(t)*migrate/3);
    MRR_hR3to4(t) = poissrnd(MRR_hR3(t)*migrate/3);
    MRR_RR3to4(t) = poissrnd(MRR_RR3(t)*migrate/3);

    % Migrants pop 4 -> pop 1:
    MHH4to1(t) = poissrnd(MHH4(t)*migrate/3);
    MHh4to1(t) = poissrnd(MHh4(t)*migrate/3);
    MHR4to1(t) = poissrnd(MHR4(t)*migrate/3);
    Mhh4to1(t) = poissrnd(Mhh4(t)*migrate/3);
    MhR4to1(t) = poissrnd(MhR4(t)*migrate/3);
    MRR4to1(t) = poissrnd(MRR4(t)*migrate/3);

    MHH_HH4to1(t) = poissrnd(MHH_HH4(t)*migrate/3);
    MHH_Hh4to1(t) = poissrnd(MHH_Hh4(t)*migrate/3);
    MHH_HR4to1(t) = poissrnd(MHH_HR4(t)*migrate/3);
    MHH_hh4to1(t) = poissrnd(MHH_hh4(t)*migrate/3);
    MHH_hR4to1(t) = poissrnd(MHH_hR4(t)*migrate/3);
    MHH_RR4to1(t) = poissrnd(MHH_RR4(t)*migrate/3);

    MHh_HH4to1(t) = poissrnd(MHh_HH4(t)*migrate/3);
    MHh_Hh4to1(t) = poissrnd(MHh_Hh4(t)*migrate/3);
    MHh_HR4to1(t) = poissrnd(MHh_HR4(t)*migrate/3);
    MHh_hh4to1(t) = poissrnd(MHh_hh4(t)*migrate/3);
    MHh_hR4to1(t) = poissrnd(MHh_hR4(t)*migrate/3);
    MHh_RR4to1(t) = poissrnd(MHh_RR4(t)*migrate/3);

    MHR_HH4to1(t) = poissrnd(MHR_HH4(t)*migrate/3);
    MHR_Hh4to1(t) = poissrnd(MHR_Hh4(t)*migrate/3);
    MHR_HR4to1(t) = poissrnd(MHR_HR4(t)*migrate/3);
    MHR_hh4to1(t) = poissrnd(MHR_hh4(t)*migrate/3);
    MHR_hR4to1(t) = poissrnd(MHR_hR4(t)*migrate/3);
    MHR_RR4to1(t) = poissrnd(MHR_RR4(t)*migrate/3);

    Mhh_HH4to1(t) = poissrnd(Mhh_HH4(t)*migrate/3);
    Mhh_Hh4to1(t) = poissrnd(Mhh_Hh4(t)*migrate/3);
    Mhh_HR4to1(t) = poissrnd(Mhh_HR4(t)*migrate/3);
    Mhh_hh4to1(t) = poissrnd(Mhh_hh4(t)*migrate/3);
    Mhh_hR4to1(t) = poissrnd(Mhh_hR4(t)*migrate/3);
    Mhh_RR4to1(t) = poissrnd(Mhh_RR4(t)*migrate/3);

    MhR_HH4to1(t) = poissrnd(MhR_HH4(t)*migrate/3);
    MhR_Hh4to1(t) = poissrnd(MhR_Hh4(t)*migrate/3);
    MhR_HR4to1(t) = poissrnd(MhR_HR4(t)*migrate/3);
    MhR_hh4to1(t) = poissrnd(MhR_hh4(t)*migrate/3);
    MhR_hR4to1(t) = poissrnd(MhR_hR4(t)*migrate/3);
    MhR_RR4to1(t) = poissrnd(MhR_RR4(t)*migrate/3);

    MRR_HH4to1(t) = poissrnd(MRR_HH4(t)*migrate/3);
    MRR_Hh4to1(t) = poissrnd(MRR_Hh4(t)*migrate/3);
    MRR_HR4to1(t) = poissrnd(MRR_HR4(t)*migrate/3);
    MRR_hh4to1(t) = poissrnd(MRR_hh4(t)*migrate/3);
    MRR_hR4to1(t) = poissrnd(MRR_hR4(t)*migrate/3);
    MRR_RR4to1(t) = poissrnd(MRR_RR4(t)*migrate/3);

    % Migrants pop 4 -> pop 2:
    MHH4to2(t) = poissrnd(MHH4(t)*migrate/3);
    MHh4to2(t) = poissrnd(MHh4(t)*migrate/3);
    MHR4to2(t) = poissrnd(MHR4(t)*migrate/3);
    Mhh4to2(t) = poissrnd(Mhh4(t)*migrate/3);
    MhR4to2(t) = poissrnd(MhR4(t)*migrate/3);
    MRR4to2(t) = poissrnd(MRR4(t)*migrate/3);

    MHH_HH4to2(t) = poissrnd(MHH_HH4(t)*migrate/3);
    MHH_Hh4to2(t) = poissrnd(MHH_Hh4(t)*migrate/3);
    MHH_HR4to2(t) = poissrnd(MHH_HR4(t)*migrate/3);
    MHH_hh4to2(t) = poissrnd(MHH_hh4(t)*migrate/3);
    MHH_hR4to2(t) = poissrnd(MHH_hR4(t)*migrate/3);
    MHH_RR4to2(t) = poissrnd(MHH_RR4(t)*migrate/3);

    MHh_HH4to2(t) = poissrnd(MHh_HH4(t)*migrate/3);
    MHh_Hh4to2(t) = poissrnd(MHh_Hh4(t)*migrate/3);
    MHh_HR4to2(t) = poissrnd(MHh_HR4(t)*migrate/3);
    MHh_hh4to2(t) = poissrnd(MHh_hh4(t)*migrate/3);
    MHh_hR4to2(t) = poissrnd(MHh_hR4(t)*migrate/3);
    MHh_RR4to2(t) = poissrnd(MHh_RR4(t)*migrate/3);

    MHR_HH4to2(t) = poissrnd(MHR_HH4(t)*migrate/3);
    MHR_Hh4to2(t) = poissrnd(MHR_Hh4(t)*migrate/3);
    MHR_HR4to2(t) = poissrnd(MHR_HR4(t)*migrate/3);
    MHR_hh4to2(t) = poissrnd(MHR_hh4(t)*migrate/3);
    MHR_hR4to2(t) = poissrnd(MHR_hR4(t)*migrate/3);
    MHR_RR4to2(t) = poissrnd(MHR_RR4(t)*migrate/3);

    Mhh_HH4to2(t) = poissrnd(Mhh_HH4(t)*migrate/3);
    Mhh_Hh4to2(t) = poissrnd(Mhh_Hh4(t)*migrate/3);
    Mhh_HR4to2(t) = poissrnd(Mhh_HR4(t)*migrate/3);
    Mhh_hh4to2(t) = poissrnd(Mhh_hh4(t)*migrate/3);
    Mhh_hR4to2(t) = poissrnd(Mhh_hR4(t)*migrate/3);
    Mhh_RR4to2(t) = poissrnd(Mhh_RR4(t)*migrate/3);

    MhR_HH4to2(t) = poissrnd(MhR_HH4(t)*migrate/3);
    MhR_Hh4to2(t) = poissrnd(MhR_Hh4(t)*migrate/3);
    MhR_HR4to2(t) = poissrnd(MhR_HR4(t)*migrate/3);
    MhR_hh4to2(t) = poissrnd(MhR_hh4(t)*migrate/3);
    MhR_hR4to2(t) = poissrnd(MhR_hR4(t)*migrate/3);
    MhR_RR4to2(t) = poissrnd(MhR_RR4(t)*migrate/3);

    MRR_HH4to2(t) = poissrnd(MRR_HH4(t)*migrate/3);
    MRR_Hh4to2(t) = poissrnd(MRR_Hh4(t)*migrate/3);
    MRR_HR4to2(t) = poissrnd(MRR_HR4(t)*migrate/3);
    MRR_hh4to2(t) = poissrnd(MRR_hh4(t)*migrate/3);
    MRR_hR4to2(t) = poissrnd(MRR_hR4(t)*migrate/3);
    MRR_RR4to2(t) = poissrnd(MRR_RR4(t)*migrate/3);

    % Migrants pop 4 -> pop 1:
    MHH4to3(t) = poissrnd(MHH4(t)*migrate/3);
    MHh4to3(t) = poissrnd(MHh4(t)*migrate/3);
    MHR4to3(t) = poissrnd(MHR4(t)*migrate/3);
    Mhh4to3(t) = poissrnd(Mhh4(t)*migrate/3);
    MhR4to3(t) = poissrnd(MhR4(t)*migrate/3);
    MRR4to3(t) = poissrnd(MRR4(t)*migrate/3);

    MHH_HH4to3(t) = poissrnd(MHH_HH4(t)*migrate/3);
    MHH_Hh4to3(t) = poissrnd(MHH_Hh4(t)*migrate/3);
    MHH_HR4to3(t) = poissrnd(MHH_HR4(t)*migrate/3);
    MHH_hh4to3(t) = poissrnd(MHH_hh4(t)*migrate/3);
    MHH_hR4to3(t) = poissrnd(MHH_hR4(t)*migrate/3);
    MHH_RR4to3(t) = poissrnd(MHH_RR4(t)*migrate/3);

    MHh_HH4to3(t) = poissrnd(MHh_HH4(t)*migrate/3);
    MHh_Hh4to3(t) = poissrnd(MHh_Hh4(t)*migrate/3);
    MHh_HR4to3(t) = poissrnd(MHh_HR4(t)*migrate/3);
    MHh_hh4to3(t) = poissrnd(MHh_hh4(t)*migrate/3);
    MHh_hR4to3(t) = poissrnd(MHh_hR4(t)*migrate/3);
    MHh_RR4to3(t) = poissrnd(MHh_RR4(t)*migrate/3);

    MHR_HH4to3(t) = poissrnd(MHR_HH4(t)*migrate/3);
    MHR_Hh4to3(t) = poissrnd(MHR_Hh4(t)*migrate/3);
    MHR_HR4to3(t) = poissrnd(MHR_HR4(t)*migrate/3);
    MHR_hh4to3(t) = poissrnd(MHR_hh4(t)*migrate/3);
    MHR_hR4to3(t) = poissrnd(MHR_hR4(t)*migrate/3);
    MHR_RR4to3(t) = poissrnd(MHR_RR4(t)*migrate/3);

    Mhh_HH4to3(t) = poissrnd(Mhh_HH4(t)*migrate/3);
    Mhh_Hh4to3(t) = poissrnd(Mhh_Hh4(t)*migrate/3);
    Mhh_HR4to3(t) = poissrnd(Mhh_HR4(t)*migrate/3);
    Mhh_hh4to3(t) = poissrnd(Mhh_hh4(t)*migrate/3);
    Mhh_hR4to3(t) = poissrnd(Mhh_hR4(t)*migrate/3);
    Mhh_RR4to3(t) = poissrnd(Mhh_RR4(t)*migrate/3);

    MhR_HH4to3(t) = poissrnd(MhR_HH4(t)*migrate/3);
    MhR_Hh4to3(t) = poissrnd(MhR_Hh4(t)*migrate/3);
    MhR_HR4to3(t) = poissrnd(MhR_HR4(t)*migrate/3);
    MhR_hh4to3(t) = poissrnd(MhR_hh4(t)*migrate/3);
    MhR_hR4to3(t) = poissrnd(MhR_hR4(t)*migrate/3);
    MhR_RR4to3(t) = poissrnd(MhR_RR4(t)*migrate/3);

    MRR_HH4to3(t) = poissrnd(MRR_HH4(t)*migrate/3);
    MRR_Hh4to3(t) = poissrnd(MRR_Hh4(t)*migrate/3);
    MRR_HR4to3(t) = poissrnd(MRR_HR4(t)*migrate/3);
    MRR_hh4to3(t) = poissrnd(MRR_hh4(t)*migrate/3);
    MRR_hR4to3(t) = poissrnd(MRR_hR4(t)*migrate/3);
    MRR_RR4to3(t) = poissrnd(MRR_RR4(t)*migrate/3);

    % Exchanging migrants:
    MHH1(t) = MHH1(t) - MHH1to2(t) - MHH1to3(t) - MHH1to4(t) + MHH2to1(t) + MHH3to1(t) + MHH4to1(t);
    MHh1(t) = MHh1(t) - MHh1to2(t) - MHh1to3(t) - MHh1to4(t) + MHh2to1(t) + MHh3to1(t) + MHh4to1(t);
    MHR1(t) = MHR1(t) - MHR1to2(t) - MHR1to3(t) - MHR1to4(t) + MHR2to1(t) + MHR3to1(t) + MHR4to1(t);
    Mhh1(t) = Mhh1(t) - Mhh1to2(t) - Mhh1to3(t) - Mhh1to4(t) + Mhh2to1(t) + Mhh3to1(t) + Mhh4to1(t);
    MhR1(t) = MhR1(t) - MhR1to2(t) - MhR1to3(t) - MhR1to4(t) + MhR2to1(t) + MhR3to1(t) + MhR4to1(t);
    MRR1(t) = MRR1(t) - MRR1to2(t) - MRR1to3(t) - MRR1to4(t) + MRR2to1(t) + MRR3to1(t) + MRR4to1(t);

    MHH_HH1(t) = MHH_HH1(t) - MHH_HH1to2(t) - MHH_HH1to3(t) - MHH_HH1to4(t) + MHH_HH2to1(t) + MHH_HH3to1(t) + MHH_HH4to1(t);
    MHH_Hh1(t) = MHH_Hh1(t) - MHH_Hh1to2(t) - MHH_Hh1to3(t) - MHH_Hh1to4(t) + MHH_Hh2to1(t) + MHH_Hh3to1(t) + MHH_Hh4to1(t);
    MHH_HR1(t) = MHH_HR1(t) - MHH_HR1to2(t) - MHH_HR1to3(t) - MHH_HR1to4(t) + MHH_HR2to1(t) + MHH_HR3to1(t) + MHH_HR4to1(t);
    MHH_hh1(t) = MHH_hh1(t) - MHH_hh1to2(t) - MHH_hh1to3(t) - MHH_hh1to4(t) + MHH_hh2to1(t) + MHH_hh3to1(t) + MHH_hh4to1(t);
    MHH_hR1(t) = MHH_hR1(t) - MHH_hR1to2(t) - MHH_hR1to3(t) - MHH_hR1to4(t) + MHH_hR2to1(t) + MHH_hR3to1(t) + MHH_hR4to1(t);
    MHH_RR1(t) = MHH_RR1(t) - MHH_RR1to2(t) - MHH_RR1to3(t) - MHH_RR1to4(t) + MHH_RR2to1(t) + MHH_RR3to1(t) + MHH_RR4to1(t);

    MHh_HH1(t) = MHh_HH1(t) - MHh_HH1to2(t) - MHh_HH1to3(t) - MHh_HH1to4(t) + MHh_HH2to1(t) + MHh_HH3to1(t) + MHh_HH4to1(t);
    MHh_Hh1(t) = MHh_Hh1(t) - MHh_Hh1to2(t) - MHh_Hh1to3(t) - MHh_Hh1to4(t) + MHh_Hh2to1(t) + MHh_Hh3to1(t) + MHh_Hh4to1(t);
    MHh_HR1(t) = MHh_HR1(t) - MHh_HR1to2(t) - MHh_HR1to3(t) - MHh_HR1to4(t) + MHh_HR2to1(t) + MHh_HR3to1(t) + MHh_HR4to1(t);
    MHh_hh1(t) = MHh_hh1(t) - MHh_hh1to2(t) - MHh_hh1to3(t) - MHh_hh1to4(t) + MHh_hh2to1(t) + MHh_hh3to1(t) + MHh_hh4to1(t);
    MHh_hR1(t) = MHh_hR1(t) - MHh_hR1to2(t) - MHh_hR1to3(t) - MHh_hR1to4(t) + MHh_hR2to1(t) + MHh_hR3to1(t) + MHh_hR4to1(t);
    MHh_RR1(t) = MHh_RR1(t) - MHh_RR1to2(t) - MHh_RR1to3(t) - MHh_RR1to4(t) + MHh_RR2to1(t) + MHh_RR3to1(t) + MHh_RR4to1(t);

    MHR_HH1(t) = MHR_HH1(t) - MHR_HH1to2(t) - MHR_HH1to3(t) - MHR_HH1to4(t) + MHR_HH2to1(t) + MHR_HH3to1(t) + MHR_HH4to1(t);
    MHR_Hh1(t) = MHR_Hh1(t) - MHR_Hh1to2(t) - MHR_Hh1to3(t) - MHR_Hh1to4(t) + MHR_Hh2to1(t) + MHR_Hh3to1(t) + MHR_Hh4to1(t);
    MHR_HR1(t) = MHR_HR1(t) - MHR_HR1to2(t) - MHR_HR1to3(t) - MHR_HR1to4(t) + MHR_HR2to1(t) + MHR_HR3to1(t) + MHR_HR4to1(t);
    MHR_hh1(t) = MHR_hh1(t) - MHR_hh1to2(t) - MHR_hh1to3(t) - MHR_hh1to4(t) + MHR_hh2to1(t) + MHR_hh3to1(t) + MHR_hh4to1(t);
    MHR_hR1(t) = MHR_hR1(t) - MHR_hR1to2(t) - MHR_hR1to3(t) - MHR_hR1to4(t) + MHR_hR2to1(t) + MHR_hR3to1(t) + MHR_hR4to1(t);
    MHR_RR1(t) = MHR_RR1(t) - MHR_RR1to2(t) - MHR_RR1to3(t) - MHR_RR1to4(t) + MHR_RR2to1(t) + MHR_RR3to1(t) + MHR_RR4to1(t);

    Mhh_HH1(t) = Mhh_HH1(t) - Mhh_HH1to2(t) - Mhh_HH1to3(t) - Mhh_HH1to4(t) + Mhh_HH2to1(t) + Mhh_HH3to1(t) + Mhh_HH4to1(t);
    Mhh_Hh1(t) = Mhh_Hh1(t) - Mhh_Hh1to2(t) - Mhh_Hh1to3(t) - Mhh_Hh1to4(t) + Mhh_Hh2to1(t) + Mhh_Hh3to1(t) + Mhh_Hh4to1(t);
    Mhh_HR1(t) = Mhh_HR1(t) - Mhh_HR1to2(t) - Mhh_HR1to3(t) - Mhh_HR1to4(t) + Mhh_HR2to1(t) + Mhh_HR3to1(t) + Mhh_HR4to1(t);
    Mhh_hh1(t) = Mhh_hh1(t) - Mhh_hh1to2(t) - Mhh_hh1to3(t) - Mhh_hh1to4(t) + Mhh_hh2to1(t) + Mhh_hh3to1(t) + Mhh_hh4to1(t);
    Mhh_hR1(t) = Mhh_hR1(t) - Mhh_hR1to2(t) - Mhh_hR1to3(t) - Mhh_hR1to4(t) + Mhh_hR2to1(t) + Mhh_hR3to1(t) + Mhh_hR4to1(t);
    Mhh_RR1(t) = Mhh_RR1(t) - Mhh_RR1to2(t) - Mhh_RR1to3(t) - Mhh_RR1to4(t) + Mhh_RR2to1(t) + Mhh_RR3to1(t) + Mhh_RR4to1(t);

    MhR_HH1(t) = MhR_HH1(t) - MhR_HH1to2(t) - MhR_HH1to3(t) - MhR_HH1to4(t) + MhR_HH2to1(t) + MhR_HH3to1(t) + MhR_HH4to1(t);
    MhR_Hh1(t) = MhR_Hh1(t) - MhR_Hh1to2(t) - MhR_Hh1to3(t) - MhR_Hh1to4(t) + MhR_Hh2to1(t) + MhR_Hh3to1(t) + MhR_Hh4to1(t);
    MhR_HR1(t) = MhR_HR1(t) - MhR_HR1to2(t) - MhR_HR1to3(t) - MhR_HR1to4(t) + MhR_HR2to1(t) + MhR_HR3to1(t) + MhR_HR4to1(t);
    MhR_hh1(t) = MhR_hh1(t) - MhR_hh1to2(t) - MhR_hh1to3(t) - MhR_hh1to4(t) + MhR_hh2to1(t) + MhR_hh3to1(t) + MhR_hh4to1(t);
    MhR_hR1(t) = MhR_hR1(t) - MhR_hR1to2(t) - MhR_hR1to3(t) - MhR_hR1to4(t) + MhR_hR2to1(t) + MhR_hR3to1(t) + MhR_hR4to1(t);
    MhR_RR1(t) = MhR_RR1(t) - MhR_RR1to2(t) - MhR_RR1to3(t) - MhR_RR1to4(t) + MhR_RR2to1(t) + MhR_RR3to1(t) + MhR_RR4to1(t);

    MRR_HH1(t) = MRR_HH1(t) - MRR_HH1to2(t) - MRR_HH1to3(t) - MRR_HH1to4(t) + MRR_HH2to1(t) + MRR_HH3to1(t) + MRR_HH4to1(t);
    MRR_Hh1(t) = MRR_Hh1(t) - MRR_Hh1to2(t) - MRR_Hh1to3(t) - MRR_Hh1to4(t) + MRR_Hh2to1(t) + MRR_Hh3to1(t) + MRR_Hh4to1(t);
    MRR_HR1(t) = MRR_HR1(t) - MRR_HR1to2(t) - MRR_HR1to3(t) - MRR_HR1to4(t) + MRR_HR2to1(t) + MRR_HR3to1(t) + MRR_HR4to1(t);
    MRR_hh1(t) = MRR_hh1(t) - MRR_hh1to2(t) - MRR_hh1to3(t) - MRR_hh1to4(t) + MRR_hh2to1(t) + MRR_hh3to1(t) + MRR_hh4to1(t);
    MRR_hR1(t) = MRR_hR1(t) - MRR_hR1to2(t) - MRR_hR1to3(t) - MRR_hR1to4(t) + MRR_hR2to1(t) + MRR_hR3to1(t) + MRR_hR4to1(t);
    MRR_RR1(t) = MRR_RR1(t) - MRR_RR1to2(t) - MRR_RR1to3(t) - MRR_RR1to4(t) + MRR_RR2to1(t) + MRR_RR3to1(t) + MRR_RR4to1(t);

    
    MHH2(t) = MHH2(t) - MHH2to1(t) - MHH2to3(t) - MHH2to4(t) + MHH1to2(t) + MHH3to2(t) + MHH4to2(t);
    MHh2(t) = MHh2(t) - MHh2to1(t) - MHh2to3(t) - MHh2to4(t) + MHh1to2(t) + MHh3to2(t) + MHh4to2(t);
    MHR2(t) = MHR2(t) - MHR2to1(t) - MHR2to3(t) - MHR2to4(t) + MHR1to2(t) + MHR3to2(t) + MHR4to2(t);
    Mhh2(t) = Mhh2(t) - Mhh2to1(t) - Mhh2to3(t) - Mhh2to4(t) + Mhh1to2(t) + Mhh3to2(t) + Mhh4to2(t);
    MhR2(t) = MhR2(t) - MhR2to1(t) - MhR2to3(t) - MhR2to4(t) + MhR1to2(t) + MhR3to2(t) + MhR4to2(t);
    MRR2(t) = MRR2(t) - MRR2to1(t) - MRR2to3(t) - MRR2to4(t) + MRR1to2(t) + MRR3to2(t) + MRR4to2(t);

    MHH_HH2(t) = MHH_HH2(t) - MHH_HH2to1(t) - MHH_HH2to3(t) - MHH_HH2to4(t) + MHH_HH1to2(t) + MHH_HH3to2(t) + MHH_HH4to2(t);
    MHH_Hh2(t) = MHH_Hh2(t) - MHH_Hh2to1(t) - MHH_Hh2to3(t) - MHH_Hh2to4(t) + MHH_Hh1to2(t) + MHH_Hh3to2(t) + MHH_Hh4to2(t);
    MHH_HR2(t) = MHH_HR2(t) - MHH_HR2to1(t) - MHH_HR2to3(t) - MHH_HR2to4(t) + MHH_HR1to2(t) + MHH_HR3to2(t) + MHH_HR4to2(t);
    MHH_hh2(t) = MHH_hh2(t) - MHH_hh2to1(t) - MHH_hh2to3(t) - MHH_hh2to4(t) + MHH_hh1to2(t) + MHH_hh3to2(t) + MHH_hh4to2(t);
    MHH_hR2(t) = MHH_hR2(t) - MHH_hR2to1(t) - MHH_hR2to3(t) - MHH_hR2to4(t) + MHH_hR1to2(t) + MHH_hR3to2(t) + MHH_hR4to2(t);
    MHH_RR2(t) = MHH_RR2(t) - MHH_RR2to1(t) - MHH_RR2to3(t) - MHH_RR2to4(t) + MHH_RR1to2(t) + MHH_RR3to2(t) + MHH_RR4to2(t);

    MHh_HH2(t) = MHh_HH2(t) - MHh_HH2to1(t) - MHh_HH2to3(t) - MHh_HH2to4(t) + MHh_HH1to2(t) + MHh_HH3to2(t) + MHh_HH4to2(t);
    MHh_Hh2(t) = MHh_Hh2(t) - MHh_Hh2to1(t) - MHh_Hh2to3(t) - MHh_Hh2to4(t) + MHh_Hh1to2(t) + MHh_Hh3to2(t) + MHh_Hh4to2(t);
    MHh_HR2(t) = MHh_HR2(t) - MHh_HR2to1(t) - MHh_HR2to3(t) - MHh_HR2to4(t) + MHh_HR1to2(t) + MHh_HR3to2(t) + MHh_HR4to2(t);
    MHh_hh2(t) = MHh_hh2(t) - MHh_hh2to1(t) - MHh_hh2to3(t) - MHh_hh2to4(t) + MHh_hh1to2(t) + MHh_hh3to2(t) + MHh_hh4to2(t);
    MHh_hR2(t) = MHh_hR2(t) - MHh_hR2to1(t) - MHh_hR2to3(t) - MHh_hR2to4(t) + MHh_hR1to2(t) + MHh_hR3to2(t) + MHh_hR4to2(t);
    MHh_RR2(t) = MHh_RR2(t) - MHh_RR2to1(t) - MHh_RR2to3(t) - MHh_RR2to4(t) + MHh_RR1to2(t) + MHh_RR3to2(t) + MHh_RR4to2(t);

    MHR_HH2(t) = MHR_HH2(t) - MHR_HH2to1(t) - MHR_HH2to3(t) - MHR_HH2to4(t) + MHR_HH1to2(t) + MHR_HH3to2(t) + MHR_HH4to2(t);
    MHR_Hh2(t) = MHR_Hh2(t) - MHR_Hh2to1(t) - MHR_Hh2to3(t) - MHR_Hh2to4(t) + MHR_Hh1to2(t) + MHR_Hh3to2(t) + MHR_Hh4to2(t);
    MHR_HR2(t) = MHR_HR2(t) - MHR_HR2to1(t) - MHR_HR2to3(t) - MHR_HR2to4(t) + MHR_HR1to2(t) + MHR_HR3to2(t) + MHR_HR4to2(t);
    MHR_hh2(t) = MHR_hh2(t) - MHR_hh2to1(t) - MHR_hh2to3(t) - MHR_hh2to4(t) + MHR_hh1to2(t) + MHR_hh3to2(t) + MHR_hh4to2(t);
    MHR_hR2(t) = MHR_hR2(t) - MHR_hR2to1(t) - MHR_hR2to3(t) - MHR_hR2to4(t) + MHR_hR1to2(t) + MHR_hR3to2(t) + MHR_hR4to2(t);
    MHR_RR2(t) = MHR_RR2(t) - MHR_RR2to1(t) - MHR_RR2to3(t) - MHR_RR2to4(t) + MHR_RR1to2(t) + MHR_RR3to2(t) + MHR_RR4to2(t);

    Mhh_HH2(t) = Mhh_HH2(t) - Mhh_HH2to1(t) - Mhh_HH2to3(t) - Mhh_HH2to4(t) + Mhh_HH1to2(t) + Mhh_HH3to2(t) + Mhh_HH4to2(t);
    Mhh_Hh2(t) = Mhh_Hh2(t) - Mhh_Hh2to1(t) - Mhh_Hh2to3(t) - Mhh_Hh2to4(t) + Mhh_Hh1to2(t) + Mhh_Hh3to2(t) + Mhh_Hh4to2(t);
    Mhh_HR2(t) = Mhh_HR2(t) - Mhh_HR2to1(t) - Mhh_HR2to3(t) - Mhh_HR2to4(t) + Mhh_HR1to2(t) + Mhh_HR3to2(t) + Mhh_HR4to2(t);
    Mhh_hh2(t) = Mhh_hh2(t) - Mhh_hh2to1(t) - Mhh_hh2to3(t) - Mhh_hh2to4(t) + Mhh_hh1to2(t) + Mhh_hh3to2(t) + Mhh_hh4to2(t);
    Mhh_hR2(t) = Mhh_hR2(t) - Mhh_hR2to1(t) - Mhh_hR2to3(t) - Mhh_hR2to4(t) + Mhh_hR1to2(t) + Mhh_hR3to2(t) + Mhh_hR4to2(t);
    Mhh_RR2(t) = Mhh_RR2(t) - Mhh_RR2to1(t) - Mhh_RR2to3(t) - Mhh_RR2to4(t) + Mhh_RR1to2(t) + Mhh_RR3to2(t) + Mhh_RR4to2(t);

    MhR_HH2(t) = MhR_HH2(t) - MhR_HH2to1(t) - MhR_HH2to3(t) - MhR_HH2to4(t) + MhR_HH1to2(t) + MhR_HH3to2(t) + MhR_HH4to2(t);
    MhR_Hh2(t) = MhR_Hh2(t) - MhR_Hh2to1(t) - MhR_Hh2to3(t) - MhR_Hh2to4(t) + MhR_Hh1to2(t) + MhR_Hh3to2(t) + MhR_Hh4to2(t);
    MhR_HR2(t) = MhR_HR2(t) - MhR_HR2to1(t) - MhR_HR2to3(t) - MhR_HR2to4(t) + MhR_HR1to2(t) + MhR_HR3to2(t) + MhR_HR4to2(t);
    MhR_hh2(t) = MhR_hh2(t) - MhR_hh2to1(t) - MhR_hh2to3(t) - MhR_hh2to4(t) + MhR_hh1to2(t) + MhR_hh3to2(t) + MhR_hh4to2(t);
    MhR_hR2(t) = MhR_hR2(t) - MhR_hR2to1(t) - MhR_hR2to3(t) - MhR_hR2to4(t) + MhR_hR1to2(t) + MhR_hR3to2(t) + MhR_hR4to2(t);
    MhR_RR2(t) = MhR_RR2(t) - MhR_RR2to1(t) - MhR_RR2to3(t) - MhR_RR2to4(t) + MhR_RR1to2(t) + MhR_RR3to2(t) + MhR_RR4to2(t);

    MRR_HH2(t) = MRR_HH2(t) - MRR_HH2to1(t) - MRR_HH2to3(t) - MRR_HH2to4(t) + MRR_HH1to2(t) + MRR_HH3to2(t) + MRR_HH4to2(t);
    MRR_Hh2(t) = MRR_Hh2(t) - MRR_Hh2to1(t) - MRR_Hh2to3(t) - MRR_Hh2to4(t) + MRR_Hh1to2(t) + MRR_Hh3to2(t) + MRR_Hh4to2(t);
    MRR_HR2(t) = MRR_HR2(t) - MRR_HR2to1(t) - MRR_HR2to3(t) - MRR_HR2to4(t) + MRR_HR1to2(t) + MRR_HR3to2(t) + MRR_HR4to2(t);
    MRR_hh2(t) = MRR_hh2(t) - MRR_hh2to1(t) - MRR_hh2to3(t) - MRR_hh2to4(t) + MRR_hh1to2(t) + MRR_hh3to2(t) + MRR_hh4to2(t);
    MRR_hR2(t) = MRR_hR2(t) - MRR_hR2to1(t) - MRR_hR2to3(t) - MRR_hR2to4(t) + MRR_hR1to2(t) + MRR_hR3to2(t) + MRR_hR4to2(t);
    MRR_RR2(t) = MRR_RR2(t) - MRR_RR2to1(t) - MRR_RR2to3(t) - MRR_RR2to4(t) + MRR_RR1to2(t) + MRR_RR3to2(t) + MRR_RR4to2(t);

    
    MHH3(t) = MHH3(t) - MHH3to1(t) - MHH3to2(t) - MHH3to4(t) + MHH1to3(t) + MHH2to3(t) + MHH4to3(t);
    MHh3(t) = MHh3(t) - MHh3to1(t) - MHh3to2(t) - MHh3to4(t) + MHh1to3(t) + MHh2to3(t) + MHh4to3(t);
    MHR3(t) = MHR3(t) - MHR3to1(t) - MHR3to2(t) - MHR3to4(t) + MHR1to3(t) + MHR2to3(t) + MHR4to3(t);
    Mhh3(t) = Mhh3(t) - Mhh3to1(t) - Mhh3to2(t) - Mhh3to4(t) + Mhh1to3(t) + Mhh2to3(t) + Mhh4to3(t);
    MhR3(t) = MhR3(t) - MhR3to1(t) - MhR3to2(t) - MhR3to4(t) + MhR1to3(t) + MhR2to3(t) + MhR4to3(t);
    MRR3(t) = MRR3(t) - MRR3to1(t) - MRR3to2(t) - MRR3to4(t) + MRR1to3(t) + MRR2to3(t) + MRR4to3(t);

    MHH_HH3(t) = MHH_HH3(t) - MHH_HH3to1(t) - MHH_HH3to2(t) - MHH_HH3to4(t) + MHH_HH1to3(t) + MHH_HH2to3(t) + MHH_HH4to3(t);
    MHH_Hh3(t) = MHH_Hh3(t) - MHH_Hh3to1(t) - MHH_Hh3to2(t) - MHH_Hh3to4(t) + MHH_Hh1to3(t) + MHH_Hh2to3(t) + MHH_Hh4to3(t);
    MHH_HR3(t) = MHH_HR3(t) - MHH_HR3to1(t) - MHH_HR3to2(t) - MHH_HR3to4(t) + MHH_HR1to3(t) + MHH_HR2to3(t) + MHH_HR4to3(t);
    MHH_hh3(t) = MHH_hh3(t) - MHH_hh3to1(t) - MHH_hh3to2(t) - MHH_hh3to4(t) + MHH_hh1to3(t) + MHH_hh2to3(t) + MHH_hh4to3(t);
    MHH_hR3(t) = MHH_hR3(t) - MHH_hR3to1(t) - MHH_hR3to2(t) - MHH_hR3to4(t) + MHH_hR1to3(t) + MHH_hR2to3(t) + MHH_hR4to3(t);
    MHH_RR3(t) = MHH_RR3(t) - MHH_RR3to1(t) - MHH_RR3to2(t) - MHH_RR3to4(t) + MHH_RR1to3(t) + MHH_RR2to3(t) + MHH_RR4to3(t);

    MHh_HH3(t) = MHh_HH3(t) - MHh_HH3to1(t) - MHh_HH3to2(t) - MHh_HH3to4(t) + MHh_HH1to3(t) + MHh_HH2to3(t) + MHh_HH4to3(t);
    MHh_Hh3(t) = MHh_Hh3(t) - MHh_Hh3to1(t) - MHh_Hh3to2(t) - MHh_Hh3to4(t) + MHh_Hh1to3(t) + MHh_Hh2to3(t) + MHh_Hh4to3(t);
    MHh_HR3(t) = MHh_HR3(t) - MHh_HR3to1(t) - MHh_HR3to2(t) - MHh_HR3to4(t) + MHh_HR1to3(t) + MHh_HR2to3(t) + MHh_HR4to3(t);
    MHh_hh3(t) = MHh_hh3(t) - MHh_hh3to1(t) - MHh_hh3to2(t) - MHh_hh3to4(t) + MHh_hh1to3(t) + MHh_hh2to3(t) + MHh_hh4to3(t);
    MHh_hR3(t) = MHh_hR3(t) - MHh_hR3to1(t) - MHh_hR3to2(t) - MHh_hR3to4(t) + MHh_hR1to3(t) + MHh_hR2to3(t) + MHh_hR4to3(t);
    MHh_RR3(t) = MHh_RR3(t) - MHh_RR3to1(t) - MHh_RR3to2(t) - MHh_RR3to4(t) + MHh_RR1to3(t) + MHh_RR2to3(t) + MHh_RR4to3(t);

    MHR_HH3(t) = MHR_HH3(t) - MHR_HH3to1(t) - MHR_HH3to2(t) - MHR_HH3to4(t) + MHR_HH1to3(t) + MHR_HH2to3(t) + MHR_HH4to3(t);
    MHR_Hh3(t) = MHR_Hh3(t) - MHR_Hh3to1(t) - MHR_Hh3to2(t) - MHR_Hh3to4(t) + MHR_Hh1to3(t) + MHR_Hh2to3(t) + MHR_Hh4to3(t);
    MHR_HR3(t) = MHR_HR3(t) - MHR_HR3to1(t) - MHR_HR3to2(t) - MHR_HR3to4(t) + MHR_HR1to3(t) + MHR_HR2to3(t) + MHR_HR4to3(t);
    MHR_hh3(t) = MHR_hh3(t) - MHR_hh3to1(t) - MHR_hh3to2(t) - MHR_hh3to4(t) + MHR_hh1to3(t) + MHR_hh2to3(t) + MHR_hh4to3(t);
    MHR_hR3(t) = MHR_hR3(t) - MHR_hR3to1(t) - MHR_hR3to2(t) - MHR_hR3to4(t) + MHR_hR1to3(t) + MHR_hR2to3(t) + MHR_hR4to3(t);
    MHR_RR3(t) = MHR_RR3(t) - MHR_RR3to1(t) - MHR_RR3to2(t) - MHR_RR3to4(t) + MHR_RR1to3(t) + MHR_RR2to3(t) + MHR_RR4to3(t);

    Mhh_HH3(t) = Mhh_HH3(t) - Mhh_HH3to1(t) - Mhh_HH1to3(t) - Mhh_HH3to4(t) + Mhh_HH1to3(t) + Mhh_HH2to3(t) + Mhh_HH4to3(t);
    Mhh_Hh3(t) = Mhh_Hh3(t) - Mhh_Hh3to1(t) - Mhh_Hh1to3(t) - Mhh_Hh3to4(t) + Mhh_Hh1to3(t) + Mhh_Hh2to3(t) + Mhh_Hh4to3(t);
    Mhh_HR3(t) = Mhh_HR3(t) - Mhh_HR3to1(t) - Mhh_HR1to3(t) - Mhh_HR3to4(t) + Mhh_HR1to3(t) + Mhh_HR2to3(t) + Mhh_HR4to3(t);
    Mhh_hh3(t) = Mhh_hh3(t) - Mhh_hh3to1(t) - Mhh_hh1to3(t) - Mhh_hh3to4(t) + Mhh_hh1to3(t) + Mhh_hh2to3(t) + Mhh_hh4to3(t);
    Mhh_hR3(t) = Mhh_hR3(t) - Mhh_hR3to1(t) - Mhh_hR1to3(t) - Mhh_hR3to4(t) + Mhh_hR1to3(t) + Mhh_hR2to3(t) + Mhh_hR4to3(t);
    Mhh_RR3(t) = Mhh_RR3(t) - Mhh_RR3to1(t) - Mhh_RR1to3(t) - Mhh_RR3to4(t) + Mhh_RR1to3(t) + Mhh_RR2to3(t) + Mhh_RR4to3(t);

    MhR_HH3(t) = MhR_HH3(t) - MhR_HH3to1(t) - MhR_HH1to3(t) - MhR_HH3to4(t) + MhR_HH1to3(t) + MhR_HH2to3(t) + MhR_HH4to3(t);
    MhR_Hh3(t) = MhR_Hh3(t) - MhR_Hh3to1(t) - MhR_Hh1to3(t) - MhR_Hh3to4(t) + MhR_Hh1to3(t) + MhR_Hh2to3(t) + MhR_Hh4to3(t);
    MhR_HR3(t) = MhR_HR3(t) - MhR_HR3to1(t) - MhR_HR1to3(t) - MhR_HR3to4(t) + MhR_HR1to3(t) + MhR_HR2to3(t) + MhR_HR4to3(t);
    MhR_hh3(t) = MhR_hh3(t) - MhR_hh3to1(t) - MhR_hh1to3(t) - MhR_hh3to4(t) + MhR_hh1to3(t) + MhR_hh2to3(t) + MhR_hh4to3(t);
    MhR_hR3(t) = MhR_hR3(t) - MhR_hR3to1(t) - MhR_hR1to3(t) - MhR_hR3to4(t) + MhR_hR1to3(t) + MhR_hR2to3(t) + MhR_hR4to3(t);
    MhR_RR3(t) = MhR_RR3(t) - MhR_RR3to1(t) - MhR_RR1to3(t) - MhR_RR3to4(t) + MhR_RR1to3(t) + MhR_RR2to3(t) + MhR_RR4to3(t);

    MRR_HH3(t) = MRR_HH3(t) - MRR_HH3to1(t) - MRR_HH1to3(t) - MRR_HH3to4(t) + MRR_HH1to3(t) + MRR_HH2to3(t) + MRR_HH4to3(t);
    MRR_Hh3(t) = MRR_Hh3(t) - MRR_Hh3to1(t) - MRR_Hh1to3(t) - MRR_Hh3to4(t) + MRR_Hh1to3(t) + MRR_Hh2to3(t) + MRR_Hh4to3(t);
    MRR_HR3(t) = MRR_HR3(t) - MRR_HR3to1(t) - MRR_HR1to3(t) - MRR_HR3to4(t) + MRR_HR1to3(t) + MRR_HR2to3(t) + MRR_HR4to3(t);
    MRR_hh3(t) = MRR_hh3(t) - MRR_hh3to1(t) - MRR_hh1to3(t) - MRR_hh3to4(t) + MRR_hh1to3(t) + MRR_hh2to3(t) + MRR_hh4to3(t);
    MRR_hR3(t) = MRR_hR3(t) - MRR_hR3to1(t) - MRR_hR1to3(t) - MRR_hR3to4(t) + MRR_hR1to3(t) + MRR_hR2to3(t) + MRR_hR4to3(t);
    MRR_RR3(t) = MRR_RR3(t) - MRR_RR3to1(t) - MRR_RR1to3(t) - MRR_RR3to4(t) + MRR_RR1to3(t) + MRR_RR2to3(t) + MRR_RR4to3(t);

 
    MHH4(t) = MHH4(t) - MHH4to1(t) - MHH4to2(t) - MHH4to3(t) + MHH1to4(t) + MHH2to4(t) + MHH3to4(t);
    MHh4(t) = MHh4(t) - MHh4to1(t) - MHh4to2(t) - MHh4to3(t) + MHh1to4(t) + MHh2to4(t) + MHh3to4(t);
    MHR4(t) = MHR4(t) - MHR4to1(t) - MHR4to2(t) - MHR4to3(t) + MHR1to4(t) + MHR2to4(t) + MHR3to4(t);
    Mhh4(t) = Mhh4(t) - Mhh4to1(t) - Mhh4to2(t) - Mhh4to3(t) + Mhh1to4(t) + Mhh2to4(t) + Mhh3to4(t);
    MhR4(t) = MhR4(t) - MhR4to1(t) - MhR4to2(t) - MhR4to3(t) + MhR1to4(t) + MhR2to4(t) + MhR3to4(t);
    MRR4(t) = MRR4(t) - MRR4to1(t) - MRR4to2(t) - MRR4to3(t) + MRR1to4(t) + MRR2to4(t) + MRR3to4(t);

    MHH_HH4(t) = MHH_HH4(t) - MHH_HH4to1(t) - MHH_HH4to2(t) - MHH_HH4to3(t) + MHH_HH1to4(t) + MHH_HH2to4(t) + MHH_HH3to4(t);
    MHH_Hh4(t) = MHH_Hh4(t) - MHH_Hh4to1(t) - MHH_Hh4to2(t) - MHH_Hh4to3(t) + MHH_Hh1to4(t) + MHH_Hh2to4(t) + MHH_Hh3to4(t);
    MHH_HR4(t) = MHH_HR4(t) - MHH_HR4to1(t) - MHH_HR4to2(t) - MHH_HR4to3(t) + MHH_HR1to4(t) + MHH_HR2to4(t) + MHH_HR3to4(t);
    MHH_hh4(t) = MHH_hh4(t) - MHH_hh4to1(t) - MHH_hh4to2(t) - MHH_hh4to3(t) + MHH_hh1to4(t) + MHH_hh2to4(t) + MHH_hh3to4(t);
    MHH_hR4(t) = MHH_hR4(t) - MHH_hR4to1(t) - MHH_hR4to2(t) - MHH_hR4to3(t) + MHH_hR1to4(t) + MHH_hR2to4(t) + MHH_hR3to4(t);
    MHH_RR4(t) = MHH_RR4(t) - MHH_RR4to1(t) - MHH_RR4to2(t) - MHH_RR4to3(t) + MHH_RR1to4(t) + MHH_RR2to4(t) + MHH_RR3to4(t);

    MHh_HH4(t) = MHh_HH4(t) - MHh_HH4to1(t) - MHh_HH4to2(t) - MHh_HH4to3(t) + MHh_HH1to4(t) + MHh_HH2to4(t) + MHh_HH3to4(t);
    MHh_Hh4(t) = MHh_Hh4(t) - MHh_Hh4to1(t) - MHh_Hh4to2(t) - MHh_Hh4to3(t) + MHh_Hh1to4(t) + MHh_Hh2to4(t) + MHh_Hh3to4(t);
    MHh_HR4(t) = MHh_HR4(t) - MHh_HR4to1(t) - MHh_HR4to2(t) - MHh_HR4to3(t) + MHh_HR1to4(t) + MHh_HR2to4(t) + MHh_HR3to4(t);
    MHh_hh4(t) = MHh_hh4(t) - MHh_hh4to1(t) - MHh_hh4to2(t) - MHh_hh4to3(t) + MHh_hh1to4(t) + MHh_hh2to4(t) + MHh_hh3to4(t);
    MHh_hR4(t) = MHh_hR4(t) - MHh_hR4to1(t) - MHh_hR4to2(t) - MHh_hR4to3(t) + MHh_hR1to4(t) + MHh_hR2to4(t) + MHh_hR3to4(t);
    MHh_RR4(t) = MHh_RR4(t) - MHh_RR4to1(t) - MHh_RR4to2(t) - MHh_RR4to3(t) + MHh_RR1to4(t) + MHh_RR2to4(t) + MHh_RR3to4(t);

    MHR_HH4(t) = MHR_HH4(t) - MHR_HH4to1(t) - MHR_HH4to2(t) - MHR_HH4to3(t) + MHR_HH1to4(t) + MHR_HH2to4(t) + MHR_HH3to4(t);
    MHR_Hh4(t) = MHR_Hh4(t) - MHR_Hh4to1(t) - MHR_Hh4to2(t) - MHR_Hh4to3(t) + MHR_Hh1to4(t) + MHR_Hh2to4(t) + MHR_Hh3to4(t);
    MHR_HR4(t) = MHR_HR4(t) - MHR_HR4to1(t) - MHR_HR4to2(t) - MHR_HR4to3(t) + MHR_HR1to4(t) + MHR_HR2to4(t) + MHR_HR3to4(t);
    MHR_hh4(t) = MHR_hh4(t) - MHR_hh4to1(t) - MHR_hh4to2(t) - MHR_hh4to3(t) + MHR_hh1to4(t) + MHR_hh2to4(t) + MHR_hh3to4(t);
    MHR_hR4(t) = MHR_hR4(t) - MHR_hR4to1(t) - MHR_hR4to2(t) - MHR_hR4to3(t) + MHR_hR1to4(t) + MHR_hR2to4(t) + MHR_hR3to4(t);
    MHR_RR4(t) = MHR_RR4(t) - MHR_RR4to1(t) - MHR_RR4to2(t) - MHR_RR4to3(t) + MHR_RR1to4(t) + MHR_RR2to4(t) + MHR_RR3to4(t);

    Mhh_HH4(t) = Mhh_HH4(t) - Mhh_HH4to1(t) - Mhh_HH4to2(t) - Mhh_HH4to3(t) + Mhh_HH1to4(t) + Mhh_HH2to4(t) + Mhh_HH3to4(t);
    Mhh_Hh4(t) = Mhh_Hh4(t) - Mhh_Hh4to1(t) - Mhh_Hh4to2(t) - Mhh_Hh4to3(t) + Mhh_Hh1to4(t) + Mhh_Hh2to4(t) + Mhh_Hh3to4(t);
    Mhh_HR4(t) = Mhh_HR4(t) - Mhh_HR4to1(t) - Mhh_HR4to2(t) - Mhh_HR4to3(t) + Mhh_HR1to4(t) + Mhh_HR2to4(t) + Mhh_HR3to4(t);
    Mhh_hh4(t) = Mhh_hh4(t) - Mhh_hh4to1(t) - Mhh_hh4to2(t) - Mhh_hh4to3(t) + Mhh_hh1to4(t) + Mhh_hh2to4(t) + Mhh_hh3to4(t);
    Mhh_hR4(t) = Mhh_hR4(t) - Mhh_hR4to1(t) - Mhh_hR4to2(t) - Mhh_hR4to3(t) + Mhh_hR1to4(t) + Mhh_hR2to4(t) + Mhh_hR3to4(t);
    Mhh_RR4(t) = Mhh_RR4(t) - Mhh_RR4to1(t) - Mhh_RR4to2(t) - Mhh_RR4to3(t) + Mhh_RR1to4(t) + Mhh_RR2to4(t) + Mhh_RR3to4(t);

    MhR_HH4(t) = MhR_HH4(t) - MhR_HH4to1(t) - MhR_HH4to2(t) - MhR_HH4to3(t) + MhR_HH1to4(t) + MhR_HH2to4(t) + MhR_HH3to4(t);
    MhR_Hh4(t) = MhR_Hh4(t) - MhR_Hh4to1(t) - MhR_Hh4to2(t) - MhR_Hh4to3(t) + MhR_Hh1to4(t) + MhR_Hh2to4(t) + MhR_Hh3to4(t);
    MhR_HR4(t) = MhR_HR4(t) - MhR_HR4to1(t) - MhR_HR4to2(t) - MhR_HR4to3(t) + MhR_HR1to4(t) + MhR_HR2to4(t) + MhR_HR3to4(t);
    MhR_hh4(t) = MhR_hh4(t) - MhR_hh4to1(t) - MhR_hh4to2(t) - MhR_hh4to3(t) + MhR_hh1to4(t) + MhR_hh2to4(t) + MhR_hh3to4(t);
    MhR_hR4(t) = MhR_hR4(t) - MhR_hR4to1(t) - MhR_hR4to2(t) - MhR_hR4to3(t) + MhR_hR1to4(t) + MhR_hR2to4(t) + MhR_hR3to4(t);
    MhR_RR4(t) = MhR_RR4(t) - MhR_RR4to1(t) - MhR_RR4to2(t) - MhR_RR4to3(t) + MhR_RR1to4(t) + MhR_RR2to4(t) + MhR_RR3to4(t);

    MRR_HH4(t) = MRR_HH4(t) - MRR_HH4to1(t) - MRR_HH4to2(t) - MRR_HH4to3(t) + MRR_HH1to4(t) + MRR_HH2to4(t) + MRR_HH3to4(t);
    MRR_Hh4(t) = MRR_Hh4(t) - MRR_Hh4to1(t) - MRR_Hh4to2(t) - MRR_Hh4to3(t) + MRR_Hh1to4(t) + MRR_Hh2to4(t) + MRR_Hh3to4(t);
    MRR_HR4(t) = MRR_HR4(t) - MRR_HR4to1(t) - MRR_HR4to2(t) - MRR_HR4to3(t) + MRR_HR1to4(t) + MRR_HR2to4(t) + MRR_HR3to4(t);
    MRR_hh4(t) = MRR_hh4(t) - MRR_hh4to1(t) - MRR_hh4to2(t) - MRR_hh4to3(t) + MRR_hh1to4(t) + MRR_hh2to4(t) + MRR_hh3to4(t);
    MRR_hR4(t) = MRR_hR4(t) - MRR_hR4to1(t) - MRR_hR4to2(t) - MRR_hR4to3(t) + MRR_hR1to4(t) + MRR_hR2to4(t) + MRR_hR3to4(t);
    MRR_RR4(t) = MRR_RR4(t) - MRR_RR4to1(t) - MRR_RR4to2(t) - MRR_RR4to3(t) + MRR_RR1to4(t) + MRR_RR2to4(t) + MRR_RR3to4(t);
 end

% Output:

t  = linspace(1,max_t,max_t);

% Vector densities:
Mm1 = MHH1 + MHh1 + MHR1 + Mhh1 + MhR1 + MRR1;
Mf1 = MHH_HH1 + MHH_Hh1+ MHH_HR1 + MHH_hh1 + MHH_hR1 + MHH_RR1 + MHh_HH1 + MHh_Hh1+ MHh_HR1 + MHh_hh1 + MHh_hR1+ MHh_RR1 + MHR_HH1 + MHR_Hh1 + MHR_HR1 + MHR_hh1 + MHR_hR1+ MHR_RR1 + Mhh_HH1 + Mhh_Hh1 + Mhh_HR1 + Mhh_hh1 + Mhh_hR1+ Mhh_RR1 + MhR_HH1 + MhR_Hh1 + MhR_HR1 + MhR_hh1 + MhR_hR1+ MhR_RR1 + MRR_HH1 + MRR_Hh1 + MRR_HR1 + MRR_hh1 + MRR_hR1+ MRR_RR1;
M1  = Mm1 + Mf1;
% Mtransgenic1 = Mm1 + Mf1 - Mhh1 - (Mhh_HH1 + Mhh_Hh1 + Mhh_HR1 + Mhh_hh1 + Mhh_hR1+ Mhh_RR1);
% Mf_transgenic1 = Mf1 - (Mhh_HH1 + Mhh_Hh1 + Mhh_HR1 + Mhh_hh1 + Mhh_hR1+ Mhh_RR1);
MhomingAllele1 = MHH1 + MHh1 + MHH_HH1 + MHH_Hh1 + MHH_HR1 + MHH_hh1 + MHH_hR1 + MHH_RR1 + MHh_HH1 + MHh_Hh1 + MHh_HR1 + MHh_hh1 + MHh_hR1 + MHh_RR1;
MresistantAllele1 = MHR1 + MhR1 + MRR1 + MHR_HH1 + MHR_Hh1 + MHR_HR1 + MHR_hh1 + MHR_hR1 + MHR_RR1 + MhR_HH1 + MhR_Hh1 + MhR_HR1 + MhR_hh1 + MhR_hR1+ MhR_RR1 + MRR_HH1 + MRR_Hh1 + MRR_HR1 + MRR_hh1 + MRR_hR1+ MRR_RR1;

Mm2 = MHH2 + MHh2 + MHR2 + Mhh2 + MhR2 + MRR2;
Mf2 = MHH_HH2 + MHH_Hh2+ MHH_HR2 + MHH_hh2 + MHH_hR2+ MHH_RR2 + MHh_HH2 + MHh_Hh2+ MHh_HR2 + MHh_hh2 + MHh_hR2+ MHh_RR2 + MHR_HH2 + MHR_Hh2 + MHR_HR2 + MHR_hh2 + MHR_hR2+ MHR_RR2 + Mhh_HH2 + Mhh_Hh2 + Mhh_HR2 + Mhh_hh2 + Mhh_hR2+ Mhh_RR2 + MhR_HH2 + MhR_Hh2 + MhR_HR2 + MhR_hh2 + MhR_hR2+ MhR_RR2 + MRR_HH2 + MRR_Hh2 + MRR_HR2 + MRR_hh2 + MRR_hR2+ MRR_RR2;
M2  = Mm2 + Mf2;
% Mtransgenic2 = Mm2 + Mf2 - Mhh2 - (Mhh_HH2 + Mhh_Hh2 + Mhh_HR2 + Mhh_hh2 + Mhh_hR2+ Mhh_RR2);
% Mf_transgenic2 = Mf2 - (Mhh_HH2 + Mhh_Hh2 + Mhh_HR2 + Mhh_hh2 + Mhh_hR2+ Mhh_RR2);
MhomingAllele2 = MHH2 + MHh2 + MHH_HH2 + MHH_Hh2 + MHH_HR2 + MHH_hh2 + MHH_hR2 + MHH_RR2 + MHh_HH2 + MHh_Hh2 + MHh_HR2 + MHh_hh2 + MHh_hR2 + MHh_RR2;
MresistantAllele2 = MHR2 + MhR2 + MRR2 + MHR_HH2 + MHR_Hh2 + MHR_HR2 + MHR_hh2 + MHR_hR2 + MHR_RR2 + MhR_HH2 + MhR_Hh2 + MhR_HR2 + MhR_hh2 + MhR_hR2 + MhR_RR2 + MRR_HH2 + MRR_Hh2 + MRR_HR2 + MRR_hh2 + MRR_hR2+ MRR_RR2;

Mm3 = MHH3 + MHh3 + MHR3 + Mhh3 + MhR3 + MRR3;
Mf3 = MHH_HH3 + MHH_Hh3+ MHH_HR3 + MHH_hh3 + MHH_hR3+ MHH_RR3 + MHh_HH3 + MHh_Hh3+ MHh_HR3 + MHh_hh3 + MHh_hR3+ MHh_RR3 + MHR_HH3 + MHR_Hh3 + MHR_HR3 + MHR_hh3 + MHR_hR3+ MHR_RR3 + Mhh_HH3 + Mhh_Hh3 + Mhh_HR3 + Mhh_hh3 + Mhh_hR3+ Mhh_RR3 + MhR_HH3 + MhR_Hh3 + MhR_HR3 + MhR_hh3 + MhR_hR3+ MhR_RR3 + MRR_HH3 + MRR_Hh3 + MRR_HR3 + MRR_hh3 + MRR_hR3+ MRR_RR3;
M3  = Mm3 + Mf3;
% Mtransgenic3 = Mm3 + Mf3 - Mhh3 - (Mhh_HH3 + Mhh_Hh3 + Mhh_HR3 + Mhh_hh3 + Mhh_hR3+ Mhh_RR3);
% Mf_transgenic3 = Mf3 - (Mhh_HH3 + Mhh_Hh3 + Mhh_HR3 + Mhh_hh3 + Mhh_hR3+ Mhh_RR3);
MhomingAllele3 = MHH3 + MHh3 + MHH_HH3 + MHH_Hh3 + MHH_HR3 + MHH_hh3 + MHH_hR3 + MHH_RR3 + MHh_HH3 + MHh_Hh3 + MHh_HR3 + MHh_hh3 + MHh_hR3 + MHh_RR3;
MresistantAllele3 = MHR3 + MhR3 + MRR3 + MHR_HH3 + MHR_Hh3 + MHR_HR3 + MHR_hh3 + MHR_hR3 + MHR_RR3 + MhR_HH3 + MhR_Hh3 + MhR_HR3 + MhR_hh3 + MhR_hR3 + MhR_RR3 + MRR_HH3 + MRR_Hh3 + MRR_HR3 + MRR_hh3 + MRR_hR3+ MRR_RR3;

Mm4 = MHH4 + MHh4 + MHR4 + Mhh4 + MhR4 + MRR4;
Mf4 = MHH_HH4 + MHH_Hh4+ MHH_HR4 + MHH_hh4 + MHH_hR4+ MHH_RR4 + MHh_HH4 + MHh_Hh4+ MHh_HR4 + MHh_hh4 + MHh_hR4+ MHh_RR4 + MHR_HH4 + MHR_Hh4 + MHR_HR4 + MHR_hh4 + MHR_hR4+ MHR_RR4 + Mhh_HH4 + Mhh_Hh4 + Mhh_HR4 + Mhh_hh4 + Mhh_hR4+ Mhh_RR4 + MhR_HH4 + MhR_Hh4 + MhR_HR4 + MhR_hh4 + MhR_hR4+ MhR_RR4 + MRR_HH4 + MRR_Hh4 + MRR_HR4 + MRR_hh4 + MRR_hR4+ MRR_RR4;
M4  = Mm4 + Mf4;
% Mtransgenic4 = Mm4 + Mf4 - Mhh4 - (Mhh_HH4 + Mhh_Hh4 + Mhh_HR4 + Mhh_hh4 + Mhh_hR4+ Mhh_RR4);
% Mf_transgenic4 = Mf4 - (Mhh_HH4 + Mhh_Hh4 + Mhh_HR4 + Mhh_hh4 + Mhh_hR4+ Mhh_RR4);
MhomingAllele4 = MHH4 + MHh4 + MHH_HH4 + MHH_Hh4 + MHH_HR4 + MHH_hh4 + MHH_hR4 + MHH_RR4 + MHh_HH4 + MHh_Hh4 + MHh_HR4 + MHh_hh4 + MHh_hR4 + MHh_RR4;
MresistantAllele4 = MHR4 + MhR4 + MRR4 + MHR_HH4 + MHR_Hh4 + MHR_HR4 + MHR_hh4 + MHR_hR4 + MHR_RR4 + MhR_HH4 + MhR_Hh4 + MhR_HR4 + MhR_hh4 + MhR_hR4 + MhR_RR4 + MRR_HH4 + MRR_Hh4 + MRR_HR4 + MRR_hh4 + MRR_hR4+ MRR_RR4;

Mm = Mm1 + Mm2 + Mm3 + Mm4;
Mf = Mf1 + Mf2 + Mf3 + Mf4;
M = M1 + M2 + M3 + M4;
MhomingAllele = MhomingAllele1 + MhomingAllele2 + MhomingAllele3 + MhomingAllele4;
MresistantAllele = MresistantAllele1 + MresistantAllele2 + MresistantAllele3 + MresistantAllele4;

finalArray = [MhomingAllele; MresistantAllele; M];
strM = 1000000000000;
idString = [int2str(floor(rho*strM)) '_' int2str(floor(e*strM)) '_' int2str(floor(s*strM)) '_' int2str(floor(M_eq*strM*4)) '_' fileID];
%csvwrite([path fileID '_head' '.csv'],[s,e,rho,M_eq]);
csvwrite([path idString '_data' '.csv'],finalArray);
end
