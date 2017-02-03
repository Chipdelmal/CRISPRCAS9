function finalArray = iterateOneSimStep(rhoValue,eValue,sValue,MeqValue,fileID,path)
 
    s       = sValue;
    e       = eValue;       % homing rate
    rho     = rhoValue;     % NHEJ rate
    M_eq    = MeqValue;     % equilibrium adult population size
    [rho,e,s,M_eq];
    
    beta     = 32;              % eggs laid per day by female mosquito
    beta_HH  = 0;               % eggs laid per day by HH female mosquito
    beta_Hh  = beta*(1 - s);    % eggs laid per day by Hh female mosquito
    beta_HR  = beta*(1 - s);    % eggs laid per day by HR female mosquito
    beta_hh  = beta;            % eggs laid per day by hh female mosquito
    beta_hR  = beta;            % eggs laid per day by hR female mosquito
    beta_RR  = beta;            % eggs laid per day by RR female mosquito
    mu_m     = 0.123;           % adult mosquito daily mortality
    mu_m_HH  = mu_m;
    mu_m_Hh  = mu_m;
    mu_m_HR  = mu_m;
    mu_m_hh  = mu_m;
    mu_m_hR  = mu_m;
    mu_m_RR  = mu_m;
    rm       = 1.096;           % daily population growth rate
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
    % migrate = 0.01/g;
    alpha   = (beta*theta0*(M_eq/2)/(Rm-1))*((1-thetaL/Rm)/(1 - (thetaL/Rm)^(1/TL)));

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
    %%%%%%%%%%%%%%%%%%%
    %finalArray = [MHH1;MHh1;MHR1;Mhh1;MhR1;MRR1;M1];
    finalArray = [MhomingAllele1; MresistantAllele1; M1];
    strM = 1000000000000;
    idString = [int2str(floor(rho*strM)) '_' int2str(floor(e*strM)) '_' int2str(floor(s*strM)) '_' int2str(floor(M_eq*strM)) '_' fileID];
    %csvwrite([path fileID '_head' '.csv'],[s,e,rho,M_eq]);
    csvwrite([path idString '_data' '.csv'],finalArray);

end
