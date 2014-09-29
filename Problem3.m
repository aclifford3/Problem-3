%Constants
    maxTime = 100; %100 ms simulation
    time = [0:0.001:maxTime];
    timeLen = length(time);
    gbar_K = 36; %mS/cm^2
    gbar_Na = 120;
    gbar_L = 0.3;
    E_K = -12; %mV
    E_Na = 115;
    E_L = 10.6;
    V_rest = -70;
    C_m = 2;
    %Step size
    stepSize = .001;
    
    %Creates pulse current of 5 microA / cm^2 and applies it for 100ms
    I_in = [zeros(1, timeLen)];
    I_in(:) = 5;
    
    %Initial conditions
    V_m = [zeros(1, timeLen)];
    V_i = 0;    
    
    a_m = 0.1*((25-V_i)/(exp((25-V_i)/10)-1));
    b_m = 4*exp(-V_i/18);
    a_n = 0.01*((10-V_i)/(exp((10-V_i)/10)-1));
    b_n = .125*exp(-V_i/80);
    a_h = 0.07*exp(-V_i/20);
    b_h = 1/(exp((30-V_i)/10)+1);
    
    
    m = a_m / (a_m + b_m);
    n = a_n / (a_n + b_n);
    h = a_h / (a_h + b_h);
    
    g_Na = [m^2*gbar_Na*h zeros(1, timeLen-1)];
    g_K = [n^4 * gbar_K zeros(1,timeLen-1)];
    g_L = [gbar_L zeros(1, timeLen-1)];
   
    %Loop
    for(i = 2:timeLen)
       
        %Update gating conditions
        a_m = 0.1*((25-V_m(i-1))/(exp((25-V_m(i-1))/10)-1));
        b_m = 4*exp(-V_m(i-1)/18);
        a_n = 0.01*((10-V_m(i-1))/(exp((10-V_m(i-1))/10)-1));
        b_n = .125*exp(-V_m(i-1)/80);
        a_h = 0.07*exp(-V_m(i-1)/20);
        b_h = 1/(exp((30-V_m(i-1))/10)+1);
        
        %Calculate differentials of gating probabilities
        dm = a_m * (1-m) - b_m * m;
        dn = a_n * (1-n) - b_n * n;
        dh = a_h * (1-h) - b_h * h;
       
        %Update channel probabilities using Euler's method
        m = m + stepSize*dm;
        n = n + stepSize*dn;
        h = h + stepSize*dh;
        
        %Calculate currents
        I_Na = m^3*gbar_Na*h*(V_m(i-1) - E_Na);
        I_K = n^4 * gbar_K * (V_m(i-1) - E_K);
        I_L = gbar_L * (V_m(i-1) - E_L);
        I_ion = I_in(i) - I_K - I_Na - I_L;
        
        %Update ion channel conductances
        g_Na(i) = m.^3*gbar_Na*h;
        g_K(i) = n^4*gbar_K;
        g_L(i) = gbar_L;
        
        %Calculate membrance potential change
        dV_m = I_ion / C_m;
  
        
        %Update membrane potential
         V_m(i) = V_m(i-1) + stepSize * dV_m;
        
        
    end
    %Shift entire system by -70 mV to enter correct biological value range
    V_m = V_m + V_rest;
    %Graphs Voltage vs. Time
    figure();
    plot(time, V_m);
    axis([0,100, -100, 100]);
    title('Voltage vs. Time');
    xlabel('Time (ms)');
    ylabel('Membrane voltage (mV)');
    
    %Graphs Ion conductance vs. time
    figure();
    plot(time, g_Na, time, g_K);
    xlabel('Time (ms)');
    ylabel('Ion conductance (mS/cm^2)');
    title('Ion conductance vs. Time');
    
