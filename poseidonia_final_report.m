% 1. Parameters struct per day
parameters.d = 10;          % depth

%Constants
parameters.resl=0.0038;       %Leaf Respiration Rate  day−1
parameters.resr= 0.0041;      %Below-ground Respiration Rate  day−1
parameters.transl= 0.005;     %Translocation Rate from Leaves  day−1
parameters.Groml=0.012;        %Leaf Maximum Specific Growth Rate Coefficient  day−1
parameters.Gromr=0.0115;      %Below-ground Specific Growth Rate Coefficient  day−1
parameters.kz=0.058;          %Water Column Light Attenuation Coefficient m−1
parameters.Io = 350;          %Surface Irradiance I0 µEm−2 s-1
parameters.kel= 225;          %Irradiance Half-saturation Constant for leaves kel µEm−2 s-1
parameters.Toptl=20;          %Optimum Temperature for Leaf Growth ◦C
parameters.Toptr=20;          %Optimum Temperature for Below-ground Growth ◦C
parameters.stl=3.6;           %Leaf Growth Temperature Dependence ◦C
parameters.str=3.6;           %Below-ground Growth Temperature Dependence str ◦C
parameters.nmin=4.28;         %Minimum Internal Nitrogen Quota mg N g−1 DW
parameters.ncrit=7.5;         %Critical Internal Nitrogen Quota mg N g−1 DW
parameters.nmax= 20;          %Maximum Internal Nitrogen Quota mg N g−1 DW
parameters.sl=750;            %Maximum Leaf Biomass g DW m−2
parameters.sr= 2500;          %Maximum Below-ground Biomass g DW m−2
parameters.ksl=500;           %Leaf Biomass Growth Dependence on Space Availability g DW m−2
parameters.ksr=300;           %Below-ground Biomass Growth Dependence on Space Availability g DW m−2
parameters.ktrans= 0.2;       %Translocation Coefficient to Above-ground Biomass 
parameters.VmlNH4=2e-4;       %Maximum Uptake Rate by leaves for NH4  g N g−1 h−1
parameters.VmlNO3=2e-4;       %Maximum Uptake Rate by leaves for NO3  g N g−1 h−1
parameters.VmrNH4=1.3e-3;     %Maximum Uptake Rate by below-ground biomass for NH4  g N g−1 h−1
parameters.klNH4=2.1e-5;      %Michaelis Constant for Water Ammonium g N l−1
parameters.klNO3= 3.01e-5;    %Michaelis Constant for Water Nitrate  g N l−1
parameters.krNH4=1.5e-4;      %Michaelis Constant for Sediment Ammonium  g N l−1

%cuadrar unidades
parameters.NH4w=100*parameters.klNH4;    
parameters.NO3w=100*parameters.klNO3;     
parameters.NH4s=100*parameters.krNH4;
parameters.conversion=240;   %upt: gN/(gN h) a mgN/(gDW day)




% 3. Set initial conditions
% Leaves biomass
L0 = 500;

% Below-ground biomass
R0 = 1600;

%Leaves and below-ground biomass internal nitrogen quota   
N0 = 10;

T = linspace(5,35,42);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Calculating values for growth limiting factors across the T range
fN_values_T = [];
fSl_values_T = [];
fSr_values_T =[];
fTr_values_T=[];
fTl_values_T=[];
fIz_values_T=[];

%fT
T_range = T; 
fTl = exp(-((T_range-parameters.Toptl)/parameters.stl).^2);
fTr = exp(-((T_range-parameters.Toptr)/parameters.str).^2);

% fIz
Iz=parameters.Io*exp(-parameters.kz*parameters.d);
fIz_values_T = Iz/(parameters.kel+Iz);
fIz_values_T = repmat(fIz_values_T, size(T_range));
fIz_values_T = fIz_values_T(:);

for i=1:length(T)
%     4. Solve ODE
      parameters.T=T(i);
      [t, Y] = ode45(@(t, Y) model(t, Y, parameters), [0 100000], [L0, R0, N0]);
      L_T(i) = Y(end,1);
      R_T(i) = Y(end,2);
      N_T(i) = Y(end,3);

    % Calculate fN for last value of N_T
    fN_T = (N_T(i) - parameters.nmin) / (parameters.ncrit - parameters.nmin);
    %fN_T = max(min(fN_T, 1), 0);
    fN_values_T = [fN_values_T; fN_T]; 

    % Calcualte dSl for last value of L_T
    fSl_T = 1 - exp((L_T(i) - parameters.sl) / parameters.ksl);
    fSl_values_T = [fSl_values_T; fSl_T]; 

    % Calculate fSr for last value of R_T
    fSr_T = 1 - exp((R_T(i) - parameters.sr) / parameters.ksr);
    fSr_values_T = [fSr_values_T; fSr_T]; 

    % Grol and Gror
    Gror(i) = 1 * fIz_values_T(i) * fTr(i) * fN_values_T(i) * fSr_values_T(i);
    Grol(i) = 1 * fIz_values_T(i) * fTl(i) * fN_values_T(i) * fSl_values_T(i);

end

%%% Leaf, Roots and Nutrients biomass
figure; 
subplot(1,2,1)
plot(T, L_T, 'color', [0 0.4470 0.7410], 'LineWidth', 2); 
hold on; 
plot(T, R_T, 'color',[0.4940 0.4940 0.4940], 'LineWidth', 2); 
xlabel('Temperature (ºC)'); 
ylabel('Biomass (g DW/m^2'); 
legend('Leaf Biomass', 'Root Biomass'); 
grid on; 

subplot(1,2,2)
plot(T, N_T, 'color', [0 0.4470 0.7410], 'LineWidth', 2); 
xlabel('Temperature (ºC)'); 
ylabel('Nitrogen concentration (mg N /mg^-1 DW)'); 
grid on; 

%%% Limiting factors plots
figure;  
plot(T, fSl_values_T, 'LineWidth', 2);
hold on;
plot(T, fN_values_T, 'LineWidth', 2);
hold on;
plot(T_range, fTl, 'LineWidth', 2);
hold on;
plot(T, fIz_values_T, 'LineWidth', 2);
hold on;
plot(T, Grol, 'LineWidth', 2);
xlabel('Temperature (ºC)'); 
%ylabel('Leaf growht rate (1/d)');
legend('fSl','fN', 'fTl','fIz', 'Grol (1/d)'); 
grid on; 

figure;  
plot(T, fSr_values_T, 'LineWidth', 2);
hold on;
plot(T, fN_values_T, 'LineWidth', 2);
hold on;
plot(T_range, fTr, 'LineWidth', 2);
hold on;
plot(T, fIz_values_T, 'LineWidth', 2);
hold on;
plot(T, Gror, 'LineWidth', 2);
xlabel('Temperature (ºC)'); 
ylabel('Root growth rate (1/d)'); 
legend('fSr','fN', 'fTr','fIz', 'Gror (1/d)'); 
grid on; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%L, R and N steady state for different temperatures 

%T_values = [15:1:20];
T_values = [20:1:25];
desired_length = 5000; 
t_interval = linspace(0, 8000, desired_length);

L_cells = cell(1, length(T_values));
R_cells = cell(1, length(T_values));
N_cells = cell(1, length(T_values));


for i = 1:length(T_values)
    parameters.T = T_values(i);
    [~, Y] = ode45(@(t, Y) model(t, Y, parameters), t_interval, [L0, R0, N0]);


    L_cells{i} = Y(:,1);
    R_cells{i} = Y(:,2);
    N_cells{i} = Y(:,3);
end

for i = 1:length(T_values)
    eval(sprintf('L_%02d = L_cells{%d};', T_values(i), i));
    eval(sprintf('R_%02d = R_cells{%d};', T_values(i), i));
    eval(sprintf('N_%02d = N_cells{%d};', T_values(i), i));
end

figure;
hold on;
for i = 1:length(T_values)
    plot(t_interval, L_cells{i}, '-', 'LineWidth', 2); 
   
end
xlabel('Time (days)');
ylabel('Leaf biomass (g DW/m^2)');
legend(arrayfun(@(x) sprintf('T = %d', x), T_values, 'UniformOutput', false));

figure;
hold on;
for i = 1:length(T_values)
    plot(t_interval, R_cells{i}, '-', 'LineWidth', 2); 
    
end
xlabel('Time (days)');
ylabel('Root biomass (g DW/m^2)');
legend(arrayfun(@(x) sprintf('T = %d', x), T_values, 'UniformOutput', false));


figure;
hold on;
for i = 1:length(T_values)
    plot(t_interval, N_cells{i}, '-', 'LineWidth', 2); 
 
end
xlabel('Time (days)');
ylabel('Nitrogen concentration (mg N /mg^-1 DW)');
legend(arrayfun(@(x) sprintf('T = %d', x), T_values, 'UniformOutput', false));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function dYdt = model(t, Y, parameters)
    d=parameters.d;          % depth
  
    %Constants
    resl=parameters.resl;       %Leaf Respiration Rate  day−1
    resr=parameters.resr;      %Below-ground Respiration Rate  day−1
    transl=parameters.transl;     %Translocation Rate from Leaves  day−1
    Groml=parameters.Groml;        %Leaf Maximum Specific Growth Rate Coefficient  day−1
    Gromr=parameters.Gromr;      %Below-ground Specific Growth Rate Coefficient  day−1
    kz=parameters.kz;          %Water Column Light Attenuation Coefficient m−1
    Io=parameters.Io;          %Surface Irradiance I0 µEm−2 s-1
    kel=parameters.kel;          %Irradiance Half-saturation Constant for leaves kel µEm−2 s-1
    Toptl=parameters.Toptl;          %Optimum Temperature for Leaf Growth ◦C
    Toptr=parameters.Toptr;          %Optimum Temperature for Below-ground Growth ◦C
    stl=parameters.stl;           %Leaf Growth Temperature Dependence ◦C
    str=parameters.str;           %Below-ground Growth Temperature Dependence str ◦C
    nmin=parameters.nmin;        %Minimum Internal Nitrogen Quota mg N g−1 DW
    ncrit=parameters.ncrit;        %Critical Internal Nitrogen Quota mg N g−1 DW
    nmax=parameters.nmax;     %Maximum Internal Nitrogen Quota mg N g−1 DW
    sl=parameters.sl;           %Maximum Leaf Biomass g DW m−2
    sr=parameters.sr;          %Maximum Below-ground Biomass g DW m−2
    ksl=parameters.ksl;            %Leaf Biomass Growth Dependence on Space Availability g DW m−2
    ksr=parameters.ksr;            %Below-ground Biomass Growth Dependence on Space Availability g DW m−2
    ktrans=parameters.ktrans;      %Translocation Coefficient to Above-ground Biomass 
    VmlNH4=parameters.VmlNH4;     %Maximum Uptake Rate by leaves for NH4  g N g−1 h−1
    VmlNO3=parameters.VmlNO3;     % Maximum Uptake Rate by leaves for NO3  g N g−1 h−1
    VmrNH4=parameters.VmrNH4;    % Maximum Uptake Rate by below-ground biomass for NH4  g N g−1 h−1
    klNH4=parameters.klNH4;    %Michaelis Constant for Water Ammonium g N l−1
    klNO3=parameters.klNO3;  %Michaelis Constant for Water Nitrate  g N l−1
    krNH4=parameters.krNH4;    %Michaelis Constant for Sediment Ammonium  g N l−1
    NH4w=parameters.NH4w;
    NH4s=parameters.NH4s;
    NO3w=parameters.NO3w;
    conversion=parameters.conversion;
    

    T=parameters.T;


    L = Y(1);
    R = Y(2);
    N= Y(3);


    % Limitations: 
    % Temperature
 
    fTl = exp(-((T-Toptl)/stl)^2);
    fTr = exp(-((T-Toptr)/str)^2);
    
    % Light intensity
    Iz=Io*exp(-kz*d);
    fIz = Iz/(kel+Iz);
    
    % Nutrientes
    fN=(N-nmin)/(ncrit-nmin);
    %fN = max(min(fN, 1), 0);

    % Space
    fSl=1-exp(-((L-sl)/ksl).^2);
    fSr=1-exp(-((R-sr)/ksr).^2);
    fSl(L > sl) = 0.1;
    fSr(R > sr) = 0.1;
    
    % Growth rates
    Grol=Groml*fIz*fTl*fN*fSl;
    Gror=Gromr*fIz*fTr*fN*fSr;

    % Translocation rates
    transr=ktrans*Gror;
    transl=1*Grol;

    % Respiration rates
    Resl=resl*fTl;
    Resr=resr*fTr;

    % Nitrogen
    fbN=(nmax-N)/(nmax-nmin);

    uptLNH4=VmlNH4*(NH4w/(NH4w+klNH4));
    uptLNO3=VmlNO3*(NO3w/(NO3w+klNO3));
    uptRNH4=VmrNH4*(NH4s/(NH4s+krNH4));

    uptL=(uptLNH4+uptLNO3);
    Upt=(uptL+uptRNH4)*fbN*conversion;

    % Calculate dP , dD and dN
    
    dLdt = (Grol-transl-Resl)*L+transr*R;
    dRdt = (Gror-transr-Resr)*R;
    dNdt= Upt-Grol*N-Gror*N;
    
    if (sum(isnan([dLdt,dRdt,dNdt]))>0)
        keyboard
    end
    dYdt = [dLdt, dRdt, dNdt]';

end
