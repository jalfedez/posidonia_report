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
parameters.nmax= 20;
parameters.sl=750;            %Maximum Leaf Biomass g DW m−2
parameters.sr= 2500;
parameters.ksl=500;             %Leaf Biomass Growth Dependence on Space Availability g DW m−2
parameters.ksr=300;             %Below-ground Biomass Growth Dependence on Space Availability g DW m−2
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
parameters.conversion=240;               %upt: gN/(gN h) a mgN/(gDW day)




% 3. Set initial conditions
% Leaves biomass
L0 = 500;

% Below-ground biomass
R0 = 1600;

%Leaves and below-ground biomass internal nitrogen quota   
N0 = 10;

T = linspace(5,35,42);

%Uncoment the variable you want to plot

%%%%%%%%%%%%%%%%%%%%%%   Groml y Gromr
% 
% Groml_values=[0.001, 0.02, 0.05, 0.12, 0.15];
% Gromr_values=[0.0115, 0.013, 0.02, 0.05, 0.06];
% 
% 
% % Preallocate arrays to store results
% L_result = zeros(length(Groml_values), length(T));
% R_result = zeros(length(Gromr_values), length(T));


% %Gromr
% figure;
% hold on;
% 
% for e = 1:length(Gromr_values)
%     parameters.Groml = 0.012;
%     parameters.Gromr = Gromr_values(e);
%     for i = 1:length(T)
%         parameters.T = T(i);
%         [t, Y] = ode45(@(t, Y) model(t, Y, parameters), [0 10000], [L0, R0, N0]);
%         L_result(e, i) = Y(end, 1);
%     end
%     plot(T, L_result(e, :), 'LineWidth', 1.5);
% end
% 
% xlabel('Temperature (ºC)', 'FontSize', 12);
% ylabel('Leaf Biomass (g DW m^{-2})', 'FontSize', 12);
% set(gca, 'FontSize', 12); 
% legend(strcat('Gromr=', num2str(Gromr_values'),"1/d"), 'FontSize', 12, 'Location', 'Best');
% 
% % (Below-Ground Biomass)
% figure;
% hold on;
% 
% for e = 1:length(Gromr_values)
%     parameters.Groml = 0.012;
%     parameters.Gromr = Gromr_values(e);
%     for i = 1:length(T)
%         parameters.T = T(i);
%         [t, Y] = ode45(@(t, Y) model(t, Y, parameters), [0 10000], [L0, R0, N0]);
%         R_result(e, i) = Y(end, 2);
%     end
%     plot(T, R_result(e, :),'LineWidth', 1.5);
% end
% 
% xlabel('Temperature (ºC)', 'FontSize', 12);
% ylabel('Below-Ground Biomass (g DW m^{-2})', 'FontSize', 12);
% set(gca, 'FontSize', 12); 
% legend(strcat('Gromr=', num2str(Gromr_values'),"1/d"), 'FontSize', 12, 'Location', 'Best'); 
% 
% 
% 
% 
% % Groml
% figure;
% hold on;
%  
% for e = 1:length(Groml_values)
%     parameters.Groml = Groml_values(e);
%     parameters.Gromr = 0.0115;
%     for i = 1:length(T)
%         parameters.T = T(i);
%         [t, Y] = ode45(@(t, Y) model(t, Y, parameters), [0 10000], [L0, R0, N0]);
%         L_result(e, i) = Y(end, 1);
%     end
%     plot(T, L_result(e, :), 'LineWidth', 1.5);
% end
% 
% xlabel('Temperature (ºC)', 'FontSize', 12);
% ylabel('Leaf Biomass (g DW m^{-2})', 'FontSize', 12);
% set(gca, 'FontSize', 12); 
% legend(strcat('Groml=', num2str(Groml_values'),"1/d"), 'FontSize', 12, 'Location', 'Best'); 
% 
% % (Below-Ground Biomass)
% figure;
% hold on;
%  
% for e = 1:length(Groml_values)
%     parameters.Groml = Groml_values(e);
%     parameters.Gromr = 0.0115;
%     for i = 1:length(T)
%         parameters.T = T(i);
%         [t, Y] = ode45(@(t, Y) model(t, Y, parameters), [0 10000], [L0, R0, N0]);
%         R_result(e, i) = Y(end, 2);
%     end
%     plot(T, R_result(e, :),'LineWidth', 1.5);
% end
% 
% xlabel('Temperature (ºC)', 'FontSize', 12);
% ylabel('Below-Ground Biomass (g DW m^{-2})', 'FontSize', 12);
% set(gca, 'FontSize', 12); 
% legend(strcat('Groml=', num2str(Groml_values'),"1/d"), 'FontSize', 12, 'Location', 'Best');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Io
Io_values = [100, 150, 250, 350, 450];

L_result = zeros(length(Io_values), length(T));
R_result = zeros(length(Io_values), length(T));

%Io
figure;
hold on;
for e = 1:length(Io_values)
    parameters.Io = Io_values(e);
    for i = 1:length(T)
        parameters.T = T(i);
        [t, Y] = ode45(@(t, Y) model(t, Y, parameters), [0 10000], [L0, R0, N0]);
        L_result(e, i) = Y(end, 1);
    end
    plot(T, L_result(e, :), 'LineWidth', 1.5);
end
xlabel('Temperature (ºC)', 'FontSize', 12);
ylabel('Leaf Biomass (g DW m^{-2})', 'FontSize', 12);
set(gca, 'FontSize', 12); 
legend(strcat('Io=', num2str(Io_values'), ' µEm^{-2} s^{-1}'), 'FontSize', 12, 'Location', 'Best'); 

%(Below-Ground Biomass)
figure;
hold on;
for e = 1:length(Io_values)
    parameters.Io = Io_values(e);
    for i = 1:length(T)
        parameters.T = T(i);
        [t, Y] = ode45(@(t, Y) model(t, Y, parameters), [0 10000], [L0, R0, N0]);
        R_result(e, i) = Y(end, 2);
    end
    plot(T, R_result(e, :), 'LineWidth', 1.5);
end
xlabel('Temperature (ºC)', 'FontSize', 12);
ylabel('Below-Ground Biomass (g DW m^{-2})', 'FontSize', 12);
set(gca, 'FontSize', 12); 
legend(strcat('Io=', num2str(Io_values'), ' µEm^{-2} s^{-1}'), 'FontSize', 12, 'Location', 'Best'); 


%%%%%%%%%%%%%%%%%%%%%%%%%%% kz
% kz_values = [0.01, 0.05, 0.1, 0.15, 0.2];
% 
% L_result = zeros(length(kz_values), length(T));
% R_result = zeros(length(kz_values), length(T));
% 
% %(Leaf Biomass)
% figure;
% hold on;
% for e = 1:length(kz_values)
%     parameters.Groml = 0.012;
%     parameters.Gromr = 0.0115;
%     parameters.kz = kz_values(e);
%     for i = 1:length(T)
%         parameters.T = T(i);
%         [t, Y] = ode45(@(t, Y) model(t, Y, parameters), [0 10000], [L0, R0, N0]);
%         L_result(e, i) = Y(end, 1);
%     end
%     plot(T, L_result(e, :), 'LineWidth', 1.5);
% end
% xlabel('Temperature (ºC)', 'FontSize', 12);
% ylabel('Leaf Biomass (g DW m^{-2})', 'FontSize', 12);
% set(gca, 'FontSize', 12); 
% legend(strcat('kz=', num2str(kz_values'), ' m^{-1}'), 'FontSize', 12, 'Location', 'Best');
% 
% % (Below-Ground Biomass)
% figure;
% hold on;
% for e = 1:length(kz_values)
%     parameters.Groml = 0.012;
%     parameters.Gromr = 0.0115;
%     parameters.kz = kz_values(e);
%     for i = 1:length(T)
%         parameters.T = T(i);
%         [t, Y] = ode45(@(t, Y) model(t, Y, parameters), [0 10000], [L0, R0, N0]);
%         R_result(e, i) = Y(end, 2);
%     end
%     plot(T, R_result(e, :), 'LineWidth', 1.5);
% end
% xlabel('Temperature (ºC)', 'FontSize', 12);
% ylabel('Below-Ground Biomass (g DW m^{-2})', 'FontSize', 12);
% set(gca, 'FontSize', 12); 
% legend(strcat('kz=', num2str(kz_values'), ' m^{-1}'), 'FontSize', 12, 'Location', 'Best'); 
% 
% 

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


    % Calcular las funciones de temperatura para las hojas y las raíces
 
    fTl = exp(-((T-Toptl)/stl)^2);
    fTr = exp(-((T-Toptr)/str)^2);
    
    % Calculate light intensity
    Iz=Io*exp(-kz*d);
    fIz = Iz/(kel+Iz);
    
    %nutrientes
    fN=(N-nmin)/(ncrit-nmin);
    %fN = max(min(fN, 1), 0);

    %space
    
    fSl=1-exp(-((L-sl)/ksl).^2);
    fSr=1-exp(-((R-sr)/ksr).^2);
    fSl(L > sl) = 0.1;
    fSr(R > sr) = 0.1;


    %fSl=1-exp((L-sl)/ksl);
    %fSr=1-exp((R-sr)/ksr);
    
    
    % 
    % if L>sl
    %    fSl=0;
    % end
    % if R>sr
    %    fSr=0;
    % end
    


    Grol=Groml*fIz*fTl*fN*fSl;
    Gror=Gromr*fIz*fTr*fN*fSr;

    transr=ktrans*Gror;
    transl=1*Grol;

    Resl=resl*fTl;
    Resr=resr*fTr;

    fbN=(nmax-N)/(nmax-nmin);

    uptLNH4=VmlNH4*(NH4w/(NH4w+klNH4));
    uptLNO3=VmlNO3*(NO3w/(NO3w+klNO3));
    uptRNH4=VmrNH4*(NH4s/(NH4s+krNH4));

    uptL=(uptLNH4+uptLNO3);
    Upt=(uptL+uptRNH4)*fbN*conversion;

    % calculate dP , dD and dN
    
    dLdt = (Grol-transl-Resl)*L+transr*R;
    dRdt = (Gror-transr-Resr)*R;
    dNdt= Upt-Grol*N-Gror*N;
    
    if (sum(isnan([dLdt,dRdt,dNdt]))>0)
        keyboard
    end
    dYdt = [dLdt, dRdt, dNdt]';

end
