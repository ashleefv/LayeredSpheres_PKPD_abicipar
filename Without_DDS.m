clc;
clear;

% Removing existing PNG files from the current directory 
pngFiles = dir('*.png');
for i = 1:length(pngFiles)
    delete(pngFiles(i).name);
end


species="human"; %can change to rabbit

figure_count = 1;
Data_time_at_target_vit_10 = []; 
Data_time_at_target_vit_50 = [];
Data_time_at_target_aq_10 = []; 
Data_time_at_target_aq_50 = [];

in_vivo_suppression_time = [];
in_vitro_suppression_time = [];

doseo = [0.01, .05, 0.1, 0.5]; %various doses to be tested in milligrams

% Time range
tinitial = 0; %days
tfinal=1000; %time (days)
time_interval = 10;
numberOfTimes = tfinal*time_interval;
t=linspace(tinitial, tfinal, numberOfTimes); %t is the time

for i = 1:length(doseo) %for loop for each dose

convf=1/34000*(10^3); %milligrams-to-micromoles conversion factor
dose=doseo(i)*convf; %convert to micromoles

time=0;
drug_dose = 0;
InitialDose=dose;
options = odeset('AbsTol',1e-10,'RelTol',1e-12); %options for ODE solver

if strcmp(species,"human")
    Ci=[1.02e-5,... %VEGFah, mcromolar
    4.37e-5,... %VEGFv, micromolar 
    0,... %DRah
    0,... %DRv
    0,... %Aah, all abicipar values are in micromoles
    InitialDose,... %Av
    0,... %Ar
    0,... %Ac
    0,]; %As  
    Vs=40495e-3; %L, human
    Vr=0.326e-3; %L, human
    Vc=0.139-3; %L, human
    Vv=4.4e-3; %L, human
    Vah=0.25e-3; %L, human
elseif strcmp(species,"rabbit")
    Ci=[1.43e-6,... %VEGFah, micromolar
    1.47e-7,... %VEGFv, micromolar 
    0,... %DRah
    0,... %DRv
    0,... %Aah, all abicipar values are in moles
    InitialDose,... %Av, micromole
    0,... %Ar
    0,... %Ac
    0,]; %As 
    Vs=2530e-3; %L, rabbit
    Vr=0.042e-3; %L,rabbit
    Vc=.0284e-3; %L, rabbit
    Vah=0.306e-3; %L, rabbit
    Vv=1.24e-3; %L, rabbit
end

soln = ode23s(@(t, y) ODEs(t, y, species, time, drug_dose), t, Ci, options);

VEGFah =deval(soln,t,1); 
VEGFv =deval(soln,t,2); 
DRah =deval(soln,t,3); 
DRv=deval(soln,t,4); 
Aah=deval(soln,t,5); 
Av=deval(soln,t,6); %6th vector has the Av values
Ar=deval(soln,t,7);
Ac=deval(soln,t,8);
As=deval(soln,t,9);

CAah=Aah./Vah; %microM
CAv=Av./Vv; %microM
CAr=Ar./Vr; %microM
CAc=Ac./Vc;%microM
CAs=As./Vs; %microM


microtonano=10^3;
microtopico=10^6;
UpperIC50=6*ones(1,length(t)); %nm
LowerIC50 =0.017*ones(1,length(t));%nm

%generates Vitreous drug concentration over time
figure(figure_count) 
semilogy(t,CAv.*microtonano,'LineWidth',2)
hold on
plot(t,UpperIC50,'m','LineWidth',2)
plot(t, LowerIC50,'b','LineWidth',2)
ylabel('Abicipar Concentration (nM)')
xlim([-5 500])
ylim([10^(-3)  10^(4)])
xlabel('Time (days)')
%legend('Vitreous Concentration','Upper IC_{50}','Lower IC_{50}')
title('Without DDS')
pbaspect([1 1 1])
axis square
set(gca, 'FontSize', 16)
exportgraphics(figure(figure_count),sprintf('drug_concentration_non_DDS%d.png', doseo(i)), 'Resolution', 300)
figure_count = figure_count+1;
hold off


% skipping first two days of data to avoid the effect of initial burst release
Index = 2*time_interval;
CAvnew = CAv(Index:end);
CAvnew = CAvnew.*microtonano;
UpperIC50 = UpperIC50(Index:end);
timenew = t(Index:end);
if (CAvnew(1)<UpperIC50(1))
index_upper_IC50 = 1;
time_at_upper_IC50 = 0;
else
index_upper_IC50 = find(CAvnew <= UpperIC50, 1);
time_at_upper_IC50 = timenew(index_upper_IC50);
end

%disp(index_upper_IC50)
%disp(time_at_upper_IC50)
in_vivo_suppression_time(i) = time_at_upper_IC50;

index_lower_IC50 = find(CAv.*microtonano <= LowerIC50, 1);
time_at_lower_IC50 = t(index_lower_IC50);

%disp(index_lower_IC50)
%disp(time_at_lower_IC50)
in_vitro_suppression_time(i) = time_at_lower_IC50;


%Calculate index of lowest value for VEGF 
[lowest_vvit, Index_vvit] = min(VEGFv);
[lowest_vaq, Index_vaq] = min(VEGFah);
Index_min = 10;

if lowest_vvit <= 0.5 * VEGFv(1) || lowest_vaq <= 0.5 * VEGFah(1)
    Index_min = min([Index_vvit, Index_vaq]);
else
    beep
    'Drug loading is not enough to reduce VEGF to 50% of its original value'
end

%disp(Index_min)
%disp(['The lowest value is:', num2str(lowest_vret), ' at index:', num2str(Index_vret)]);
%disp(['The lowest value is:', num2str(lowest_vvit), ' at index:', num2str(Index_vvit)]);
%disp(['The lowest value is:', num2str(lowest_vaq), ' at index:', num2str(Index_vaq)]);

%Calculates 10% Free VEGF Suppression Time
editedC_vvit = VEGFv(Index_min:end);
editedC_vaq = VEGFah(Index_min:end);
editedt = t(Index_min:end);

% Target concentrations for 10% suppression
target_concentration_vit_10 = 0.1 * VEGFv(1);
target_concentration_aq_10 = 0.1 * VEGFah(1);

% Find the times for 10% suppression
index_vit_10 = find(editedC_vvit >= target_concentration_vit_10, 1);
time_at_target_vit_10 = editedt(index_vit_10);
index_aq_10 = find(editedC_vaq >= target_concentration_aq_10, 1);
time_at_target_aq_10 = editedt(index_aq_10);
%fprintf('With DDS 10 percent VEGF suppression for the retina chamber is: %.2f\n With DDS 10 percent VEGF suppression for the vitreous chamber is: %.2f\n With DDS 10 percent VEGF suppression for the aqueous chamber is: %.2f\n', time_at_target_ret_10, time_at_target_vit_10, time_at_target_aq_10);

% Target concentrations for 50% suppression
target_concentration_vit_50 = 0.5 * VEGFv(1);
target_concentration_aq_50 = 0.5 * VEGFah(1);

% Find the times for 50% suppression
index_vit_50 = find(editedC_vvit >= target_concentration_vit_50, 1);
time_at_target_vit_50 = editedt(index_vit_50);
index_aq_50 = find(editedC_vaq >= target_concentration_aq_50, 1);
time_at_target_aq_50 = editedt(index_aq_50);
%fprintf('With DDS 50 percent VEGF suppression for the retina chamber is: %.2f\n With DDS 50 percent VEGF suppression for the vitreous chamber is: %.2f\n With DDS 50 percent VEGF suppression for the aqueous chamber is: %.2f\n', time_at_target_ret_50, time_at_target_vit_50, time_at_target_aq_50);

Data_time_at_target_vit_10(i) = time_at_target_vit_10; 
Data_time_at_target_vit_50(i) = time_at_target_vit_50;
Data_time_at_target_aq_10(i) = time_at_target_aq_10; 
Data_time_at_target_aq_50(i) = time_at_target_aq_50;


VEGFdoseaq(i, :) = VEGFah;
VEGFdosev(i, :) = VEGFv;
Drugaq(i, :) = (CAah).*microtonano;
Drugc(i,:) = (CAc).*microtonano;
Drugr(i,:) = (CAr).*microtonano;
Drugs(i,:) = (CAs).*microtonano;
Drugv(i,:) = (CAv).*microtonano;

end

% bar plot of 50% VEGF suppression
barWidth = 1;
figure(figure_count);
hold on
hBar1 = bar([Data_time_at_target_vit_50',Data_time_at_target_aq_50'], 'grouped', 'LineWidth', 1.5); 
xticks(1:length(doseo)); 
xticklabels({'0.01', '.05', '0.1', '0.5'});
ylabel('VEGF Suppression Time (Days)');
xlabel('Drug Amount (mg)');
set(hBar1, 'BarWidth', barWidth);
ylim([0  500])
%title(sprintf('Dose %0.2f mg',doseo(j)));
set(gca, 'FontSize', 14)
legend({'Vitreous 50%', 'Aqueous 50%'}, 'FontSize',14, 'Location', 'northwest');
pbaspect([1 1 1])
axis square
box on
exportgraphics(figure(figure_count),sprintf('bar_50_percent_VEGF.png'), 'Resolution', 300)
figure_count = figure_count+1;
hold off

% bar plot of 50% VEGF suppression
barWidth = 1;
figure(figure_count);
hold on
hBar1 = bar([Data_time_at_target_vit_10',Data_time_at_target_aq_10'], 'grouped', 'LineWidth', 1.5); 
xticks(1:length(doseo)); 
xticklabels({'0.01', '.05', '0.1', '0.5'});
ylabel('VEGF Suppression Time (Days)');
xlabel('Drug Amount (mg)');
set(hBar1, 'BarWidth', barWidth);
%title(sprintf('Dose %0.2f mg',doseo(j)));
ylim([0  500])
set(gca, 'FontSize', 14)
legend({'Vitreous 10%', 'Aqueous 10%'}, 'FontSize',14, 'Location', 'northwest');
pbaspect([1 1 1])
axis square
box on
exportgraphics(figure(figure_count),sprintf('bar_10_percent_VEGF.png'), 'Resolution', 300)
figure_count = figure_count+1;
hold off


% bar plot IC50
barWidth = 1;
bar_graph_data = [in_vivo_suppression_time; in_vitro_suppression_time];
% Dose response
figure(figure_count);
hold on
hBar1 = bar(bar_graph_data', 'grouped', 'LineWidth', 1.5); 
xticks(1:length(doseo)); 
xticklabels({'0.01', '.05', '0.1', '0.5'}); 
ylabel('Time (days)');
xlabel('Drug Amount (mg)');
ylim([0  500])
set(hBar1, 'BarWidth', barWidth);
title('Without DDS (Vitreous)');
set(gca, 'FontSize', 14)
legend({'Upper IC_{50}', 'Lower IC_{50}'}, 'FontSize',14, 'Location', 'northwest'); % Legend for rows
pbaspect([1 1 1])
axis square
box on
exportgraphics(figure(figure_count),sprintf('bar_dose_response_non_DDS.png'), 'Resolution', 300)
figure_count = figure_count+1;
hold off


%VEGF Concentration aqueous with different doses
figure(figure_count)
set(gca, 'FontSize', 14)
hold on
box on
for j = 1:length(doseo)
    plot(t, (VEGFdoseaq(j, :)).*microtopico, 'LineWidth', 2)
end
%legend('Dose = 0.01 mg', 'Dose = 0.05 mg', 'Dose = 0.1 mg', 'Dose = 0.5 mg', 'FontSize',12, 'Location','southeast');
legend('Dose 0.01 mg', 'Dose 0.05 mg', 'Dose 0.1 mg', 'Dose 0.5 mg', 'FontSize',12, 'Location','southeast');
xlabel('Time (days)');
ylabel('VEGF Concentration (pM)');
%xlim([-5 500])
title('VEGF Aqueous');
%fontsize(figure(3), 17, "points")
pbaspect([1 1 1])
axis square
exportgraphics(figure(figure_count),sprintf('VEGFaq_non_DDS.png'), 'Resolution', 300)
figure_count = figure_count+1;
hold off


%VEGF Concentration vitreous with different doses
figure(figure_count)
set(gca, 'FontSize', 14)
hold on
box on
for j = 1:length(doseo)
    plot(t, (VEGFdosev(j, :)).*microtopico, 'LineWidth', 2)
end
%legend('Dose = 0.01mg', 'Dose = 0.05mg', 'Dose = 0.1mg', 'Dose = 0.5mg', 'Location','southeast');
legend('Dose 0.01 mg', 'Dose 0.05 mg', 'Dose 0.1 mg', 'Dose 0.5 mg', 'FontSize',12, 'Location','southeast');
xlabel('Time (days)');
ylabel('VEGF Concentration (pM)');
%xlim([-5 500])
title('VEGF Vitreous');
%fontsize(figure(4), 17, "points")
pbaspect([1 1 1])
axis square
exportgraphics(figure(figure_count),sprintf('VEGFv_non_DDS.png'), 'Resolution', 300)
figure_count = figure_count+1;
hold off

save("Without_DDS.mat","VEGFdoseaq","VEGFdosev","Drugaq","Drugc","Drugr","Drugs","Drugv","Data_time_at_target_vit_10","Data_time_at_target_vit_50", "Data_time_at_target_aq_10", "Data_time_at_target_aq_50", "in_vivo_suppression_time","in_vitro_suppression_time")


  function derivVector=ODEs(t,y,species,time,drug_dose) %must be (denominator, numerators)
   
   % Note: this entire section depends on the constants and variables and
   % differential equations and is customized for each problem.
   
    % Set parameters also known as rate constants. Use comments to indicate
    % the units for each after giving the values
 if strcmp(species,"human")
    kd=0.911e-6; %Human, micromolar
    kvrfactor=1;
    kvr=11.5*kvrfactor; %human, day^-1
    kva=0.0256; %human, day^-1
    krv=609; %human, day^-1
    krc=25700; %human, day^-1
    krs=4620; %human, day^-1
    kcr=35500; %human, day^-1
    kcs=6460; %human, day^-1
    kas=0.402; %human, day^-1
    ksc=18.8; %human, day^-1
    ks=0.0563; %human, day^-1
    kint=9.25; %human, day^-1. %This parameter was changed from 55.4 to 9.25
    Vs=40495e-3; %L, human
    Vr=0.326e-3; %L, human
    Vc=0.139-3; %L, human
    Vv=4.4e-3; %L, human
    Vah=0.25e-3; %L, human
    baselineV=4.37e-5; %micromolar
    baselineAH=1.02e-5; %micromolar, end of human parameters

 elseif strcmp(species, "rabbit")
    kvr=3.37; %day^-1
    kva=.00217;
    krv=385;
    krc=34300;
    krs=869;
    kcr=30400;
    kcs=466;
    kas=.253;
    ksc=3.56;
    ks=.155;
    kint=404;
    baselineV=1.47e-7; %microMolar
    baselineAH=1.43e-6; %microMolar
    kd=4.34e-6; %microM
    Vv=1.24e-3; %L
    Vah=.306e-3; %L
 end
 
 kdeg=log(2)/2.46 * 24; %day^-1
 kon=8 * 86400; %1/(microM*day)
 koff=kd*kon;
 ksynV=kdeg*baselineV; %micromol/L*day (microM/day)
 ksynAH=kdeg*baselineAH; %,micromol/L*day (microM/day)
 
    VEGFah = y(1);
    VEGFv = y(2);
    DRah =y(3);
    DRv=y(4);
    Aah=y(5);
    Av=y(6);
    Ar=y(7);
    Ac=y(8);
    As=y(9);

    %7 ODEs
    dVEGFahdt=ksynAH-(kdeg*VEGFah)+ (koff*DRah) - (kon*VEGFah*Aah/Vah) ; %microM/day
    dVEGFvdt=ksynV-(kdeg*VEGFv)-(kon*VEGFv*Av/Vv) + (koff*DRv); %microM/day
    dDRahdt=(kon*VEGFah*Aah/Vah)-DRah*(koff+kint); %microM/day
    dDRvdt=(kon*VEGFv*Av/Vv)-DRv*(koff+kint); %microM/day
    dAahdt= -kon*VEGFah*Aah + koff*DRah*Vah -kas*Aah + kva*Av; %micromol/day
    DDSTerm=0; 
    dAvdt=-kon*VEGFv*Av + koff*DRv*Vv + krv*Ar - kvr*Av - kva*Av + DDSTerm; %micromol/day
    dArdt=kvr*Av - krv*Ar + kcr*Ac - krs*Ar - krc*Ar; %micromol/day
    dAcdt=krc*Ar - kcr*Ac - kcs*Ac + ksc*As; %micromol/day
    dAsdt=kas*Aah + krs*Ar - ks*As + kcs*Ac -ksc*As; %micromol/day

    derivVector = [dVEGFahdt,dVEGFvdt,dDRahdt,dDRvdt,dAahdt, dAvdt, dArdt, dAcdt, dAsdt]'; %The ' is very important as it properly shapes the vector for the ODE solver
    end