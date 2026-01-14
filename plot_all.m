clc;
clear;


% Time range
tinitial = 0; %days
tfinal=1000; %time (days)
time_interval = 10;
numberOfTimes = tfinal*time_interval;
t=linspace(tinitial, tfinal, numberOfTimes); %t is the time

% conversion_factor
microtonano=10^3;
microtopico=10^6;

doseo = [0.01, .05, 0.1, 0.5]; %various doses to be tested in milligrams


%%%%%%%%%%%%%%%% figure 2 %%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(2);
figname = 'figure2';
figure_count = 1;

load ("Without_DDS.mat")

for i = 1:length(doseo) %for loop for each dose

subplot(2,4,figure_count) 
% Vitreous abicipar concentration
UpperIC50=6*ones(1,length(t)); %nm
LowerIC50 =0.017*ones(1,length(t));%nm

%generates Vitreous drug concentration over time
semilogy(t,Drugv(i, :),'LineWidth',2)
hold on
plot(t,UpperIC50,'m','LineWidth',2)
plot(t, LowerIC50,'b','LineWidth',2)
ylabel('Abicipar Concentration (nM)')
xlim([-5 500])
ylim([10^(-3)  10^(4)])
xlabel('Time (days)')
if figure_count == 1
legend('Vitreous Concentration','Upper IC_{50}','Lower IC_{50}','FontSize',14, 'Location', 'northeast')  %uncomment to add legend
end
title('Without DDS')
pbaspect([1 1 1])
axis square
set(gca, 'FontSize', 16)
%exportgraphics(figure(figure_count),sprintf('drug_concentration_non_DDS%d.png', doseo(i)), 'Resolution', 300)
figure_count = figure_count+1;
hold off
end

load ("DDS_doses.mat")
% DDS geometry
DDS_geometry = "Chitosan_PCL";
%DDS_geometry = "Chitosan";
%DDS_geometry = "PCL";


radius_scale = [1, 1, 1, 1];
thickness_scale = [1, 1, 1, 1];

for i = 1:length(doseo) %for loop for each dose

subplot(2,4,figure_count) 
% Vitreous abicipar concentration
UpperIC50=6*ones(1,length(t)); %nm
LowerIC50 =0.017*ones(1,length(t));%nm
semilogy(t,Drugv(i, :),'LineWidth',2)
hold on
plot(t,UpperIC50,'m','LineWidth',2)
plot(t, LowerIC50,'b','LineWidth',2)
ylabel('Abicipar Concentration (nM)')
xlim([-5 500])
ylim([10^(-3)  10^(4)])
xlabel('Time (days)')
%legend('Vitreous Concentration','UpperIC50','LowerIC50')  %uncomment to add legend
title('With DDS')
pbaspect([1 1 1])
axis square
%fontsize(figure(2), 20, "points")
set(gca, 'FontSize', 16)
%exportgraphics(figure(figure_count),sprintf('drug_concentration%d.png', doseo(i)), 'Resolution', 300)
figure_count = figure_count+1;
hold off
end

labelstring = {'A', 'B', 'C', 'D','E', 'F', 'G', 'H'};
for v = 1:8
    subplot(2,4,v)
    hold on
    text(-0.2, 1.1, labelstring(v)', 'Units', 'normalized', 'FontWeight', 'bold','FontSize', 20)
end

widthInches = 20;
heightInches = 10;
run('ScriptForExportingImages.m')


%%%%%%%%%%%%%%%% figure 3 %%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(3);
figname = 'figure3';
figure_count = 1;

load ("Without_DDS.mat")

% bar plot IC50
barWidth = 1;
bar_graph_data = [in_vivo_suppression_time; in_vitro_suppression_time];
% Dose response
subplot(2,5,figure_count)
hold on
hBar1 = bar(bar_graph_data', 'grouped', 'LineWidth', 1.5); 
xticks(1:length(doseo)); 
xticklabels({'0.01', '.05', '0.1', '0.5'}); 
ylabel('PD Suppression Time (days)');
xlabel('Drug Amount (mg)');
ylim([0  500])
set(hBar1, 'BarWidth', barWidth);
title('Without DDS (Vitreous)');
set(gca, 'FontSize', 14)
legend({'Upper IC_{50}', 'Lower IC_{50}'}, 'FontSize',14, 'Location', 'northwest'); % Legend for rows
pbaspect([1 1 1])
axis square
box on
%exportgraphics(figure(figure_count),sprintf('bar_dose_response_non_DDS.png'), 'Resolution', 300)
figure_count = figure_count+1;
hold off

% bar plot of 10% VEGF suppression
barWidth = 1;
subplot(2,5,figure_count)
hold on
hBar1 = bar([Data_time_at_target_vit_10',Data_time_at_target_aq_10'], 'grouped', 'LineWidth', 1.5); 
xticks(1:length(doseo)); 
xticklabels({'0.01', '.05', '0.1', '0.5'});
ylabel('PD Suppression Time (Days)');
xlabel('Drug Amount (mg)');
set(hBar1, 'BarWidth', barWidth);
%title(sprintf('Dose %0.2f mg',doseo(j)));
ylim([0  500])
set(gca, 'FontSize', 14)
legend({'Vitreous 10%', 'Aqueous 10%'}, 'FontSize',14, 'Location', 'northwest');
pbaspect([1 1 1])
axis square
box on
%exportgraphics(figure(figure_count),sprintf('bar_10_percent_VEGF.png'), 'Resolution', 300)
figure_count = figure_count+1;
hold off

% bar plot of 50% VEGF suppression
barWidth = 1;
subplot(2,5,figure_count)
hold on
hBar1 = bar([Data_time_at_target_vit_50',Data_time_at_target_aq_50'], 'grouped', 'LineWidth', 1.5); 
xticks(1:length(doseo)); 
xticklabels({'0.01', '.05', '0.1', '0.5'});
ylabel('PD Suppression Time (Days)');
xlabel('Drug Amount (mg)');
set(hBar1, 'BarWidth', barWidth);
ylim([0  500])
%title(sprintf('Dose %0.2f mg',doseo(j)));
set(gca, 'FontSize', 14)
legend({'Vitreous 50%', 'Aqueous 50%'}, 'FontSize',14, 'Location', 'northwest');
pbaspect([1 1 1])
axis square
box on
%exportgraphics(figure(figure_count),sprintf('bar_50_percent_VEGF.png'), 'Resolution', 300)
figure_count = figure_count+1;
hold off


%VEGF Concentration vitreous with different doses
subplot(2,5,figure_count)
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
%exportgraphics(figure(figure_count),sprintf('VEGFv_non_DDS.png'), 'Resolution', 300)
figure_count = figure_count+1;
hold off

%VEGF Concentration aqueous with different doses
subplot(2,5,figure_count)
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
%exportgraphics(figure(figure_count),sprintf('VEGFaq_non_DDS.png'), 'Resolution', 300)
figure_count = figure_count+1;
hold off


load ("DDS_doses.mat")

% bar plot IC50
barWidth = 1;
bar_graph_data = [in_vivo_suppression_time; in_vitro_suppression_time];
% Dose response
subplot(2,5,figure_count)
hold on
hBar1 = bar(bar_graph_data', 'grouped', 'LineWidth', 1.5);
ylim([0  500])
xticks(1:length(doseo)); 
xticklabels({'0.01', '.05', '0.1', '0.5'}); 
ylabel('PD Suppression Time (days)');
xlabel('Drug Amount (mg)');
set(hBar1, 'BarWidth', barWidth);
title('With DDS (Vitreous)');
set(gca, 'FontSize', 14)
legend({'Upper IC_{50}', 'Lower IC_{50}'}, 'FontSize',14, 'Location', 'northwest'); % Legend for rows
pbaspect([1 1 1])
axis square
box on
%exportgraphics(figure(figure_count),sprintf('bar_dose_response.png'), 'Resolution', 300)
figure_count = figure_count+1;
hold off

% bar plot of 10% VEGF suppression
barWidth = 1;
subplot(2,5,figure_count)
hold on
hBar1 = bar([Data_time_at_target_vit_10',Data_time_at_target_aq_10'], 'grouped', 'LineWidth', 1.5); 
xticks(1:length(doseo)); 
xticklabels({'0.01', '.05', '0.1', '0.5'});
ylabel('PD Suppression Time (Days)');
xlabel('Drug Amount (mg)');
set(hBar1, 'BarWidth', barWidth);
ylim([0  500])
%title(sprintf('Dose %0.2f mg',doseo(j)));
set(gca, 'FontSize', 14)
legend({'Vitreous 10%', 'Aqueous 10%'}, 'FontSize',14, 'Location', 'northwest');
pbaspect([1 1 1])
axis square
box on
%exportgraphics(figure(figure_count),sprintf('bar_10_percent_VEGF.png'), 'Resolution', 300)
figure_count = figure_count+1;
hold off

% bar plot of 50% VEGF suppression
barWidth = 1;
subplot(2,5,figure_count)
hold on
hBar1 = bar([Data_time_at_target_vit_50',Data_time_at_target_aq_50'], 'grouped', 'LineWidth', 1.5); 
xticks(1:length(doseo)); 
xticklabels({'0.01', '.05', '0.1', '0.5'});
ylabel('PD Suppression Time (Days)');
xlabel('Drug Amount (mg)');
set(hBar1, 'BarWidth', barWidth);
ylim([0  500])
%title(sprintf('Dose %0.2f mg',doseo(j)));
set(gca, 'FontSize', 14)
legend({'Vitreous 50%', 'Aqueous 50%'}, 'FontSize',14, 'Location', 'northwest');
pbaspect([1 1 1])
axis square
box on
%exportgraphics(figure(figure_count),sprintf('bar_50_percent_VEGF.png'), 'Resolution', 300)
figure_count = figure_count+1;
hold off

%VEGF Concentration vitreous with different doses
subplot(2,5,figure_count)
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
%exportgraphics(figure(figure_count),sprintf('VEGFv.png'), 'Resolution', 300)
figure_count = figure_count+1;
hold off

%VEGF Concentration aqueous with different doses
subplot(2,5,figure_count)
set(gca, 'FontSize', 14)
hold on
box on
for j = 1:length(doseo)
    plot(t, (VEGFdoseaq(j, :)).*microtopico, 'LineWidth', 2)
end
%xline(in_vivo_suppression_time(j),'k','LineWidth',2)
%xline(in_vitro_suppression_time(j),'k','LineWidth',2)
%legend('Dose = 0.01 mg', 'Dose = 0.05 mg', 'Dose = 0.1 mg', 'Dose = 0.5 mg', 'FontSize',12, 'Location','southeast');
legend('Dose 0.01 mg', 'Dose 0.05 mg', 'Dose 0.1 mg', 'Dose 0.5 mg', 'FontSize',12, 'Location','southeast');
xlabel('Time (days)');
ylabel('VEGF Concentration (pM)');
%xlim([-5 500])
title('VEGF Aqueous');
%fontsize(figure(3), 17, "points")
pbaspect([1 1 1])
axis square
figure_count = figure_count+1;
hold off

labelstring = {'A', 'B', 'C', 'D','E', 'F', 'G', 'H', 'I', 'J'};
for v = 1:10
    subplot(2,5,v)
    hold on
    text(-0.2, 1.1, labelstring(v)', 'Units', 'normalized', 'FontWeight', 'bold','FontSize', 20)
end

widthInches = 25;
heightInches = 10;
run('ScriptForExportingImages.m')



%%%%%%%%%%%%%%%% figure 4 %%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(4);
figname = 'figure4';
figure_count = 1;

load ("DDS_doses.mat")

for i = 1:length(doseo) %for loop for each dose
% plot_DDS_drug_release_dynamics
subplot(2,4,figure_count)
set(gca, 'FontSize', 14)
hold on
box on
plot(drug_release_time(i,2:end), drug_release_profile(i,2:end), 'LineWidth', 2)
xlim([-2 200])
xlabel('Time (days)')
ylabel('DDS Drug release rate (mg/day)');
legend(sprintf('Dose %0.2f mg', doseo(i)), 'Location', 'northeast')
%fontsize(figure(3), 17, "points")
pbaspect([1 1 1])
axis square
%exportgraphics(figure(figure_count),sprintf('drug_release_DDS%d.png', doseo(i)), 'Resolution', 300)
figure_count = figure_count+1;
hold off
end

%Drug Concentration aqueous with different doses
subplot(2,4,figure_count)
set(gca, 'FontSize', 14)
box on
hold on
for j = 1:length(doseo)
    semilogy(t, Drugaq(j, :), 'LineWidth', 2)
end
legend('Dose 0.01 mg', 'Dose 0.05 mg', 'Dose 0.1 mg', 'Dose 0.5 mg', 'FontSize',12, 'Location','northeast');
xlabel('Time (days)');
ylabel('Abicipar Concentration (nM)');
xlim([-5 500])
title('Aqueous');
%fontsize(figure(5), 17, "points")
pbaspect([1 1 1])
axis square
%exportgraphics(figure(figure_count),sprintf('Drugaq.png'), 'Resolution', 300)
figure_count = figure_count+1;
hold off

%Drug Concentration retina with different doses
subplot(2,4,figure_count)
set(gca, 'FontSize', 14)
hold on
box on
for j = 1:length(doseo)
    plot(t, Drugr(j, :), 'LineWidth', 2)
end
%legend('Dose = 0.01mg', 'Dose = 0.05mg', 'Dose = 0.1mg', 'Dose = 0.5mg', 'Location','southeast');
xlabel('Time (days)');
ylabel('Abicipar Concentration (nM)');
xlim([-5 500])
title('Retina');
%fontsize(figure(7), 17, "points")
pbaspect([1 1 1])
axis square
%exportgraphics(figure(figure_count),sprintf('drugretina.png'), 'Resolution', 300)
figure_count = figure_count+1;
hold off

%Drug Concentration choroid with different doses
subplot(2,4,figure_count)
set(gca, 'FontSize', 14)
hold on
box on
for j = 1:length(doseo)
    plot(t, Drugc(j, :), 'LineWidth', 2)
end
%legend('Dose = 0.01mg', 'Dose = 0.05mg', 'Dose = 0.1mg', 'Dose = 0.5mg', 'Location','northeast');
xlabel('Time (days)');
ylabel('Abicipar Concentration (nM)');
xlim([-5 500])
title('Choroid');
%fontsize(figure(6), 17, "points")
pbaspect([1 1 1])
axis square
%exportgraphics(figure(figure_count),sprintf('drugchoroid.png'), 'Resolution', 300)
figure_count = figure_count+1;
hold off

%Drug Concentration serum with different doses
subplot(2,4,figure_count)
set(gca, 'FontSize', 14)
hold on
box on
for j = 1:length(doseo)
    plot(t, Drugs(j, :), 'LineWidth', 2)
end
%legend('Dose = 0.01mg', 'Dose = 0.05mg', 'Dose = 0.1mg', 'Dose = 0.5mg', 'Location','southeast');
xlabel('Time (days)');
ylabel('Abicipar Concentration (nM)');
xlim([-5 500])
title('Serum');
%fontsize(figure(8), 17, "points")
pbaspect([1 1 1])
axis square
%exportgraphics(figure(figure_count),sprintf('Drugserum.png'), 'Resolution', 300)
figure_count = figure_count+1;
hold off

labelstring = {'A', 'B', 'C', 'D','E', 'F', 'G', 'H'};
for v = 1:8
    subplot(2,4,v)
    hold on
    text(-0.2, 1.1, labelstring(v)', 'Units', 'normalized', 'FontWeight', 'bold','FontSize', 20)
end

widthInches = 20;
heightInches = 10;
run('ScriptForExportingImages.m')



%%%%%%%%%%%%%%%% figure 5 %%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(5);
figname = 'figure5';
figure_count = 1;


% chitosan single layer

doseo = [0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1]; %various doses to be tested in milligrams
radius_scale = [0.001, 1, 2, 3, 5, 10, 20];
thickness_scale = [1, 1, 1, 1, 1, 1, 1];

load ("chitosan_single.mat")

% bar plot chitosan single layer
barWidth = 1;
bar_graph_data = [in_vivo_suppression_time; in_vitro_suppression_time];

%chitosan single layer
subplot(2,3,figure_count)
hold on
hBar1 = bar(bar_graph_data', 'grouped', 'LineWidth', 1.5); 
xticks(1:length(radius_scale)); 
xticklabels({'0.001R_{C}', 'R_{C}', '2R_{C}', '3R_{C}','5R_{C}', '10R_{C}', '20R_{C}'}); 
ylabel('Time (days)');
set(hBar1, 'BarWidth', barWidth);
ylim([0  1000])
title('Chitosan Single Layer Varied R_{C}');
set(gca, 'FontSize', 12)
legend({'Upper IC_{50}', 'Lower IC_{50}'}, 'FontSize',14, 'Location', 'northwest'); % Legend for rows
pbaspect([1 1 1])
figure_count = figure_count+1;
axis square
%exportgraphics(figure(figure_count),sprintf('barsingleC.png'), 'Resolution', 300)
hold off

% bar plot of 10% VEGF suppression
barWidth = 1;
subplot(2,3,figure_count)
hold on
hBar1 = bar([Data_time_at_target_vit_10',Data_time_at_target_aq_10'], 'grouped', 'LineWidth', 1.5); 
xticks(1:length(radius_scale)); 
xticklabels({'0.001R_{C}', 'R_{C}', '2R_{C}', '3R_{C}','5R_{C}', '10R_{C}', '20R_{C}'});
ylabel('VEGF Suppression Time (Days)');
set(hBar1, 'BarWidth', barWidth);
ylim([0  1000])
title('Chitosan Single Layer Varied R_{C}');
set(gca, 'FontSize', 12)
legend({'Vitreous 10%', 'Aqueous 10%'}, 'FontSize',14, 'Location', 'northwest');
pbaspect([1 1 1])
axis square
box on
%exportgraphics(figure(figure_count),sprintf('barsingleC_10_percent_VEGF.png'), 'Resolution', 300)
figure_count = figure_count+1;
hold off

% bar plot of 50% VEGF suppression
barWidth = 1;
subplot(2,3,figure_count)
hold on
hBar1 = bar([Data_time_at_target_vit_50',Data_time_at_target_aq_50'], 'grouped', 'LineWidth', 1.5); 
xticks(1:length(radius_scale)); 
xticklabels({'0.001R_{C}', 'R_{C}', '2R_{C}', '3R_{C}','5R_{C}', '10R_{C}', '20R_{C}'});
ylabel('VEGF Suppression Time (Days)');
set(hBar1, 'BarWidth', barWidth);
ylim([0  1000])
title('Chitosan Single Layer Varied R_{C}');
set(gca, 'FontSize', 12)
legend({'Vitreous 50%', 'Aqueous 50%'}, 'FontSize',14, 'Location', 'northwest');
pbaspect([1 1 1])
axis square
box on
%exportgraphics(figure(figure_count),sprintf('barsingleC_50_percent_VEGF.png'), 'Resolution', 300)
figure_count = figure_count+1;
hold off


% PCL single
doseo = [0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1]; %various doses to be tested in milligrams
radius_scale = [0.001, 1, 10, 25, 50, 75, 100];
thickness_scale = [1, 1, 1, 1, 1, 1, 1];

load ("PCL_single.mat")

% bar plot
barWidth = 1;
bar_graph_data = [in_vivo_suppression_time; in_vitro_suppression_time];

%PCL single layer
subplot(2,3,figure_count)
hold on
hBar1 = bar(bar_graph_data', 'grouped', 'LineWidth', 1.5); 
xticks(1:length(radius_scale)); 
xticklabels({'0.001R_{P}', 'R_{P}', '10R_{P}', '25R_{P}', '50R_{P}', '75R_{P}', '100R_{P}'}); 
ylabel('Time (days)');
set(hBar1, 'BarWidth', barWidth);
ylim([0  1000])
title('PCL Single Layer Varied R_{P}');
set(gca, 'FontSize', 12)
legend({'Upper IC_{50}', 'Lower IC_{50}'}, 'FontSize',14, 'Location', 'northwest'); % Legend for rows
pbaspect([1 1 1])
axis square
%exportgraphics(figure(figure_count),sprintf('barsinglePCL.png'), 'Resolution', 300)
figure_count = figure_count+1;
hold off

% bar plot of 10% VEGF suppression
barWidth = 1;
subplot(2,3,figure_count)
hold on
hBar1 = bar([Data_time_at_target_vit_10',Data_time_at_target_aq_10'], 'grouped', 'LineWidth', 1.5); 
xticks(1:length(radius_scale)); 
xticklabels({'0.001R_{P}', 'R_{P}', '10R_{P}', '25R_{P}', '50R_{P}', '75R_{P}', '100R_{P}'}); 
ylabel('VEGF Suppression Time (Days)');
set(hBar1, 'BarWidth', barWidth);
ylim([0  1000])
title('PCL Single Layer Varied R_{P}');
set(gca, 'FontSize', 12)
legend({'Vitreous 10%', 'Aqueous 10%'}, 'FontSize',14, 'Location', 'northwest');
pbaspect([1 1 1])
axis square
box on
%exportgraphics(figure(figure_count),sprintf('barsinglePCL_10_percent_VEGF.png'), 'Resolution', 300)
figure_count = figure_count+1;
hold off

% bar plot of 50% VEGF suppression
barWidth = 1;
subplot(2,3,figure_count)
hold on
hBar1 = bar([Data_time_at_target_vit_50',Data_time_at_target_aq_50'], 'grouped', 'LineWidth', 1.5); 
xticks(1:length(radius_scale)); 
xticklabels({'0.001R_{P}', 'R_{P}', '10R_{P}', '25R_{P}', '50R_{P}', '75R_{P}', '100R_{P}'}); 
ylabel('VEGF Suppression Time (Days)');
set(hBar1, 'BarWidth', barWidth);
ylim([0  1000])
title('PCL Single Layer Varied R_{P}');
set(gca, 'FontSize', 12)
legend({'Vitreous 50%', 'Aqueous 50%'}, 'FontSize',14, 'Location', 'northwest');
pbaspect([1 1 1])
axis square
box on
%exportgraphics(figure(figure_count),sprintf('barsinglePCL_50_percent_VEGF.png'), 'Resolution', 300)
figure_count = figure_count+1;
hold off

labelstring = {'A', 'B', 'C', 'D','E', 'F'};
for v = 1:6
    subplot(2,3,v)
    hold on
    text(-0.2, 1.1, labelstring(v)', 'Units', 'normalized', 'FontWeight', 'bold','FontSize', 20)
end

widthInches = 15;
heightInches = 10;
run('ScriptForExportingImages.m')


%%%%%%%%%%%%%%%% figure 6 %%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(6);
figname = 'figure6';
figure_count = 1;


% chitosan single layer

doseo = [0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1]; %various doses to be tested in milligrams
radius_scale = [0.001, 1, 2, 3, 5, 10, 20];
thickness_scale = [1, 1, 1, 1, 1, 1, 1];

load ("chitosan_single.mat")


for i = 1:length(doseo) %for loop for each dose

% plot_DDS_drug_release_dynamics
subplot(2,7,figure_count)
hold on
box on
plot(drug_release_time(i,2:end), drug_release_profile(i,2:end), 'LineWidth', 2)
xlim([-2 200])
xlabel('Time (days)')
ylabel('DDS Drug release rate (mg/day)');
legend(sprintf('%0.3fR_C', radius_scale(i)), 'Location', 'northeast')
%fontsize(figure(3), 17, "points")
set(gca, 'FontSize', 14)
pbaspect([1 1 1])
axis square
%exportgraphics(figure(figure_count),sprintf('drug_release_DDS%d.png', radius_scale(i)), 'Resolution', 300)
figure_count = figure_count+1;
hold off
end


for i = 1:length(doseo) %for loop for each dose
% Vitreous abicipar concentration
UpperIC50=6*ones(1,length(t)); %nm
LowerIC50 =0.017*ones(1,length(t));%nm

subplot(2,7,figure_count)
semilogy(t,Drugv(i, :),'LineWidth',2)
hold on
plot(t,UpperIC50,'m','LineWidth',2)
plot(t, LowerIC50,'b','LineWidth',2)
ylabel('Abicipar Concentration (nM)')
%xlim([-2 200])
ylim([10^(-3)  10^(4)])
xlabel('Time (days)')
%legend('Vitreous Concentration','UpperIC50','LowerIC50')  %uncomment to add legend
%title('Vitreous with DDS')
pbaspect([1 1 1])
axis square
%fontsize(figure(2), 20, "points")
set(gca, 'FontSize', 14)
%exportgraphics(figure(figure_count),sprintf('drug_concentration%d.png', radius_scale(i)), 'Resolution', 300)
figure_count = figure_count+1;
hold off
end



labelstring = {'A', 'B', 'C', 'D','E', 'F','G', 'H', 'I', 'J','K', 'L','M','N'};
for v = 1:14
    subplot(2,7,v)
    hold on
    text(-0.2, 1.1, labelstring(v)', 'Units', 'normalized', 'FontWeight', 'bold','FontSize', 20)
end

widthInches = 35;
heightInches = 10;
run('ScriptForExportingImages.m')


%%%%%%%%%%%%%%%% figure 7 %%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(7);
figname = 'figure7';
figure_count = 1;


% PCL single layer

doseo = [0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1]; %various doses to be tested in milligrams
radius_scale = [0.001, 1, 10, 25, 50, 75, 100];
thickness_scale = [1, 1, 1, 1, 1, 1, 1];

load ("PCL_single.mat")



for i = 1:length(doseo) %for loop for each dose

% plot_DDS_drug_release_dynamics
subplot(2,7,figure_count)
set(gca, 'FontSize', 14)
hold on
box on
plot(drug_release_time(i,2:end), drug_release_profile(i,2:end), 'LineWidth', 2)
xlim([-2 200])
xlabel('Time (days)')
ylabel('DDS Drug release rate (mg/day)');
legend(sprintf('%0.3fR_P', radius_scale(i)), 'Location', 'northeast')
%fontsize(figure(3), 17, "points")
pbaspect([1 1 1])
axis square
%exportgraphics(figure(figure_count),sprintf('drug_release_DDS%d.png', radius_scale(i)), 'Resolution', 300)
figure_count = figure_count+1;
hold off
end


for i = 1:length(doseo) %for loop for each dose
% Vitreous abicipar concentration
UpperIC50=6*ones(1,length(t)); %nm
LowerIC50 =0.017*ones(1,length(t));%nm


subplot(2,7,figure_count)
semilogy(t,Drugv(i, :),'LineWidth',2)
hold on
plot(t,UpperIC50,'m','LineWidth',2)
plot(t, LowerIC50,'b','LineWidth',2)
ylabel('Abicipar Concentration (nM)')
%xlim([-2 200])
ylim([10^(-3)  10^(4)])
xlabel('Time (days)')
%legend('Vitreous Concentration','UpperIC50','LowerIC50')  %uncomment to add legend
%title('Vitreous with DDS')
pbaspect([1 1 1])
axis square
%fontsize(figure(2), 20, "points")
set(gca, 'FontSize', 14)
%exportgraphics(figure(figure_count),sprintf('drug_concentration%d.png', radius_scale(i)), 'Resolution', 300)
figure_count = figure_count+1;
hold off
end

labelstring = {'A', 'B', 'C', 'D','E', 'F','G', 'H', 'I', 'J','K', 'L','M','N'};
for v = 1:14
    subplot(2,7,v)
    hold on
    text(-0.2, 1.1, labelstring(v)', 'Units', 'normalized', 'FontWeight', 'bold','FontSize', 20)
end

widthInches = 35;
heightInches = 10;
run('ScriptForExportingImages.m')


%%%%%%%%%%%%%%%% figure 8 %%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(8);
figname = 'figure8';
figure_count = 1;


% chitosan bi layer

doseo = [0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1]; %various doses to be tested in milligrams
radius_scale = [0.01, 1, 1.5, 2, 3, 5, 10];
thickness_scale = [1, 1, 1, 1, 1, 1, 1];

load ("bi_layer_chitosan.mat")

% bar plot
barWidth = 1;
bar_graph_data = [in_vivo_suppression_time; in_vitro_suppression_time];

%chitosan double layer
subplot(3,3,figure_count);
hold on
box on
hBar1 = bar(bar_graph_data', 'grouped', 'LineWidth', 1.5); 
xticks(1:length(radius_scale)); 
xticklabels({'0.01R_{C}', 'R_{C}', '1.5R_{C}', '2R_{C}', '3R_{C}', '5R_{C}', '10R_{C}'}); 
ylabel('Time (days)');
set(hBar1, 'BarWidth', barWidth);
ylim([0  1000])
title('Bi-Layered Varied R_{C}, Constant \DeltaR');
set(gca, 'FontSize', 8)
legend({'Upper IC_{50}', 'Lower IC_{50}'}, 'FontSize',8, 'Location', 'northwest'); % Legend for rows
pbaspect([1 1 1])
axis square
%exportgraphics(figure(figure_count),sprintf('barbilayerC.png'), 'Resolution', 300)
figure_count = figure_count+1;
hold off

% bar plot of 10% VEGF suppression
barWidth = 1;
subplot(3,3,figure_count);
hold on
hBar1 = bar([Data_time_at_target_vit_10',Data_time_at_target_aq_10'], 'grouped', 'LineWidth', 1.5); 
xticks(1:length(radius_scale)); 
xticklabels({'0.01R_{C}', 'R_{C}', '1.5R_{C}', '2R_{C}', '3R_{C}', '5R_{C}', '10R_{C}'}); 
ylabel('VEGF Suppression Time (Days)');
set(hBar1, 'BarWidth', barWidth);
ylim([0  1000])
title('Bi-Layered Varied R_{C}, Constant \DeltaR');
set(gca, 'FontSize', 8)
legend({'Vitreous 10%', 'Aqueous 10%'}, 'FontSize',8, 'Location', 'northwest');
pbaspect([1 1 1])
axis square
box on
%exportgraphics(figure(figure_count),sprintf('barbilayerC_10_percent_VEGF.png'), 'Resolution', 300)
figure_count = figure_count+1;
hold off

% bar plot of 50% VEGF suppression
barWidth = 1;
subplot(3,3,figure_count);
hold on
hBar1 = bar([Data_time_at_target_vit_50',Data_time_at_target_aq_50'], 'grouped', 'LineWidth', 1.5); 
xticks(1:length(radius_scale)); 
xticklabels({'0.01R_{C}', 'R_{C}', '1.5R_{C}', '2R_{C}', '3R_{C}', '5R_{C}', '10R_{C}'}); 
ylabel('VEGF Suppression Time (Days)');
set(hBar1, 'BarWidth', barWidth);
ylim([0  1000])
title('Bi-Layered Varied R_{C}, Constant \DeltaR');
set(gca, 'FontSize', 8)
legend({'Vitreous 50%', 'Aqueous 50%'}, 'FontSize',8, 'Location', 'northwest');
pbaspect([1 1 1])
axis square
box on
%exportgraphics(figure(figure_count),sprintf('barbilayerC_50_percent_VEGF.png'), 'Resolution', 300)
figure_count = figure_count+1;
hold off


% PCL bi layer
doseo = [0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1]; %various doses to be tested in milligrams
radius_scale = [1, 1, 1, 1, 1, 1, 1];
thickness_scale = [0.01, 0.1, 1, 10, 20, 25, 30];

load ("bi_layered_PCL.mat")

% bar plot
barWidth = 1;
bar_graph_data = [in_vivo_suppression_time; in_vitro_suppression_time];

%PCL double layer
subplot(3,3,figure_count)
hold on
box on
hBar1 = bar(bar_graph_data', 'grouped', 'LineWidth', 1.5); 
xticks(1:length(thickness_scale)); 
xticklabels({'0.01\DeltaR', '0.1\DeltaR', '\DeltaR', '10\DeltaR', '20\DeltaR', '25\DeltaR', '30\DeltaR'}); 
ylabel('Time (days)');
set(hBar1, 'BarWidth', barWidth);
ylim([0  1000])
title('Bi-Layered Varied \DeltaR, Constant R_{C}');
set(gca, 'FontSize', 8)
legend({'Upper IC_{50}', 'Lower IC_{50}'}, 'FontSize',8, 'Location', 'northwest');
pbaspect([1 1 1])
axis square
%exportgraphics(figure(figure_count),sprintf('barbilayerPCL.png'), 'Resolution', 300)
figure_count = figure_count+1;
hold off


% bar plot of 10% VEGF suppression
barWidth = 1;
subplot(3,3,figure_count)
hold on
hBar1 = bar([Data_time_at_target_vit_10',Data_time_at_target_aq_10'], 'grouped', 'LineWidth', 1.5); 
xticks(1:length(thickness_scale)); 
xticklabels({'0.01\DeltaR', '0.1\DeltaR', '\DeltaR', '10\DeltaR', '20\DeltaR', '25\DeltaR', '30\DeltaR'}); 
ylabel('VEGF Suppression Time (Days)');
set(hBar1, 'BarWidth', barWidth);
ylim([0  1000])
title('Bi-Layered Varied \DeltaR, Constant R_{C}');
set(gca, 'FontSize', 8)
legend({'Vitreous 10%', 'Aqueous 10%'}, 'FontSize',8, 'Location', 'northwest');
pbaspect([1 1 1])
axis square
box on
%exportgraphics(figure(figure_count),sprintf('barbilayerC_10_percent_VEGF.png'), 'Resolution', 300)
figure_count = figure_count+1;
hold off

% bar plot of 50% VEGF suppression
barWidth = 1;
subplot(3,3,figure_count)
hold on
hBar1 = bar([Data_time_at_target_vit_50',Data_time_at_target_aq_50'], 'grouped', 'LineWidth', 1.5); 
xticks(1:length(thickness_scale)); 
xticklabels({'0.01\DeltaR', '0.1\DeltaR', '\DeltaR', '10\DeltaR', '20\DeltaR', '25\DeltaR', '30\DeltaR'}); 
ylabel('VEGF Suppression Time (Days)');
set(hBar1, 'BarWidth', barWidth);
ylim([0  1000])
title('Bi-Layered Varied \DeltaR, Constant R_{C}');
set(gca, 'FontSize', 8)
legend({'Vitreous 50%', 'Aqueous 50%'}, 'FontSize',8, 'Location', 'northwest');
pbaspect([1 1 1])
axis square
box on
%exportgraphics(figure(figure_count),sprintf('barbilayerC_50_percent_VEGF.png'), 'Resolution', 300)
figure_count = figure_count+1;
hold off

% chitosan pcl bi layer

doseo = [0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1]; %various doses to be tested in milligrams
radius_scale = [0.1, 0.5, 1.5, 2, 3, 4, 5];
thickness_scale = [10, 2, 3, 5, 6, 2, 0.5];

load ("bi_layer_changing_both.mat")


% bar plot
barWidth = 1;
bar_graph_data = [in_vivo_suppression_time; in_vitro_suppression_time];
% changing both
subplot(3,3,figure_count)
hold on
box on
hBar1 = bar(bar_graph_data', 'grouped', 'LineWidth', 1.5); 
xticks(1:length(thickness_scale)); 
xticklabels({'0.1R_{C},10\DeltaR', '0.5R_{C},2\DeltaR', '1.5R_{C},3\DeltaR', '2R_{C},5\DeltaR', '3R_{C},6\DeltaR', '4R_{C},2\DeltaR', '5R_{C},0.5\DeltaR'}); 
ylabel('Time (days)');
set(hBar1, 'BarWidth', barWidth);
ylim([0  1000])
title('Bi-Layered Varied R_{C} and \DeltaR');
set(gca, 'FontSize', 8)
legend({'Upper IC_{50}', 'Lower IC_{50}'}, 'FontSize',8, 'Location', 'northwest');
pbaspect([1 1 1])
axis square
%hold off
%exportgraphics(figure(figure_count),sprintf('changingboth.png'), 'Resolution', 300)
figure_count = figure_count+1;
hold off

% bar plot of 10% VEGF suppression
barWidth = 1;
subplot(3,3,figure_count)
hold on
hBar1 = bar([Data_time_at_target_vit_10',Data_time_at_target_aq_10'], 'grouped', 'LineWidth', 1.5); 
xticks(1:length(thickness_scale)); 
xticklabels({'0.1R_{C},10\DeltaR', '0.5R_{C},2\DeltaR', '1.5R_{C},3\DeltaR', '2R_{C},5\DeltaR', '3R_{C},6\DeltaR', '4R_{C},2\DeltaR', '5R_{C},0.5\DeltaR'}); 
ylabel('VEGF Suppression Time (Days)');
set(hBar1, 'BarWidth', barWidth);
ylim([0  1000])
title('Bi-Layered Varied R_{C} and \DeltaR');
set(gca, 'FontSize', 8)
legend({'Vitreous 10%', 'Aqueous 10%'}, 'FontSize',8, 'Location', 'northwest');
pbaspect([1 1 1])
axis square
box on
%exportgraphics(figure(figure_count),sprintf('barbilayerC_10_percent_VEGF.png'), 'Resolution', 300)
figure_count = figure_count+1;
hold off

% bar plot of 50% VEGF suppression
barWidth = 1;
subplot(3,3,figure_count)
hold on
hBar1 = bar([Data_time_at_target_vit_50',Data_time_at_target_aq_50'], 'grouped', 'LineWidth', 1.5); 
xticks(1:length(thickness_scale)); 
xticklabels({'0.1R_{C},10\DeltaR', '0.5R_{C},2\DeltaR', '1.5R_{C},3\DeltaR', '2R_{C},5\DeltaR', '3R_{C},6\DeltaR', '4R_{C},2\DeltaR', '5R_{C},0.5\DeltaR'}); 
ylabel('VEGF Suppression Time (Days)');
set(hBar1, 'BarWidth', barWidth);
ylim([0  1000])
title('Bi-Layered Varied R_{C} and \DeltaR');
set(gca, 'FontSize', 8)
legend({'Vitreous 50%', 'Aqueous 50%'}, 'FontSize',8, 'Location', 'northwest');
pbaspect([1 1 1])
axis square
box on
%exportgraphics(figure(figure_count),sprintf('barbilayerC_50_percent_VEGF.png'), 'Resolution', 300)
figure_count = figure_count+1;
hold off

labelstring = {'A', 'B', 'C', 'D','E', 'F','G','H', 'I'};
for v = 1:9
    subplot(3,3,v)
    hold on
    text(-0.2, 1.1, labelstring(v)', 'Units', 'normalized', 'FontWeight', 'bold','FontSize', 20)
end

widthInches = 15;
heightInches = 15;
run('ScriptForExportingImages.m')








%%%%%%%%%%%%%%%% figure 9 %%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(9);
figname = 'figure9';
figure_count = 1;


% Chitosan bi-layer

doseo = [0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1]; %various doses to be tested in milligrams
radius_scale = [0.01, 1, 1.5, 2, 3, 5, 10];
thickness_scale = [1, 1, 1, 1, 1, 1, 1];
load ("bi_layer_chitosan.mat")


for i = 1:length(doseo) %for loop for each dose

% plot_DDS_drug_release_dynamics
subplot(2,7,figure_count)
set(gca, 'FontSize', 14)
hold on
box on
plot(drug_release_time(i,2:end), drug_release_profile(i,2:end), 'LineWidth', 2)
xlim([-2 200])
xlabel('Time (days)')
ylabel('DDS Drug release rate (mg/day)');
%fontsize(figure(3), 17, "points")
legend(sprintf('%0.2fR_C', radius_scale(i)), 'Location', 'northeast')
pbaspect([1 1 1])
axis square
%exportgraphics(figure(figure_count),sprintf('drug_release_DDS%d.png', radius_scale(i)), 'Resolution', 300)
figure_count = figure_count+1;
hold off
end

for i = 1:length(doseo) %for loop for each dose
% Vitreous abicipar concentration
UpperIC50=6*ones(1,length(t)); %nm
LowerIC50 =0.017*ones(1,length(t));%nm


subplot(2,7,figure_count)
semilogy(t,Drugv(i, :),'LineWidth',2)
hold on
plot(t,UpperIC50,'m','LineWidth',2)
plot(t, LowerIC50,'b','LineWidth',2)
ylabel('Abicipar Concentration (nM)')
%xlim([-2 200])
ylim([10^(-3)  10^(4)])
xlabel('Time (days)')
%legend('Vitreous Concentration','UpperIC50','LowerIC50')  %uncomment to add legend
%title('Vitreous with DDS')
pbaspect([1 1 1])
axis square
set(gca, 'FontSize', 14)
%fontsize(figure(2), 20, "points")
%exportgraphics(figure(figure_count),sprintf('drug_concentration%d.png', radius_scale(i)), 'Resolution', 300)
figure_count = figure_count+1;
hold off
end


labelstring = {'A', 'B', 'C', 'D','E', 'F','G', 'H', 'I', 'J','K', 'L','M','N'};
for v = 1:14
    subplot(2,7,v)
    hold on
    text(-0.2, 1.1, labelstring(v)', 'Units', 'normalized', 'FontWeight', 'bold','FontSize', 20)
end

widthInches = 35;
heightInches = 10;
run('ScriptForExportingImages.m')


%%%%%%%%%%%%%%%% figure 10 %%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(10);
figname = 'figure10';
figure_count = 1;


% PCL bi-layer

doseo = [0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1]; %various doses to be tested in milligrams
radius_scale = [1, 1, 1, 1, 1, 1, 1];
thickness_scale = [0.01, 0.1, 1, 10, 20, 25, 30];
load ("bi_layered_PCL.mat")


for i = 1:length(doseo) %for loop for each dose

% plot_DDS_drug_release_dynamics
subplot(2,7,figure_count)
set(gca, 'FontSize', 14)
hold on
box on
plot(drug_release_time(i,2:end), drug_release_profile(i,2:end), 'LineWidth', 2)
xlim([-2 200])
xlabel('Time (days)')
ylabel('DDS Drug release rate (mg/day)');
legend(sprintf('%0.2fΔR', thickness_scale(i)), 'Location', 'northeast')
%fontsize(figure(3), 17, "points")
pbaspect([1 1 1])
axis square
%exportgraphics(figure(figure_count),sprintf('drug_release_DDS%d.png', radius_scale(i)), 'Resolution', 300)
%exportgraphics(figure(figure_count),sprintf('drug_release_DDS%d.png', thickness_scale(i)), 'Resolution', 300)
figure_count = figure_count+1;
hold off
end

for i = 1:length(doseo) %for loop for each dose
% Vitreous abicipar concentration
UpperIC50=6*ones(1,length(t)); %nm
LowerIC50 =0.017*ones(1,length(t));%nm


subplot(2,7,figure_count)
semilogy(t,Drugv(i, :),'LineWidth',2)
hold on
plot(t,UpperIC50,'m','LineWidth',2)
plot(t, LowerIC50,'b','LineWidth',2)
ylabel('Abicipar Concentration (nM)')
%xlim([-2 200])
ylim([10^(-3)  10^(4)])
xlabel('Time (days)')
%legend('Vitreous Concentration','UpperIC50','LowerIC50')  %uncomment to add legend
%title('Vitreous with DDS')
pbaspect([1 1 1])
axis square
set(gca, 'FontSize', 14)
%fontsize(figure(2), 20, "points")
%exportgraphics(figure(figure_count),sprintf('drug_concentration%d.png', radius_scale(i)), 'Resolution', 300)
%exportgraphics(figure(figure_count),sprintf('drug_concentration%d.png', thickness_scale(i)), 'Resolution', 300)
figure_count = figure_count+1;
hold off
end


labelstring = {'A', 'B', 'C', 'D','E', 'F','G', 'H', 'I', 'J','K', 'L','M','N'};
for v = 1:14
    subplot(2,7,v)
    hold on
    text(-0.2, 1.1, labelstring(v)', 'Units', 'normalized', 'FontWeight', 'bold','FontSize', 20)
end

widthInches = 35;
heightInches = 10;
run('ScriptForExportingImages.m')

%%%%%%%%%%%%%%%% figure 11 %%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(11);
figname = 'figure11';
figure_count = 1;


% Chitosan PCL bi-layer

doseo = [0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1]; %various doses to be tested in milligrams
radius_scale = [0.1, 0.5, 1.5, 2, 3, 4, 5];
thickness_scale = [10, 2, 3, 5, 6, 2, 0.5];
load ("bi_layer_changing_both.mat")

for i = 1:length(doseo) %for loop for each dose

% plot_DDS_drug_release_dynamics
subplot(2,7,figure_count)
set(gca, 'FontSize', 14)
hold on
box on
plot(drug_release_time(i,2:end), drug_release_profile(i,2:end), 'LineWidth', 2)
xlim([-2 200])
xlabel('Time (days)')
ylabel('DDS Drug release rate (mg/day)');
legend(sprintf('%0.2fR_C, %0.2fΔR', radius_scale(i),thickness_scale(i)), 'Location', 'northeast')
%fontsize(figure(3), 17, "points")
pbaspect([1 1 1])
axis square
%exportgraphics(figure(figure_count),sprintf('drug_release_DDS%d.png', radius_scale(i)), 'Resolution', 300)
%exportgraphics(figure(figure_count),sprintf('drug_release_DDS%d.png', radius_scale(i)), 'Resolution', 300)
figure_count = figure_count+1;
hold off
end

for i = 1:length(doseo) %for loop for each dose
% Vitreous abicipar concentration
UpperIC50=6*ones(1,length(t)); %nm
LowerIC50 =0.017*ones(1,length(t));%nm


subplot(2,7,figure_count)
semilogy(t,Drugv(i, :),'LineWidth',2)
hold on
plot(t,UpperIC50,'m','LineWidth',2)
plot(t, LowerIC50,'b','LineWidth',2)
ylabel('Abicipar Concentration (nM)')
%xlim([-2 200])
ylim([10^(-3)  10^(4)])
xlabel('Time (days)')
%legend('Vitreous Concentration','UpperIC50','LowerIC50')  %uncomment to add legend
%title('Vitreous with DDS')
pbaspect([1 1 1])
axis square
set(gca, 'FontSize', 14)
%fontsize(figure(2), 20, "points")
%exportgraphics(figure(figure_count),sprintf('drug_concentration%d.png', radius_scale(i)), 'Resolution', 300)
%exportgraphics(figure(figure_count),sprintf('drug_concentration%d.png', radius_scale(i)), 'Resolution', 300)
figure_count = figure_count+1;
hold off
end


labelstring = {'A', 'B', 'C', 'D','E', 'F','G', 'H', 'I', 'J','K', 'L','M','N'};
for v = 1:14
    subplot(2,7,v)
    hold on
    text(-0.2, 1.1, labelstring(v)', 'Units', 'normalized', 'FontWeight', 'bold','FontSize', 20)
end

widthInches = 35;
heightInches = 10;
run('ScriptForExportingImages.m')