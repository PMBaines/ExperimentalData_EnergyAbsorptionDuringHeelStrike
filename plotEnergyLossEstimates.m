%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot figures presented in article                                        %
% -------------------------------------------------------------------------%
% Experimental data and processing files corresponsing to article:         %
% Experimental estimation of energy absorption during heel strike in human %
% barefoot walking.                                                        %
% by Patricia M. Baines*, A.L. Schwab* and A.J. van Soest^                 %
% * Delft University of Technology                                         %
% ^ Vrije Universiteit Amsterdam                                           %
%                                                                          %
% Corresponding author: p.m.baines@tudelft.nl                              %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% close all
% clc
% clear all

%% load: choose to run getEnergyLossEstimates or load saved mat file

getEnergyLossEstimates
% load('AllData.mat')

%% Default Figure Settings
set(0,'DefaultAxesLineStyleOrder','-|-.|--|:','DefaultLineLineWidth',2.5,'DefaultAxesFontsize',13,'DefaultAxesGridLineStyle','--')

% set(0,'DefaultFigureVisible','off')
set(0,'DefaultTextInterpreter','latex')
set(0,'DefaultLegendInterpreter','latex')
set(0,'DefaultTextFontname','Times New Roman')
set(0,'DefaultAxesFontName','Times New Roman')
set(0,'DefaultAxesLabelFontSizeMultiplier',1)
set(0,'DefaultFigureColor','w')

%% Choose subject and run to display

%choose subject and run
plot_subject = 6; % choose subject article
plot_run = Data.subject(plot_subject).n_normal(2); % choose trial article
% plot_subject = 7; % choose other subject
% plot_run = Data.subject(plot_subject).n_normal(2); % choose other trial

% pick data for readability
D = Data.subject(plot_subject).n(plot_run); 
t_kp = D.t_kp_adjusted;
t_opto = D.t_opto_filtered_adjusted;

toggle_save = 'off' %'export' %%% set to export if you want to save the plots, default 'off'
screen_shift = 0; % add a screenshift if you want to project on an external screen

folder = 'C:\Users'; % choose folder you want to save to


if strcmp(toggle_save,'on')
    disp('Warning: Automatic saving of plots is on!')
end

%% Add age and BMI
Age = [29,29,27,35,34,32,27,23,21,21,33,20];
mean(Age);
std(Age);
height = [1.68,1.80,1.75,1.68,1.76,1.74,1.98,1.89,1.77,1.75,1.86,1.64];
for index_subject = Data.subjects
    bodymass(index_subject-1) = Data.subject(index_subject).bodyweight/Data.g;
    BMI(index_subject-1) = bodymass(index_subject-1)/(height(index_subject-1))^2;
end
mean(BMI);
std(BMI);
[BMI_sort, BMI_index] = sort(BMI);
mean(bodymass);
std(bodymass);

%% Figure 3 GRF
figure(100+plot_subject)
close


COPcolor1 = [0.5,0.5,0.5];
lgrey = [0.8,0.8,0.8];

figure(100+plot_subject)

plot(t_kp,(D.Fx-D.Fx(end))/Data.subject(plot_subject).bodyweight*100,'-.')
hold on
plot(t_kp,(D.Fy-D.Fy(end))/Data.subject(plot_subject).bodyweight*100,'--')
plot(t_kp,D.Fz/Data.subject(plot_subject).bodyweight*100,'-')
temp_ylim = get(gca,'YLim');

h_legend = legend({'$F_X$','$F_Y$','$F_Z$'},'Location','best');
%         set(h_legend,'Position' ,[0.1512 0.7431 0.1210 0.1632])
h_legend.Box = 'off';

drawnow

plot([0.08 0.04 +0.04 -0.005 -0.005 +0.04 0.04 ],[10 10 70 70 -20 -20 10 ],'k:','LineWidth',1.5);%,'Color',COPcolor1)
text(0.09,10,'shown in Figure 4','FontSize',11.7)
plot([-0.2,0.9],[0,0],'--','LineWidth',0.5,'Color',lgrey)

%%% line through origin
plot([0,0],temp_ylim+[-10 10],'--','LineWidth',0.5,'Color',lgrey)
set(h_legend,'Position' ,[0.7632 0.6731 0.1210 0.1780])

% grid on
xlim([-0.10 +0.85])
ylim([temp_ylim(1)-10 temp_ylim(2)+10])
ylabel('Ground reaction force (in \% bodyweight)')
xlabel('Time (in s)')

%%% save figure to eps file
h_figure = gcf;
set(h_figure,'Position' ,[263   413   634   385] )

% eps fix
if strcmp(toggle_save,'on')
    set(h_figure,'Units','points')
    set(h_figure,'PaperUnits','points')
    temp_size = get(h_figure,'Position');
    temp_size = temp_size(3:4);
    set(h_figure,'PaperSize',temp_size)
    set(h_figure,'PaperPosition',[0,0,temp_size(1),temp_size(2)])
    file_name = strcat(folder,'\Fig3.pdf');
    
    %             file_name = strcat(folder,'\Force_s',num2str(i-1),'.eps');
    saveas(h_figure,file_name,'pdf')
end

if strcmp(toggle_save,'export')
    export_fig(strcat(folder,'\Fig3.pdf'))
end

Data.subject(plot_subject).bodyweight;

%% Figure 4 GRF and foot ankle deformation

COPcolor1 = [0.5,0.5,0.5];

figure(651+plot_subject)
close

figure(651+plot_subject)

subplot(211)
plot(t_kp*1000,(D.Fx-D.Fx(end))/Data.subject(plot_subject).bodyweight*100,'-.')
hold on
plot(t_kp*1000,(D.Fy-D.Fy(end))/Data.subject(plot_subject).bodyweight*100,'--')
plot(t_kp*1000,D.Fz/Data.subject(plot_subject).bodyweight*100,'-')
h_figure = gcf;
temp_ylim = [-20 80];
plot([t_kp(D.t_kp_index_0)*1000,t_kp(D.t_kp_index_0)*1000,t_kp(D.t_kp_index_0)*1000],[temp_ylim(1)-10,D.Fz( D.t_kp_index_0)/Data.subject(plot_subject).bodyweight*100,temp_ylim(2)+10],':','LineWidth',1.5,'Color',COPcolor1)
plot([t_kp(D.t_kp_index_peak)*1000,t_kp(D.t_kp_index_peak)*1000,t_kp(D.t_kp_index_peak)*1000],[temp_ylim(1)-10,D.Fz( D.t_kp_index_peak)/Data.subject(plot_subject).bodyweight*100,temp_ylim(2)+10],':','LineWidth',1.5,'Color',COPcolor1)
plot([t_kp(D.t_kp_index_zd0_COP)*1000,t_kp(D.t_kp_index_zd0_COP)*1000,t_kp(D.t_kp_index_zd0_COP)*1000],[temp_ylim(1)-10,D.Fz( D.t_kp_index_zd0_COP)/Data.subject(plot_subject).bodyweight*100,temp_ylim(2)+10],':','LineWidth',1.5,'Color',COPcolor1)
plot([-10,50],[0,0],'--','LineWidth',0.5,'Color',lgrey)

ylabel('Ground reaction force (in \% bodyweight)')
xlim([-0.005*1000 +0.04*1000])
ylim(temp_ylim-[0 10])
h_legend = legend({'$F_X$','$F_Y$','$F_Z$'},'Interpreter','Latex','Location','best');
        set(h_legend,'Position' ,[0.2747 0.8468 0.0980 0.0700])
h_legend.Box = 'off';
ticks_x = get(gca,'XTick');
ticks_xlabel = get(gca,'XTickLabel');
h_temp1 = gca;
[hx,hy] = format_ticks(gca,{'$t_0$','$t_p$','$t_e$'},[],[0,t_kp(D.t_kp_index_peak)*1000,t_opto(D.t_opto_filtered_index_zd0_COP)*1000],[],[],[],-10);

h_temp1.XTick = ticks_x;
h_temp1.Position = h_temp1.Position + [0 -0.05 0 0.09];



subplot(212)

plot(t_opto(D.t_opto_filtered_index_0:end)*1000,(D.x_COP_filtered(D.t_opto_filtered_index_0:end)-D.x_COP_filtered(D.t_opto_filtered_index_0))*1000,'-.')
hold on

plot(t_opto(D.t_opto_filtered_index_0:end)*1000,(D.y_COP_filtered(D.t_opto_filtered_index_0:end)-D.y_COP_filtered(D.t_opto_filtered_index_0))*1000,'--')
plot(t_opto(D.t_opto_filtered_index_0:end)*1000,(D.z_COP_filtered(D.t_opto_filtered_index_0:end)-D.z_COP_filtered(D.t_opto_filtered_index_0))*1000,'-')

xlim([-0.005*1000 +0.04*1000])

temp_ylim = get(gca,'YLim');

plot([t_kp(D.t_kp_index_0)*1000,t_kp(D.t_kp_index_0)*1000,t_kp(D.t_kp_index_0)*1000],[temp_ylim(1)-10,D.Fz( D.t_kp_index_0)/Data.subject(plot_subject).bodyweight*100,temp_ylim(2)+10],':','LineWidth',1.5,'Color',COPcolor1)
plot([t_kp(D.t_kp_index_peak)*1000,t_kp(D.t_kp_index_peak)*1000,t_kp(D.t_kp_index_peak)*1000],[temp_ylim(1)-10,D.Fz( D.t_kp_index_peak)/Data.subject(plot_subject).bodyweight*100,temp_ylim(2)+10],':','LineWidth',1.5,'Color',COPcolor1)

plot([t_opto(D.t_opto_filtered_index_zd0_COP)*1000,t_opto(D.t_opto_filtered_index_zd0_COP)*1000,t_opto(D.t_opto_filtered_index_zd0_COP)*1000],[temp_ylim(1)-10,D.y_COP_filtered(D.t_opto_filtered_index_zd0_COP),temp_ylim(2)+10],':','LineWidth',1.5,'Color',COPcolor1)
plot([-10,50],[0,0],'--','LineWidth',0.5,'Color',lgrey)

ylim([-15 25])

drawnow

ticks_x = get(gca,'XTick');
ticks_xlabel = get(gca,'XTickLabel');

% grid on
ylabel('Foot-ankle deformation (in mm)')
xlabel('Time (in ms)')
h_legend2 = legend({'$S_X$','$S_Y$','$S_Z$'},'interpreter','latex','Location','best');
h_legend2.Position = [0.2747 0.3730 0.0980 0.0700];
h_legend2.Box = 'off';

h_temp = gca;
[hx,hy] = format_ticks(gca,{'$t_0$','$t_p$','$t_e$'},[],[0,t_kp(D.t_kp_index_peak)*1000,t_opto(D.t_opto_filtered_index_zd0_COP)*1000],[],[],[],-1.08);

h_temp.XTick = ticks_x;
h_temp.XTickLabel = ticks_xlabel;

h_temp.Position = h_temp.Position + [0 -0.05 0 0.09];
%%% save figure to eps file

h_figure = gcf;
set(h_figure,'Position' ,[365.00 225.00 632.00 713.00] )

% eps fix
if strcmp(toggle_save,'on')
    set(h_figure,'Units','points')
    set(h_figure,'PaperUnits','points')
    temp_size = get(h_figure,'Position');
    temp_size = temp_size(3:4);
    set(h_figure,'PaperSize',temp_size)
    set(h_figure,'PaperPosition',[0,0,temp_size(1),temp_size(2)])
    
    file_name = strcat(folder,'\Results2a','.pdf');
    saveas(h_figure,file_name,'pdf')
end

if strcmp(toggle_save,'export')
    export_fig(strcat(folder,'\Fig4.pdf'))
end

%% Figure 5 COP

D = Data.subject(plot_subject).n(plot_run);

COPcolor1 = [0.6,0.6,0.6];
COPcolor2 = [0.8,0.8,0.8];

heelcolor1 = [182/255, 1, 87/255];
% heelcolor2 = [126/255, 47/255, 142/255]; % purple
heelcolor2 = [226/255, 161/255, 96/255];

begin_index = 1;

figure(20+plot_subject)
% clf

figure(20+plot_subject)

subplot(2,2,2)
hold off
% plot3(D.x_p(1:begin_index)*1E3,D.y_p(1:begin_index)*1E3,D.z_p(1:begin_index)*1E3,'.','Color',COPcolor2,'Linewidth',0.5,'MarkerSize',10)
% hold on

% plot3(D.x_p(begin_index)*1E3,D.y_p(begin_index)*1E3,D.z_p(begin_index)*1E3,'r.')
% plot3(D.x_p(end)*1E3,D.y_p(end)*1E3,D.z_p(end)*1E3,'r.')

plot3(D.x_COP_p*1E3,D.y_COP_p*1E3,D.z_COP_p*1E3,'.','Color',heelcolor1,'MarkerSize',36)
hold on
plot3(D.x_p(begin_index:end)*1E3,D.y_p(begin_index:end)*1E3,D.z_p(begin_index:end)*1E3,'.','Color',COPcolor1,'Linewidth',0.5,'MarkerSize',10)
plot3(D.x_COP_p*1E3+10,D.y_COP_p*1E3-10,D.z_COP_p*1E3-10,'.','Color',heelcolor2,'MarkerSize',25)


xlabel('x (in mm)')
ylabel('y (in mm)')
zlabel('z (in mm)')
axis equal
axislimit = 15;
xlim([D.x_COP_p*1E3-axislimit D.x_COP_p*1E3+axislimit])
ylim([D.y_COP_p*1E3-axislimit D.y_COP_p*1E3+axislimit])
zlim([D.z_COP_p*1E3-axislimit D.z_COP_p*1E3+axislimit])
grid on

subplot(2,2,3)
hold off
% plot(D.y_p(1:begin_index)*1E3,D.z_p(1:begin_index)*1E3,'.','Color',COPcolor2,'Linewidth',0.5,'MarkerSize',10)
% hold on
% plot(D.y_p*1E3,D.z_p*1E3,'.-','Color',COPcolor1,'Linewidth',0.5)
% plot(D.y_p(begin_index)*1E3,D.z_p(begin_index)*1E3,'r.')
% plot(D.y_p(end)*1E3,D.z_p(end)*1E3,'r.')

plot(D.y_COP_p*1E3,D.z_COP_p*1E3,'.','Color',heelcolor1,'MarkerSize',36)
hold on
plot(D.y_p(begin_index:end)*1E3,D.z_p(begin_index:end)*1E3,'.','Color',COPcolor1,'Linewidth',0.5,'MarkerSize',10)
plot(D.y_COP_p*1E3-10,D.z_COP_p*1E3-10,'.','Color',heelcolor2,'MarkerSize',25)

xlabel('y (in mm)')
ylabel('z (in mm)')

% grid on

axis equal
axislimit = 13;
xlim([D.y_COP_p*1E3-axislimit D.y_COP_p*1E3+axislimit])
ylim([D.z_COP_p*1E3-axislimit D.z_COP_p*1E3+axislimit])
set(gca,'XDir','reverse')

h_legend = legend({'COP','$H$','$Q$'},'Location','South');
h_legend.Box = 'off';
% h_legend.Position = [0.2279 0.2296 0.1396 0.1295];
h_legend.Position = [0.8358 0.8319 0.1396 0.1296];

subplot(2,2,1)
hold off
% plot(D.x_p(1:begin_index)*1E3,D.y_p(1:begin_index)*1E3,'.','Color',COPcolor2,'Linewidth',0.5,'MarkerSize',10)
% hold on
% plot(D.x_p(begin_index)*1E3,D.y_p(begin_index)*1E3,'r.')
% plot(D.x_p(end)*1E3,D.y_p(end)*1E3,'r.')

plot(D.x_COP_p*1E3,D.y_COP_p*1E3,'.','Color',heelcolor1,'MarkerSize',36)
hold on
plot(D.x_p(begin_index:end)*1E3,D.y_p(begin_index:end)*1E3,'.','Color',COPcolor1,'Linewidth',0.5,'MarkerSize',10)
plot(D.x_COP_p*1E3+10,D.y_COP_p*1E3-10,'.','Color',heelcolor2,'MarkerSize',25)


xlabel('x (in mm)')
ylabel('y (in mm)')

% grid on
axis equal
xlim([D.x_COP_p*1E3-axislimit D.x_COP_p*1E3+axislimit])
ylim([D.y_COP_p*1E3-axislimit D.y_COP_p*1E3+axislimit])

subplot(2,2,4)
hold off
% plot(D.x_p(1:begin_index)*1E3,D.y_p(1:begin_index)*1E3,'.','Color',COPcolor2,'Linewidth',0.5,'MarkerSize',10)
% hold on
% plot(D.x_p(begin_index)*1E3,D.y_p(begin_index)*1E3,'r.')
% plot(D.x_p(end)*1E3,D.y_p(end)*1E3,'r.')

plot(D.x_COP_p*1E3,D.z_COP_p*1E3,'.','Color',heelcolor1,'MarkerSize',36)
hold on
plot(D.x_p(begin_index:end)*1E3,D.z_p(begin_index:end)*1E3,'.','Color',COPcolor1,'Linewidth',0.5,'MarkerSize',10)
plot(D.x_COP_p*1E3+10,D.z_COP_p*1E3-10,'.','Color',heelcolor2,'MarkerSize',25)


xlabel('x (in mm)')
ylabel('z (in mm)')

% grid on
axis equal
xlim([D.x_COP_p*1E3-axislimit D.x_COP_p*1E3+axislimit])
ylim([D.z_COP_p*1E3-axislimit D.z_COP_p*1E3+axislimit])

% h_axis = gca;

set(gcf,'Position',[325 312 645 485])

if strcmp(toggle_save,'export')
    export_fig(strcat(folder,'\Fig5.pdf'))
end

%% Figure 6 Scatter foot-ankle deformation - energy loss, work Force subject specific
dgrey = [0.25,0.25,0.25];

figure(6)
close
figure(6)


matrix_Work_zd0_plot_normal = -Matrix.Work_zd0.*Matrix.n_normal;
matrix_heel_deformation_plot_normal = Matrix.heel_deformation.*Matrix.n_normal;


matrix_Work_zd0_plot = matrix_Work_zd0_plot_normal(Data.subjects,:);
matrix_heel_deformation_plot = matrix_heel_deformation_plot_normal(Data.subjects,:);
[~, BMI_index] = sort(bodymass);

cmap = parula(13);
vector_hd_pr_plot = [8.8,9.8,9.7,8.5,10.3,10.3,10,11.8,10.8,11.8,11.5,12.1,12.8,12.1]'*1E-3;
vector_E_pr_plot = [0.65,0.85,0.95,1.1,1.3,1.45,1.55,1.85,1.95,2.45,2.55,3.3,4.2,5.3]';


h1 = plot(vector_E_pr_plot(1),vector_hd_pr_plot(1)*1000);
hold on
h1.HandleVisibility = 'off';

color_plot = ['>';'.';'x';'*';'o';'p';'+';'^';'h';'d';'v';'s';'<';'.'];
legendstring = {};
colorindex = 0;
for index_subject = BMI_index
    colorindex = colorindex + 1;
    if colorindex == 2
        mark = 200;
    else
        mark = 36;
    end
    figure(6)
    scatter(matrix_Work_zd0_plot(index_subject,matrix_Work_zd0_plot(index_subject,:)~=0),matrix_heel_deformation_plot(index_subject,matrix_Work_zd0_plot(index_subject,:)~=0)*1000,mark,cmap(colorindex,:),color_plot(colorindex,:),'Linewidth',1.5)
    hold on
    legendstring = [legendstring, strcat('$M_{',num2str(index_subject),'} = \ $',num2str(bodymass(index_subject),'%3.0f'),'\ kg')];
end
hold on
scatter(vector_E_pr_plot,vector_hd_pr_plot*1000,60,dgrey,'s','filled')
scatter(vector_E_pr_plot,vector_hd_pr_plot*1000,60,'w','s','filled','visible','off')
h_legend = legend([legendstring,'De Clercq','et al. (1994)'],'interpreter','latex','Location','EastOutside');

% h_legend = legend([legendstring,'De Clercq et al. (1994)'],'interpreter','latex','Location','East');
h_legend.Box = 'off';
xlabel('Energy loss per heel strike (in J)')
ylabel({'Maximum foot-ankle','compression (in mm)'})
xlim([0 6])

h_figure = gcf;
set(h_figure,'Position' ,[444 439 820 348] )

% eps fix

if strcmp(toggle_save,'on')
    
    set(h_figure,'Units','points')
    set(h_figure,'PaperUnits','points')
    temp_size = get(h_figure,'Position');
    temp_size = temp_size(3:4);
    set(h_figure,'PaperSize',temp_size)
    set(h_figure,'PaperPosition',[0,0,temp_size(1),temp_size(2)])
    %
    file_name = strcat(folder,'\Results3','.pdf');
    saveas(h_figure,file_name,'pdf')
end

if strcmp(toggle_save,'export')
    export_fig(strcat(folder,'\Fig6.pdf'))
end

%% Figure 7 Heel point acceleration

figure(700+plot_subject)
close

figure(700+plot_subject)
ax1 = axes;
plot(t_opto*1000,D.zdd_COP_filtered,':','Color',[0.9290 0.6940 0.1250])
hold on
xlim([-0.005*1000 +0.04*1000])
h_figure = gcf;
temp_ylim = get(gca,'YLim');
plot([t_kp(D.t_kp_index_0)*1000,t_kp(D.t_kp_index_0)*1000,t_kp(D.t_kp_index_0)*1000],[temp_ylim(1)-10,D.Fz( D.t_kp_index_0)/Data.subject(plot_subject).bodyweight*100,temp_ylim(2)+10],':','LineWidth',1.5,'Color',COPcolor1)
plot([t_kp(D.t_kp_index_peak)*1000,t_kp(D.t_kp_index_peak)*1000,t_kp(D.t_kp_index_peak)*1000],[temp_ylim(1)-10,D.Fz( D.t_kp_index_peak)/Data.subject(plot_subject).bodyweight*100,temp_ylim(2)+10],':','LineWidth',1.5,'Color',COPcolor1)
plot([t_kp(D.t_kp_index_zd0_COP)*1000,t_kp(D.t_kp_index_zd0_COP)*1000,t_kp(D.t_kp_index_zd0_COP)*1000],[temp_ylim(1)-10,D.Fz( D.t_kp_index_zd0_COP)/Data.subject(plot_subject).bodyweight*100,temp_ylim(2)+10],':','LineWidth',1.5,'Color',COPcolor1)
% plot([-10,50],[0,0],'--','LineWidth',0.5,'Color',lgrey)
plot([-10,50],[-9.81,-9.81],':','LineWidth',1.5,'Color',COPcolor1)
ylabel('Heel point acceleration (in m/s$^2$)')
xlim([-0.005*1000 +0.04*1000])
% ylim(temp_ylim)
ylim([-20 70])
ticks_x = get(gca,'XTick');
ticks_xlabel = get(gca,'XTickLabel');
ticks_y = get(gca,'YTick');
ticks_ylabel = get(gca,'YTickLabel');
h_temp1 = gca;
[~,~] = format_ticks(gca,[],{'$g$'},[],[-9.81],[],[],[]);
h_temp1.XTick = ticks_x;
h_temp1.YTick = ticks_y(1:2:end);
h_temp1.YTickLabel = ticks_ylabel(1:2:end);


ax2 = axes;
plot(t_opto*1000,D.zdd_COP_filtered,'-','Color',[0.9290 0.6940 0.1250])
hold on
ax2.YAxisLocation = 'right';
plot(t_kp*1000,D.Fz/Data.subject(plot_subject).bodyweight*100,':','Color',[0.9290 0.6940 0.1250])
h_figure = gcf;
temp_ylim = [-20 80];
plot([-10,50],[-9.81,-9.81],':','LineWidth',1.5,'Color',COPcolor1)
plot([t_kp(D.t_kp_index_0)*1000,t_kp(D.t_kp_index_0)*1000,t_kp(D.t_kp_index_0)*1000],[temp_ylim(1)-10,D.Fz( D.t_kp_index_0)/Data.subject(plot_subject).bodyweight*100,temp_ylim(2)+10],':','LineWidth',1.5,'Color',COPcolor1)
plot([t_kp(D.t_kp_index_peak)*1000,t_kp(D.t_kp_index_peak)*1000,t_kp(D.t_kp_index_peak)*1000],[temp_ylim(1)-10,D.Fz( D.t_kp_index_peak)/Data.subject(plot_subject).bodyweight*100,temp_ylim(2)+10],':','LineWidth',1.5,'Color',COPcolor1)
plot([t_kp(D.t_kp_index_zd0_COP)*1000,t_kp(D.t_kp_index_zd0_COP)*1000,t_kp(D.t_kp_index_zd0_COP)*1000],[temp_ylim(1)-10,D.Fz( D.t_kp_index_zd0_COP)/Data.subject(plot_subject).bodyweight*100,temp_ylim(2)+10],':','LineWidth',1.5,'Color',COPcolor1)
plot([-10,50],[0,0],'--','LineWidth',0.5,'Color',lgrey)

ylabel('Ground reaction force (in \% bodyweight)')
xlim([-0.005*1000 +0.04*1000])
ylim(temp_ylim-[0 10])
h_legend = legend({'$\ddot{Z}_H$','$F_Z$'},'Interpreter','Latex','Location','best');
set(h_legend,'Position' ,[0.2747 0.7364 0.1221 0.1161])
h_legend.Box = 'off';
ticks_x = get(gca,'XTick');
ticks_xlabel = get(gca,'XTickLabel');
h_temp1 = gca;
[~,~] = format_ticks(gca,{'$t_0$','$t_p$','$t_e$'},[],[0,t_kp(D.t_kp_index_peak)*1000,t_opto(D.t_opto_filtered_index_zd0_COP)*1000],[],[],[],-1.075);

h_temp1.XTick = ticks_x;
h_temp1.XTickLabel = ticks_xlabel;
h_temp1.YTick = ticks_y(1:2:end);
h_temp1.YTickLabel = ticks_ylabel(1:2:end);
xlabel('Time (in ms)')


%%% save figure to eps file

h_figure = gcf;
set(h_figure,'Position' ,[263   413   634   385] )
ax1.Position = ax2.Position;

% eps fix
if strcmp(toggle_save,'on')
    set(h_figure,'Units','points')
    set(h_figure,'PaperUnits','points')
    temp_size = get(h_figure,'Position');
    temp_size = temp_size(3:4);
    set(h_figure,'PaperSize',temp_size)
    set(h_figure,'PaperPosition',[0,0,temp_size(1),temp_size(2)])
    
    file_name = strcat(folder,'\Results7','.pdf');
    saveas(h_figure,file_name,'pdf')
end

if strcmp(toggle_save,'export')
    export_fig(strcat(folder,'\Fig7.pdf'))
end

%% Figure 8 Instantaneous effective foot mass

% M_eff_COP_instant

figure(653+plot_subject)
close

figure(653+plot_subject)
plot(t_kp(D.t_kp_index_0:D.t_kp_index_zd0_COP)*1000,D.M_eff_COP_instant(1:length((D.t_kp_index_0:D.t_kp_index_zd0_COP)))/(Data.subject(plot_subject).bodyweight/Data.g)*100,'-','Color',[0.9290 0.6940 0.1250])
hold on
h_figure = gcf;
h_axes = gca
temp_ylim = h_axes.YLim;
plot([t_kp(D.t_kp_index_0)*1000,t_kp(D.t_kp_index_0)*1000,t_kp(D.t_kp_index_0)*1000],[temp_ylim(1)-10,D.Fz( D.t_kp_index_0)/Data.subject(plot_subject).bodyweight*100,temp_ylim(2)+10],':','LineWidth',1.5,'Color',COPcolor1)
plot([t_kp(D.t_kp_index_peak)*1000,t_kp(D.t_kp_index_peak)*1000,t_kp(D.t_kp_index_peak)*1000],[temp_ylim(1)-10,D.Fz( D.t_kp_index_peak)/Data.subject(plot_subject).bodyweight*100,temp_ylim(2)+10],':','LineWidth',1.5,'Color',COPcolor1)
plot([t_kp(D.t_kp_index_zd0_COP)*1000,t_kp(D.t_kp_index_zd0_COP)*1000,t_kp(D.t_kp_index_zd0_COP)*1000],[temp_ylim(1)-10,D.Fz( D.t_kp_index_zd0_COP)/Data.subject(plot_subject).bodyweight*100,temp_ylim(2)+10],':','LineWidth',1.5,'Color',COPcolor1)
plot([-10,50],[0,0],'--','LineWidth',0.5,'Color',lgrey)

ylabel('Effective foot mass (in \% bodymass)')
xlim([-0.005*1000 +0.04*1000])
ylim(temp_ylim)
h_legend = legend({'$M_{eff}(t)$'},'Interpreter','Latex','Location','best');
set(h_legend,'Position' ,[0.2747 0.7364 0.1221 0.1161])

h_legend.Box = 'off';
ticks_x = get(gca,'XTick');
ticks_xlabel = get(gca,'XTickLabel');
h_temp1 = gca;
% [hx,hy] = format_ticks(gca,{'$t_0$','$t_p$','$t_e$'},[],[0,t_kp(D.t_kp_index_peak)*1000,t_opto(D.t_opto_filtered_index_zd0_COP)*1000],[],[],[],-10);
xlabel('Time (in ms)')

h_axes2 = axes;
xlim([-0.005*1000 +0.04*1000])
[~,~] = format_ticks(gca,{'$t_0$','$t_p$','$t_e$'},[{'.'}],[0,t_kp(D.t_kp_index_peak)*1000,t_opto(D.t_opto_filtered_index_zd0_COP)*1000],[0],[],[],-1.075);

h_temp1.XTick = ticks_x;
h_temp1.Position = h_temp1.Position + [0 -0.05 0 0.09];

h_axes2.Position = [0.1300    0.1273    0.7750    0.7977];
h_axes.Position = h_axes2.Position;
axes(h_axes)

%%% save figure to eps file

h_figure = gcf;
set(h_figure,'Position' ,[263   413   634   385] )

% eps fix
if strcmp(toggle_save,'on')
    set(h_figure,'Units','points')
    set(h_figure,'PaperUnits','points')
    temp_size = get(h_figure,'Position');
    temp_size = temp_size(3:4);
    set(h_figure,'PaperSize',temp_size)
    set(h_figure,'PaperPosition',[0,0,temp_size(1),temp_size(2)])
    
    file_name = strcat(folder,'\Fig8','.pdf');
    saveas(h_figure,file_name,'pdf')
end

if strcmp(toggle_save,'export')
    export_fig(strcat(folder,'\Fig8.pdf'))
end

%% Table 1
[mean_total,mean_min,mean_max,within_std,between_std]=function_GetStats(Matrix.v_average,Matrix.n_normal,Data.subjects);
v_stats = [mean_total,mean_min,mean_max,between_std,within_std]

[mean_total,mean_min,mean_max,within_std,between_std]=function_GetStats(Matrix.v_average,Matrix.n_slow,Data.subjects);
v_stats_slow = [mean_total,mean_min,mean_max,between_std,within_std]

[mean_total,mean_min,mean_max,within_std,between_std]=function_GetStats(Matrix.v_average./Matrix.step_length,Matrix.n_normal,Data.subjects);
step_stats = [mean_total,mean_min,mean_max,between_std,within_std]

[mean_total,mean_min,mean_max,within_std,between_std]=function_GetStats(Matrix.step_length,Matrix.n_normal,Data.subjects);
steplength_stats = [mean_total,mean_min,mean_max,between_std,within_std]

[mean_total,mean_min,mean_max,within_std,between_std]=function_GetStats(Matrix.heel_deformation_x_peak,Matrix.n_normal,Data.subjects);
heeldeformx_peak_stats = [mean_total,mean_min,mean_max,between_std,within_std]*1000
[mean_total,mean_min,mean_max,within_std,between_std]=function_GetStats(Matrix.heel_deformation_y_peak,Matrix.n_normal,Data.subjects);
heeldeformy_peak_stats = [mean_total,mean_min,mean_max,between_std,within_std]*1000
[mean_total,mean_min,mean_max,within_std,between_std]=function_GetStats(Matrix.heel_deformation_z_peak,Matrix.n_normal,Data.subjects);
heeldeformz_peak_stats = [mean_total,mean_min,mean_max,between_std,within_std]*1000
[mean_total,mean_min,mean_max,within_std,between_std]=function_GetStats(Matrix.total_heel_deformation_peak,Matrix.n_normal,Data.subjects);
total_heeldeform_peak_stats = [mean_total,mean_min,mean_max,between_std,within_std]*1000

[mean_total,mean_min,mean_max,within_std,between_std]=function_GetStats(Matrix.heel_deformation_x,Matrix.n_normal,Data.subjects);
heeldeformx_stats = [mean_total,mean_min,mean_max,between_std,within_std]*1000
[mean_total,mean_min,mean_max,within_std,between_std]=function_GetStats(Matrix.heel_deformation_y,Matrix.n_normal,Data.subjects);
heeldeformy_stats = [mean_total,mean_min,mean_max,between_std,within_std]*1000
[mean_total,mean_min,mean_max,within_std,between_std]=function_GetStats(Matrix.heel_deformation_z,Matrix.n_normal,Data.subjects);
heeldeformz_stats = [mean_total,mean_min,mean_max,between_std,within_std]*1000
[mean_total,mean_min,mean_max,within_std,between_std]=function_GetStats(Matrix.total_heel_deformation,Matrix.n_normal,Data.subjects);
total_heeldeform_stats = [mean_total,mean_min,mean_max,between_std,within_std]*1000

[mean_total,mean_min,mean_max,within_std,between_std]=function_GetStats(Matrix.zd_0_adjusted,Matrix.n_normal,Data.subjects);
zd0_stats = [mean_total,mean_min,mean_max,between_std,within_std]

[mean_total,mean_min,mean_max,within_std,between_std]=function_GetStats(Matrix.zd_peak_adjusted,Matrix.n_normal,Data.subjects);
zdpeak_stats = [mean_total,mean_min,mean_max,between_std,within_std]

[mean_total,mean_min,mean_max,within_std,between_std]=function_GetStats(Matrix.zd_peak_adjusted-Matrix.zd_0_adjusted,Matrix.n_normal,Data.subjects);
zddif_stats = [mean_total,mean_min,mean_max,between_std,within_std]

[mean_total,mean_min,mean_max,within_std,between_std]=function_GetStats(Matrix.t_peak_adjusted,Matrix.n_normal,Data.subjects);
tpeak_stats = [mean_total,mean_min,mean_max,between_std,within_std]*1000

[mean_total,mean_min,mean_max,within_std,between_std]=function_GetStats(Matrix.t_zd0_adjusted,Matrix.n_normal,Data.subjects);
tzd0_stats = [mean_total,mean_min,mean_max,between_std,within_std]*1000

[mean_total,mean_min,mean_max,within_std,between_std]=function_GetStats(Matrix.t_zd0_adjusted-Matrix.t_peak_adjusted,Matrix.n_normal,Data.subjects);
tdiff_stats = [mean_total,mean_min,mean_max,between_std,within_std]*1000

[mean_total,mean_min,mean_max,within_std,between_std]=function_GetStats(Matrix.M_eff_zd0_COP_percentageMB,Matrix.n_normal,Data.subjects);
Meff_stats = [mean_total,mean_min,mean_max,between_std,within_std]

[mean_total,mean_min,mean_max,within_std,between_std]=function_GetStats(Matrix.M_eff_peak_COP_percentageMB,Matrix.n_normal,Data.subjects);
Meff_peak_stats = [mean_total,mean_min,mean_max,between_std,within_std]

[mean_total,mean_min,mean_max,within_std,between_std]=function_GetStats(Matrix.M_eff_Chi_COP_percentageMB,Matrix.n_normal,Data.subjects);
MeffChi_stats = [mean_total,mean_min,mean_max,between_std,within_std]

[mean_total,mean_min,mean_max,within_std,between_std,EMeff_means]=function_GetStats(Matrix.E_M_eff_zd0,Matrix.n_normal,Data.subjects);
EMeff_stats = [mean_total,mean_min,mean_max,between_std,within_std]

[mean_total,mean_min,mean_max,within_std,between_std,EMeffPeak_means]=function_GetStats(Matrix.E_M_eff_peak,Matrix.n_normal,Data.subjects);
EMeffPeak_stats = [mean_total,mean_min,mean_max,between_std,within_std]

[mean_total,mean_min,mean_max,within_std,between_std,EMeffChi_means]=function_GetStats(Matrix.E_M_eff_Chi,Matrix.n_normal,Data.subjects);
EMeffChi_stats = [mean_total,mean_min,mean_max,between_std,within_std]

[mean_total,mean_min,mean_max,within_std,between_std]=function_GetStats(Matrix.E_metabolic_Kuo,Matrix.n_normal,Data.subjects);
E_metabolic_Kuo_stats = [mean_total,mean_min,mean_max,between_std,within_std]

[mean_total,mean_min,mean_max,within_std,between_std,Workz_peak_means]=function_GetStats(Matrix.Work_z_peak,Matrix.n_normal,Data.subjects);
Workz_peak_stats = [mean_total,mean_min,mean_max,between_std,within_std]
[mean_total,mean_min,mean_max,within_std,between_std,Workx_peak_means]=function_GetStats(Matrix.Work_x_peak,Matrix.n_normal,Data.subjects);
Workx_peak_stats = [mean_total,mean_min,mean_max,between_std,within_std]
[mean_total,mean_min,mean_max,within_std,between_std,Worky_peak_means]=function_GetStats(Matrix.Work_y_peak,Matrix.n_normal,Data.subjects);
Worky_peak_stats = [mean_total,mean_min,mean_max,between_std,within_std]

[mean_total,mean_min,mean_max,within_std,between_std,Work_peak_means]=function_GetStats(Matrix.Work_peak,Matrix.n_normal,Data.subjects);
Work_peak_stats = [mean_total,mean_min,mean_max,between_std,within_std]

[mean_total,mean_min,mean_max,within_std,between_std,Workz_means]=function_GetStats(Matrix.Work_zd0,Matrix.n_normal,Data.subjects);
Workz_stats = [mean_total,mean_min,mean_max,between_std,within_std]
[mean_total,mean_min,mean_max,within_std,between_std,Workx_means]=function_GetStats(Matrix.Work_x,Matrix.n_normal,Data.subjects);
Workx_stats = [mean_total,mean_min,mean_max,between_std,within_std]
[mean_total,mean_min,mean_max,within_std,between_std,Worky_means]=function_GetStats(Matrix.Work_y,Matrix.n_normal,Data.subjects);
Worky_stats = [mean_total,mean_min,mean_max,between_std,within_std]

[mean_total,mean_min,mean_max,within_std,between_std,Work_means]=function_GetStats(Matrix.Work,Matrix.n_normal,Data.subjects);
Work_stats = [mean_total,mean_min,mean_max,between_std,within_std]

[mean_total,mean_min,mean_max,within_std,between_std]=function_GetStats(Matrix.Work,Matrix.n_slow,Data.subjects);
Work_stats_slow = [mean_total,mean_min,mean_max,between_std,within_std]

[mean_total,mean_min,mean_max,within_std,between_std]=function_GetStats(Matrix.Work-Matrix.E_M_eff_zd0,Matrix.n_normal,Data.subjects);
Ediff_stats = [mean_total,mean_min,mean_max,between_std,within_std]

%% Extra variables
-Work_stats(1)/(0.2*E_metabolic_Kuo_stats(1))
-Work_stats_slow(1)/(0.2*E_metabolic_Kuo_stats(1))

(tzd0_stats(1)/1000)/(1/step_stats(1))

%% Statistics

%%% Kolmogorov-Smirnov tests
[hW,pW] = kstest(Work_means)
[hWz,pWz] = kstest(Workz_means)
[hE,pE] = kstest(EMeff_means)
[hp,pp] = kstest(EMeffPeak_means)
[hW,pW] = kstest(EMeffChi_means)

%%% => sets assumed to have non-normal distribution: 
%%% Friedman and signed rank tests
% [h,p,ci,stats] = ttest(Workz_means,EMeff_means)
[p12,h12,stat12] = signrank(Work_means,Workz_means,'method','approximate');%
[p13,h13,stat13] = signrank(Work_means,EMeff_means,'method','approximate');
[p14,h14,stat14] = signrank(Work_means,EMeffPeak_means,'method','approximate');
[p15,h15,stat15] = signrank(Work_means,EMeffChi_means,'method','approximate');
[p23,h23,stat23] = signrank(Workz_means,EMeff_means,'method','approximate');
[p24,h24,stat24] = signrank(Workz_means,EMeffPeak_means,'method','approximate');
[p25,h25,stat25] = signrank(Workz_means,EMeffChi_means,'method','approximate');
[p34,h34,stat34] = signrank(EMeff_means,EMeffPeak_means,'method','approximate');
[p35,h35,stat35] = signrank(EMeff_means,EMeffChi_means,'method','approximate');
[p45,h45,stat45] = signrank(EMeffPeak_means,EMeffChi_means,'method','approximate');


[p,tbl,stats] = friedman([Work_means' Workz_means' EMeff_means' EMeffPeak_means' EMeffChi_means'])
[p,tbl,stats] = friedman([Work_means' Workz_means' EMeff_means' EMeffChi_means'])
[p12,p13,p14,p15,p23,p24,p25,p34,p35,p45]
[h12,h13,h14,h15,h23,h24,h25,h34,h35,h45]

siglevel = 0.05/5

%% Figure X scatter t_e t_p subject specific
%%% extra figure, not included in manuscript

[mean_total,mean_min,mean_max,within_std,between_std]=function_GetStats(Matrix.t_zd0_adjusted./Matrix.t_peak_adjusted,Matrix.n_normal,Data.subjects);
tpeak_scale_stats = [mean_total,mean_min,mean_max,between_std,within_std];

[mean_total,mean_min,mean_max,within_std,between_std]=function_GetStats(Matrix.t_zd0_adjusted-Matrix.t_peak_adjusted,Matrix.n_normal,Data.subjects);
diff_tzd0_tpeak_stats = [mean_total,mean_min,mean_max,between_std,within_std]*1000;

[mean_total,mean_min,mean_max,within_std,between_std]=function_GetStats(Matrix.t_zd0_adjusted./Matrix.t_peak_adjusted,Matrix.n_slow,Data.subjects);
tpeak_scale_stats = [mean_total,mean_min,mean_max,between_std,within_std];

[mean_total,mean_min,mean_max,within_std,between_std]=function_GetStats(Matrix.t_zd0_adjusted-Matrix.t_peak_adjusted,Matrix.n_slow,Data.subjects);
diff_tzd0_tpeak_stats = [mean_total,mean_min,mean_max,between_std,within_std]*1000;

dgrey = [0.25,0.25,0.25];
lgrey = [0.75,0.75,0.75];

figure(1513)
close
figure(1513)


% clf
matrix_Work_zd0_plot_normal = Matrix.t_zd0_adjusted.*Matrix.n_normal;
matrix_heel_deformation_plot_normal = Matrix.t_peak_adjusted.*Matrix.n_normal;

matrix_Work_zd0_plot = matrix_Work_zd0_plot_normal(Data.subjects,:);
matrix_heel_deformation_plot = matrix_heel_deformation_plot_normal(Data.subjects,:);
[BMI_sort, BMI_index] = sort(bodymass);

cmap = parula(13);

%%%
t_e_vec = matrix_Work_zd0_plot(isfinite(matrix_Work_zd0_plot(:)));
t_p_vec = matrix_heel_deformation_plot(isfinite(matrix_heel_deformation_plot(:)));
[FO, G, O] = fit(t_e_vec,t_p_vec,'poly1');
%%%
color_plot = ['>';'.';'x';'*';'o';'p';'+';'^';'h';'d';'v';'s';'<';'.'];
legendstring = {};
colorindex = 0;
for index_subject = BMI_index
    colorindex = colorindex + 1;
    if colorindex == 2
        mark = 200;
    else
        mark = 36;
    end
    figure(1513)
    scatter(matrix_heel_deformation_plot(index_subject,matrix_Work_zd0_plot(index_subject,:)~=0).*1000,matrix_Work_zd0_plot(index_subject,matrix_Work_zd0_plot(index_subject,:)~=0).*1000,mark,cmap(colorindex,:),color_plot(colorindex,:),'Linewidth',1.5)
    hold on
    legendstring = [legendstring, strcat('$M_{',num2str(index_subject),'} = \ $',num2str(bodymass(index_subject),'%3.0f'),'\ kg')];
end
hold on
plot([0,60/1.70],[0,60],'--','Color',dgrey)
h_legend = legend([legendstring,'line corresponding to $t_e$ = 1.70 $t_p$'],'interpreter','latex','Location','EastOutside');
h_legend.Box = 'off';
xlabel('Time to impact force peak $t_p$ (in ms)')
ylabel('Time to end of heel strike $t_e$ (in ms)')
axis equal

xlim([0 30])
ylim([0 60])

h_figure = gcf;
set(h_figure,'Position' ,[137         242        1000         348] )


if strcmp(toggle_save,'export')
    export_fig(strcat(folder,'\Results4.pdf'))
end

if strcmp(toggle_save,'on')
    
    set(h_figure,'Units','points')
    set(h_figure,'PaperUnits','points')
    temp_size = get(h_figure,'Position');
    temp_size = temp_size(3:4);
    set(h_figure,'PaperSize',temp_size)
    set(h_figure,'PaperPosition',[0,0,temp_size(1),temp_size(2)])
    %
    file_name = strcat(folder,'\Results4','.pdf');
    saveas(h_figure,file_name,'pdf')
end


