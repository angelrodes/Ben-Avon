clear
close all hidden

%% Core location 57.109150, -3.462608
% Converted to British National Grid (m)
x_core=311510;
y_core=802971+9.7; % see correction in Google Earth image


%% define X and Y with a polar mesh
disp(['Building DEM...'])
dd = [0.001,0.01,0.05:0.01:1,1.1:0.1:3,3.2:0.2:5,5.5:0.5:10,11:1:20,25:5:50,60:10:100,200:100:1000]; % space radial distances to get a good resolution of shielding
azaz = (0:.01:0.99)*2*pi; % azimuths
[az,dist] = meshgrid(azaz,dd);
xq = x_core+dist.*cos(az);
yq = y_core+dist.*sin(az);
disp(['Total points:' num2str(numel(xq)) ' (similar to a ' num2str(ceil(numel(xq)^0.5)) 'x' num2str(ceil(numel(xq)^0.5)) ' grid)'])

%% Import terrain-5-dtm-NJ10SW.xyz
% Load XYZ data from the file
data = importdata('terrain-5-dtm-NJ10SW.xyz');
% Extract x, y, and z coordinates
x = data(:, 1);
y = data(:, 2);
z = data(:, 3);
zq = griddata(x, y, z, xq, yq, 'natural'); % Interpolate the data onto the grid
z_core=griddata(xq,yq,zq,x_core,y_core);

%% Field measurements E-W (made in cm and converted to m)
x_relative=1/100*[0  -40 -41 -40-200 -40-300                  1    110  111     110+170 110+171     110+170+40  110+170+40+100                  ];  
z_relative=1/100*[0 0   -85 -85     -85-100/tan(25*2*pi/360) -123 -123 -123-60 -123-60 -123-60+120 -123-60+120 -123-60+120-100/tan(25*2*pi/360)];
x_field=x_core+x_relative;
y_field=y_core*ones(size(x_relative));
z_field=z_core+z_relative;

% incorporate this data
x2=[x;x_field'];
y2=[y;y_field'];
z2=[z;z_field'];
% zq = griddata(x2, y2, z2, xq, yq, 'natural'); % Interpolate the data onto the grid

%% S-N measurements on picture 20190806-1141-DSC_0136.jpg (in m)
% Digitalised using g3data (Frantz, 2000)
% https://manpages.ubuntu.com/manpages/jammy/man1/g3data.1.html
data = importdata('20190806-1141-DSC_0136.JPG.dat');
z_pic = data(:, 1)+z_core-max(data(:, 1)); % core at the top
y_pic_temp = data(:, 2);
y_pic=(y_pic_temp-min(y_pic_temp))*1.2+802955; % correct picture deformation (1.2) and place
x_pic=x_core*ones(size(y_pic)); % same longitude as core

% incorporate this data (cropping last points to avoid artifact)
x3=[x2;x_pic(1:37)];
y3=[y2;y_pic(1:37)];
z3=[z2;z_pic(1:37)];
zq = griddata(x3, y3, z3, xq, yq, 'natural'); % Interpolate the data onto the grid


%% Identify tors

% Calculate suface without tors (rough)
areas_with_tors= ( xq>311490 & xq< 311540 & yq>802960 & yq<803020 ) |( yq<802940 & xq>311450 & xq<311520);
z_no_tors=griddata(xq(~areas_with_tors),yq(~areas_with_tors),zq(~areas_with_tors),xq,yq,'natural');

% Select only areas with > 0.5 m of diference (fine)
areas_with_tors=zq-z_no_tors>0.5;
z_no_tors(~areas_with_tors)=zq(~areas_with_tors);

disp('Done.')

%% Samples and cosmognic concentrations
sample_names=[{'CNG19-02'},{'CNG19-01A'},{'CNG19-01B'},{'CNG19-01C'},{'CNG19-01D'},{'CNG19-01E'},{'CNG19-01F'}];
z_core=griddata(xq,yq,zq,x_core,y_core); % fine (re)place core at the surface
z_samples=z_core-1/100*[2 143.5 170 202.75 241 346 447];
x_samples=x_core*ones(size(z_samples));
y_samples=y_core*ones(size(z_samples));
C10=[1717234    387190   252335   157529   103638   45781   93634 ];
dC10=[42732    9776   7014   4491   3780   1730   3076 ];
C26=[17812505    3074041   1633063   993319   981233   448317   925777 ];
dC26=[633818    110093   67911   45210   51715   28841   58737 ];
C21=[48673916    92122045   71151798   64087845   40799371   22051193   30817597 ];
dC21=[1123580    1025012   688710   884258   893349   707819   829855 ];

% consider C21 concentration of the second lowest sample as average inheritance
% C21=C21-C21(end-1);
% disp('Minimum cosmogenic [^{21}Ne] considered (excess respect lowest sample)')

%% Surface production rates

% CRONUS input (https://hess.ess.washington.edu/math/v3/v3_age_in.html):
% CNG19-02  57.10915  -3.462608   1015 std  0 2.6  1 0  2019 ;
% CNG19-02 Be-10 quartz 1.7172E+06 4.2732E+04  07KNSTD ;
% CNG19-02 Al-26 quartz 1.7813E+07 6.3382E+05 KNSTD ;
% CNG19-02 Ne-21 quartz 4.8674E+07 1.1236E+06 CRONUS-A 3.38E+08 ; 

% CRONUS LSDn outputs:
% wrapper: 3.0.2; get_age: 3.0.2; muons: 1A, alpha = 1; validate: validate_v3_input.m - 3.0; consts: 2020-08-26 
% Nuclide       Age (yr) Interr (yr) Exterr (yr)
% Be-10 (qtz) 	154248 	3990 	10268 	
% Al-26 (qtz) 	244433 	9831 	25715
% Ne-21 (qtz) 	989864 	22850 	65587 	

% Using https://github.com/angelrodes/average_cosmogenic_production_rate_calculator
P10=[11.46492   0.02464   0.00821   0.06961];
P26=[80.9226   0.177   0.059   0.8233];
P21=[47.393   1.335   0.445   0];

% attenuationg lengths
L=[160 850 5000 500]; % g/cm2. See Rodés, Á. (2021) https://doi.org/10.3390/geosciences11090362

% decay contants
l10=-log(0.5)./1.387e6; % 1/years
l26=-log(0.5)./0.705e6; % 1/years
earth_age=4.543e9; % years
l21=-log(0.5)./(1000*earth_age); % see https://angelrodes.wordpress.com/2022/01/25/are-he-3-and-ne-21-stable-isotopes/

%% Erosion models
disp(['Computing erosion models...'])
%lowerings=1/100*[0,10:10:100,120:20:250,300:50:500,6e2:1e2:10e2,20e2:10e2:100e2,1000e2,10000e2]; % m
lowerings=1/100*[0,10:10:100,120:20:250,300:50:500,6e2:1e2:50e2,55e2:5e2:100e2,1000e2,10000e2]; % m
% Model A: inherted relief
% Model B: exhumed tor
% Model C: lateral erosion
zq_modelA=zeros([size(zq),numel(lowerings)]);
zq_modelB=zeros([size(zq),numel(lowerings)]);
zq_modelC=zeros([size(zq),numel(lowerings)]);
original_surface_modelC=max(zq,z_no_tors+z_core-griddata(xq,yq,z_no_tors,x_core,y_core));
for step=1:numel(lowerings)
    low=lowerings(step);
    zq_modelA(:,:,step)=zq+low;
    zq_modelB(:,:,step)=max(z_no_tors+low,zq);
    captured=(xq-x_core).^2+(yq-y_core).^2>(low)^2;
    zq_modelC(:,:,step)=original_surface_modelC.*(captured==0)+zq.*(captured==1);
end
disp('Done.')


%% Map shielding factors

for model_letter=[{'A'},{'B'},{'C'}]
    disp(['Computing shielding of model ' model_letter{:} '...'])
    eval(['DEM_model=zq_model' model_letter{:} ';'])
    % Shielding parameters
    rho=2.6; % Density of granite
    Shielding_table=zeros(numel(lowerings),numel(z_samples),numel(L));
    for lowering_i=1:numel(lowerings)
        DEM=DEM_model(:,:,lowering_i);
        for L_i=1:numel(L)
            for sample_i=1:numel(z_samples)
                sample_elv=z_samples(sample_i);
                z_dem=DEM-sample_elv;
                angles=atan((DEM-sample_elv)./dist);
                depths=(z_dem.^2+dist.^2).^0.5.*(angles>0); % compute no shielding for negative angles
                wehight1=sin(abs(angles)).^2.3; % weight for radiation angle
                angles_reference=linspace(atan(0/1),atan(1/0),100);
                [counts,angles_ref] = hist(abs(angles(:)),angles_reference);
                wehight2=1./interp1(angles_ref,counts,abs(angles),'nearest'); % weight for angle distribution
                wehights=wehight1.*wehight2;
                attenuation_factors=(1-exp(-depths*100*rho/L(L_i)));
                shielding=1-sum(attenuation_factors.*wehights)/sum(wehights);  
                Shielding_table(lowering_i,sample_i,L_i)=shielding;
            end
        end
    end
    eval(['Shielding_table_model' model_letter{:} '=Shielding_table;'])
    disp('Done.')
end

%% Create and fit models A, B, and C

% Discretise time
earth_age=4.543e9; % years
start_time=[0,logspace(2,log10(earth_age),1000)]; % end of time section
end_time=[logspace(2,log10(earth_age),1000),earth_age]; % start of time section
time=(start_time+end_time)./2; % center of time section
delta_time=end_time-start_time; % length of time section

% Discretise exhumation rates
exrates=[0.001,logspace(-1,4,1000),100000]/1e6; % m/year

% Run accumulation models
for model_letter=[{'A'},{'B'},{'C'}]
    disp(['Concentrations model ' model_letter{:} '...'])
    eval(['Shielding_table=Shielding_table_model' model_letter{:} ';'])
    C10_model=zeros(numel(exrates),numel(z_samples));
    C26_model=zeros(numel(exrates),numel(z_samples));
    C21_model=zeros(numel(exrates),numel(z_samples));
    for exrates_i=1:numel(exrates)
        % Shielding_table=zeros(numel(lowerings),numel(z_samples),numel(L));
        lowerings_model=time*exrates(exrates_i);
        indexes_lowerings=interp1(lowerings,lowerings_model,'linear','extrap');
        indexes_lowerings=max(1,min(numel(lowerings),indexes_lowerings)); % use nearest neighbour for extrapolation to avoid artifacts
        shieldings_model=min(1,max(0,interp1(Shielding_table,indexes_lowerings,'linear'))); % dimensions: time,sample_i,L-i
        
        % Accumulation
        C10_model(exrates_i,:)=sum(sum(permute(P10,[1 3 2])./l10.*shieldings_model.*(1-exp(-l10.*permute(delta_time,[2 1]))).*exp(-l10*permute(time,[2 1])),1),3);
        C26_model(exrates_i,:)=sum(sum(permute(P26,[1 3 2])./l26.*shieldings_model.*(1-exp(-l26.*permute(delta_time,[2 1]))).*exp(-l26*permute(time,[2 1])),1),3);
        C21_model(exrates_i,:)=sum(sum(permute(P21,[1 3 2])./l21.*shieldings_model.*(1-exp(-l21.*permute(delta_time,[2 1]))).*exp(-l21*permute(time,[2 1])),1),3);
    end
    eval(['C10_model' model_letter{:} '=C10_model;'])
    eval(['C26_model' model_letter{:} '=C26_model;'])
    eval(['C21_model' model_letter{:} '=C21_model;'])
    disp('Done.')
end

% Fit exhumation rates
plot_count=0;
for model_letter=[{'A'},{'B'},{'C'}]
    disp(['Fitting model ' model_letter{:} '...'])
    eval(['C10_model=C10_model' model_letter{:} ';'])
    eval(['C26_model=C26_model' model_letter{:} ';'])
    eval(['C21_model=C21_model' model_letter{:} ';'])
    chisq=sum(((C10-C10_model)./dC10).^2+((C26-C26_model)./dC26).^2+((C21-C21_model)./dC21).^2,2); % chi squared
    best_model=find(chisq==min(chisq),1,'first');
    fitting_exhumation_rate=exrates(best_model);
    fitting_chisq=min(chisq);
    fitting_C10=C10_model(best_model,:);
    fitting_C26=C26_model(best_model,:);
    fitting_C21=C21_model(best_model,:);
    eval(['Model' model_letter{:} '_exrate= fitting_exhumation_rate ;'])
    eval(['Model' model_letter{:} '_chisq= fitting_chisq ;'])
    eval(['Model' model_letter{:} '_C10= fitting_C10 ;'])
    eval(['Model' model_letter{:} '_C26= fitting_C26 ;'])
    eval(['Model' model_letter{:} '_C21= fitting_C21 ;'])
    
    disp(['Model' model_letter{:} ': exrate=' num2str(fitting_exhumation_rate*1e6) ' ; chisq=' num2str(fitting_chisq)])
    
    % plot models
    if plot_count==0
        h=figure('units','normalized','outerposition',[0 0 1 1]);
        yplot=z_samples-max(z_samples);
    end
    plot_count=plot_count+1;
    subplot(3,3,plot_count)
    xlabel('[^{10}Be]')
    set(gca, 'XAxisLocation', 'top')
    set(gca,'YTickLabel',[]);
    hold on
    ylim([min(yplot) max(yplot)])
    plot(fitting_C10,yplot,'-b','LineWidth',3)
    plot(C10-dC10,yplot,'<r','MarkerSize',3)
    plot(C10+dC10,yplot,'>r','MarkerSize',3)
    plot(C10,yplot,'pr','MarkerFaceColor','r','MarkerSize',10)
    box on
    plot_count=plot_count+1;
    subplot(3,3,plot_count)
    xlabel('[^{26}Al]')
    set(gca, 'XAxisLocation', 'top')
    set(gca,'YTickLabel',[]);
    hold on
    ylim([min(yplot) max(yplot)])
    plot(fitting_C26,yplot,'-b','LineWidth',3)
    plot(C26-dC26,yplot,'<r','MarkerSize',3)
    plot(C26+dC26,yplot,'>r','MarkerSize',3)
    plot(C26,yplot,'pr','MarkerFaceColor','r','MarkerSize',10)
    box on
    ylim([min(yplot) max(yplot)])
    title(['Model ' model_letter{:} ' ;  Ex. rate = ' num2str(fitting_exhumation_rate*1e6,3) ' m/Ma'])
    plot_count=plot_count+1;
    subplot(3,3,plot_count)
    % xlabel('[^{21}Ne] (not fitted)')
    xlabel('[^{21}Ne]')
    set(gca, 'XAxisLocation', 'top')
    set(gca,'YTickLabel',[]);
    hold on
    ylim([min(yplot) max(yplot)])
    plot(fitting_C21,yplot,'-b','LineWidth',3)
    plot(C21-dC21,yplot,'<r','MarkerSize',3)
    plot(C21+dC21,yplot,'>r','MarkerSize',3)
    plot(C21,yplot,'pr','MarkerFaceColor','r','MarkerSize',10)
    box on
    
end
saveas(h,['./figs/' 'results_erosion_models'],'jpg')

%% Minimum age model
disp(['Fitting no-erosion model...'])
Shielding_table=Shielding_table_modelA(1,:,:);
P10_current=sum(permute(P10,[1 3 2]).*Shielding_table,3);
P26_current=sum(permute(P26,[1 3 2]).*Shielding_table,3);
P21_current=sum(permute(P21,[1 3 2]).*Shielding_table,3);
C10_model=P10_current/l10.*(1-exp(-l10*time'));
C26_model=P26_current/l26.*(1-exp(-l26*time'));
C21_model=P21_current/l21.*(1-exp(-l21*time'));
chisq=sum(((C10-C10_model)./dC10).^2+((C26-C26_model)./dC26).^2+((C21-C21_model)./dC21).^2,2); % chi squared
best_model=find(chisq==min(chisq),1,'first');
fitting_age=time(best_model);
fitting_chisq=min(chisq);
fitting_C10=C10_model(best_model,:);
fitting_C26=C26_model(best_model,:);
fitting_C21=C21_model(best_model,:);
disp(['Minimum age model ;  T = ' num2str(fitting_age/1e3,3) ' ka' ' ; chisq=' num2str(fitting_chisq)])
% plot model
plot_count=0;
h=figure('units','normalized','outerposition',[0 0 1 1]);
yplot=z_samples-max(z_samples);
plot_count=plot_count+1;
subplot(1,3,plot_count)
xlabel('[^{10}Be]')
set(gca, 'XAxisLocation', 'top')
set(gca,'YTickLabel',[]);
hold on
ylim([min(yplot) max(yplot)])
plot(fitting_C10,yplot,'-b','LineWidth',3)
plot(C10-dC10,yplot,'<r','MarkerSize',3)
plot(C10+dC10,yplot,'>r','MarkerSize',3)
plot(C10,yplot,'pr','MarkerFaceColor','r','MarkerSize',10)
box on
plot_count=plot_count+1;
subplot(1,3,plot_count)
xlabel('[^{26}Al]')
set(gca, 'XAxisLocation', 'top')
set(gca,'YTickLabel',[]);
hold on
ylim([min(yplot) max(yplot)])
plot(fitting_C26,yplot,'-b','LineWidth',3)
plot(C26-dC26,yplot,'<r','MarkerSize',3)
plot(C26+dC26,yplot,'>r','MarkerSize',3)
plot(C26,yplot,'pr','MarkerFaceColor','r','MarkerSize',10)
box on
ylim([min(yplot) max(yplot)])
title(['Minimum age model ;  T = ' num2str(fitting_age/1e3,3) ' ka'])
plot_count=plot_count+1;
subplot(1,3,plot_count)
% xlabel('[^{21}Ne] (not fitted)')
xlabel('[^{21}Ne]')
set(gca, 'XAxisLocation', 'top')
set(gca,'YTickLabel',[]);
hold on
ylim([min(yplot) max(yplot)])
plot(fitting_C21,yplot,'-b','LineWidth',3)
plot(C21-dC21,yplot,'<r','MarkerSize',3)
plot(C21+dC21,yplot,'>r','MarkerSize',3)
plot(C21,yplot,'pr','MarkerFaceColor','r','MarkerSize',10)
box on
saveas(h,['./figs/' 'results_no_erosion'],'jpg')

%% Glacial model
% This model considers the tor to be shielded from cosmic radiation when
% covered by a glacier. It also considers uplift, assuming that the tor was
% ice-free for the last galciations (d18O threshold> max(d18O)).

% Glacial data from Rodés, Á. (2021) https://doi.org/10.3390/geosciences11090362
make_climatecurves
load('climatecurves.mat')
T=climatecurves.age;
dT=diff([climatecurves.age,4543e6]);
D=climatecurves.d18O;

% build d18-deglaciation age reference
step=0;
for n=2:numel(T)
    if D(n)>max(D(1:n-1))
        step=step+1;
        deglaciation_ref.d18(step)=D(n);
        deglaciation_ref.age(step)=T(n);
    end
end

% get current production rates
Shielding_table=Shielding_table_modelA(1,:,:);
P10_current=sum(permute(P10,[1 3 2]).*Shielding_table,3);
P26_current=sum(permute(P26,[1 3 2]).*Shielding_table,3);
P21_current=sum(permute(P21,[1 3 2]).*Shielding_table,3);

disp(['Fitting glacial model...'])

deglaciation_d18_test=unique(climatecurves.d18O); % discrtise thresholds

% Make input matrices 
[sample_index,test_index,climate_index]=ndgrid(1:numel(z_samples),1:numel(deglaciation_d18_test),1:numel(climatecurves.age));
% Dimensions:
% 1 samples: P_current
% 2 test_model
% 3 climate curves
% Rearrange data with these dimensions:
P10_current=P10_current(sample_index); P26_current=P26_current(sample_index); P21_current=P21_current(sample_index);
T=T(climate_index); dT=dT(climate_index); D=D(climate_index);
deglaciation_d18_test=deglaciation_d18_test(test_index);
% calculate glaciated times
glaciated = (D>deglaciation_d18_test);

% Get concentrations
C10_model=cumsum(P10_current.*~glaciated/l10.*(1-exp(-l10.*dT)).*exp(-l10.*T),3);
C26_model=cumsum(P26_current.*~glaciated/l26.*(1-exp(-l26.*dT)).*exp(-l26.*T),3);
C21_model=cumsum(P21_current.*~glaciated/l21.*(1-exp(-l21.*dT)).*exp(-l21.*T),3);
% Dimensions:
% 1 samples: P_current
% 2 test_model
chisq=sum(((C10'-C10_model)./dC10').^2+((C26'-C26_model)./dC26').^2+((C21'-C21_model)./dC21').^2,1);
bestmodel=find(chisq==min(chisq(:)),1,'first');
[xi, yi, zi] = ind2sub(size(chisq), bestmodel);
fitting_d18=deglaciation_d18_test(xi, yi, zi);
fitting_age=T(xi, yi, zi);
fitting_chisq=min(chisq(:));
fitting_C10=C10_model(:,yi, zi);
fitting_C26=C26_model(:,yi, zi);
fitting_C21=C21_model(:,yi, zi);


disp(['Glacial model ;  D18 threshold = ' num2str(fitting_d18,3) ' ; First exposure = ' num2str(fitting_age/1e3,3) ' ka ; chisq=' num2str(fitting_chisq)])



%% Plot glacial model results
plot_count=0;
h=figure('units','normalized','outerposition',[0 0 1 1]);
yplot=z_samples-max(z_samples);
plot_count=plot_count+1;
subplot(2,3,plot_count)
xlabel('[^{10}Be]')
set(gca, 'XAxisLocation', 'top')
set(gca,'YTickLabel',[]);
hold on
ylim([min(yplot) max(yplot)])
plot(fitting_C10,yplot,'-b','LineWidth',3)
plot(C10-dC10,yplot,'<r','MarkerSize',3)
plot(C10+dC10,yplot,'>r','MarkerSize',3)
plot(C10,yplot,'pr','MarkerFaceColor','r','MarkerSize',10)
box on
plot_count=plot_count+1;
subplot(2,3,plot_count)
xlabel('[^{26}Al]')
set(gca, 'XAxisLocation', 'top')
set(gca,'YTickLabel',[]);
hold on
ylim([min(yplot) max(yplot)])
plot(fitting_C26,yplot,'-b','LineWidth',3)
plot(C26-dC26,yplot,'<r','MarkerSize',3)
plot(C26+dC26,yplot,'>r','MarkerSize',3)
plot(C26,yplot,'pr','MarkerFaceColor','r','MarkerSize',10)
box on
ylim([min(yplot) max(yplot)])
title(['Glacial model'])
plot_count=plot_count+1;
subplot(2,3,plot_count)
% xlabel('[^{21}Ne] (not fitted)')
xlabel('[^{21}Ne]')
set(gca, 'XAxisLocation', 'top')
set(gca,'YTickLabel',[]);
hold on
ylim([min(yplot) max(yplot)])
plot(fitting_C21,yplot,'-b','LineWidth',3)
plot(C21-dC21,yplot,'<r','MarkerSize',3)
plot(C21+dC21,yplot,'>r','MarkerSize',3)
plot(C21,yplot,'pr','MarkerFaceColor','r','MarkerSize',10)
box on
plot_count=plot_count+1;
subplot(2,3,[plot_count plot_count+2])
hold on
plot(climatecurves.age/1e3,climatecurves.d18O,'-','Color','k','LineWidth',2)
plot([0,fitting_age]/1e3,[fitting_d18,fitting_d18],'-','Color','b','LineWidth',1)
plot(fitting_age/1e3,fitting_d18,'p','Color','b','MarkerFaceColor','b','MarkerSize',10)
ylabel('\delta^{18}O')
xlabel('Age (ka)')
set(gca, 'Xdir', 'reverse')
maxage=1.2*fitting_age;
xlim([0 maxage/1e3])
ylim([min(min(fitting_d18,climatecurves.d18O(climatecurves.age<maxage)))-0.5 max(fitting_d18,max(climatecurves.d18O))+0.5])
set (gca, 'Clipping', 'on');
box on
grid on
title(['Best fit: \delta^{18}O_{0} = ' num2str(fitting_d18,'%-5.3f') ' ; t_{0} = ' num2str(fitting_age/1e3,'%-5.0f') ' ka'])

saveas(h,['./figs/' 'results_glacial_model'],'jpg')

%% Discarded models

return %--------------------------------------------STOP------------------

% %% Minimum age model at depth
% low_index=33;
% disp(['Fitting no-erosion model at depth ' num2str(lowerings(low_index))])
% Shielding_table=Shielding_table_modelB(low_index,:,:);
% P10_current=sum(permute(P10,[1 3 2]).*Shielding_table,3);
% P26_current=sum(permute(P26,[1 3 2]).*Shielding_table,3);
% P21_current=sum(permute(P21,[1 3 2]).*Shielding_table,3);
% C10_model=P10_current/l10.*(1-exp(-l10*time'));
% C26_model=P26_current/l26.*(1-exp(-l26*time'));
% C21_model=P21_current/l21.*(1-exp(-l21*time'));
% chisq=sum(((C10-C10_model)./dC10).^2+((C26-C26_model)./dC26).^2+((C21-C21_model)./dC21).^2,2); % chi squared
% best_model=find(chisq==min(chisq),1,'first');
% fitting_age=time(best_model);
% fitting_chisq=min(chisq);
% fitting_C10=C10_model(best_model,:);
% fitting_C26=C26_model(best_model,:);
% fitting_C21=C21_model(best_model,:);
% disp(['Minimum age model at depth ' num2str(lowerings(low_index)) ' ;  T = ' num2str(fitting_age/1e3,3) ' ka' ' ; chisq=' num2str(fitting_chisq)])
% % plot model
% plot_count=0;
% h=figure('units','normalized','outerposition',[0 0 1 1]);
% yplot=z_samples-max(z_samples);
% plot_count=plot_count+1;
% subplot(1,3,plot_count)
% xlabel('[^{10}Be]')
% set(gca, 'XAxisLocation', 'top')
% set(gca,'YTickLabel',[]);
% hold on
% ylim([min(yplot) max(yplot)])
% plot(fitting_C10,yplot,'-b','LineWidth',3)
% plot(C10-dC10,yplot,'<r','MarkerSize',3)
% plot(C10+dC10,yplot,'>r','MarkerSize',3)
% plot(C10,yplot,'pr','MarkerFaceColor','r','MarkerSize',10)
% box on
% plot_count=plot_count+1;
% subplot(1,3,plot_count)
% xlabel('[^{26}Al]')
% set(gca, 'XAxisLocation', 'top')
% set(gca,'YTickLabel',[]);
% hold on
% ylim([min(yplot) max(yplot)])
% plot(fitting_C26,yplot,'-b','LineWidth',3)
% plot(C26-dC26,yplot,'<r','MarkerSize',3)
% plot(C26+dC26,yplot,'>r','MarkerSize',3)
% plot(C26,yplot,'pr','MarkerFaceColor','r','MarkerSize',10)
% box on
% ylim([min(yplot) max(yplot)])
% title(['Minimum age model at depth ' num2str(lowerings(low_index)) ' m ;  T = ' num2str(fitting_age/1e3,3) ' ka'])
% plot_count=plot_count+1;
% subplot(1,3,plot_count)
% % xlabel('[^{21}Ne] (not fitted)')
% xlabel('[^{21}Ne]')
% set(gca, 'XAxisLocation', 'top')
% set(gca,'YTickLabel',[]);
% hold on
% ylim([min(yplot) max(yplot)])
% plot(fitting_C21,yplot,'-b','LineWidth',3)
% plot(C21-dC21,yplot,'<r','MarkerSize',3)
% plot(C21+dC21,yplot,'>r','MarkerSize',3)
% plot(C21,yplot,'pr','MarkerFaceColor','r','MarkerSize',10)
% box on


%% Stepped erosion model
% This model assumes that current tor shape is very recent and most of the
% radiation recieved by the samples ocurred qhen the surface was a few
% meters above the tor top. We can emulate this by assuming constant
% exposure for DEMs generated in model B.

% Shielding_table_modelB; % dimensions lowerings,samples,L
% zs_model=[0.5:0.5:30,32:2:50,55:5:100,110:10:1000];
% Shielding_table_thismodel=interp1(lowerings',Shielding_table_modelB,zs_model');
zs_model=lowerings;
Shielding_table_thismodel=Shielding_table_modelB;
C10_model=sum(permute(P10,[1 3 2]).*Shielding_table_thismodel,3)/l10;
C26_model=sum(permute(P26,[1 3 2]).*Shielding_table_thismodel,3)/l26;
C21_model=sum(permute(P21,[1 3 2]).*Shielding_table_thismodel,3)/l21;
% chisq=sum(((C10-C10_model)./dC10).^2+((C26-C26_model)./dC26).^2+((C21-C21_model)./dC21).^2,2); % chi squared
chisq=sum(((C10-C10_model)./dC10).^2+((C26-C26_model)./dC26).^2,2); % chi squared 10 and 26 only
% chisq=sum(((C21-C21_model)./dC21).^2,2); % chi squared 21 only
best_model=find(chisq==min(chisq),1,'first');
fitting_depth=zs_model(best_model);
fitting_chisq=min(chisq);
fitting_C10=C10_model(best_model,:);
fitting_C26=C26_model(best_model,:);
fitting_C21=C21_model(best_model,:);
figure; plot(zs_model,chisq,'.-'); set(gca, 'YScale', 'log')

% plot model
plot_count=0;
h=figure('units','normalized','outerposition',[0 0 1 1]);
yplot=z_samples-max(z_samples);
plot_count=plot_count+1;
subplot(1,3,plot_count)
xlabel('[^{10}Be]')
set(gca, 'XAxisLocation', 'top')
set(gca,'YTickLabel',[]);
hold on
ylim([min(yplot) max(yplot)])
plot(fitting_C10,yplot,'-b','LineWidth',3)
plot(C10-dC10,yplot,'<r','MarkerSize',3)
plot(C10+dC10,yplot,'>r','MarkerSize',3)
plot(C10,yplot,'pr','MarkerFaceColor','r','MarkerSize',10)
box on
plot_count=plot_count+1;
subplot(1,3,plot_count)
xlabel('[^{26}Al]')
set(gca, 'XAxisLocation', 'top')
set(gca,'YTickLabel',[]);
hold on
ylim([min(yplot) max(yplot)])
plot(fitting_C26,yplot,'-b','LineWidth',3)
plot(C26-dC26,yplot,'<r','MarkerSize',3)
plot(C26+dC26,yplot,'>r','MarkerSize',3)
plot(C26,yplot,'pr','MarkerFaceColor','r','MarkerSize',10)
box on
ylim([min(yplot) max(yplot)])
title(['Depth model ;  Z = ' num2str(fitting_depth,3) ' m'])
plot_count=plot_count+1;
subplot(1,3,plot_count)
% xlabel('[^{21}Ne] (not fitted)')
xlabel('[^{21}Ne]')
set(gca, 'XAxisLocation', 'top')
set(gca,'YTickLabel',[]);
hold on
ylim([min(yplot) max(yplot)])
plot(fitting_C21,yplot,'-b','LineWidth',3)
plot(C21-dC21,yplot,'<r','MarkerSize',3)
plot(C21+dC21,yplot,'>r','MarkerSize',3)
plot(C21,yplot,'pr','MarkerFaceColor','r','MarkerSize',10)
box on

% Surface C21 and C10 for 2.2 Ma of exposure followed by 4 Ma of burial
% 11e7/sum(P21)
% sum(P10)/l10*(1-exp(-l10*11e7/sum(P21)))*exp(-l10*4e6)

return %--------------------------------------------STOP------------------

%% Other figures %%

%% DEM
h=figure('Name','Contours','units','normalized','outerposition',[0 0 1 1])
hold on
% levels=floor(min(z)):1:ceil(max(z)); % every 0.5 m
levels=975:0.5:1031;
caxis([min(levels) max(levels)])
contourf(xq, yq, zq,levels,'ShowText','off');
colormap(bone)
contour(xq, yq,areas_with_tors,[-1,0.5],'EdgeColor','r','LineWidth',3)
plot(x_core,y_core,'pr','MarkerFaceColor','r','MarkerSize',25)
xlabel('Lon.');
ylabel('Lat.');
title('DEM, core location, and tor areas');
grid on; box on;
axis equal;
% zoom to tor
distance=(2*100^2)^0.5/2;
xlim([x_core-distance x_core+distance])
ylim([y_core-distance y_core+distance])
saveas(h,['./figs/' 'DEM'],'jpg')

% figure('Name','Tors only')
% hold on
% levels=0:0.5:15; % every 0.5 m
% contour(xq, yq, zq-z_no_tors,levels);
% plot(x_core,y_core,'*r')
% xlabel('Lon.');
% ylabel('Lat.');
% title('Contour Plot of tors');
% grid on; box on;
% axis equal;

%% 3-D DEM
figure('Name','DEM')
hold on
surf(xq, yq, zq,'FaceColor','none');
plot3(x_samples,y_samples,z_samples,'*r')
xlabel('Lon.');
ylabel('Lat.');
zlabel('Elv.')
grid on; box on;
axis equal;
view(-37.5+180,30)
% distance=(2*100^2)^0.5/2;
distance=50;
xlim([x_core-distance x_core+distance])
ylim([y_core-distance y_core+distance])
zlim([975 1050])

%% Shielding method
h=figure('Name','Shielding method');
hold on
distance=1;
zsample=z_samples(3);
dirs=unique(az(:))';
dir=dirs(2);
az_i=dir;
xplot=dist(az==az_i)*(2*(az_i<mean(az(:)))-1);
zplot=zq(az==az_i);
az_i=dir+pi;
xplot=[flipud(xplot),;dist(az==az_i)*(2*(az_i<mean(az(:)))-1)];
zplot=[flipud(zplot);zq(az==az_i)];
yplotlimits=[zsample-distance/10 zplot(round(numel(zplot)/2))+distance/3];

% plot cosmic rays
for n=1:numel(xplot)
    zray=[max(yplotlimits) zsample];
    xray=interp1([zplot(n),zsample],[xplot(n),0],zray,'linear','extrap');
    plot(xray,zray,'-','Color',[0.9 0.9 0.9]) % one ray
end

% plot one ray
selectray=round(numel(zplot)/2)+45;
zray=[zplot(selectray) zsample];
xray=[xplot(selectray),0];
plot(xray,zray,'-b','LineWidth',2) % selected ray
text(mean(xray),mean(zray),'d','Color','b','VerticalAlignment','bottom') % distance

plot([-distance/20 0],zsample*[1 1],'-b') % angle
text(-distance/20,zsample,'\phi','Color','b','VerticalAlignment','bottom','HorizontalAlignment','right') % angle


plot(xplot,zplot,'o-k','LineWidth',1,'MarkerFaceColor','k','MarkerSize',3) % surface
text(xplot(round(numel(zplot)/2)),zplot(round(numel(zplot)/2)),' Surface','VerticalAlignment','bottom')
text(0,zsample,'   Sample','Color','r') % sample
plot(0,zsample,'pr','MarkerFaceColor','r','MarkerSize',10) % sample

xlim([-distance*0.7 distance*0.2])
ylim(yplotlimits)

title('Shielding calculations')
xlabel('Distance to core (m)')
ylabel('Elv.')
set (gca, 'Clipping', 'on');
box on

saveas(h,['./figs/' 'Shielding_method'],'jpg')


%% Shielding profile
h=figure('Name','Shielding profile','units','normalized','outerposition',[0 0 1 1])
subplot(1,2,1)
hold on
distance=10;
% xplot=x_core-distance:0.005:x_core+distance;
% yplot=y_core*ones(size(xplot));
% zplot=griddata(xq,yq,zq,xplot,yplot);
% plot(xplot,zplot,'-k','LineWidth',2)
for az_i=unique(az(:))'
    xplot=dist(az==az_i)*(2*(az_i<mean(az(:)))-1);
    zplot=zq(az==az_i);
    plot(xplot,zplot,'-k','LineWidth',1)
end
plot(x_samples.*0,z_samples,'pr','MarkerFaceColor','r')
for n=1:numel(z_samples)
    text(x_samples(n).*0,z_samples(n),[' ' sample_names{n}],'Color','r')
end
title('Radial cross sections')
xlabel('Distance to core (m)')
ylabel('Elv.')
xlim([-10 10])
ylim([min(z_samples)-1 max(z_samples)+1.5])
box on

subplot(1,2,2)
hold on
yplot=z_samples;
Shielding_table=Shielding_table_modelA(1,:,:);
% P10_current=sum(permute(P10,[1 3 2]).*Shielding_table,3);
% P26_current=sum(permute(P26,[1 3 2]).*Shielding_table,3);
% P21_current=sum(permute(P21,[1 3 2]).*Shielding_table,3);
% plot(P10_current/sum(P10),yplot,':*')
% plot(P26_current/sum(P26),yplot,':*')
% plot(P21_current/sum(P21),yplot,':*')
plot(Shielding_table(:,:,1),yplot,'--pr','MarkerFaceColor','r')
for n=2:size(Shielding_table,3)
    plot(Shielding_table(:,:,n),yplot,'--pb','MarkerFaceColor','b')
end
legend('Spallation','Muons','AutoUpdate','off','Location','southeast')
title('Shielding profiles (P_{i}/P_{i,0})')
% xlabel('P/P_{0}')
ylabel('Elv.')
ylim([min(z_samples)-1 max(z_samples)+1.5])
set(gca, 'XAxisLocation', 'top')
set(gca, 'XScale', 'log')
box on
grid on
saveas(h,['./figs/' 'Shielding_profile_calculations'],'jpg')


%% Sources of z values
h=figure('Name','Raw elevation data','units','normalized','outerposition',[0 0 1 1])
subplot(2,2,[1 3]); hold on
distance=100;
x_values=x_core-distance:5:x_core+distance;
y_values=y_core-distance:5:y_core+distance;
[x1,y1] = meshgrid(x_values,y_values);
z1 = griddata(x, y, z, x1, y1, 'natural'); % Interpolate the data onto the grid
levels=floor(min(z1(:)/10))*10:2:ceil(max(z1(:)/10))*10;
contourf(x1, y1, z1,levels,'ShowText','on')
% plot(x_pic,y_pic,'-g','LineWidth',2)
% plot(x_field,y_field,'-b','LineWidth',2)
plot(x_core,y_core,'pr','MarkerFaceColor','r','MarkerSize',10)
axis equal; box on
title('5-dtm (digimap.edina.ac.uk)')    
xlabel('Long.')
ylabel('Lat.')
subplot(2,2,2); hold on
[~,order]=sort(x_field);
plot(x_field(order),z_field(order),'.-b','LineWidth',2)
plot(x_samples,z_samples,'pr','MarkerFaceColor','r')
axis equal; box on
title('Field measurements')
xlabel('Long.')
ylabel('Elv.')
subplot(2,2,4); hold on
plot(y_pic,z_pic,'.-g','LineWidth',2)
plot(y_samples,z_samples,'pr','MarkerFaceColor','r')
axis equal; box on
title('Picture measurements')
xlabel('Lat.')
ylabel('Elv.')
saveas(h,['./figs/' 'Elevation_data_sources'],'jpg')

%% Plot erosion models 
% Combine generated frames with: 
% convert -delay 10 -loop 0 *.jpg "$(pwd | awk -F'/' '{print $NF}').gif"

% plot Model A
frame=0;
for step=numel(lowerings):-1:1
    frame=frame+1;

    zplot=zq_modelA(:,:,step);
    
    h=figure;
    hold on
    surf(xq, yq, zq,'FaceColor','none','EdgeColor',[0.8 0.8 0.8]);
    surf(xq, yq, zplot,'FaceColor','none');
    plot3(x_samples,y_samples,z_samples,'*r')
    xlabel('Lon.');
    ylabel('Lat.');
    zlabel('Elv.')
    title(['Lowering=' num2str(lowerings(step),'%-5.2f') 'm'])
    grid on; box on;
    axis equal;
    view(-37.5+180,30)
    % distance=(2*100^2)^0.5/2;
    distance=50;
    xlim([x_core-distance x_core+distance])
    ylim([y_core-distance y_core+distance])
    zlim([975 1050])
    drawnow
    % exportgraphics(gcf,'ModelB.gif','Append',true);
    saveas(h,['./figs/ModelA/' num2str(1000+frame)],'jpg')
    close(h)
end

% plot Model B
frame=0;
for step=numel(lowerings):-1:1
    frame=frame+1;

    zplot=zq_modelB(:,:,step);
    
    h=figure;
    hold on
    surf(xq, yq, zq,'FaceColor','none','EdgeColor',[0.8 0.8 0.8]);
    surf(xq, yq, zplot,'FaceColor','none');
    plot3(x_samples,y_samples,z_samples,'*r')
    xlabel('Lon.');
    ylabel('Lat.');
    zlabel('Elv.')
    title(['Lowering=' num2str(lowerings(step),'%-5.2f') 'm'])
    grid on; box on;
    axis equal;
    view(-37.5+180,30)
    % distance=(2*100^2)^0.5/2;
    distance=50;
    xlim([x_core-distance x_core+distance])
    ylim([y_core-distance y_core+distance])
    zlim([975 1050])
    
    % exportgraphics(gcf,'ModelB.gif','Append',true);
    saveas(h,['./figs/ModelB/' num2str(1000+frame)],'jpg')
    close(h)
end

% plot Model C
frame=0;
for step=numel(lowerings):-1:1
    frame=frame+1;

    zplot=zq_modelC(:,:,step);
    
    h=figure;
    hold on
    surf(xq, yq, zq,'FaceColor','none','EdgeColor',[0.8 0.8 0.8]);
    surf(xq, yq, zplot,'FaceColor','none');
    plot3(x_samples,y_samples,z_samples,'*r')
    xlabel('Lon.');
    ylabel('Lat.');
    zlabel('Elv.')
    title(['Horz.lowering=' num2str(lowerings(step),'%-5.2f') 'm'])
    grid on; box on;
    axis equal;
    view(-37.5+180,30)
    % distance=(2*100^2)^0.5/2;
    distance=50;
    xlim([x_core-distance x_core+distance])
    ylim([y_core-distance y_core+distance])
    zlim([975 1050])
    
    % exportgraphics(gcf,'ModelB.gif','Append',true);
    saveas(h,['./figs/ModelC/' num2str(1000+frame)],'jpg')
    close(h)
end

%% Shielding models plot
h=figure('Name','Shielding models','units','normalized','outerposition',[0 0 1 1]);
row=0;
plot_number=0;
for model_letter=[{'A'},{'B'},{'C'}]
    row=row+1;
    disp(['Plotting shielding of model ' model_letter{:} '...'])
    eval(['Shielding_table=Shielding_table_model' model_letter{:} ';'])
    for column=1:numel(L)
        plot_number=plot_number+1;
        subplot(3,numel(L),plot_number); hold on
        if plot_number==1 
            ylabel('Surf. lowerings (m)')
        elseif plot_number==5
            ylabel('Base lowerings (m)')
        elseif plot_number==9
            ylabel('Horiz. lowerings (m)')
        end
        if plot_number>8
            xlabel('Sample elevation (m)')
        end
        thisplot=Shielding_table(:,:,column);
        title(['Model ' model_letter{:} ' \Lambda=' num2str(L(column))])
        [sample_elevations,lowerings_model] = meshgrid(z_samples,lowerings);
        levels=-0.01:0.005:2;
        contourf(sample_elevations,lowerings_model,thisplot,levels, 'LineColor', 'none');
%         colormap(flipud(hot))colormap(hot)
        levels=[0.1,0.2,0.5,0.9];
        contour(sample_elevations,lowerings_model,thisplot,levels, 'LineColor', 'k','ShowText','on');
        set(gca, 'YScale', 'log')
        ylim([lowerings(2) lowerings(end-1)])
        box on
        if max(thisplot(:))>1
            disp(['Model ' model_letter{:} ' \Lambda=' num2str(L(column)) 'max>1 !!!'])
        end
    end
end
saveas(h,['./figs/' 'Shielding_models'],'jpg')