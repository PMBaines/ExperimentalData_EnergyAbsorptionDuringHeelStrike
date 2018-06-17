%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Process raw data to get energy loss estimates and other results          %
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

% clear all
% clc

%% load: choose to run getRawDataStruture or load saved mat file 

getRawDataStruture
% load('RawData.mat')

%% Voltage-force gain
% for Force Plate in z direction per subject

for i = Data.subjects
    counter = 0;
    gain_total = zeros(length(Data.subject(i).n_25),1);
    for j = Data.subject(i).n_25
        counter = counter + 1;
        
        %%% uncomment to plot all the voltage data in 1 plot to compare and check results
%         figure(1)
%         plot(Data.subject(i).n(j).Vz)
%         hold on
        
        gain_total(counter) = Data.w_25kg/mean(Data.subject(i).n(j).Vz);
        
    end
%     gain_total
%     gaindif(i)=gain_total(2)-gain_total(1)
    Data.subject(i).gain_z = mean(gain_total); % [N/V]
end
Data.g = 9.81; % [m/s^2] gravitaional acceleration

% figure(1)
% title('check calibration measurement')
% xlabel('sample')
% ylabel('voltage [V]')

%%% Calibration data checks out:
% gain according to manual specs = 131.5789
% highest gain found = 133.4219 [N/V]
% lowest gain found = 132.5754
% difference:  0.8465
% highest difference within subject measurement: 0.1874

%% Force and acceleration data
%%% Create force data from voltage FP and velocity and acceleration from position data

for i = Data.subjects
    temp_markers = median(Data.subject(i).n(Data.subject(i).n_bodyweight(1)).z(:,2)) - median(Data.subject(i).n(Data.subject(i).n_bodyweight(1)).z(:,1));
    for j = Data.subject(i).n_total
        Data.subject(i).n(j).Fx = Data.subject(i).n(j).Vx * Data.gain_xy;
        Data.subject(i).n(j).Fy = -Data.subject(i).n(j).Vy * Data.gain_xy;
        Data.subject(i).n(j).Fz = Data.subject(i).n(j).Vz * Data.subject(i).gain_z;
        Data.subject(i).n(j).Fz1 = Data.subject(i).n(j).kp(:,5) * Data.subject(i).gain_z;
        Data.subject(i).n(j).Fz2 = Data.subject(i).n(j).kp(:,6) * Data.subject(i).gain_z;
        Data.subject(i).n(j).Fz3 = Data.subject(i).n(j).kp(:,7) * Data.subject(i).gain_z;        
        Data.subject(i).n(j).Fz4 = Data.subject(i).n(j).kp(:,8) * Data.subject(i).gain_z;
        
        Data.subject(i).n(j).t_kp = (1/Data.f_kp:1/Data.f_kp:length(Data.subject(i).n(j).Vz)/Data.f_kp)';
        Data.subject(i).n(j).t_opto = (1/Data.f_opto:1/Data.f_opto:length(Data.subject(i).n(j).z)/Data.f_opto)';
        
        if temp_markers < 0
            Data.subject(i).n(j).x = [Data.subject(i).n(j).x(:,2) Data.subject(i).n(j).x(:,1)];
            Data.subject(i).n(j).y = [Data.subject(i).n(j).y(:,2) Data.subject(i).n(j).y(:,1)];
            Data.subject(i).n(j).z = [Data.subject(i).n(j).z(:,2) Data.subject(i).n(j).z(:,1)];
        end
        
        
        Data.subject(i).n(j).zd(:,1) = gradient(Data.subject(i).n(j).z(:,1),1/Data.f_opto);
        Data.subject(i).n(j).yd(:,1) = gradient(Data.subject(i).n(j).y(:,1),1/Data.f_opto);
        Data.subject(i).n(j).xd(:,1) = gradient(Data.subject(i).n(j).x(:,1),1/Data.f_opto);
        
        Data.subject(i).n(j).zd(:,2) = gradient(Data.subject(i).n(j).z(:,2),1/Data.f_opto);
        Data.subject(i).n(j).yd(:,2) = gradient(Data.subject(i).n(j).y(:,2),1/Data.f_opto);
        Data.subject(i).n(j).xd(:,2) = gradient(Data.subject(i).n(j).x(:,2),1/Data.f_opto);
        
        Data.subject(i).n(j).zdd(:,1) = gradient(Data.subject(i).n(j).zd(:,1),1/Data.f_opto);
        Data.subject(i).n(j).ydd(:,1) = gradient(Data.subject(i).n(j).yd(:,1),1/Data.f_opto);
        Data.subject(i).n(j).xdd(:,1) = gradient(Data.subject(i).n(j).xd(:,1),1/Data.f_opto);
        
        Data.subject(i).n(j).zdd(:,2) = gradient(Data.subject(i).n(j).zd(:,2),1/Data.f_opto);
        Data.subject(i).n(j).ydd(:,2) = gradient(Data.subject(i).n(j).yd(:,2),1/Data.f_opto);
        Data.subject(i).n(j).xdd(:,2) = gradient(Data.subject(i).n(j).xd(:,2),1/Data.f_opto);
        
%         Data.subject(i).n(j).v_yz = sqrt(Data.subject(i).n(j).zd.^2 + Data.subject(i).n(j).yd.^2);
%         Data.subject(i).n(j).v = sqrt(Data.subject(i).n(j).zd.^2 + Data.subject(i).n(j).yd.^2 + Data.subject(i).n(j).xd.^2);
        
    end
end

%% Body weight

for i = Data.subjects
    counter = 0;
    weight_total = zeros(length(Data.subject(i).n_bodyweight),1);
    for j = Data.subject(i).n_bodyweight
        counter = counter+1;
        l = length(Data.subject(i).n(j).Fz);
        weight_total(counter) = mean(Data.subject(i).n(j).Fz(round(l/3):round(2*l/3)));
    end
    % weight_total
    weight = mean(weight_total); % body weight [N]
    std_weight = std(weight_total); % standard deviation of body weight calculation
    Data.subject(i).bodyweight = weight;
end

%% Drift compensation
% Adjust offset of force data to compensate for drift
% Needed because amplifier is not reseted before each walking test

% n_walking = 9;
clc
close all
countx = 0
for i = Data.subjects
    for j = [Data.subject(i).n_normal Data.subject(i).n_slow]
        index_offset = 50;
        if std(Data.subject(i).n(j).Fx(250:250+index_offset)) > 1
            % standard deviation should be less than 0.4, otherwise it could
            % indicate force data used for offset calculation is not only noise
            disp('Warning: index used for calculating offset x due to drift is not valid for:')
            disp(strcat('subject',num2str(i),';test',num2str(j),';std',num2str(std(Data.subject(i).n(j).Fx(1:index_offset)))))
            countx = countx+1;
        end
        
        if std(Data.subject(i).n(j).Fy(1:index_offset)) > 1
            % standard deviation should be less than 0.4, otherwise it could
            % indicate force data used for offset calculation is not only noise
            disp('Warning: index used for calculating offset y due to drift is not valid for:')
            disp(strcat('subject',num2str(i),';test',num2str(j),';std',num2str(std(Data.subject(i).n(j).Fy(1:index_offset)))))
        end
        
        if std(Data.subject(i).n(j).Fz(1:index_offset)) > 1
            % standard deviation should be less than 0.4, otherwise it could
            % indicate force data used for offset calculation is not only noise
            disp('Warning: index used for calculating offset z due to drift is not valid for:')
            disp(strcat('subject',num2str(i),';test',num2str(j),';std',num2str(std(Data.subject(i).n(j).Fz(1:index_offset)))))
        end
        
        offsetx = median(Data.subject(i).n(j).Fx(250:250+index_offset));
        Data.subject(i).n(j).Fx = Data.subject(i).n(j).Fx-offsetx;
        
        offsety = median(Data.subject(i).n(j).Fy(1:index_offset));
        Data.subject(i).n(j).Fy = Data.subject(i).n(j).Fy-offsety;
        
        offsetz = median(Data.subject(i).n(j).Fz(1:index_offset));
        Data.subject(i).n(j).Fz = Data.subject(i).n(j).Fz-offsetz;
        Data.subject(i).n(j).Fz1 = Data.subject(i).n(j).Fz1-offsetz/4;
        Data.subject(i).n(j).Fz2 = Data.subject(i).n(j).Fz2-offsetz/4;
        Data.subject(i).n(j).Fz3 = Data.subject(i).n(j).Fz3-offsetz/4;
        Data.subject(i).n(j).Fz4 = Data.subject(i).n(j).Fz4-offsetz/4;
        Data.subject(i).n(j).Fz_separate = [Data.subject(i).n(j).Fz1,Data.subject(i).n(j).Fz2,Data.subject(i).n(j).Fz3,Data.subject(i).n(j).Fz4];
        
    end
end
countx
%%% Use for tests, with displayed warning, visual chosen offset indices

close all

i = 3
j = 11
offset_std = std(Data.subject(i).n(j).Fz(15:65))
offsetz = mean(Data.subject(i).n(j).Fz(15:65));
Data.subject(i).n(j).Fz = Data.subject(i).n(j).Fz-offsetz;
Data.subject(i).n(j).Fz1 = Data.subject(i).n(j).Fz1-offsetz/4;
Data.subject(i).n(j).Fz2 = Data.subject(i).n(j).Fz2-offsetz/4;
Data.subject(i).n(j).Fz3 = Data.subject(i).n(j).Fz3-offsetz/4;
Data.subject(i).n(j).Fz4 = Data.subject(i).n(j).Fz4-offsetz/4;
Data.subject(i).n(j).Fz_separate = [Data.subject(i).n(j).Fz1,Data.subject(i).n(j).Fz2,Data.subject(i).n(j).Fz3,Data.subject(i).n(j).Fz4];


i = 8
j = 19
offset_std = std(Data.subject(i).n(j).Fz(1:20))
offsetz = mean(Data.subject(i).n(j).Fz(1:20));
Data.subject(i).n(j).Fz = Data.subject(i).n(j).Fz-offsetz;
Data.subject(i).n(j).Fz1 = Data.subject(i).n(j).Fz1-offsetz/4;
Data.subject(i).n(j).Fz2 = Data.subject(i).n(j).Fz2-offsetz/4;
Data.subject(i).n(j).Fz3 = Data.subject(i).n(j).Fz3-offsetz/4;
Data.subject(i).n(j).Fz4 = Data.subject(i).n(j).Fz4-offsetz/4;
Data.subject(i).n(j).Fz_separate = [Data.subject(i).n(j).Fz1,Data.subject(i).n(j).Fz2,Data.subject(i).n(j).Fz3,Data.subject(i).n(j).Fz4];
offsetx = mean(Data.subject(i).n(j).Fx(10:40));
offsetx_std = std(Data.subject(i).n(j).Fz(10:40))
Data.subject(i).n(j).Fx = Data.subject(i).n(j).Fx-offsetx;

i = 8
j = 19
offsetx = mean(Data.subject(i).n(j).Fx(10:40));
offsetx_std = std(Data.subject(i).n(j).Fz(10:40))
Data.subject(i).n(j).Fx = Data.subject(i).n(j).Fx-offsetx;

i = 11
j = 5
offsetx = mean(Data.subject(i).n(j).Fx(20:40));
offsetx_std = std(Data.subject(i).n(j).Fz(20:40))
Data.subject(i).n(j).Fx = Data.subject(i).n(j).Fx-offsetx;

i = 2
j = 30
offsetx = mean(Data.subject(i).n(j).Fx(200:250));
offsetx_std = std(Data.subject(i).n(j).Fz(200:250))
Data.subject(i).n(j).Fx = Data.subject(i).n(j).Fx-offsetx;


i = 2
j = 32
offsetx = mean(Data.subject(i).n(j).Fx(140:180));
offsetx_std = std(Data.subject(i).n(j).Fz(140:180))
Data.subject(i).n(j).Fx = Data.subject(i).n(j).Fx-offsetx;

%% impact peak
%%% Find time corresponding to impact in force data
clc
close all

for i = Data.subjects
    for j = [Data.subject(i).n_normal Data.subject(i).n_slow]
        
        % find time derivative of force, used later to find peak values in force
        dFz = gradient(Data.subject(i).n(j).Fz);
        
        % find arbitrary position of first force increase during impact
        temp1 = find(Data.subject(i).n(j).Fz > 100);
        index_slope = temp1(1);
        
        % find beginning of slope using last below zero value force before slope        
        temp2 = find(Data.subject(i).n(j).Fz(1:index_slope)<0);
        t0_index(i) = temp2(end);
        
        % find force peak using first zero derivative after slope
        temp3 = find(dFz(index_slope:end)<0);
        tpeak_index(i) = temp3(1)+(index_slope-1)-1;
              
        Data.subject(i).n(j).t_period = Data.subject(i).n(j).t_kp(tpeak_index(i))-Data.subject(i).n(j).t_kp(t0_index(i));
        Data.subject(i).n(j).t_0 = Data.subject(i).n(j).t_kp(t0_index(i));
        Data.subject(i).n(j).t_peak = Data.subject(i).n(j).t_kp(tpeak_index(i))-Data.subject(i).n(j).t_kp(t0_index(i));
        
        Data.subject(i).n(j).t_kp_index_0 = t0_index(i);
        Data.subject(i).n(j).t_kp_index_peak = tpeak_index(i);
        
        % find beginning of slope at 25% peak force        
        temp4 = find(Data.subject(i).n(j).Fz(1:Data.subject(i).n(j).t_kp_index_peak)>(0.25*Data.subject(i).n(j).Fz(Data.subject(i).n(j).t_kp_index_peak)));
        Data.subject(i).n(j).t_kp_index_25percentpeak = temp4(1);
        
%         Data.subject(i).n(j).impulse_peak = trapz(Data.subject(i).n(j).Fz(t0_index(i):tpeak_index(i)),Data.subject(i).n(j).t_kp(t0_index(i):tpeak_index(i)));
        
        Data.subject(i).n(j).t_opto_index_0 = round(t0_index(i)*Data.f_opto/Data.f_kp);
        Data.subject(i).n(j).t_opto_index_peak = round(tpeak_index(i)*Data.f_opto/Data.f_kp);
        
%         Data.subject(i).n(j).t_kp = (1/Data.f_kp:1/Data.f_kp:length(Data.subject(i).n(j).Vz)/Data.f_kp)';
%         Data.subject(i).n(j).t_kp = Data.subject(i).n(j).t_kp - Data.subject(i).n(j).t_0;
%         Data.subject(i).n(j).t_opto = (1/Data.f_opto:1/Data.f_opto:length(Data.subject(i).n(j).z)/Data.f_opto)';
%         Data.subject(i).n(j).t_opto = Data.subject(i).n(j).t_opto - Data.subject(i).n(j).t_0;
    end
end

%% Marker distance
% Check distance between markers
% needs to be constant for rigid body assumption
% look for 'shock wave' pattern due to impact

clc
close all

for i = Data.subjects
    temp_distance = [];
    temp_h1 = [];
    temp_h2 = [];
    for j = [Data.subject(i).n_bodyweight] %Data.subject(i).n_total
        
        Data.subject(i).n(j).marker_vector = [Data.subject(i).n(j).x(:,2)-Data.subject(i).n(j).x(:,1) Data.subject(i).n(j).y(:,2)-Data.subject(i).n(j).y(:,1) Data.subject(i).n(j).z(:,2)-Data.subject(i).n(j).z(:,1)];
        Data.subject(i).n(j).marker_distance = sqrt(Data.subject(i).n(j).marker_vector(:,1).^2+Data.subject(i).n(j).marker_vector(:,2).^2+Data.subject(i).n(j).marker_vector(:,3).^2);
        
%         figure(100+i)
%         subplot(211)
%         plot(Data.subject(i).n(j).t_opto,Data.subject(i).n(j).marker_distance)
%         hold on
        
        temp_distance = [temp_distance, median(Data.subject(i).n(j).marker_distance(isfinite(Data.subject(i).n(j).marker_distance)))];
        temp_h1 = [temp_h1, median(Data.subject(i).n(j).z(isfinite(Data.subject(i).n(j).z(:,1)),1))];
        temp_h2 = [temp_h2, median(Data.subject(i).n(j).z(isfinite(Data.subject(i).n(j).z(:,2)),2))];
    end
    
    Data.subject(i).height_marker12 = median(temp_distance);
    Data.subject(i).height_marker1 = median(temp_h1);
    Data.subject(i).height_marker2 = median(temp_h2);
    
    for j = [Data.subject(i).n_normal Data.subject(i).n_slow]%Data.subject(i).n_total
        
        Data.subject(i).n(j).marker_vector = [Data.subject(i).n(j).x(:,2)-Data.subject(i).n(j).x(:,1) Data.subject(i).n(j).y(:,2)-Data.subject(i).n(j).y(:,1) Data.subject(i).n(j).z(:,2)-Data.subject(i).n(j).z(:,1)];
        Data.subject(i).n(j).marker_distance = sqrt(Data.subject(i).n(j).marker_vector(:,1).^2+Data.subject(i).n(j).marker_vector(:,2).^2+Data.subject(i).n(j).marker_vector(:,3).^2);
    end
end

%% COP

% position of the four force plate sensors in coordinate system of optotrak
Data.sensor1 = [0.080,0.100+0.400];
Data.sensor2 = [0.080+0.240,0.100+0.400];
Data.sensor3 = [0.080+0.240,0.100];
Data.sensor4 = [0.080,0.100];
Data.size_kp = [0.400,0.600];


for i = Data.subjects
    for j = [Data.subject(i).n_normal Data.subject(i).n_slow]
        COP_kpx = ((Data.subject(i).n(j).Fz2 + Data.subject(i).n(j).Fz3 - Data.subject(i).n(j).Fz1 - Data.subject(i).n(j).Fz4) * (Data.sensor2(1)-Data.sensor1(1))/2)./Data.subject(i).n(j).Fz;
        COP_kpy = -((Data.subject(i).n(j).Fz1 + Data.subject(i).n(j).Fz2 - Data.subject(i).n(j).Fz3 - Data.subject(i).n(j).Fz4) * (Data.sensor1(2)-Data.sensor4(2))/2)./Data.subject(i).n(j).Fz;
        Data.subject(i).n(j).x_COP_kp = Data.size_kp(1)/2 + COP_kpx;
        Data.subject(i).n(j).y_COP_kp = Data.size_kp(2)/2 - COP_kpy;
    end
end

%% heel
%%% Calculate heel and other lower leg data from marker positions

for i = Data.subjects
    for j = [Data.subject(i).n_normal Data.subject(i).n_slow]
        
        % rotations of lower leg, z rotation is assumed zero
        Data.subject(i).n(j).psi = asin((Data.subject(i).n(j).x(:,2)-Data.subject(i).n(j).x(:,1))./Data.subject(i).n(j).marker_distance); % rotation about y axis
        Data.subject(i).n(j).phi = sign(Data.subject(i).n(j).y(:,1)-Data.subject(i).n(j).y(:,2)) .* acos((Data.subject(i).n(j).z(:,2)-Data.subject(i).n(j).z(:,1))./(Data.subject(i).n(j).marker_distance.*cos(Data.subject(i).n(j).psi))); % rotation about the x axis
        
        % in order to calculate the local lower leg position of COP, we
        % need to interpolate the opto data to match the kp data
        
        index_opto_interp1 = Data.subject(i).n(j).t_opto_index_0:Data.subject(i).n(j).t_opto_index_peak+(Data.subject(i).n(j).t_opto_index_peak-Data.subject(i).n(j).t_opto_index_0);
        index_kp_interp1 = Data.subject(i).n(j).t_kp_index_25percentpeak:(Data.subject(i).n(j).t_kp_index_peak+Data.subject(i).n(j).t_kp_index_peak-Data.subject(i).n(j).t_kp_index_25percentpeak);%+round((Data.subject(i).n(j).t_kp_index_peak-Data.subject(i).n(j).t_kp_index_0)/2);
        Data.subject(i).n(j).x_interp1 = [interp1(Data.subject(i).n(j).t_opto(index_opto_interp1),Data.subject(i).n(j).x(index_opto_interp1,1),Data.subject(i).n(j).t_kp(index_kp_interp1)) interp1(Data.subject(i).n(j).t_opto(index_opto_interp1),Data.subject(i).n(j).x(index_opto_interp1,2),Data.subject(i).n(j).t_kp(index_kp_interp1))];
        Data.subject(i).n(j).y_interp1 = [interp1(Data.subject(i).n(j).t_opto(index_opto_interp1),Data.subject(i).n(j).y(index_opto_interp1,1),Data.subject(i).n(j).t_kp(index_kp_interp1)) interp1(Data.subject(i).n(j).t_opto(index_opto_interp1),Data.subject(i).n(j).y(index_opto_interp1,2),Data.subject(i).n(j).t_kp(index_kp_interp1))];
        Data.subject(i).n(j).z_interp1 = [interp1(Data.subject(i).n(j).t_opto(index_opto_interp1),Data.subject(i).n(j).z(index_opto_interp1,1),Data.subject(i).n(j).t_kp(index_kp_interp1)) interp1(Data.subject(i).n(j).t_opto(index_opto_interp1),Data.subject(i).n(j).z(index_opto_interp1,2),Data.subject(i).n(j).t_kp(index_kp_interp1))];
        Data.subject(i).n(j).psi_interp1 = interp1(Data.subject(i).n(j).t_opto(index_opto_interp1),Data.subject(i).n(j).psi(index_opto_interp1),Data.subject(i).n(j).t_kp(index_kp_interp1));
        Data.subject(i).n(j).phi_interp1 = interp1(Data.subject(i).n(j).t_opto(index_opto_interp1),Data.subject(i).n(j).phi(index_opto_interp1),Data.subject(i).n(j).t_kp(index_kp_interp1));
        
        % calculate local position COP on lower leg using rotation matrices
        % and a coordinate system originated in marker 1, z axis through
        % marker 1 and 2, y axis stays in y-z plane of global opto system
        Data.subject(i).n(j).x_COP_p = mean(Data.subject(i).n(j).z_interp1(:,1).*cos(Data.subject(i).n(j).phi_interp1).*sin(Data.subject(i).n(j).psi_interp1) - cos(Data.subject(i).n(j).psi_interp1).*(Data.subject(i).n(j).x_interp1(:,1)-Data.subject(i).n(j).x_COP_kp(index_kp_interp1)) - sin(Data.subject(i).n(j).psi_interp1).*sin(Data.subject(i).n(j).phi_interp1).*(Data.subject(i).n(j).y_interp1(:,1)-Data.subject(i).n(j).y_COP_kp(index_kp_interp1)));
        Data.subject(i).n(j).y_COP_p = mean(- Data.subject(i).n(j).z_interp1(:,1).*sin(Data.subject(i).n(j).phi_interp1) - cos(Data.subject(i).n(j).phi_interp1).*(Data.subject(i).n(j).y_interp1(:,1) - Data.subject(i).n(j).y_COP_kp(index_kp_interp1)));
        Data.subject(i).n(j).z_COP_p = mean(cos(Data.subject(i).n(j).psi_interp1).*sin(Data.subject(i).n(j).phi_interp1).*(Data.subject(i).n(j).y_interp1(:,1) - Data.subject(i).n(j).y_COP_kp(index_kp_interp1)) - Data.subject(i).n(j).z_interp1(:,1).*cos(Data.subject(i).n(j).phi_interp1).*cos(Data.subject(i).n(j).psi_interp1) - sin(Data.subject(i).n(j).psi_interp1).*(Data.subject(i).n(j).x_interp1(:,1)-Data.subject(i).n(j).x_COP_kp(index_kp_interp1)));
        
        Data.subject(i).n(j).x_p = Data.subject(i).n(j).z_interp1(:,1).*cos(Data.subject(i).n(j).phi_interp1).*sin(Data.subject(i).n(j).psi_interp1) - cos(Data.subject(i).n(j).psi_interp1).*(Data.subject(i).n(j).x_interp1(:,1)-Data.subject(i).n(j).x_COP_kp(index_kp_interp1)) - sin(Data.subject(i).n(j).psi_interp1).*sin(Data.subject(i).n(j).phi_interp1).*(Data.subject(i).n(j).y_interp1(:,1)-Data.subject(i).n(j).y_COP_kp(index_kp_interp1));
        Data.subject(i).n(j).y_p = - Data.subject(i).n(j).z_interp1(:,1).*sin(Data.subject(i).n(j).phi_interp1) - cos(Data.subject(i).n(j).phi_interp1).*(Data.subject(i).n(j).y_interp1(:,1) - Data.subject(i).n(j).y_COP_kp(index_kp_interp1));
        Data.subject(i).n(j).z_p = cos(Data.subject(i).n(j).psi_interp1).*sin(Data.subject(i).n(j).phi_interp1).*(Data.subject(i).n(j).y_interp1(:,1) - Data.subject(i).n(j).y_COP_kp(index_kp_interp1)) - Data.subject(i).n(j).z_interp1(:,1).*cos(Data.subject(i).n(j).phi_interp1).*cos(Data.subject(i).n(j).psi_interp1) - sin(Data.subject(i).n(j).psi_interp1).*(Data.subject(i).n(j).x_interp1(:,1)-Data.subject(i).n(j).x_COP_kp(index_kp_interp1));
        
        Data.subject(i).n(j).x_COP_Q = Data.subject(i).n(j).x_COP_p + 0.01;
        Data.subject(i).n(j).y_COP_Q = Data.subject(i).n(j).y_COP_p - 0.01;
        Data.subject(i).n(j).z_COP_Q = Data.subject(i).n(j).z_COP_p - 0.01;
       
        
        % calculate the global position of the COP prime using rotation
        % matrices and angles phi and psi
        
        %%%% CORRECT HEEL POINT H, this is called point COP throughout de code, since we started by investigated multiple possible heel points using various methods, this was not presented in the article
        Data.subject(i).n(j).x_COP = Data.subject(i).n(j).x(:,1) + Data.subject(i).n(j).x_COP_p.*cos(Data.subject(i).n(j).psi) + Data.subject(i).n(j).z_COP_p.*sin(Data.subject(i).n(j).psi);
        Data.subject(i).n(j).y_COP = Data.subject(i).n(j).y(:,1) + Data.subject(i).n(j).y_COP_p.*cos(Data.subject(i).n(j).phi) - Data.subject(i).n(j).z_COP_p.*cos(Data.subject(i).n(j).psi).*sin(Data.subject(i).n(j).phi) + Data.subject(i).n(j).x_COP_p.*sin(Data.subject(i).n(j).phi).*sin(Data.subject(i).n(j).psi);
        Data.subject(i).n(j).z_COP = Data.subject(i).n(j).z(:,1) + Data.subject(i).n(j).y_COP_p.*sin(Data.subject(i).n(j).phi) + Data.subject(i).n(j).z_COP_p.*cos(Data.subject(i).n(j).phi).*cos(Data.subject(i).n(j).psi) - Data.subject(i).n(j).x_COP_p.*cos(Data.subject(i).n(j).phi).*sin(Data.subject(i).n(j).psi);
        
        %%%% Use Q HEEL POINT - UNCOMMENT FOR PARAMETER SENSITIVITY ANALYSIS
%         Data.subject(i).n(j).x_COP = Data.subject(i).n(j).x(:,1) + Data.subject(i).n(j).x_COP_Q.*cos(Data.subject(i).n(j).psi) + Data.subject(i).n(j).z_COP_Q.*sin(Data.subject(i).n(j).psi);
%         Data.subject(i).n(j).y_COP = Data.subject(i).n(j).y(:,1) + Data.subject(i).n(j).y_COP_Q.*cos(Data.subject(i).n(j).phi) - Data.subject(i).n(j).z_COP_Q.*cos(Data.subject(i).n(j).psi).*sin(Data.subject(i).n(j).phi) + Data.subject(i).n(j).x_COP_Q.*sin(Data.subject(i).n(j).phi).*sin(Data.subject(i).n(j).psi);
%         Data.subject(i).n(j).z_COP = Data.subject(i).n(j).z(:,1) + Data.subject(i).n(j).y_COP_Q.*sin(Data.subject(i).n(j).phi) + Data.subject(i).n(j).z_COP_Q.*cos(Data.subject(i).n(j).phi).*cos(Data.subject(i).n(j).psi) - Data.subject(i).n(j).x_COP_Q.*cos(Data.subject(i).n(j).phi).*sin(Data.subject(i).n(j).psi);
                
    end
end

%% filter 
%%% filter position data
% close all

for i = Data.subjects
    for j = [Data.subject(i).n_normal Data.subject(i).n_slow]
        
        % Find index of data surrounding impact    
        index_0 = Data.subject(i).n(j).t_opto_index_0;
        temp_before = find(isnan(Data.subject(i).n(j).z_COP(1:index_0)));
        if isempty(temp_before)
            index_first = 1;
        else
            index_first = (temp_before(end)+1);
        end
        temp_after = (find(isnan(Data.subject(i).n(j).z_COP(index_0:end))))+(index_0-1);
        data_index = index_first:(temp_after(1)-1);
       
        % Filter this index of data using a forward-backward butterworth
        % filter with cutoff frequency = 0.5 * 0.5 * samplerate = 323.3750 Hz
        % 1293.5*0.5*0.3 = 194.0250
        % 1293.5*0.5*0.16 = 103.4800
        % forward-backward filtering compensates for delay caused by filter
        [B,A] = butter(2,0.15); % according to article; vectors B (numerator) and A (denominator).The cutoff frequency Wn must be 0.0 < Wn < 1.0, with 1.0 corresponding to half the sample rate
%         [B,A] = butter(2,0.3); % testing cut-off frequency; vectors B (numerator) and A (denominator).The cutoff frequency Wn must be 0.0 < Wn < 1.0, with 1.0 corresponding to half the sample rate
        Data.subject(i).n(j).x_filtered = [filtfilt(B, A, Data.subject(i).n(j).x(data_index,1)) filtfilt(B, A, Data.subject(i).n(j).x(data_index,2))];
        Data.subject(i).n(j).y_filtered = [filtfilt(B, A, Data.subject(i).n(j).y(data_index,1)) filtfilt(B, A, Data.subject(i).n(j).y(data_index,2))];
        Data.subject(i).n(j).z_filtered = [filtfilt(B, A, Data.subject(i).n(j).z(data_index,1)) filtfilt(B, A, Data.subject(i).n(j).z(data_index,2))];
        
        Data.subject(i).n(j).z_COP_filtered = filtfilt(B, A, Data.subject(i).n(j).z_COP(data_index));
        Data.subject(i).n(j).y_COP_filtered = filtfilt(B, A, Data.subject(i).n(j).y_COP(data_index));
        Data.subject(i).n(j).x_COP_filtered = filtfilt(B, A, Data.subject(i).n(j).x_COP(data_index));

        %%% UNCOMMENT TO CHECK INFLUENCE OF FILTER ON RESULT (CHECKED: negliable, would not influence values stated in article)
%         Data.subject(i).n(j).z_COP_filtered = Data.subject(i).n(j).z_COP(data_index);
%         Data.subject(i).n(j).y_COP_filtered = Data.subject(i).n(j).y_COP(data_index);
%         Data.subject(i).n(j).x_COP_filtered = Data.subject(i).n(j).x_COP(data_index);
        
        Data.subject(i).n(j).t_opto_filtered = Data.subject(i).n(j).t_opto(data_index);
        Data.subject(i).n(j).index_filtered = data_index;
        Data.subject(i).n(j).t_opto_filtered_index_0 = Data.subject(i).n(j).t_opto_index_0-(data_index(1)-1);
        Data.subject(i).n(j).t_opto_filtered_index_peak = Data.subject(i).n(j).t_opto_index_peak-(data_index(1)-1);
        
        % get velocity and acceleration from position 
        Data.subject(i).n(j).zd_filtered(:,1) = gradient(Data.subject(i).n(j).z_filtered(:,1),1/Data.f_opto);
        Data.subject(i).n(j).yd_filtered(:,1) = gradient(Data.subject(i).n(j).y_filtered(:,1),1/Data.f_opto);
        Data.subject(i).n(j).xd_filtered(:,1) = gradient(Data.subject(i).n(j).x_filtered(:,1),1/Data.f_opto);
        
        Data.subject(i).n(j).zd_filtered(:,2) = gradient(Data.subject(i).n(j).z_filtered(:,2),1/Data.f_opto);
        Data.subject(i).n(j).yd_filtered(:,2) = gradient(Data.subject(i).n(j).y_filtered(:,2),1/Data.f_opto);
        Data.subject(i).n(j).xd_filtered(:,2) = gradient(Data.subject(i).n(j).x_filtered(:,2),1/Data.f_opto);
        
        Data.subject(i).n(j).zdd_filtered(:,1) = gradient(Data.subject(i).n(j).zd_filtered(:,1),1/Data.f_opto);
        Data.subject(i).n(j).ydd_filtered(:,1) = gradient(Data.subject(i).n(j).yd_filtered(:,1),1/Data.f_opto);
        Data.subject(i).n(j).xdd_filtered(:,1) = gradient(Data.subject(i).n(j).xd_filtered(:,1),1/Data.f_opto);
        
        Data.subject(i).n(j).zdd_filtered(:,2) = gradient(Data.subject(i).n(j).zd_filtered(:,2),1/Data.f_opto);
        Data.subject(i).n(j).ydd_filtered(:,2) = gradient(Data.subject(i).n(j).yd_filtered(:,2),1/Data.f_opto);
        Data.subject(i).n(j).xdd_filtered(:,2) = gradient(Data.subject(i).n(j).xd_filtered(:,2),1/Data.f_opto);
  
        Data.subject(i).n(j).zd_COP_filtered = gradient(Data.subject(i).n(j).z_COP_filtered,1/Data.f_opto);
        Data.subject(i).n(j).yd_COP_filtered = gradient(Data.subject(i).n(j).y_COP_filtered,1/Data.f_opto);
        Data.subject(i).n(j).xd_COP_filtered = gradient(Data.subject(i).n(j).x_COP_filtered,1/Data.f_opto);
        Data.subject(i).n(j).zdd_COP_filtered = gradient(Data.subject(i).n(j).zd_COP_filtered,1/Data.f_opto);
        Data.subject(i).n(j).ydd_COP_filtered = gradient(Data.subject(i).n(j).yd_COP_filtered,1/Data.f_opto);
        Data.subject(i).n(j).xdd_COP_filtered = gradient(Data.subject(i).n(j).xd_COP_filtered,1/Data.f_opto);
    end
end

%% zd0
%%% Find time where heel z velocity is zero after peak impact, and calculate impulse using this point

count_heel = 0;
count_COP = 0;

for i = Data.subjects
    for j = [Data.subject(i).n_normal Data.subject(i).n_slow]
                        
        %%% calculations for COP point
        temp_zd0_COP = find(Data.subject(i).n(j).zd_COP_filtered(Data.subject(i).n(j).t_opto_filtered_index_peak:Data.subject(i).n(j).t_opto_filtered_index_peak+2*(Data.subject(i).n(j).t_opto_filtered_index_peak-Data.subject(i).n(j).t_opto_filtered_index_0))>=0);
        if isempty(temp_zd0_COP)
            count_COP = count_COP+1;
            temp_zd0_COP = find(Data.subject(i).n(j).zdd_COP_filtered(Data.subject(i).n(j).t_opto_filtered_index_peak:Data.subject(i).n(j).t_opto_filtered_index_peak+2*(Data.subject(i).n(j).t_opto_filtered_index_peak-Data.subject(i).n(j).t_opto_filtered_index_0))<=0);
            if isempty(temp_zd0_COP)
                disp('Error: cannot find zd_COP = 0 or zdd_COP=0')
            end
        end
        Data.subject(i).n(j).t_opto_filtered_index_zd0_COP = temp_zd0_COP(1) + (Data.subject(i).n(j).t_opto_filtered_index_peak-1);
        Data.subject(i).n(j).t_opto_index_zd0_COP = Data.subject(i).n(j).t_opto_filtered_index_zd0_COP + (Data.subject(i).n(j).index_filtered(1)-1);
        Data.subject(i).n(j).t_kp_index_zd0_COP = round(Data.subject(i).n(j).t_opto_index_zd0_COP*Data.f_kp/Data.f_opto);
    end
end

%% step time

Matrix.non_complete_position = ones(13,32);
for i = Data.subjects
    for j = [Data.subject(i).n_normal Data.subject(i).n_slow]
        
        % find time that the right foot is lifted
        Data.subject(i).n(j).t_kp_index_lift = Data.subject(i).n(j).t_kp_index_0 + find(Data.subject(i).n(j).Fz(Data.subject(i).n(j).t_kp_index_0+1:end)<0,1);
        if isempty(Data.subject(i).n(j).t_kp_index_lift) % not found in 3 cases
            Matrix.non_complete_force(i,j) = 1;
            Data.subject(i).n(j).t_kp_index_lift = length(Data.subject(i).n(j).Fz);
        end
        Data.subject(i).n(j).t_opto_index_lift = round(Data.subject(i).n(j).t_kp_index_lift*Data.f_opto/Data.f_kp);
        
        minimum = max(min(Data.subject(i).n(j).z_COP(1:Data.subject(i).n(j).t_opto_index_lift)),min(Data.subject(i).n(j).z_COP(Data.subject(i).n(j).t_opto_index_lift:end)));
        maximum = min(max(Data.subject(i).n(j).z_COP(1:Data.subject(i).n(j).t_opto_index_lift)),max(Data.subject(i).n(j).z_COP(Data.subject(i).n(j).t_opto_index_lift:end)));
        if minimum>maximum
%             disp(strcat('Error: minimum > maximum for i=',num2str(i),', j=',num2str(j),', minimum=',num2str(minimum),', maximum=',num2str(maximum)))
            Matrix.non_complete_position(i,j) = NaN; % this error happens if beginning of swing is not captured, and therefore not 2 swings are seen in data
            Data.subject(i).n(j).step_index_1 = 1;
            Data.subject(i).n(j).step_index_end = find(Data.subject(i).n(j).z_COP>0,1,'last');
        else
            target = (maximum-minimum)/2+minimum;
            Data.subject(i).n(j).step_index_1 = find(Data.subject(i).n(j).z_COP(1:Data.subject(i).n(j).t_opto_index_lift)<target,1);
            Data.subject(i).n(j).step_index_end = Data.subject(i).n(j).t_opto_index_lift - 1 + find(Data.subject(i).n(j).z_COP(Data.subject(i).n(j).t_opto_index_lift:end)<target,1);
            if isempty(Data.subject(i).n(j).step_index_end)==1
                Matrix.non_complete_position(i,j) = NaN;
                Data.subject(i).n(j).step_index_1 = 1;
                Data.subject(i).n(j).step_index_end = find(Data.subject(i).n(j).z_COP>0,1,'last');
            end
        end
        Data.subject(i).n(j).step_time = (Data.subject(i).n(j).t_opto(Data.subject(i).n(j).step_index_end)-Data.subject(i).n(j).t_opto(Data.subject(i).n(j).step_index_1))/2;
        Data.subject(i).n(j).step_length = (Data.subject(i).n(j).y(Data.subject(i).n(j).step_index_end,1)-Data.subject(i).n(j).y(Data.subject(i).n(j).step_index_1,1))/2;
        Data.subject(i).n(j).step_index_opto = round((Data.subject(i).n(j).step_index_end - Data.subject(i).n(j).step_index_1)/2);
        Data.subject(i).n(j).step_index_kp = round((Data.subject(i).n(j).step_index_end - Data.subject(i).n(j).step_index_1)/2*Data.f_kp/Data.f_opto);
        Data.subject(i).n(j).v_average =  Data.subject(i).n(j).step_length/ Data.subject(i).n(j).step_time;
        Matrix.v_average(i,j) = Data.subject(i).n(j).v_average;
    end
end

%% Effective mass
clc
Data.g = 9.81; % [m/s^2] gravitaional acceleration

% M_{eff} = \frac{\int_{t_0}^{t_{end}} F_{pad}\ \text{d}t}{\Delta v + g\ \Delta t}
for i = Data.subjects
    for j = [Data.subject(i).n_normal Data.subject(i).n_slow]
        
        Data.subject(i).n(j).t_kp_adjusted = (1/Data.f_kp:1/Data.f_kp:length(Data.subject(i).n(j).Vz)/Data.f_kp)'-Data.subject(i).n(j).t_0;
        Data.subject(i).n(j).t_opto_adjusted = (1/Data.f_opto:1/Data.f_opto:length(Data.subject(i).n(j).z)/Data.f_opto)'-Data.subject(i).n(j).t_0;
        Data.subject(i).n(j).t_opto_filtered_adjusted = Data.subject(i).n(j).t_opto_filtered - Data.subject(i).n(j).t_0;
        
        Data.subject(i).n(j).impulse_zd0_COP = trapz(Data.subject(i).n(j).Fz(Data.subject(i).n(j).t_kp_index_0:Data.subject(i).n(j).t_kp_index_zd0_COP))/Data.f_kp;
        Data.subject(i).n(j).impulse_peak = trapz(Data.subject(i).n(j).Fz(Data.subject(i).n(j).t_kp_index_0:Data.subject(i).n(j).t_kp_index_peak))/Data.f_kp;
        
        Data.subject(i).n(j).impulse_y_zd0_COP = trapz(Data.subject(i).n(j).Fy(Data.subject(i).n(j).t_kp_index_0:Data.subject(i).n(j).t_kp_index_zd0_COP))/Data.f_kp;
        Data.subject(i).n(j).impulse_y_peak = trapz(Data.subject(i).n(j).Fy(Data.subject(i).n(j).t_kp_index_0:Data.subject(i).n(j).t_kp_index_peak))/Data.f_kp;

        Data.subject(i).n(j).M_eff_zd0_COP = Data.subject(i).n(j).impulse_zd0_COP/((Data.subject(i).n(j).zd_COP_filtered(Data.subject(i).n(j).t_opto_filtered_index_zd0_COP)-Data.subject(i).n(j).zd_COP_filtered(Data.subject(i).n(j).t_opto_filtered_index_0))+Data.g*(Data.subject(i).n(j).t_opto(Data.subject(i).n(j).t_opto_index_zd0_COP)-Data.subject(i).n(j).t_0));
        Data.subject(i).n(j).M_eff_zd0_COP_percentageMB = Data.subject(i).n(j).M_eff_zd0_COP/(Data.subject(i).bodyweight/Data.g)*100;
        Matrix.M_eff_zd0_COP_percentageMB(i,j) = Data.subject(i).n(j).M_eff_zd0_COP_percentageMB;
        
        Data.subject(i).n(j).zdd_COP_filtered_interp1 = interp1(Data.subject(i).n(j).t_opto_filtered(Data.subject(i).n(j).t_opto_filtered_index_0-5:Data.subject(i).n(j).t_opto_filtered_index_zd0_COP+100),Data.subject(i).n(j).zdd_COP_filtered(Data.subject(i).n(j).t_opto_filtered_index_0-5:Data.subject(i).n(j).t_opto_filtered_index_zd0_COP+100),Data.subject(i).n(j).t_kp(Data.subject(i).n(j).t_kp_index_0:Data.subject(i).n(j).t_kp_index_zd0_COP+50));
        Data.subject(i).n(j).M_eff_COP_instant = Data.subject(i).n(j).Fz(Data.subject(i).n(j).t_kp_index_0:Data.subject(i).n(j).t_kp_index_zd0_COP+50)./(Data.subject(i).n(j).zdd_COP_filtered_interp1+Data.g);
        
        Data.subject(i).n(j).M_eff_COP_instant_max = max(Data.subject(i).n(j).M_eff_COP_instant);
        Data.subject(i).n(j).M_eff_COP_instant_max_percentageMB = Data.subject(i).n(j).M_eff_COP_instant_max/(Data.subject(i).bodyweight/Data.g)*100;
        Matrix.M_eff_COP_instant_max_percentageMB(i,j) = Data.subject(i).n(j).M_eff_COP_instant_max_percentageMB;
        
        Data.subject(i).n(j).M_eff_peak_COP = Data.subject(i).n(j).impulse_peak/((Data.subject(i).n(j).zd_COP_filtered(Data.subject(i).n(j).t_opto_filtered_index_peak)-Data.subject(i).n(j).zd_COP_filtered(Data.subject(i).n(j).t_opto_filtered_index_0))+Data.g*(Data.subject(i).n(j).t_opto(Data.subject(i).n(j).t_opto_index_peak)-Data.subject(i).n(j).t_0));
        Data.subject(i).n(j).M_eff_peak_COP_percentageMB = Data.subject(i).n(j).M_eff_peak_COP/(Data.subject(i).bodyweight/Data.g)*100;
        Matrix.M_eff_peak_COP_percentageMB(i,j) = Data.subject(i).n(j).M_eff_peak_COP_percentageMB;

        Data.subject(i).n(j).M_eff_Chi_COP = Data.subject(i).n(j).impulse_peak/((-Data.subject(i).n(j).zd_COP_filtered(Data.subject(i).n(j).t_opto_filtered_index_0))+Data.g*(Data.subject(i).n(j).t_opto(Data.subject(i).n(j).t_opto_index_peak)-Data.subject(i).n(j).t_0));
        Data.subject(i).n(j).M_eff_Chi_COP_percentageMB = Data.subject(i).n(j).M_eff_Chi_COP/(Data.subject(i).bodyweight/Data.g)*100;
        Matrix.M_eff_Chi_COP_percentageMB(i,j) = Data.subject(i).n(j).M_eff_Chi_COP_percentageMB;
        
    end
end

%% Matrices for plot assistance

Matrix.n_normal = NaN(13,32);
Matrix.n_slow = NaN(13,32);
for i = Data.subjects
    for j = [Data.subject(i).n_normal]
        Matrix.n_normal(i,j) = 1;
    end 
    for j = [Data.subject(i).n_slow]
        Matrix.n_slow(i,j) = 1;
    end
    for j = [Data.subject(i).n_normal,Data.subject(i).n_slow]
        Matrix.n_all(i,j) = 1;
    end
end

%% Energy loss Meff

for i = Data.subjects
    for j = [Data.subject(i).n_normal Data.subject(i).n_slow]
         
        Data.subject(i).n(j).E_M_eff_peak = 0.5 * Data.subject(i).n(j).M_eff_peak_COP * (Data.subject(i).n(j).zd_COP_filtered(Data.subject(i).n(j).t_opto_filtered_index_peak)^2 - Data.subject(i).n(j).zd_COP_filtered(Data.subject(i).n(j).t_opto_filtered_index_0)^2) + Data.subject(i).n(j).M_eff_peak_COP * Data.g * (Data.subject(i).n(j).z_COP_filtered(Data.subject(i).n(j).t_opto_filtered_index_peak) - Data.subject(i).n(j).z_COP_filtered(Data.subject(i).n(j).t_opto_filtered_index_0));
        Matrix.E_M_eff_peak(i,j) = Data.subject(i).n(j).E_M_eff_peak;

        Data.subject(i).n(j).E_M_eff_Chi = 0.5 * Data.subject(i).n(j).M_eff_Chi_COP * (- Data.subject(i).n(j).zd_COP_filtered(Data.subject(i).n(j).t_opto_filtered_index_0)^2) + Data.subject(i).n(j).M_eff_Chi_COP * Data.g * (Data.subject(i).n(j).z_COP_filtered(Data.subject(i).n(j).t_opto_filtered_index_peak) - Data.subject(i).n(j).z_COP_filtered(Data.subject(i).n(j).t_opto_filtered_index_0));
        Matrix.E_M_eff_Chi(i,j) = Data.subject(i).n(j).E_M_eff_Chi;
        
        Data.subject(i).n(j).E_M_eff_zd0 = 0.5 * Data.subject(i).n(j).M_eff_zd0_COP * (Data.subject(i).n(j).zd_COP_filtered(Data.subject(i).n(j).t_opto_filtered_index_zd0_COP)^2 - Data.subject(i).n(j).zd_COP_filtered(Data.subject(i).n(j).t_opto_filtered_index_0)^2) + Data.subject(i).n(j).M_eff_zd0_COP * Data.g * (Data.subject(i).n(j).z_COP_filtered(Data.subject(i).n(j).t_opto_filtered_index_zd0_COP) - Data.subject(i).n(j).z_COP_filtered(Data.subject(i).n(j).t_opto_filtered_index_0));
        Matrix.E_M_eff_zd0(i,j) = Data.subject(i).n(j).E_M_eff_zd0;

        Data.subject(i).n(j).E_M_eff_peak_zd0 = 0.5 * Data.subject(i).n(j).M_eff_peak_COP * (Data.subject(i).n(j).zd_COP_filtered(Data.subject(i).n(j).t_opto_filtered_index_zd0_COP)^2 - Data.subject(i).n(j).zd_COP_filtered(Data.subject(i).n(j).t_opto_filtered_index_0)^2) + Data.subject(i).n(j).M_eff_peak_COP * Data.g * (Data.subject(i).n(j).z_COP_filtered(Data.subject(i).n(j).t_opto_filtered_index_zd0_COP) - Data.subject(i).n(j).z_COP_filtered(Data.subject(i).n(j).t_opto_filtered_index_0));
        Matrix.E_M_eff_peak_zd0(i,j) = Data.subject(i).n(j).E_M_eff_peak_zd0;
        
        Data.subject(i).n(j).E_metabolic = (32 + 0.005 * (Data.subject(i).n(j).v_average*60)^2) * (Data.subject(i).n(j).step_time/60) * (Data.subject(i).bodyweight/Data.g) * 4.184;
        Matrix.E_metabolic(i,j) = Data.subject(i).n(j).E_metabolic;
        
        Data.subject(i).n(j).E_metabolic_Kuo = 0.2 * Data.subject(i).bodyweight * Data.subject(i).n(j).step_length;
        Matrix.E_metabolic_Kuo(i,j) = Data.subject(i).n(j).E_metabolic_Kuo;
    end
end
 
%% Work
 
%  lineair inpterpolation: YI = interp1(X,Y,XI)
%  Z = trapz(X,Y) computes the integral of Y with respect to X
%   Q = quad(FUN,A,B) tries to approximate the integral of scalar-valued function FUN from A to B to within an error of 1.e-6 using recursive adaptive Simpson quadrature

for i = Data.subjects
    for j = [Data.subject(i).n_normal Data.subject(i).n_slow]
        Data.subject(i).n(j).z_COP_filtered_interp1_peak = interp1(Data.subject(i).n(j).t_opto_filtered(Data.subject(i).n(j).t_opto_filtered_index_0-5:Data.subject(i).n(j).t_opto_filtered_index_peak+5),Data.subject(i).n(j).z_COP_filtered(Data.subject(i).n(j).t_opto_filtered_index_0-5:Data.subject(i).n(j).t_opto_filtered_index_peak+5),Data.subject(i).n(j).t_kp(Data.subject(i).n(j).t_kp_index_0:Data.subject(i).n(j).t_kp_index_peak));
        Data.subject(i).n(j).Work_z_peak = trapz(Data.subject(i).n(j).z_COP_filtered_interp1_peak,Data.subject(i).n(j).Fz(Data.subject(i).n(j).t_kp_index_0:Data.subject(i).n(j).t_kp_index_peak));
        Data.subject(i).n(j).x_COP_filtered_interp1_peak = interp1(Data.subject(i).n(j).t_opto_filtered(Data.subject(i).n(j).t_opto_filtered_index_0-5:Data.subject(i).n(j).t_opto_filtered_index_peak+5),Data.subject(i).n(j).x_COP_filtered(Data.subject(i).n(j).t_opto_filtered_index_0-5:Data.subject(i).n(j).t_opto_filtered_index_peak+5),Data.subject(i).n(j).t_kp(Data.subject(i).n(j).t_kp_index_0:Data.subject(i).n(j).t_kp_index_peak));
        Data.subject(i).n(j).Work_x_peak = trapz(Data.subject(i).n(j).x_COP_filtered_interp1_peak,Data.subject(i).n(j).Fx(Data.subject(i).n(j).t_kp_index_0:Data.subject(i).n(j).t_kp_index_peak));
        Data.subject(i).n(j).y_COP_filtered_interp1_peak = interp1(Data.subject(i).n(j).t_opto_filtered(Data.subject(i).n(j).t_opto_filtered_index_0-5:Data.subject(i).n(j).t_opto_filtered_index_peak+5),Data.subject(i).n(j).y_COP_filtered(Data.subject(i).n(j).t_opto_filtered_index_0-5:Data.subject(i).n(j).t_opto_filtered_index_peak+5),Data.subject(i).n(j).t_kp(Data.subject(i).n(j).t_kp_index_0:Data.subject(i).n(j).t_kp_index_peak));
        Data.subject(i).n(j).Work_y_peak = trapz(Data.subject(i).n(j).y_COP_filtered_interp1_peak,Data.subject(i).n(j).Fy(Data.subject(i).n(j).t_kp_index_0:Data.subject(i).n(j).t_kp_index_peak));
        Data.subject(i).n(j).Work_peak = Data.subject(i).n(j).Work_x_peak+Data.subject(i).n(j).Work_y_peak+Data.subject(i).n(j).Work_z_peak;

        Matrix.Work_x_peak(i,j) = Data.subject(i).n(j).Work_x_peak;
        Matrix.Work_y_peak(i,j) = Data.subject(i).n(j).Work_y_peak;
        Matrix.Work_z_peak(i,j) = Data.subject(i).n(j).Work_z_peak;
        Matrix.Work_peak(i,j) = Data.subject(i).n(j).Work_peak;
        
        Data.subject(i).n(j).z_COP_filtered_interp1_zd0 = interp1(Data.subject(i).n(j).t_opto_filtered(Data.subject(i).n(j).t_opto_filtered_index_0-5:Data.subject(i).n(j).t_opto_filtered_index_zd0_COP+5),Data.subject(i).n(j).z_COP_filtered(Data.subject(i).n(j).t_opto_filtered_index_0-5:Data.subject(i).n(j).t_opto_filtered_index_zd0_COP+5),Data.subject(i).n(j).t_kp(Data.subject(i).n(j).t_kp_index_0:Data.subject(i).n(j).t_kp_index_zd0_COP));
        Data.subject(i).n(j).Work_zd0 = trapz(Data.subject(i).n(j).z_COP_filtered_interp1_zd0,Data.subject(i).n(j).Fz(Data.subject(i).n(j).t_kp_index_0:Data.subject(i).n(j).t_kp_index_zd0_COP));
        
        Data.subject(i).n(j).x_COP_filtered_interp1_zd0 = interp1(Data.subject(i).n(j).t_opto_filtered(Data.subject(i).n(j).t_opto_filtered_index_0-5:Data.subject(i).n(j).t_opto_filtered_index_zd0_COP+5),Data.subject(i).n(j).x_COP_filtered(Data.subject(i).n(j).t_opto_filtered_index_0-5:Data.subject(i).n(j).t_opto_filtered_index_zd0_COP+5),Data.subject(i).n(j).t_kp(Data.subject(i).n(j).t_kp_index_0:Data.subject(i).n(j).t_kp_index_zd0_COP));
        Data.subject(i).n(j).Work_x = trapz(Data.subject(i).n(j).x_COP_filtered_interp1_zd0,Data.subject(i).n(j).Fx(Data.subject(i).n(j).t_kp_index_0:Data.subject(i).n(j).t_kp_index_zd0_COP));
        Data.subject(i).n(j).y_COP_filtered_interp1_zd0 = interp1(Data.subject(i).n(j).t_opto_filtered(Data.subject(i).n(j).t_opto_filtered_index_0-5:Data.subject(i).n(j).t_opto_filtered_index_zd0_COP+5),Data.subject(i).n(j).y_COP_filtered(Data.subject(i).n(j).t_opto_filtered_index_0-5:Data.subject(i).n(j).t_opto_filtered_index_zd0_COP+5),Data.subject(i).n(j).t_kp(Data.subject(i).n(j).t_kp_index_0:Data.subject(i).n(j).t_kp_index_zd0_COP));
        Data.subject(i).n(j).Work_y = trapz(Data.subject(i).n(j).y_COP_filtered_interp1_zd0,Data.subject(i).n(j).Fy(Data.subject(i).n(j).t_kp_index_0:Data.subject(i).n(j).t_kp_index_zd0_COP));
        Data.subject(i).n(j).Work = Data.subject(i).n(j).Work_x+Data.subject(i).n(j).Work_y+Data.subject(i).n(j).Work_zd0;
        
        Matrix.Work_zd0(i,j) = Data.subject(i).n(j).Work_zd0;
        Matrix.Work_x(i,j) = Data.subject(i).n(j).Work_x;
        Matrix.Work_y(i,j) = Data.subject(i).n(j).Work_y;
        Matrix.Work(i,j) = Data.subject(i).n(j).Work;
      
    end
end

%% Put plot data in matrix

for i = Data.subjects
    for j = [Data.subject(i).n_normal Data.subject(i).n_slow]
       % impact duration
       Matrix.t_peak_adjusted(i,j) = Data.subject(i).n(j).t_kp_adjusted(Data.subject(i).n(j).t_kp_index_peak);
       Matrix.t_zd0_adjusted(i,j) = Data.subject(i).n(j).t_opto_adjusted(Data.subject(i).n(j).t_opto_index_zd0_COP);
       
       % heel velocity
       Matrix.zd_0_adjusted(i,j) = Data.subject(i).n(j).zd_COP_filtered(Data.subject(i).n(j).t_opto_filtered_index_0);
       Matrix.zd_peak_adjusted(i,j) = Data.subject(i).n(j).zd_COP_filtered(Data.subject(i).n(j).t_opto_filtered_index_peak);
       
       % max deformation
       Data.subject(i).n(j).heel_deformation = Data.subject(i).n(j).z_COP_filtered(Data.subject(i).n(j).t_opto_filtered_index_0)-Data.subject(i).n(j).z_COP_filtered(Data.subject(i).n(j).t_opto_filtered_index_zd0_COP);
       Matrix.heel_deformation(i,j) = Data.subject(i).n(j).heel_deformation;
       
       % max deformation components
       Data.subject(i).n(j).heel_deformation_x = Data.subject(i).n(j).x_COP_filtered(Data.subject(i).n(j).t_opto_filtered_index_zd0_COP)-Data.subject(i).n(j).x_COP_filtered(Data.subject(i).n(j).t_opto_filtered_index_0);
       Matrix.heel_deformation_x(i,j) = Data.subject(i).n(j).heel_deformation_x;
       Matrix.heel_deformation_x_peak(i,j) = Data.subject(i).n(j).x_COP_filtered(Data.subject(i).n(j).t_opto_filtered_index_peak)-Data.subject(i).n(j).x_COP_filtered(Data.subject(i).n(j).t_opto_filtered_index_0);
       Data.subject(i).n(j).heel_deformation_y = Data.subject(i).n(j).y_COP_filtered(Data.subject(i).n(j).t_opto_filtered_index_zd0_COP)-Data.subject(i).n(j).y_COP_filtered(Data.subject(i).n(j).t_opto_filtered_index_0);
       Matrix.heel_deformation_y(i,j) = Data.subject(i).n(j).heel_deformation_y;
       Matrix.heel_deformation_y_peak(i,j) = Data.subject(i).n(j).y_COP_filtered(Data.subject(i).n(j).t_opto_filtered_index_peak)-Data.subject(i).n(j).y_COP_filtered(Data.subject(i).n(j).t_opto_filtered_index_0);
       Data.subject(i).n(j).heel_deformation_z = Data.subject(i).n(j).z_COP_filtered(Data.subject(i).n(j).t_opto_filtered_index_zd0_COP)-Data.subject(i).n(j).z_COP_filtered(Data.subject(i).n(j).t_opto_filtered_index_0);
       Matrix.heel_deformation_z(i,j) = Data.subject(i).n(j).heel_deformation_z;
       Matrix.heel_deformation_z_peak(i,j) = Data.subject(i).n(j).z_COP_filtered(Data.subject(i).n(j).t_opto_filtered_index_peak)-Data.subject(i).n(j).z_COP_filtered(Data.subject(i).n(j).t_opto_filtered_index_0);
       % max deformation
       Data.subject(i).n(j).total_heel_deformation = norm([Data.subject(i).n(j).heel_deformation_x,Data.subject(i).n(j).heel_deformation_y,Data.subject(i).n(j).heel_deformation_z]);
       Matrix.total_heel_deformation(i,j) = Data.subject(i).n(j).total_heel_deformation;
       Matrix.total_heel_deformation_peak(i,j) = norm([Matrix.heel_deformation_x_peak(i,j),Matrix.heel_deformation_y_peak(i,j),Matrix.heel_deformation_z_peak(i,j)]);
       
       % step length
       Matrix.step_length(i,j) = Data.subject(i).n(j).step_length;
       
    end
    % body mass
    Matrix.BodyMass(i) = Data.subject(i).bodyweight/Data.g;
end

%% Optional: uncomment to save structure

% save('AllData.mat','Data','Matrix','-v7.3')


