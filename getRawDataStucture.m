%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Load and structure experimental data                                     %
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

% clc
% clear all

%% Trial numbers and heights

Data.f_kp = 3000; % [Hz] sampling frequency FP (Force Plate - Ground Reaction Force Data) 
Data.f_opto = 1293.5; % [Hz] sampling frequency opto (Optotrak - Motion Capture Data) 

Data.subjects = 2:13; % [-] subject 1 pilot data is not included since experimental protocol was different

Data.subject(2).height_knee = 0.485; % [m] approximate height knee joint from ground, purely indicative, measured with tape measure
Data.subject(2).height_hip = 0.88; % [m] approximate height hip joint from ground, purely indicative, measured with tape measure
Data.subject(2).height_navel = 1.035; % [m] approximate height navel from ground, purely indicative, measured with tape measure
Data.subject(2).n_total = 1:34; % [-] total trials numbers
Data.subject(2).n_25 = [2,3]; % [-] trial numbers for force plate calibration data
Data.subject(2).n_bodyweight = [4,5,33,34]; % [-] trial numbers for static measurement on forceplate, usefull for static marker location during stance and subject weight
Data.subject(2).n_normal = [7:11,21:24,26]; % [-] trial numbers for heel strike on forceplate during normal forward walking speed 
Data.subject(2).n_slow = [12,13,16,18,19,28:32]; % [-] EXTRA! trial numbers for heel strike on forceplate during slow forward walking speed, these are NOT used for the publication, however could be valuable if there is an interest in walking speed dependence

Data.subject(3).height_knee = 0.515;
Data.subject(3).height_hip = 0.95;
Data.subject(3).height_navel = 1.09;
Data.subject(3).n_total = [1:9,11:31];
Data.subject(3).n_25 = [30,31];
Data.subject(3).n_bodyweight = [1,2,28,29];
Data.subject(3).n_normal = [3:7,18:22];
Data.subject(3).n_slow = [9,11:13,16,23:27];

Data.subject(4).height_knee = 0.515;
Data.subject(4).height_hip = 0.90;
Data.subject(4).height_navel = 1.02;
Data.subject(4).n_total = 1:30;
Data.subject(4).n_25 = [28,30];
Data.subject(4).n_bodyweight = [1,2,26,27];
Data.subject(4).n_normal = [3,5,6,9,10,16:20];
Data.subject(4).n_slow = [11:15,21:25];

Data.subject(5).height_knee = 0.495;
Data.subject(5).height_hip = 0.865;
Data.subject(5).height_navel = 1.055;
Data.subject(5).n_total = 1:27;
Data.subject(5).n_25 = [26,27];
Data.subject(5).n_bodyweight = [1,2,24,25];
Data.subject(5).n_normal = [4:8,14:18];
Data.subject(5).n_slow = [9:13,19:23];

Data.subject(6).height_knee = 0.515;
Data.subject(6).height_hip = 0.905;
Data.subject(6).height_navel = 1.095;
Data.subject(6).n_total = 1:31;
Data.subject(6).n_25 = [30,31];
Data.subject(6).n_bodyweight = [1,2,28,29]; % 2? maybe not reliable, not static stance?
Data.subject(6).n_normal = [4,6:9,17:21];
Data.subject(6).n_slow = [11:15,23:27];
Data.subject(6).n_check = [5]; % check for effect of treadmill (other experiment in lab) in forceplate data
Data.subject(7).height_knee = 0.5;
Data.subject(7).height_hip = 0.91;
Data.subject(7).height_navel = 1.045;
Data.subject(7).n_total = 1:27;
Data.subject(7).n_25 = [26,27];
Data.subject(7).n_bodyweight = [1,2,24,25];
Data.subject(7).n_normal = [3:7,13:17];
Data.subject(7).n_slow = [8:12,19:23];

Data.subject(8).height_knee = 0.575;
Data.subject(8).height_hip = 1.03;
Data.subject(8).height_navel = 1.22;
Data.subject(8).n_total = 1:28;
Data.subject(8).n_25 = [27,28];
Data.subject(8).n_bodyweight = [1,2,25,26];
Data.subject(8).n_normal = [3:7,13:16,18];
Data.subject(8).n_slow = [8:12,19:21,23,24];

Data.subject(9).height_knee = 0.56;
Data.subject(9).height_hip = 0.965;
Data.subject(9).height_navel = 1.13;
Data.subject(9).n_total = 1:36;
Data.subject(9).n_25 = [35,36];
Data.subject(9).n_bodyweight = [2,3,33,34];
Data.subject(9).n_normal = [4:8,21:25];
Data.subject(9).n_slow = [9,15,16,19,20,26:28,30,32];
Data.subject(9).n_check = [1]; % check for effect of treadmill (other experiment in lab) in forceplate data

Data.subject(10).height_knee = 0.505;
Data.subject(10).height_hip = 0.94;
Data.subject(10).height_navel = 1.08;
Data.subject(10).n_total = 1:26;
Data.subject(10).n_25 = [25,26];
Data.subject(10).n_bodyweight = [1,2,23,24];
Data.subject(10).n_normal = [3,5:7,13:17];
Data.subject(10).n_slow = [8:12,18:22];
% camerabeelden: 3 ~ 8440 - 22 ~ 8459

Data.subject(11).height_knee = 0.495;
Data.subject(11).height_hip = 0.93;
Data.subject(11).height_navel = 0.965; 
Data.subject(11).n_total = 1:27;
Data.subject(11).n_25 = [26,27];
Data.subject(11).n_bodyweight = [1,2,24,25];
Data.subject(11).n_normal = [3:7,13:17];
Data.subject(11).n_slow = [8:12,18:21,23];
% camerabeelden: 4 x 5 tests

Data.subject(12).height_knee = 0.54;
Data.subject(12).height_hip = 0.965;
Data.subject(12).height_navel = 1.14;
Data.subject(12).n_total = 1:30;
Data.subject(12).n_25 = [29,30];
Data.subject(12).n_bodyweight = [1,2,27,28];
Data.subject(12).n_normal = [3:6,8,16:19]; % 1 test missing in second run!
Data.subject(12).n_slow = [9:11,13,14,20,22,24:26];
% camerabeelden: 4 x 5 tests

Data.subject(13).height_knee = 0.465;
Data.subject(13).height_hip = 0.873;
Data.subject(13).height_navel = 0.98;
Data.subject(13).n_total = 1:28;
Data.subject(13).n_25 = [27,28];
Data.subject(13).n_bodyweight = [1,2,25,26];
Data.subject(13).n_normal = [3:7,14,15,17:19];
Data.subject(13).n_slow = [8:12,20:24];
% camerabeelden: 4 x 5 tests

%% load Force plate Data and OptoTrak Data into Data structure

for i = Data.subjects
    for j = Data.subject(i).n_total
        
        % create string with correct number of zeros such that file name always has 6 numbers
        n_string = num2str(j);
        n_zeros = 6-length(n_string);
        zero_string = '';
        for temp = 1:n_zeros;
            zero_string = strcat(zero_string,'0');
        end
        
        % create strings for the correct data file to load into Data structure
        kp_string = strcat(num2str(i),'/TN',zero_string,n_string,'.gai');
        opto_string = strcat(num2str(i),'/TN',zero_string,n_string,'.ndf');
        
        % load raw Data into Data structure
        Data.subject(i).n(j).kp = load(kp_string);
        Data.subject(i).n(j).Vx = sum(Data.subject(i).n(j).kp(:,1:2),2); % sum channels in x direction
        Data.subject(i).n(j).Vy = sum(Data.subject(i).n(j).kp(:,3:4),2); % sum channels in y direction
        Data.subject(i).n(j).Vz = sum(Data.subject(i).n(j).kp(:,5:8),2); % sum channels in z direction
        Data.subject(i).n(j).Trig = Data.subject(i).n(j).kp(:,9);
        [Data.subject(i).n(j).x,Data.subject(i).n(j).y,Data.subject(i).n(j).z] = readndf(opto_string);
        Data.subject(i).n(j).x = Data.subject(i).n(j).x*1E-3;Data.subject(i).n(j).y = Data.subject(i).n(j).y*1E-3;Data.subject(i).n(j).z = Data.subject(i).n(j).z*1E-3;% change into SI units, meters
    end
end

%% Calculate Force gain

% Manual force plate:
% Emfinlichkeitsstreuung
% Fxy = -8 pC/N;
% Fz = -3.8 pC/N;

% Range of amplifyer:
% Range1 = 1000 pC;
% Range2 = 5000 pC;
% Range3 = 10000 pC;
% Range4 = 50000 pC;

% Theoretical gain in Range 2
Data.gain_xy = 5000/8/10; % [N/V]
Data.gain_z = 5000/3.8/10 % [N/V]

% Waardes gemeten met gebruiken verschillende gewichten
% value_znull = mean(Data(?).Vz(1:4000));
% value_z25(1) = mean(Data(1).Vz);
% value_z25(2) = mean(Data(2).Vz);
% value_z25(3) = mean(Data(47).Vz);
% value_z5 = mean(Data(?).Vz);
% value_x5 = mean(Data(?).Vx);
% value_y5 = mean(Data(?).Vy);

% Data.w_5kg = 4.9868*9.81;
Data.w_25kg = 25.01*9.81; % Calibration weight used during test

% gain of Forceplate in z-direction
% gain_z25 = mean(w_25kg./value_z25)
% gain_z5 = w_5kg/value_z5
% gain_x5 = w_5kg/value_x5
% gain_y5 = w_5kg/value_y5
% 
% gain_xy =
% 
%    62.5000
% 
% 
% gain_z =
% 
%   131.5789
% 
% 

gain_z25 =  133.1067

%% Optional: uncomment to save structure

% save('RawData.mat','Data')


