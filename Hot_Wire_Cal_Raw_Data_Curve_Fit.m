close all
clear
clc

%Postion 2
% Adds Data and code path
addpath 'C:\Lvm Data Files - HWA';
addpath 'C:\Lvm Data Files - HWA\71550003 SY hc';


fg_tiff_list = dir(strcat('C:\Lvm Data Files - HWA\71550003 SY hc\**.lvm')); %Generates alphabetic list of data files
 
%Sorts through multiple .lvm files in a folder
for j=1:length(fg_tiff_list)
    file_name{j} = fg_tiff_list(j).name;
    rawdata= lvm_import(file_name{j}); %imports raw lvm data into a string from single file 
    delta_Pa(j,1)= strread(file_name{j})/10000; %reads the change of pressure from the file name
    ts_ma(j,1)= mean((rawdata.Segment1.data).^2); %Voltage mean average squared
end

%Calculating theroitcal laminar velocity
rho= 1.225; %kg/m^3
ua= sqrt(2*delta_Pa*1000/rho); %Bernoulli's Equaiton
%Curves fits data and theroitcal values, King's Law Fit
fa= fit(ua,ts_ma,'power2');
ca= coeffvalues(fa);
Qa= ((ts_ma-ca(3))/ca(1)).^(1/ca(2)); %Velocity Calibration from Buresti and Talamelli (1991)

figure(1)
plot(ts_ma,Qa)
hold on
scatter(ts_ma,ua)

% 
% %% Postion 4
% % Adds Data and code path
% 
% addpath 'C:\Lvm Data Files - HWA\71550003 SY kb';
% 
% 
% fg_tiff_list = dir(strcat('C:\Lvm Data Files - HWA\71550003 SY kb\**.lvm')); %Generates alphabetic list of data files
%  
% %Sorts through multiple .lvm files in a folder
% for j=1:length(fg_tiff_list)
%     file_name{j} = fg_tiff_list(j).name;
%     rawdata= lvm_import(file_name{j}); %imports raw lvm data into a string from single file 
%     delta_Pb(j,1)= strread(file_name{j})/10000; %reads the change of pressure from the file name
%     ts_mb(j,1)= mean((rawdata.Segment1.data).^2); %Voltage mean average squared
% end
% 
% %Calculating theroitcal laminar velocity
% rho= 1.225; %kg/m^3
% ub= sqrt(2*delta_Pb*1000/rho); %Bernoulli's Equaiton
% %Curves fits data and theroitcal values, King's Law Fit
% fb= fit(ub,ts_mb,'power2');
% cb= coeffvalues(fb);
% Qb= ((ts_mb-cb(3))/cb(1)).^(1/cb(2)); %Velocity Calibration from Buresti and Talamelli (1991)
% 
% figure(2)
% plot(ts_mb,Qb)
% hold on
% scatter(ts_mb,ub)

% 
% %% Postion 5
% % Adds Data and code path
% 
% addpath 'C:\Lvm Data Files - HWA\71550003 SY hc';
% 
% 
% fg_tiff_list = dir(strcat('C:\Lvm Data Files - HWA\71550003 SY hc\**.lvm')); %Generates alphabetic list of data files
%  
% %Sorts through multiple .lvm files in a folder
% for j=1:length(fg_tiff_list)
%     file_name{j} = fg_tiff_list(j).name;
%     rawdata= lvm_import(file_name{j}); %imports raw lvm data into a string from single file 
%     delta_Pc(j,1)= strread(file_name{j})/10000; %reads the change of pressure from the file name
%     ts_mc(j,1)= mean((rawdata.Segment1.data).^2); %Voltage mean average squared
% end
% 
% %Calculating theroitcal laminar velocity
% rho= 1.225; %kg/m^3
% uc= sqrt(2*delta_Pc*1000/rho); %Bernoulli's Equaiton
% %Curves fits data and theroitcal values, King's Law Fit
% fc= fit(uc,ts_mc,'power2');
% cc= coeffvalues(fc);
% Qc= ((ts_mc-cc(3))/cc(1)).^(1/cc(2)); %Velocity Calibration from Buresti and Talamelli (1991)
% 
% figure(3)
% plot(ts_mc,Qc)
% hold on
% scatter(ts_mc,uc)

% %% Postion 6
% % Adds Data and code path
% 
% addpath 'C:\Lvm Data Files - HWA\Cal_Postion_6_+45';
% 
% 
% fg_tiff_list = dir(strcat('C:\Lvm Data Files - HWA\Cal_Postion_6_+45\**.lvm')); %Generates alphabetic list of data files
%  
% %Sorts through multiple .lvm files in a folder
% for j=1:length(fg_tiff_list)
%     file_name{j} = fg_tiff_list(j).name;
%     rawdata= lvm_import(file_name{j}); %imports raw lvm data into a string from single file 
%     delta_Pd(j,1)= strread(file_name{j})/100; %reads the change of pressure from the file name
%     ts_md(j,1)= mean((rawdata.Segment1.data).^2); %Voltage mean average squared
% end
% 
% %Calculating theroitcal laminar velocity
% rho= 1.225; %kg/m^3
% ud= sqrt(2*delta_Pd*1000/rho); %Bernoulli's Equaiton
% %Curves fits data and theroitcal values, King's Law Fit
% fd= fit(ud,ts_md,'power2');
% cd= coeffvalues(fd);
% Qd= ((ts_md-cd(3))/cd(1)).^(1/cd(2)); %Velocity Calibration from Buresti and Talamelli (1991)
% 
% figure(3)
% plot(ts_md,Qd)
% hold on
% scatter(ts_md,ud)
