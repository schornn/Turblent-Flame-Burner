clear
clc
close all
%Raw Data to Velocities and Constants
% Adds Data and code path
addpath 'C:\Lvm Data Files - HWA';
addpath 'C:\Lvm Data Files - HWA\10000 HTI 3D';

fg_tiff_list = dir(strcat('C:\Lvm Data Files - HWA\10000 HTI 3D\**.lvm')); %Generates alphabetic list of data files
 
%Sorts through multiple .lvm files in a folder
for j=1:length(fg_tiff_list)
    file_name{j} = fg_tiff_list(j).name;
    rawdata= lvm_import(file_name{j}); %imports raw lvm data into a string from single file 
    ts_m(j,:)= rawdata.Segment1.data.^2; %Voltage mean average squared
end

%Curve fits from laminar flow, Converting Bridge Voltage to Effective Velocity
Qa_h= ((ts_m(7,:)+7.1395)/5.219).^(1/.5013);
Qb_h= ((ts_m(9,:)+4.1908)/3.0664).^(1/.5856);
Qc_h= ((ts_m(11,:)+2.4597)/3.8922).^(1/.4623);
Qa_k= ((ts_m(8,:)+7.4743)/5.5094).^(1/.4893);
Qb_k= ((ts_m(10,:)+3.6202)/2.7321).^(1/.6147);
Qc_k= ((ts_m(12,:)+6.3567)/4.6927).^(1/.5194);
Ve1= ((ts_m(1,:)+7.4743)/5.5094).^(1/.4893);
Ve2= ((ts_m(2,:)+2.8685)/2.2023).^(1/.4151);
Ve3= ((ts_m(3,:)+7.1395)/5.219).^(1/.5013);
Ve4= ((ts_m(4,:)+6.3241)/4.5928).^(1/.5283);
Ve5= ((ts_m(5,:)+8.2607)/6.067).^(1/.4814);
Ve6= ((ts_m(6,:)+7.7387)/5.6603).^(1/.4985);
Qa_h_bar= mean(Qa_h);
Qb_h_bar= mean(Qb_h);
Qc_h_bar= mean(Qc_h);
Qa_k_bar= mean(Qa_k);
Qb_k_bar= mean(Qb_k);
Qc_k_bar= mean(Qc_k);
Ve1_bar= mean((mean(Ve1)-Ve1).^2);
Ve2_bar= mean((mean(Ve2)-Ve2).^2);
Ve3_bar= mean((mean(Ve3)-Ve3).^2);
Ve4_bar= mean((mean(Ve4)-Ve4).^2);
Ve5_bar= mean((mean(Ve5)-Ve5).^2);
Ve6_bar= mean((mean(Ve6)-Ve6).^2);

%% Inputs
vbare= [Ve1_bar Ve2_bar Ve3_bar Ve4_bar Ve5_bar Ve6_bar]; %Six time averaged flucating effective velocity terms taken at the 6 different postions 

%% Constants
gamma= [0 180 90 270 45 315]; %The six postions at different roll angles with respect to the probe
alpha= 45; %Angle of the normal vector to the flow
h= .5*((Qa_h_bar+Qc_h_bar)/Qb_h_bar)^2-1; 
k= .5*((Qa_k_bar+Qc_k_bar)/Qb_k_bar)^2-1; 
A1= cosd(alpha)^2+k^2*sind(alpha)^2;
B1= (A1)^.5;
C1= A1;
D1= A1;
D2= sind(alpha)^2+k^2*cosd(alpha)^2;
D3= (1-k^2)*sind(2*alpha);
D4= D3^2/(4*D1);
D5= (D2-D4)/(4*D1);
D6= h^2/(4*D1);

%% Equations for solving for Mean Air Flow Velocitie components based on measurments
Vbar= (mean(Ve1)-mean(Ve2))/(D3/B1); %Mean Velocity in the Y direction with probe as refrence, with mesaruements takes at postions 1 and 2
Wbar= (mean(Ve4)-mean(Ve3))/(D3/B1); %Mean Veloctiy in the Z diretion with probe as refrence, with mesaruements takes at postions 3 and 4
syms x
Ubarsym = solve(x^2-x*1/(2*B1)*(mean(Ve1)+mean(Ve2)) + D5*Vbar^2+D6*Wbar^2 == 0, 'PrincipalValue', false); %Solving the mean Velocity in the X direction while ignoring turblence effects
Ubarvalue= double(Ubarsym);
if Ubarvalue(1) > abs(Vbar) && Ubarvalue(1) > abs(Wbar)
    Ubar= Ubarvalue(1);
elseif Ubarvalue(2) > abs(Vbar) && Ubarvalue(2) > abs(Wbar)
    Ubar= Ubarvalue(2);
end

%%Main Equation
R= 1;
Count= 1;
while R > .01

for i=1:length(gamma)
    % Coefficients
    A2= (sind(alpha)^2+k^2*cosd(alpha)^2)*cosd(gamma(1,i))^2+h^2*sind(gamma(1,i))^2;
    A3= (sind(alpha)^2+k^2*cosd(alpha)^2)*sind(gamma(1,i))^2+h^2*cosd(gamma(1,i))^2;
    A4= (1-k^2)*cosd(gamma(1,i))*sind(2*alpha);
    A5= -(1-k^2)*sind(gamma(1,i))*sind(2*alpha);
    A6= -(sind(alpha)^2+k^2*cosd(alpha)^2-h^2)*sind(2*gamma(1,i));
    B2= A4/(2*B1);
    B3= A5/(2*B1);
    B4= (A2-B2^2)/(2*B1); 
    B5= (A3-B3^2)/(2*B1);
    B6= (A6-2*B2*B3)/(2*B1);
    C2= B2^2;
    C3= B3^2;
    C4= A4;
    C5= A5;
    C6= 2*B2*B3;
    C7=2*B1*B4;
    C8= 2*B1*B6;
    C9= 2*B1*B5;
    C10= 2*B2*B4;
    C11= 2*B2*B6+4*B3*B4;
    C12= 2*B2*B6;
    C13= 2*B3*B6;
    C14= 2*B3*B6 + 4*B2*B5;
    C15= 2*B3*B5;
    % ubar^2 vbar^2 wbar^2 uvbar uwbar uwbar
    C(i,:)= [C1 (C2+2*C10*Vbar/Ubar+C12*Wbar/Ubar) (C3+C13*Vbar/Ubar+2*C15*Wbar/Ubar) (C4+2*C7*Vbar/Ubar+C8*Wbar/Ubar) (C5 + C8*Vbar/Ubar+2*C9*Wbar/Ubar) (C6+C11*Vbar/Ubar+C14*Wbar/Ubar)];
end

u= inv(C)*vbare'; % Final values of ubar^2 vbar^2 wbar^2 uvbar uwbar uwbar
Ubar_corrected= 1/(2*B1)*(mean(Ve1)+mean(Ve2))-D5*((Vbar^2+u(2,1))/Ubar)-D6*((Wbar^2+u(3,1))/Ubar);
R= abs(Ubar-Ubar_corrected);
Ubar= Ubar_corrected;
Count= Count+1;
end





