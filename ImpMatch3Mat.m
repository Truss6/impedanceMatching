%%  Shock Φysics Impedence Matching
clc;clear;close all;figure();TRformat;axis square;  %   init & format plot
%%   Preset Mats
%   from meyers
Cu=Hug(8930,3940,1.489,'Copper',[0.875,0.5,0.25]);
PMMA=Hug(1190,2600,1.52,'PMMA',[0.75,0.5,.75]);
Concrete=Hug2(2340,551,4.52,2235,1.745,'Concrete',[0.75,0.75,0.25]);
Al=Hug(2700,5220,1.37,'Al',[0.7,0.75,0.7]);
Fe=Hug(7850,3570,1.92,'Fe',[0.75,0.5,0.25]);
PE=Hug(820,2900,1.48,'PE',[0.25,0.25,0.5]);
Void=Hug(0,0,0,'Void',[0.1,0.1,0.1]);
U=Hug(18950,2490,2.20,'U',[0.5,0.5,0.25]);
W=Hug(19220,4030,1.24,'W',[0.25,0.25,0.251]);
% Sugar=Hug(2740,5528,1.343,'Sugar',[0.5,1,0.3]);
Sugar=Hug(1580,3040,2.05,'Sugar',[0.5,1,0.3]);
SiO2=Hug(2200,4000,0.01,'Quartz',[0.7,0.7,1]);
WC=Hug(15560,4930,1.39,'WC',[0.8,0.4,0.4]);
Ta=Hug(16690,3414,1.2,'Ta',[0.2,0.4,0.5]);
LiF=Hug(2638,5150,1.35,'LiF',[0.75,0.75,1]);
AlO=Hug(3895,11190,0.97,'Alumina',[0.5,0.5,0.5]);
Nylon=Hug(1140,2290,1.63,'Nylon',[0.5,0.5,0.5]);
Water=Hugf(998,@(u) 1483+10999*log(1+u/5190),'Water',[0.2,0.3,1]); %from forbes
%%  INIT program & Hugoniot params
% Set flyer target, and window materials. Set Backer=Void if no window is
% used
Flyer = PMMA;
Target = Cu;
Backer = Void;

Flyer_Velocity=600; % Set Flyer Velocity in m/s

u0=0;   %   init target velocity (keep 0)

n_solves=5; % set number of reverb solves
%%  Repeted Solver
for i=1:n_solves
    %   Solve Fe-Al interaction
    [u1,P1]=hugSolver(Flyer,Target,Flyer_Velocity,u0);
        %   solve Flyer@Flyer_V & Target@u0
    u(2*i-1)=u1;    P(2*i-1)=P1;    r(2*i-1)=1/Target.v(u1-u0);
        %   Store req'd data
    
    u_prime=u0+2*(u1-u0);   
        %calculate base of Right-hand Al P-u hugoniot
    
    [u2,P2]=hugSolver(Target,Backer,u_prime,0);
        %   sove Target@u_prime & Backer@0
    u(2*i)=u2;      P(2*i)=P2;    r(2*i)=1/Target.v(u_prime-u2);
        %   "
    
    u0=2*u2-u_prime;    
        %calculate base of Left-hand Al P-u hugoniot
end
ylim([min(P)*(min(P)<0),1.2*max(P)]/1e9);
% xlim([min([0,u]),max(u)]*1.20);
% xlim([0,500]);
legend([Flyer.name,' Flyer'],[Target.name,' Target']);%,'','',[Backer.name,' Backer']);
xlabel('Particle Velocity [u_p] (m/s)');
ylabel('Pressure [P] (GPa)');
title("Impedence Matching");
%   Output data in spec'd units
up=u/1e0
Pressure=P/1e9
% Density=r/1e0
% if false
%     figure(2);TRformat;hold on;axis square;
%     plot((Target.r0)./Density,Pressure,'ko');
%     title('PV Hugoniot');
%     xlabel('Volumetric Strain')
% end
%%  hugsolver
function [u_out,P_out]=hugSolver(Mat1,Mat2,u1,u2)
    %   hug solver will solve RH: Mat1@u1 & LH: Mat2&u2
    P1=@(u) Mat1.P(u1-u);   %   RH P-u
    P2=@(u) Mat2.P(u-u2);   %   LH P-u
    
    %   Plot P-u hugs
    fplot(@(u) P1(u)/1e9,range([u2,u1]),'color',Mat1.color,'linewidth',2)
    fplot(@(u) P2(u)/1e9,range([u2,u1]),'color',Mat2.color,'linewidth',2)
    
    %   make a guess for solver init
    u_guess=(u1+u2)/2;
    
    %   solve for u
    u_out=fzero(@(u) P1(u)-P2(u),u_guess);
    
    %   calc P
    P_out=P1(u_out);
    
    %   plot match
    plot(u_out,P_out/1e9,'ko');
end
function [m]=range(u)
    m=[min(u),max(u)*1.00];  %   !!!
end