% SHELLSIM MAIN FILE Last updated: 12/17/2019

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                            Description
% This file introduces the parameters and state variables that are used to
% simulate Pacific oyster growth. The temperature, chlorophyl, particulate
% organic matter, and particulate organic carbon estimates used in this
% estimate were collected by [INSERT WHO WE SHOULD CREDIT] in [AND WHERE IT
% WAS COLLECTED HERE], Oregon. The model calls the function
% shellsimtransition.m. The function sets up the structure fo the energetic
% model and simulates shell and tissue growth.
%
%                         *State Variables*
% SE: Shell energy (endogenous) TISE: soft tissue energy (endogenous) DSTW:
% Dry soft tissue weight (endogenous) DSW: Dry shell weight (endogenous)
% TEMP: Temperature (exogenous) CHL: Chlorophyl (exogenous) POC:
% Particulate organice carbon (exogenous) POM: Particulate organic
% matter(exogenous) OA: Ocean acidification (exogenous)
% 
%                         *Internal Switches*
% OA: Ocean Acidification SPAWN: weather we all spawning to happen
%
%                         *Other notes*
% The step size in this simulator is a single day. The mesh.horizon
% parameter specifies the number of days for which you want to run the
% simulation.
%
% Reference: Hawkins, A. J. S., Pascoe, P. L., Parry, H., Brinsley, M.,
% Black, K. D.,
%   McGonigle, C., ... & O'Loan, B. (2013). Shellsim: a generic model of
%   growth and environmental effects validated across contrasting habitats
%   in bivalve shellfish. Journal of Shellfish Research, 32(2), 237-253.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Housekeeping
dbstop if error
clear all; close all; clc;

% location = 'D:\Users\klingd\Dropbox';
location = 'C:\Users\vasant\Dropbox\';


% BEGIN SIMULATION Call the file that has the exogenous variables. The
% first column should include the exogenous time series data for TEMP, the
% second for CHL, third for POM, and the fourth for POC. 

% if strcmp(location,'C:\Users\vasant\Dropbox\')
EXOGENOUS = textread(strcat(location,'\oashellfisheconm\timeseries data\exogenoustime.txt'),'','delimiter',' ');
% else
	% EXOGENOUS = textread(strcat(location,'\oashellfisheconm\timeseries
	% data\exogenoustime.txt'),'','delimiter',' ');
% end

TEMP = EXOGENOUS(:,1); 
CHL = EXOGENOUS(:,2);
POM = EXOGENOUS(:,3);
POC = EXOGENOUS(:,4);


% The following section specifies parameters for the model. We use the
% parameter structure here. It is simplified in the transition function
% file. 

% parameter for EREM
p.EREM_a = 0.632;
p.EREM_b = 0.086;
p.EREM_c = 1000;
p.EREM_d = 100;
p.EREM_e = 4.187;
p.phytoplankton = 23.5;
p.EREM_f = 8.25;
p.EREM_g = 21.24;
p.EREM_h = -2.79;
p.EREM_j = -0.174;
p.EREM_k = 20.48;

%parameters for SELORG
p.SELORG_chl_a = 50;
p.SELORG_chl_b = 12;
p.SELORG_a = 0.38; 
p.SELORG_b = 1000;

%parameters for NIRSELORG; there are two sets of parameters what should we
%do here? page 7
p.NIRS_a = -0.33;
p.NIRS_b = 4.11;
p.NIRS_c = -0.62; 

%parameters for NIRREMORG
p.NIRM_a = 8.21;
p.NIRM_b = 0.34;

%parameters for NEA
p.NEA_a = 23.5;
p.NEA_b = 0.15;
p.NEA_c = 0.82;
p.NEA_d = 24;

%parameters for MHL
p.MHL_a = 4.005;
p.MHL_b = -0.72;
p.MHL_c = 24;

%we is the dry soft tisse weight (experimental data, depends on species
%etc. does george have data for this?) ws (sw sometimes) is the standard
%weight (1g;
%https://pdfs.semanticscholar.org/60a8/95b39f9458fe493948940a7d283612cbc14d.pdf)
p.WS = 1 ; %this is definitely just 1
p.WE = .25; % We need experimental data here 


% Clearance rate function parameter
p.CR_a = 0.320;
p.CR_b = 0.323;
p.CR_c = -0.011;
p.CRcost = 15; 

p.CR15 = p.CR_a + p.CR_b*p.CRcost + p.CR_c*p.CRcost^2;

% Chlorophyll threshold for NIRSELORG
p.CHLbarN=0.01;

p.MNEA = 1350;
p.MTA = 0.76;
p.THL_a = 0.23;

%parameters for ON
p.ON_a = 10;
p.ON_b = 200;

% Parameters for THL
p.EL_a = 14.06;
p.EL_b = 16;
p.EL_c = 14;
p.EL_d = 1000;

% Parameters for NEB
p.NEB_a = 0.02428;

%there is a nirselorg discrepency in the code. need to match what is in the
%text to what is in the code table
 
%STATE VARIABLE PARAMETERS
p.DSW_a =1000;
p.DSTW_a = 23.5;
p.DSTW_b = 1000;

%other constants/unit
p.OrgSh = 0.0035 ; %(frac)
p.ECS = p.OrgSh*46; %(J/total dry mg) energy content of shell
p.WCS = 0.189;
p.WCT = 0.914;
p.SCW = 1.115;
p.a = 2.767; % C. gigas SL conversion
p.b = 0.327; % C. gigas SL conversion
p.SLM = 5;
p.TTS = 19;
p.PSTL = 0.44;


% Initial conditions
SE(1)= 0.0056; 
DSW(1) = 0.035;
TISE(1) = 0.0786;
DSTW (1)= 0.005;

SEnext = 0; 
DSWnext = 0;
TISEnext = 0;
DSTWnext = 0;

% Number of steps in simulation
mesh.horizon = 1000;

% State variable vectors that will be populated with each iteration
SE = [SE(1);zeros(mesh.horizon-1,1)]; 
DSTW = [DSTW(1);zeros(mesh.horizon-1,1)];
DSW = [DSW(1);zeros(mesh.horizon-1,1)]; 
TISE = [TISE(1);zeros(mesh.horizon-1,1)]; 



for j=1:mesh.horizon-1
% At the end of each round, the new value of the state function was the
% ending value from the previous step. This is true for each step, except
% the first step. In the first time step, the state variables are
% equivalent to the initial conditions

SE(j)= SEnext;
TISE(j) = TISEnext;
DSW(j) = DSWnext;
DSTW(j) = DSTWnext;
SE(1)= 0.0056; 
DSW(1) = 0.035;
TISE(1) = 0.0786;
DSTW (1)= 0.005;

% Calling the transition function. In doing this we will populate a vector
% for each state variable and we will have the final x(j+1) step value. 
[DSWnext, DSTWnext, SEnext, TISEnext] = shellsimtransition(p, DSW(j), DSTW(j), SE(j), TISE(j), TEMP(j), POM(j), POC(j), CHL(j));

end

