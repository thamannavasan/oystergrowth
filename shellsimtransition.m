% SHELLSIM TRANSITION FILE
% Last updated: 12/17/2019 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The outputs of this function are shell energy, soft tissue, dry shell
% weight, and dry soft tissue weight in the next period. The inputs are
% shell energy and soft tissue energy. These endogenous state variables are
% determined by a number of exogenous parameters and exogenous variables.
% Together TISE and SE will give us the total fresh weight.
%
%
%                   *State Variables*
% SE: Shell energy (endogenous) 
% TISE: soft tissue energy (endogenous) 
% DSTW: Dry soft tissue weight (endogenous) 
% DSW: Dry shell weight (endogenous)
% TEMP: Temperature (exogenous) 
% CHL: Chlorophyl (exogenous) 
% POC: Particulate organice carbon (exogenous) 
% POM: Particulate organic matter(exogenous) 
% OA: Ocean acidification (exogenous)
% 
%                   *Internal Switches*
% OA: Ocean Acidification 
% SPAWN: weather we all spawning to happen
%
%                   *Other notes and variabels of note*
% TG: energy allocation to total soft tissues
% SG: energy allocation to shell growth
% Note that SE and TISE are functions of the exogenous parameters and of SE
% and TISE at the start of the period. However, DSW and DSTW are also
% functions of SE and TISE as well as DSW and DSTW. 
% 
% Reference:
% Hawkins, A. J. S., Pascoe, P. L., Parry, H., Brinsley, M., Black, K. D.,
%   McGonigle, C., ... & O'Loan, B. (2013). Shellsim: a generic model of
%   growth and environmental effects validated across contrasting habitats
%   in bivalve shellfish. Journal of Shellfish Research, 32(2), 237-253.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [DSWnext, DSTWnext, SEnext, TISEnext] = shellsimtransition(p, DSW, DSTW, SE, TISE,TEMP, POM, POC, CHL,DSWnext, DSTWnext, SEnext, TISEnext);

%%%%%%% Unpacking parameters from the "p" structure

% parameter for erem
EREM_a = p.EREM_a;
EREM_b = p.EREM_b;
EREM_c = p.EREM_c;
EREM_d = p.EREM_d;
EREM_e = p.EREM_e;
phytoplankton = p.phytoplankton;
EREM_f = p.EREM_f;
EREM_g = p.EREM_g;
EREM_h = p.EREM_h;
EREM_j = p.EREM_j;
EREM_k = p.EREM_k;


%parameters for selorg 
SELORG_chl_a = p.SELORG_chl_a ;
SELORG_chl_b = p.SELORG_chl_b ;
SELORG_a = p.SELORG_a;
SELORG_b = p.SELORG_b;


%parameters for NIRSELORG; there are two sets of parameters what should we do here? page 7
NIRS_a = p.NIRS_a;
NIRS_b = p.NIRS_b;
NIRS_c = p.NIRS_c; 

%parameters for NIRREMORG
NIRM_a = p.NIRM_a;
NIRM_b = p.NIRM_b;

%parameters for NEA 
NEA_a = p.NEA_a;
NEA_b = p.NEA_b;
NEA_c  = p.NEA_c;
NEA_d = p.NEA_d;

%parameters for MHL
MHL_a  = p.MHL_a;
MHL_b = p.MHL_b;
MHL_c = p.MHL_c;

%we is the dry soft tisse weight (experimental data, depends on species etc. does george have data for this?)
%ws (sw sometimes) is the standard weight (1g; https://pdfs.semanticscholar.org/60a8/95b39f9458fe493948940a7d283612cbc14d.pdf)
WS = p.WS; %this is definitely just 1
WE = p.WE; % We need experimental data here 

%clearance rate function parameter
CR_a = p.CR_a;
CR_b = p.CR_b;
CR_c = p.CR_c;
CRcost = p.CRcost; 

CR15 = CR_a + CR_b*CRcost + CR_c*CRcost^2;


% Chlorophyll threshold for NIRSELORG
CHLbarN = p.CHLbarN;

MNEA = p.MNEA;
MTA = p.MTA;
THL_a = p.THL_a;

%parameters for ON
ON_a = p.ON_a;
ON_b = p.ON_b;

%parameters for THL 
EL_a = p.EL_a;
EL_b = p.EL_b;
EL_c = p.EL_c;
EL_d = p.EL_d;

%parameter for NEB 
NEB_a = p.NEB_a;

%STATE VARIABLE PARAMETERS
DSW_a = p.DSW_a;
DSTW_a = p.DSTW_a ;
DSTW_b = p.DSTW_b ;

%unpack other parameters 
OrgSh = p.OrgSh;
ECS = p.ECS ;
WCS = p.WCS;
WCT = p.WCT;
SCW = p.SCW ;
a = p.a ;
b = p.b ;
SLM = p.SLM;
TTS = p.TTS;
PSTL = p.PSTL;

%Here are the options for OA and spawning. Turn options.SPAWN = 0 if triploid, options.OA = 1 when we want to add OA effects to the model
options.SPAWN = 1;
options.OA = 0;

% the following functions make use of the exogenous variable at time step
% j. The model then runs through a process where the energy is processed by
% the organism, it then grows as a result. 

	if CHL >0 & POC == 0 & POM== 0
	 SELORG = ((CHL/SELORG_b)*SELORG_chl_a)/SELORG_a;
		elseif (CHL >0 & POC == 0 & POM > 0 | CHL >0 & POC > 0 & POM > 0)
	 SELORG = ((CHL/SELORG_b)*SELORG_chl_b)/SELORG_a;
		else
	 SELORG = 0;
	end
 
	REMORG = POM-SELORG; %if we don't have POM and POC data this will be negative
	
	if (CHL >0 & POC > 0 & POM > 0) %no option for no pom or poc data?
	 EREM = ((POM.*((EREM_a + (EREM_b*(POC./(POM*EREM_c))*EREM_d))*EREM_e)) - (SELORG*phytoplankton))./REMORG;
		elseif (CHL >0 & POC == 0 & POM > 0)
	 EREM = (EREM_f + (EREM_g*(1 - exp(EREM_h*SELORG))))+(EREM_j*REMORG);
		elseif (CHL >0 & POC > 0 & POM == 0)
	 EREM = EREM_k;
		else
	 EREM = 0;
    end
    
% clearance rate (this is only for one of the species...is this the correct
% one?)
	CR = CR_a+CR_b*TEMP+CR_c*TEMP.^2;
% Temperature effect on feeding (temp is a function of this so we define it
% here)m
	TEF = CR/CR15;

% Temperature effects on maintenance heat losses; 
	TEM = exp(CR_a*TEMP)/exp(CR_a*CRcost);
	
	if (CHL <0.01)
	 NIRSELORG = 0; %discrepency between code and text
		else
		NIRSELORG = (NIRS_b*SELORG).*TEF.*(WE./WS).^NIRS_c;
	end
	
	NIRREMORG = NIRM_a.*(1-exp(-NIRM_b.*REMORG)).*TEF*(WE/WS)^NIRS_c;
	
	NEA = ((NIRSELORG.*NEA_a)+(NIRREMORG*NEA_b*EREM))*NEA_c*NEA_d;
	
	MHL = MHL_a.*TEF.*(WE/WS).^(MHL_b).*MHL_c;
	
	THL = MHL+THL_a*NEA;
	
	ON = ON_a+((ON_b-ON_a)./MNEA).*NEA;
	
	EL = (((THL./EL_a)/EL_b)/ON).*EL_c.*EL_d;
	
	NEB = NEA - THL -(EL.*NEB_a);
    
    COND = TISE./(TISE+SE);
	
	if (COND>=MTA)
		SG=(1-MTA).*NEB;
	else
		SG=0;
	end
	
	if (COND<MTA)
		TG = NEB;
    elseif (COND>=MTA)
		TG = MTA.*NEB;
	else
        TG = 0;
    end
    
	SL = a*DSW.^b;
    

	if (SL>=SLM & TEMP >= TTS & COND>=0.95*MTA)
		if (options.SPAWN == 1)
			SPAWN = 23.5*PSTL*DSTW;
		else
			SPAWN = 0;
		end
	end


%	if 
%		if (options.OA ==1)
%		else
%			OA = 0;
%		end
%	end



 
% at the end, the resulting ending value of each state variable is a
% function of the state variable in the previous step and total growth (TG)
	SEnext = SE + SG;
	TISEnext = TISE  + (TG);
	DSWnext = DSW+(SE./(DSW_a*ECS));
	DSTWnext = DSTW + (TISE/(DSTW_a*DSTW_b));
    
    

	
end
