%% Phantom Recipe
% Ashley Laughney
% November 3rd, 2010

%% INPUT
% totalVolume = total phantom volume
% desired_uM_Hb  = micromolar concentration of blood (15uM typical breast)
% desired muS' = desired reduced scattering coefficient (inverse mm)
% wavelength_muSp = wavelength given muS' specified at
%% OUTPUT
% Ingredients for mixture A, B, C
% B gives IL volume based on desired muS' using Kienle paper and also based upon empirical
% observation (18g IL per 480 mL volume gives mus'~ 1 inverse)
% 1) Put mixture A into the microwave until boils for 15 to 20 seconds
% 2) Quickly put mixture A on stirplate and add mixture B slowly, monitor temperature
% 3) When the temperature drops to 38d C add mixture C (the blood)
% 4) Let the final mixture stir for ~2min and pour into mold
%%
% function [A B C] = phantom_calculator(totalVolume,desired_uM_Hb, desired_muSp, wavelength_muSp)
function [A B C] = phantom_calculator(V,desired_uM_Hb, hematocrit, IL)

%%%%%%%% USER INPUT%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% V = 480; % mL, total volume
% desired_uM_Hb = 15; % desired uM Hb
% hematocrit = 149.32; % g/L, measure with Hemocue
% desired_muSp = 1.2; % inverse mm


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute mL IL based upon desired mus'
% Parameters for muSp for 20% Intralipid
% constants derived from Michels, Foschum paper (2008)
% lambda = wavelength_muSp; %nm
% y0 = 8.261E+1;
% a = -1.288E-1;
% b = 6.093E-5;
% muSp_calc = y0 + a*lambda + b*lambda.^2; % inverse mm
% xIL = 0.2*desired_muSp/muSp_calc; % percent of IL required for desired muS'
% mL_ILc = xIL*V;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Determine needed volume of blood
uM_Hb = 10^6*(hematocrit/64500);
mL_Hb = desired_uM_Hb*V/uM_Hb;

% Container A
A.mL_PBS = V/2;
A.g_Agarose = V/100;

% Container B
B.mL_IL = IL*V/480;  % based on Mike's experience with 20%IL
B.mL_PBS = (V/2) - (B.mL_IL+mL_Hb);
% B.mL_IL_calc = mL_ILc;
% B.mL_PBS_calc = (V/2)-(B.mL_IL_calc+mL_Hb);

% Container C (syringe)
C.mL_Hb = mL_Hb;

end