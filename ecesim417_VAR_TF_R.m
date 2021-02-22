close all; clear all; clc;
%         Given constants
% k_b = 1.38e-23;   %Boltzmann constant (J/K)
k_b = 8.62e-5;      %Boltzmann constant (eV/K) this has the q
q = 1.602e-19;      %Electron charge (C)
T = 300;            %Room Temperature (K)
e_0 = 8.854e-14;    %Permittivity of vacuum (F/cm)
e_si = 11.8*e_0;    %Permittivity of Si (F/cm)
n_i = 1.0e10;       %Intrinsic carrier concentration (cm^-3)
W_e = 2.5e-5;       %Emitter thickness (cm)
W_b = 2.5e-5;       %Base thickness (cm)
W_c = 5e-6;         %Collector thickness (cm)
%         Additional constants
Vth = 0.0259;       %Thermal Voltage (V)
Thk_e = 0.25e-4;    %Emitter thickness (cm)
Thk_b = 0.25e-4;    %Base thickness (cm)
Thk_c = 0.5e-4;     %Collector thickness (cm)
A_e = Thk_e^2;      %BJT surface area (cm^2)
A_b = Thk_b^2;      %Base surface area (cm^2)
A_c = Thk_c^2;      %Collector surface area (cm^2)
t_0 = 5e-7;         %Initial minority carrier lifetime (s)
N_0 = 5e16;         %Doping density (cm^-3)
Va_e = 0;           %Emitter applied bias (V)
Va_b = -0.75;       %Base applied bias (V)
Va_c = -1.5;        %Collector applied bias (V)
V_ce = -2:0.001:0;  %Collector-Emitter voltage sweep (V)
                    
% Part I                    
% Finding I_s Transport Saturation Current (both junctions R)
%         Concentrations
N_e = 1e19;         %Emitter concentration (cm^-3)
N_b = 1e17;         %Base concentration (cm^-3)
N_c = 2e16;         %Collector concentration (cm^-3)
%         Additional constants2
ni2 = (1e10)^2;     %intrinsic carrier concentration Si
% n_e0 = ni2/N_e;   %minority electrons in emitter
p_b0 = ni2/N_b;     %minority holes in base
n_c0 = ni2/N_c;     %minority electrons in collector
V_be = [0, -0.66, -0.68, -0.70, -0.72, -0.74]; 
% Bandgap narrowing
delE_g = 0.009 * (log(N_e/1e17)+sqrt( (log(N_e/1e17))^2+0.5));
ni_prime = sqrt(ni2*exp(delE_g / (k_b*T)));
% for the emitter base
n_e0 = (ni_prime^2) / N_e; %minority electrons in emitter, adjusted for BGN

%         Base-Emitter voltage sweep (V)
%         Built-in Voltages
V_JE = Vth*log((N_e*N_b)/(ni_prime^2)); %(eV)
%         Emitter built-in voltage
V_JC = Vth*log((N_c*N_b)/ni2); %(eV)
%         Collector built-in voltage

%         Mobilities
Mup_e = 54.3 + (407 / (1 + N_e / 2.67e17)); %(V/s)
% majority hole in emitter
mun_e = 88 + (1252 / (1 + N_e / 1.432e17)); %(V/s)
% minority electron in emitter
Mun_b = 88 + (1252 / (1 + N_b / 1.432e17)); %(V/s)
% majority electron in base
mup_b = 54.3 + (407 / (1 + N_b / 2.67e17)); %(V/s)
% minority hole in base
Mup_c = 54.3 + (407 / (1 + N_c / 2.67e17)); %(V/s)
% majority hole in collector
mun_c = 88 + (1252 / (1 + N_c / 1.432e17)); %(V/s)
% minority electron in collector

%         Einstein relation for Diffusivities
D_e = Mup_e * Vth; %(cm^2/s)
% diffusivity of majority holes in emitter
d_e = mun_e * Vth; %(cm^2/s)
% df of min electrons in emitter
D_b = Mun_b * Vth; %(cm^2/s)
% df of maj electrons in base
d_b = mup_b * Vth; %(cm^2/s)
% df of min hole in base
D_c = Mup_c * Vth; %(cm^2/s)
% df of maj hole in collector
d_c = mun_c * Vth; %(cm^2/s)
% df of min electron in collector

%         Minority carrier lifetimes
t_e = t_0 / (1 + N_e / N_0); %(s) 
t_b = t_0 / (1 + N_b / N_0); %(s)
t_c = t_0 / (1 + N_c / N_0); %(s)
%         Diffusion length 
L_e = sqrt(D_e * t_e); 
% df length of majority holes in emitter
l_e = sqrt(d_e * t_e);
% df length of min electrons in emitter
L_b = sqrt(D_b * t_b);
% df length of maj electrons in base
l_b = sqrt(d_b * t_b); 
% df length of min hole in base
L_c = sqrt(D_c * t_c);
% df length of maj hole in collector
l_c = sqrt(d_c * t_c);
% df length of min electron in collector


% Base width modulation
% measurements for effective base width
x_ne = sqrt( ((2*e_si)/q)* (N_e / (N_b*(N_e+N_b)))*(V_JE - Va_e));
x_nc = sqrt( ((2*e_si)/q)* (N_c / (N_b*(N_c+N_b)))*(V_JC - Va_c));
W_eff = W_b - (x_ne + x_nc);

x_pe = sqrt( ((2*e_si)/q)* (N_b / (N_e*(N_e+N_b)))*(V_JE - Va_e));

x_pc = sqrt( ((2*e_si)/q)* (N_b / (N_c*(N_c+N_b)))*(V_JC - Va_c));

% lengths for resistances
L_e = Thk_e - x_pe;
L_b = W_eff;
L_c = Thk_c - x_pc;

%         Collector current 
I_s = q*(ni_prime^2)*A_e*(D_b / (N_b * W_eff));
% Jessica used this one


% I_s = (q*A_c) * (((D_c/L_c) * N_c) + ((D_b/L_b) * N_b)) * (exp((q*Va_c/k_b*T))-1);
%         Forward and Reverse gain
%         Constants for a
a_21 = (q*A_b)*(((d_b * p_b0) / l_b)*(1/(sinh(W_eff / l_b))));
% (q*A_b)*(((d_b*p_b0)/l_b)*(1/(sinh(W_eff/l_b))));
    a_12 = a_21;
a_11 = (q*A_b)*(((d_e*n_e0)/l_e)+(((d_b*p_b0)/l_b)*(coth(W_eff/l_b))));
a_22 = (q*A_b)*(((d_c*n_c0)/l_c)+((d_b*p_b0)/l_b)*(coth(W_eff/l_b)));
B_f = a_21 / (a_11 - a_21);
B_r = a_12 / (a_22 - a_12);

% I_s = q*A_e*((D_b*ni_prime^2) / (L_b * N_c) + (D_c * ni2) / (L_c*N_b));
Ia_c = I_s*exp(Va_e/(k_b*T)) - (I_s*(1 + (1/B_r))*exp(Va_c/(k_b*T)))
% check with the beta values 
% I_s = q*A_c*( (D_b*ni2) / (L_b * p_b0) + (D_c * ni2) / (L_c * n_c0))
% I_cbias1 = I_s*exp( (Va_e / k_b *T) - 1)
I_cbias = I_s * exp((q*Va_c / (k_b*T)) -1)

% I_cbias = (a_21*(exp((-Va_e/(k_b*T)))-1)) - a_22*(exp((-Va_c/(k_b*T)))-1)

% Part II
for i = 1:6
    for j = 1:length(V_ce);
V_bc(j) = V_be(i) - V_ce(j);
I_c(i,j) = (a_21*(exp((-1*V_be(i)/(k_b*T)))-1)) - a_22*(exp((-1*V_bc(j)/(k_b*T)))-1);
    end
end
% plot(V_ce, I_c);
% legend('VBE=0', 'VBE=-0.66', 'VBE=-0.68', 'VBE=-0.70', 'VBE=-0.72', 'VBE=-0.74');

% Part III
% Parasitic resistances (assuming bulk dominated)

%         Resistivities
rho_e = 1 / (q * ((Mup_e * N_e) + (mun_e * n_e0))); %(Ohms)

rho_b = 1 / (q * ((Mun_b * N_b) + (mup_b * p_b0))); %(Ohms)

rho_c = 1 / (q * ((Mup_c * N_c) + (mun_c * n_c0))); %(Ohms)
%         Resistances
R_e = (rho_e * L_e) / A_e
R_b = (rho_b * L_b) / A_b
R_c = (rho_c * L_c) / A_c
% Part IV
% Emitter and Collector Capacitances

% C_JE = (e_si*A_e) / x_ne;
% C_JC = (e_si*A_c) / x_nc;
C_JE0 = sqrt((e_si*q*N_b)/(2*V_JE));
C_JC0 = sqrt((e_si*q*N_b)/(2*V_JC));
% 
% C_JE0 = sqrt((e_si*q*N_b)/(2*V_JE))*A_e;
% C_JC0 = sqrt((e_si*q*N_b)/(2*V_JC))*A_c;
% 
C_JE = C_JE0 / sqrt( 1 - (Va_e / V_JE));
C_JC = C_JC0 / sqrt( 1 - (Va_c / V_JC));

% Part V
% Transit times
% TF = ((Vth / I_cbias) * C_JC) + ((W_eff^2) / (2 * D_c))
% TR = B_f*(TF / B_r);

TF = ((Vth * C_JC) / I_cbias) + (W_eff^2 / (2*D_c));
TR = B_f*(TF / B_r);
% TR = ((Vth * C_JE) / I_cbias) + (W_eff^2 / (2*D_e));
% Part IV
% Early voltages 
VAF = (q*N_b*W_eff*A_c)/C_JC0;
% Forward early voltage (CB)
VAR = (q*N_b*W_eff*A_e)/C_JE0;
% Reverse early voltage (BE)

% VAF = (q*N_b*W_eff)/C_JC
% % Forward early voltage (CB)
% VAR = (q*N_b*W_eff)/C_JE
% % Reverse early voltage (BE)
% Part IIV
% Knee currents
v_sat = 1e7; %(cm^2)
I_kf = q*N_c*v_sat*A_c;
% Forward knee current
I_kr = q*N_e*v_sat*A_e;
% Reverse knee current
% Part IIIV
% Small-signal model parameters
g_m = I_cbias / Vth;
% transconductance
g_pi = g_m / B_f;
% input conductance
g_o = Vth * g_m;
% output conductance
g_mu = g_o / B_f;
% reverse feedback conductance