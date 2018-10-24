% Mattia Poinelli
% JPL, June 2018
% Summary of LEFM results and elaboration of features' characteristics table

clear all, close all, clc
format short
%% LENGTH
run length_crack.m

Length_km = LENGTH';
NAME = struct2cell(S);
NUMBER =1:26; NUMBER = NUMBER';
clear w A a_ang angular_distance b_ang Bulk_Modulus C d DONE e E END epsilon eta F g G H h_d H_p h_s i l_d l_s lambda LENGTH Mass mu n Name Norm_ve nu q R rho_a rho_i rho_s rho_w ROT S T_europa T_ns T_ns_years Tang_ve TOU TS a

%% DEPTH & HEIGHT
Max_depth  = [105 75 75 80 70 50 95 100 70 35 80 45 20 85 80 95 90 60 70 75 10 100 5 80 60 70]';
Max_height = [1280 935 920 980 845 620 1130 1235 830 485 965 590 365 1025 980 1145 1085 740 860 920 305 1250 245 980 740 875]';

%% HOR PROP
run horizontal_propagation.m
clear A a_ang angular_distance b b_ang Bulk_Modulus c C check COMP  Cycles_to_build_filtered d d_Normal_stress d_Shear_stress d_sigma_phi d_sigma_theta d_STRAIN_R d_T d_T_R d_tau Delta depth T dF dV e E END eta F g G H h_d H_p h_s i j K_I l_d l_s lambda Lambda LEGEND Mass Mean_depth mu n N_S_feat Norm_ve nu  phi phi_rad q R raster rho_a rho_i rho_s rho_w ROT S Shear_stress sigma_phi sigma_theta START STRAIN_R STRING T_europa T_ns T_ns_years T_R Tang_ve tau theta theta_rad time TIME TOU TS V w Z

%% TABLE WRITING
Max_op_width_m= Max_opening_width';
Max_prop_rate_m_s = (Max_prop_rate)';
Mean_prop_rate_m_s = Mean_prop_rate'; Mean_prop_rate_m_hr = Mean_prop_rate_m_s.*3600;
Max_op_rate_mm_yrs = Max_opening_rate.*1000.*(60*60*24*365);Max_op_rate_mm_yrs = Max_op_rate_mm_yrs';
Cycles_to_devel = ceil(Cycles_to_build)';
Perc_to_devel = Perc_completion';

Max_op_width_m     (Max_opening_width < 0) = 0;
Max_op_rate_mm_yrs(Max_opening_width < 0) = 0;
Max_op_rate_mm_yrs(Max_op_rate_mm_yrs < 0) = 0;
Max_prop_rate_m_s (Max_opening_width < 0) = 0;
Mean_prop_rate_m_s (Max_opening_width < 0) = 0;
Cycles_to_devel    (Max_opening_width < 0) = 0;
Perc_to_devel      (Max_opening_width < 0) = 0;

T = table(Length_km,Max_depth,Max_height,Max_op_width_m,Max_op_rate_mm_yrs,Max_prop_rate_m_s,Mean_prop_rate_m_hr,Cycles_to_devel,Perc_to_devel,'RowNames',NAME(5,:))
%writetable(T,'Europa_results.xls','WriteRowNames',true)
