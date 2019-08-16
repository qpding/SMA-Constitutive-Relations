% Liang and Rogers Model
%-----------------------------------------------------%
% -- Introduction:
% Numerical simulation of constitutive relations for one-dimentional SMA
% based on Liang and Rogers model
%
% -- Reference:
% One-Dimensional Thermomechanical Constitutive Relations for Shape Memory Materials
%
% -- Time:
% Aug. 15th, 2019
%-----------------------------------------------------%

close all;
clear;
clc;
FontSize = 12;

% extract parameters
[ coeffDic, TDic, RDic ] = loadParameters('L_and_R_Article.xml');
D     = coeffDic.('YoungsModulus');
THETA = coeffDic.('ThermoelasticTensor');
OMEGA = coeffDic.('TransformationTensor');
M_f   = TDic.('MartensiteFinish');
M_s   = TDic.('MartensiteStart');
A_s   = TDic.('AusteniteStart');
A_f   = TDic.('AusteniteFinish');
C_A   = RDic.('C_A');
C_M   = RDic.('C_M');
epsilon_L = -OMEGA/D;
a_A       = pi/(A_f - A_s);
a_M       = pi/(M_s - M_f);
b_A       = -a_A/C_A;
b_M       = -a_M/C_M;

% FREE RECOVERY
% Figure 19. Recovery strain vs. temperature of free recovery
[ curveCounts, N ] = deal( 3, 1000 );
epsilon_r          = zeros(N, curveCounts);
% initial condition (residual strain)
epsilon_res        = linspace(epsilon_L*0.1, epsilon_L, curveCounts);
T                  = linspace(-26, -6, N)';
f = @(t)1.*(t<=A_s) + cos(a_A*(t-A_s)).*(A_s<t & t<A_f) + (-1).*(t >= A_f);
for i = 1:curveCounts
    epsilon_r(:, i) = epsilon_res(:, i)...
                      - ( THETA*(T-A_s) + OMEGA/2 * epsilon_res(:, i)/epsilon_L * (f(T)-1) ) / D;
end

figure(1);
FigureHandle = subplot(2,2,1);
hold on;
box on;
p1 = plot(T, epsilon_r(:, 1), '-',  'Color', [0/255 115/255 174/255], 'LineWidth', 1.5);
p2 = plot(T, epsilon_r(:, 2), '--', 'Color', [115/255 0/255 174/255], 'LineWidth', 1.5);
p3 = plot(T, epsilon_r(:, 3), ':',  'Color', [174/255 115/255 0/255], 'LineWidth', 1.5);
legend([p3 p2 p1], {['\epsilon _{res} = ' num2str(epsilon_res(:,3)*100) '%'], ...
                    ['\epsilon _{res} = ' num2str(epsilon_res(:,2)*100) '%'], ...
                    ['\epsilon _{res} = ' num2str(epsilon_res(:,1)*100) '%']}, ...
       'Box', 'off', ...
       'Orientation', 'vertical', ...
       'Location', 'northeastoutside');
text(0, 0, 'Heating');%-30, 50, 
xlabel('Temperature (Deg. C)','FontName','Times New Roman','FontSize',FontSize);
ylabel('Recovery Strain','FontName','Times New Roman','FontSize',FontSize);
title('Figure 19. Recovery strain vs. temperature of free recovery', 'FontName','Times New Roman','FontSize',FontSize);

% RESTRAINED RECOVERY
% Heating: Figure 16. Recovery stress vs. temperature of restrained SMA wire
% Cooling: Figure 17. Predicted hysteresis of a copper based SMA with 0.2% initial strain
[ curveCounts, N ] = deal( 3, 1000 );
sigma_r_h          = zeros(N, curveCounts);
sigma_r_c          = zeros(N, curveCounts);
epsilon_res        = [0.2, 0.6, 1.0] * 0.01;
T                  = linspace(-90, 50, N)';
options            = optimoptions('fsolve','Display','none');

for i = 1:curveCounts
       % Heating
       xi_0       = epsilon_res(:, i) / epsilon_L;
       T_M        = M_f + acos(2*(xi_0-0.5))/a_M;
       A_s_m      = (C_A*A_s-THETA*T_M) / (C_A-THETA);
       sigma_As_r = THETA * (A_s_m-T_M);
       A_f_m      = (a_A*A_s - b_A*sigma_As_r + b_A*OMEGA*xi_0 + b_A*THETA*A_s_m + pi) ...
                  / (a_A + b_A*THETA);
       sigma_Af_r = THETA * (A_f_m-A_s_m) - OMEGA*xi_0 + sigma_As_r;
       ce = @(s, t)THETA*(t-A_s_m) + OMEGA/2*xi_0*(cos(a_A*(t-A_s)+b_A*s)-1) ...
                                   + sigma_As_r - s;

%        r = @(t) fsolve(@(s) ce(s, t), [sigma_As_r sigma_Af_r]);
%        relation = @(t)(t<=T_M)             .* 0 ...
%                      +(t>T_M & t<=A_s_m)   .* (THETA*(t-T_M)) ...
%                      +(t>A_s_m & t<=A_f_m) .* r(t) ...
%                      +(t>A_f_m)            .* (THETA*(t-A_f_m) + sigma_Af_r);
%     sigma_r_h(:, i) = relation(T);
    for j = 1:N
       if T(j)<=T_M
              sigma_r_h(j, i) = 0;
       elseif T(j)>T_M && T(j)<=A_s_m
              sigma_r_h(j, i) = (THETA*(T(j)-T_M));
       elseif T(j)>A_s_m && T(j)<=A_f_m
              sigma_r_h(j, i) = fsolve(@(s) ce(s, T(j)), sigma_r_h(j-1, i), options);
       else
              sigma_r_h(j, i) = (THETA*(T(j)-A_f_m) + sigma_Af_r);
       end  
    end

    % Cooling from temperature above A_f_m
       sigma_c_r  = sigma_r_h(N, i);
       xi_c       = 0;
       T_c        = T(N);
       M_s_m      = (a_M*M_f - b_M*sigma_c_r + b_M*THETA*T_c + pi) / (a_M + b_M*THETA);
       sigma_Ms_r = sigma_c_r + THETA * (M_s_m-T_c);
       M_f_m      = (a_M*M_f - b_M*sigma_Ms_r - b_M*OMEGA*(1-xi_c) + b_M*THETA*M_s_m) ...
                     / (a_M + b_M*THETA);
       sigma_Mf_r = THETA * (M_f_m-M_s_m) + OMEGA*(1-xi_c) + sigma_Ms_r;
       ce = @(s, t)THETA*(t-M_s_m) + OMEGA/2*(1-xi_c)*(cos(a_M*(t-M_f)+b_M*s)+1) ...
                                   + sigma_Ms_r - s;

    for j = 1:N
       if T(j)<=M_f_m
              sigma_r_c(j, i) = THETA*(T(j)-M_f_m) + OMEGA*(1-xi_c) + sigma_Mf_r;
       elseif T(j)>M_f_m && T(j)<=M_s_m
              sigma_r_c(j, i) = fsolve(@(s) ce(s, T(j)), sigma_r_c(j, i), options);
       else
              sigma_r_c(j, i) = (THETA*(T(j)-A_f_m) + sigma_Af_r);
       end  
    end
end

subplot(2, 2, 2);
hold on;
box on;
p1 = plot(T, sigma_r_h(:, 1), '-',  'Color', [0/255 115/255 174/255], 'LineWidth', 1.5);
p2 = plot(T, sigma_r_h(:, 2), '--', 'Color', [115/255 0/255 174/255], 'LineWidth', 1.5);
p3 = plot(T, sigma_r_h(:, 3), ':',  'Color', [174/255 115/255 0/255], 'LineWidth', 1.5);
% p4 = plot(T, sigma_r_c(:, 1), '-',  'Color', [0/255 115/255 174/255], 'LineWidth', 1.5);
p5 = plot(T, sigma_r_c(:, 2), '--', 'Color', [115/255 0/255 174/255], 'LineWidth', 1.5);
% p6 = plot(T, sigma_r_c(:, 3), ':',  'Color', [174/255 115/255 0/255], 'LineWidth', 1.5);
legend([p3 p2 p1], {['\epsilon _{res} = ' num2str(epsilon_res(:,3)*100) '%'], ...
                    ['\epsilon _{res} = ' num2str(epsilon_res(:,2)*100) '%'], ...
                    ['\epsilon _{res} = ' num2str(epsilon_res(:,1)*100) '%']}, ...
       'Box', 'off', ...
       'Orientation', 'vertical', ...
       'Location', 'northeastoutside');
xlabel('Temperature (Deg. C)','FontName','Times New Roman','FontSize',FontSize);
ylabel('Recovery Stress (MPa)','FontName','Times New Roman','FontSize',FontSize);
title('Figure 16. Recovery stress vs. temperature of restrained SMA wire', 'FontName','Times New Roman','FontSize',FontSize);


residualStrainIndex = 2;

for i = residualStrainIndex

end

set(gcf, 'Position', [9,49,1300,600]);

%--------------------End of code----------------------%
%-----------------------------------------------------%
% text('Position', [-10, 6], 'String', ['A_s = ' num2str(A_s) ', A_f = ' num2str(A_f)]);
% xlim([2 3.7]);
% ylim([-2.5 0.5]);
% set(gca,'XTick',[]);
% set(gca,'YTick',[]);
