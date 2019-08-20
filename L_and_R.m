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

%---------------Code for debug control----------------%
% userDEBUG = true;
timerDEBUG = true;
%----------------End of debug control-----------------%

FontSize   = 12;
figureRows = 2;
figureCols = 3;
figureIndex= 1;
[ curveCounts, N ] = deal( 3, 1000 );

if exist('timerDEBUG', 'var')
    tic
end

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

% STRESS AND STRAIN RELATIONS
% Figure 14. Prediction of stress-strain relation for a copper based SMA.
sigma   = linspace(0, 30, N)';
epsilon = zeros(N, 2);
T = [-20 -27];

sigma_linFunc = @(t) C_M * (t - M_s);
sigma_endOfM  = @(t) C_M * (t - M_f);
xi_A  = 0;
for i = 1:2
    for j = 1:N
       if sigma(j) <= sigma_linFunc(T(i))
           epsilon(j, i) = sigma(j) / D;
       elseif sigma(j) <= sigma_endOfM(T(i))
           xi = (1-xi_A)/2 * cos(a_M*(T(i)-M_f)+b_M*sigma(j)) + (1+xi_A)/2;
           epsilon(j, i) = (sigma(j) - OMEGA*xi) / D;
       else
           epsilon(j, i) = (sigma(j) - OMEGA*1) / D;
       end
    end
end

figure(1);
subplot(figureRows, figureCols, figureIndex);
figureIndex  = figureIndex + 1;
hold on;
box on;
p1 = plot(epsilon(:, 1), sigma, '-',  'Color', [0/255 115/255 174/255], 'LineWidth', 1.5);
p2 = plot(epsilon(:, 2), sigma, '--', 'Color', [115/255 0/255 174/255], 'LineWidth', 1.5);
legend([p1 p2], {['T = ' num2str(T(1))], ...
                 ['T = ' num2str(T(2))]}, ...
       'Box', 'off', ...
       'Orientation', 'vertical', ...
       'Location', 'northwest');
xlabel('Strain','FontName','Times New Roman','FontSize',FontSize);
ylabel('Stress (MPa)','FontName','Times New Roman','FontSize',FontSize);
title({'Figure 14. Prediction of stress-strain', 'relation for a copper based SMA.'}, ...
      'FontName', 'Times New Roman','FontSize',FontSize);

% FREE RECOVERY
% Figure 19. Recovery strain vs. temperature of free recovery
epsilon_r          = zeros(N, curveCounts);
% initial condition (residual strain)
epsilon_res        = linspace(epsilon_L*0.1, epsilon_L, curveCounts);
T                  = linspace(-26, -6, N)';
f = @(t)1.*(t<=A_s) + cos(a_A*(t-A_s)).*(A_s<t & t<A_f) + (-1).*(t >= A_f);
for i = 1:curveCounts
    epsilon_r(:, i) = epsilon_res(:, i)...
                      - ( THETA*(T-A_s) + OMEGA/2 * epsilon_res(:, i)/epsilon_L * (f(T)-1) ) / D;
end

FigureHandle = subplot(figureRows,figureCols, figureIndex);
figureIndex  = figureIndex + 1;
hold on;
box on;
p1 = plot(T, epsilon_r(:, 1), '-',  'Color', [0/255 115/255 174/255], 'LineWidth', 1.5);
p2 = plot(T, epsilon_r(:, 2), '--', 'Color', [115/255 0/255 174/255], 'LineWidth', 1.5);
p3 = plot(T, epsilon_r(:, 3), ':',  'Color', [174/255 115/255 0/255], 'LineWidth', 1.5);
legend([p3 p2 p1], {['\epsilon_{res} = ' num2str(epsilon_res(:,3)*100) '%'], ...
                    ['\epsilon_{res} = ' num2str(epsilon_res(:,2)*100) '%'], ...
                    ['\epsilon_{res} = ' num2str(epsilon_res(:,1)*100) '%']}, ...
       'Box', 'off', ...
       'Orientation', 'vertical', ...
       'Location', 'northeast');
text(-10, 4e-3, 'A_{s}=-25');
xlabel('Temperature (Deg. C)','FontName','Times New Roman','FontSize',FontSize);
ylabel('Recovery Strain','FontName','Times New Roman','FontSize',FontSize);
title({'Figure 19. Recovery strain vs. temperature','of free recovery'}, ...
       'FontName','Times New Roman','FontSize',FontSize);

% RESTRAINED RECOVERY
% Heating: Figure 16. Recovery stress vs. temperature of restrained SMA wire
% Cooling: Figure 17. Predicted hysteresis of a copper based SMA with 0.2% initial strain
[ curveCounts, N ] = deal( 3, 1000 );
sigma_r_h          = zeros(N, curveCounts);
sigma_r_c          = zeros(N, curveCounts);
epsilon_res        = [0.2, 0.6, 1.0] * 0.01;
T                  = linspace(-90, 50, N)';
options            = optimoptions('fsolve','Display','none');
L                  = zeros(1, curveCounts);

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
       T_c        = T(end);
       M_s_m      = (a_M*M_f - b_M*sigma_c_r + b_M*THETA*T_c + pi) / (a_M + b_M*THETA);
       sigma_Ms_r = sigma_c_r + THETA * (M_s_m-T_c);
       M_f_m      = (a_M*M_f - b_M*sigma_Ms_r - b_M*OMEGA*(1-xi_c) + b_M*THETA*M_s_m) ...
                     / (a_M + b_M*THETA);
       sigma_Mf_r = THETA * (M_f_m-M_s_m) + OMEGA*(1-xi_c) + sigma_Ms_r;
       ce = @(s, t)THETA*(t-M_s_m) + OMEGA/2*(1-xi_c)*(cos(a_M*(t-M_f)+b_M*s)+1) ...
                                   + sigma_Ms_r - s;

% Code for debug
    if exist('userDEBUG', 'var')
       sigmaFuncValue = zeros(N, 1);
       sigmaTemp      = linspace(0, 30, N);
       for index = 1:N
              sigmaFuncValue(index) = ce(sigmaTemp(index), T(floor(N/2)));
       end
       % subplot(2, 2, 3);
       % plot(sigmaTemp, sigmaFuncValue);
    end
% End of debug

% Code for debug
    if exist('userDEBUG', 'var')
       temperatureDic  = struct('M_f', M_f, 'M_s', M_s, 'A_s', A_s, 'A_f', A_f);
       temperatureMDic = struct('M_f_m', M_f_m, 'M_s_m', M_s_m, 'A_s_m', A_s_m, 'A_f_m', A_f_m);
       display(temperatureDic);
       display(temperatureMDic);
    end
% End of debug

    % Note: calculation of \sigma_r_c should be in a decremental direction,
    %       because 'ce' has multiple roots in [0, sigma_Ms_r]
    for j = N:-1:1
       if T(j) < M_s
              break;
       end
       if T(j)<=M_f_m
           sigma_r_c(j, i) = THETA*(T(j)-M_f_m) + OMEGA*(1-xi_c) + sigma_Mf_r;
       elseif T(j)>M_f_m && T(j)<=M_s_m
           sigma_r_c(j, i) = fsolve(@(s) ce(s, T(j)), sigma_r_c(j+1, i), options);
       else
           sigma_r_c(j, i) = (THETA*(T(j)-A_f_m) + sigma_Af_r);
       end
    end
    L(i) = j;
end

subplot(figureRows,figureCols, figureIndex);
figureIndex  = figureIndex + 1;
hold on;
box on;
p1 = plot(T, sigma_r_h(:, 1), '-',  'Color', [0/255 115/255 174/255], 'LineWidth', 1.5);
p2 = plot(T, sigma_r_h(:, 2), '--', 'Color', [115/255 0/255 174/255], 'LineWidth', 1.5);
p3 = plot(T, sigma_r_h(:, 3), ':',  'Color', [174/255 115/255 0/255], 'LineWidth', 1.5);
% p4 = plot(T(L(1):end), sigma_r_c(L(1):end, 1), '-',  'Color', [0/255 115/255 174/255], 'LineWidth', 1.5);
p5 = plot(T(L(2):end), sigma_r_c(L(2):end, 2), '--', 'Color', [115/255 0/255 174/255], 'LineWidth', 1.5);
% p6 = plot(T(L(3):end), sigma_r_c(L(3):end, 3), ':',  'Color', [174/255 115/255 0/255], 'LineWidth', 1.5);
legend([p3 p2 p1], {['\epsilon_{res} = ' num2str(epsilon_res(:,3)*100) '%'], ...
                    ['\epsilon_{res} = ' num2str(epsilon_res(:,2)*100) '%'], ...
                    ['\epsilon_{res} = ' num2str(epsilon_res(:,1)*100) '%']}, ...
       'Box', 'off', ...
       'Orientation', 'vertical', ...
       'Location', 'northwest');
xlabel('Temperature (Deg. C)','FontName','Times New Roman','FontSize',FontSize);
ylabel('Recovery Stress (MPa)','FontName','Times New Roman','FontSize',FontSize);
title({'Figure 16. Recovery stress vs. temperature', 'of restrained SMA wire'}, ...
       'FontName','Times New Roman','FontSize',FontSize);
ylim([-10 90]);

% Martensite fraction
% Figure 5. Martensite fraction vs. temperature
% M to A transformation, i.e. heating
T          = linspace(M_f-10, A_f+10, N)';
xi_M = 1;
xiFunc_M2A = @(t, s) (t<=A_s)        .* xi_M + ...
                     (t>A_s & t<A_f) .* xi_M / 2 ...
                                     .* (cos(a_A*(t-A_s)+b_A*s) + 1) + ...
                     (t>=A_f)        .* 0;
xi_M2A = xiFunc_M2A(T, 0);

% A to M transformation, i.e. cooling
xi_A = 0;
xiFunc_A2M = @(t, s) (t<=M_f)        .* 1.0 + ...
                     (t>M_f & t<M_s) .* ((1-xi_A) / 2 ...
                                     .*  cos(a_M*(t-M_f)+b_M*s) + (1+xi_A)/2) + ...
                     (t>=M_s)        .* xi_A;
xi_A2M = xiFunc_A2M(T, 0);

subplot(figureRows,figureCols, figureIndex);
figureIndex  = figureIndex + 1;
hold on;
box on;
p1 = plot(T, xi_M2A, '--', 'Color', [0/255 115/255 174/255], 'LineWidth', 1.5);
p2 = plot(T, xi_A2M, '-',  'Color', [115/255 0/255 174/255], 'LineWidth', 1.5);
legend([p1 p2], {'heating', 'cooling'}, ...
       'Box', 'off', ...
       'Orientation', 'vertical', ...
       'Location', 'northeast');
xlabel('Temperature (Deg. C)','FontName','Times New Roman','FontSize',FontSize);
ylabel('Martensite Fraction','FontName','Times New Roman','FontSize',FontSize);
title('Figure 5. Martensite fraction vs. temperature', 'FontName','Times New Roman','FontSize',FontSize);
ylim([-0.2 1.2]);

% Figure 15. Martensite fraction vs. stress, A to M
% A to M transformation, i.e. cooling
sigma  = linspace(0, 30, N)';
xi_A2M = zeros(N, 2);
sigma_linFunc = @(t) C_M * (t - M_s);
sigma_endOfM  = @(t) C_M * (t - M_f);
xi_A  = 0;
xiFunc_A2M = @(t, s) (s<=sigma_linFunc(t))                    .* xi_A + ...
                     (s>sigma_linFunc(t) & s<sigma_endOfM(t)) .* ((1-xi_A) / 2 ...
                                                              .*  cos(a_M*(t-M_f)+b_M*s) + (1+xi_A)/2) + ...
                     (s>=sigma_endOfM(t))                     .* 1;

xi_A2M(:, 1) = xiFunc_A2M(-27, sigma);
xi_A2M(:, 2) = xiFunc_A2M(-20, sigma);

subplot(figureRows,figureCols, figureIndex);
figureIndex  = figureIndex + 1;
hold on;
box on;
p1 = plot(sigma, xi_A2M(:, 1), '--', 'Color', [0/255 115/255 174/255], 'LineWidth', 1.5);
p2 = plot(sigma, xi_A2M(:, 2), '-',  'Color', [115/255 0/255 174/255], 'LineWidth', 1.5);
legend([p1 p2], {'T = -27', 'T = -20'}, ...
       'Box', 'off', ...
       'Orientation', 'vertical', ...
       'Location', 'southeast');
xlabel('Stress (MPa)','FontName','Times New Roman','FontSize',FontSize);
ylabel('Martensite Fraction','FontName','Times New Roman','FontSize',FontSize);
title('Figure 15. Martensite fraction vs. stress', 'FontName','Times New Roman','FontSize',FontSize);
ylim([-0.2 1.2]);

% Figure 18. \xi vs. \sigma^{r} and T for a restrained SMA wire with 0.2% initial strain
% M to A transformation, i.e. heating (recovery)
T      = linspace(-90, 50, N)';
sigma  = linspace(0, 30, N)';
xi_M2A = zeros(N, 1);
sigma_linFunc = @(t) C_A * (t - A_f);
sigma_endOfA  = @(t) C_A * (t - A_s);
xi_M  = 1;
xiFunc_M2A = @(t, s) (s<=sigma_linFunc(t))                    .* xi_M + ...
                     (s>sigma_linFunc(t) & s<sigma_endOfA(t)) .* (xi_M / 2 ...
                                                              .*  (cos(a_A*(t-A_s)+b_A*s) + 1)) + ...
                     (s>=sigma_endOfA(t))                     .* 0;

% xi_M2A(:, 1) = xiFunc_M2A(-27, sigma);
% xi_M2A(:, 2) = xiFunc_M2A(-20, sigma);
for i = 1:N
    xi_M2A(i, 1) = xiFunc_M2A(T(i), sigma);
end

subplot(figureRows,figureCols, figureIndex);
figureIndex  = figureIndex + 1;
hold on;
box on;
p1 = plot3(T, sigma, xi_M2A, '-', 'Color', [0/255 115/255 174/255], 'LineWidth', 1.5);
% legend([p1 p2], {'T = -27', 'T = -20'}, ...
%        'Box', 'off', ...
%        'Orientation', 'vertical', ...
%        'Location', 'southeast');
xlabel('Stress (MPa)','FontName','Times New Roman','FontSize',FontSize);
ylabel('Martensite Fraction','FontName','Times New Roman','FontSize',FontSize);
title({'Figure 18. \xi vs. \sigma^{r} and T', 'for a restrained SMA wire with 0.2% initial strain'}, ...
       'FontName','Times New Roman','FontSize',FontSize);
% ylim([-0.2 1.2]);

set(gcf, 'Position', [9,49,1300,600]); 

if exist('timerDEBUG', 'var')
    toc
end

%--------------------End of code----------------------%
%-----------------------------------------------------%
% text('Position', [-10, 6], 'String', ['A_s = ' num2str(A_s) ', A_f = ' num2str(A_f)]);
% xlim([2 3.7]);
% ylim([-2.5 0.5]);
% set(gca,'XTick',[]);
% set(gca,'YTick',[]);


