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

hold on;
box on;
p1 = plot(T, epsilon_r(:, 1), '-',  'Color', [0/255 115/255 174/255], 'LineWidth', 1.5);
p2 = plot(T, epsilon_r(:, 2), '--', 'Color', [115/255 0/255 174/255], 'LineWidth', 1.5);
p3 = plot(T, epsilon_r(:, 3), ':',  'Color', [174/255 115/255 0/255], 'LineWidth', 1.5);
legend([p1 p2 p3], {['\epsilon _{res} = ' num2str(epsilon_res(:,1))*100 '%'], ...
                    ['\epsilon _{res} = ' num2str(epsilon_res(:,2))*100 '%'], ...
                    ['\epsilon _{res} = ' num2str(epsilon_res(:,3))*100 '%']}, ...
       'Box', 'off');
% text('Position', [-10, 6], 'String', ['A_s = ' num2str(A_s) ', A_f = ' num2str(A_f)]);
xlabel('Temperature (Deg. C)','FontName','Times New Roman','FontSize',10);
ylabel('Recovery Strain','FontName','Times New Roman','FontSize',10);
% title('Figure 19. Recovery strain vs. temperature of free recovery', 'position', [-15, -4],'FontName','Times New Roman','FontSize',10);

% RESTRAINED RECOVERY
% Figure 16. Recovery stress vs. temperature of restrained SMA wire
[ curveCounts, N ] = deal( 3, 1000 );
sigma_r            = zeros(N, curveCounts);
epsilon_res        = [0.2, 0.6, 1.0] * 0.01;
T                  = linspace(-26, -6, N)';

for i = 1:curveCounts
       xi_0       = epsilon_res(:, i) / epsilon_L;
       T_M        = M_f + acos(2*(xi_0-0.5))/a_M;
       A_s_m      = (C_A*A_s-THETA*T_M) / (C_A-THETA);
       sigma_As_r = THETA * (A_s_m-T_M);
       A_f_m      = (a_A*A_s - b_A*sigma_As_r + b_A*OMEGA*xi_0 + b_A*THETA*A_s_m + pi) ...
                  / (a_A + b_A*THETA);
       sigma_Af_r = THETA * (A_f_m-A_s_m) - OMEGA*epsilon_res(:, i)/epsilon_L + sigma_As_r;
       ce = @(s, t)THETA*(t-T_M) + OMEGA/2*epsilon_res(:, i)/epsilon_L ...
                                   *(cos(a_A*(t-A_s)+b_A*s)-1) ...
                                 + sigma_As_r - s;
       
       figure(2);
       tempVec = zeros(N, 1);
       sigmaVec = linspace(sigma_As_r, sigma_Af_r, N);
       for temp = 1:N
              tempVec(temp, 1) = THETA*(-20-T_M) + OMEGA/2*epsilon_res(:, i)/epsilon_L ...
              *(cos(a_A*(-20-A_s)+b_A*sigmaVec(temp))-1) ...
            + sigma_As_r - temp;
              plot(sigmaVec, tempVec);
            hold on;
       end

       fTest = @(s) ce(s, -20);
%        r = @(t) fsolve(@(s) ce(s, t), [sigma_As_r sigma_Af_r]);
%        relation = @(t)(t<=T_M)             .* 0 ...
%                      +(t>T_M & t<=A_s_m)   .* (THETA*(t-T_M)) ...
%                      +(t>A_s_m & t<=A_f_m) .* r(t) ...
%                      +(t>A_f_m)            .* (THETA*(t-A_f_m) + sigma_Af_r);
%     sigma_r(:, i) = relation(T);
    for j = 1:N
       if T(j)<=T_M
              sigma_r(j, i) = 0;
       elseif T(j)>T_M && T(j)<=A_s_m
              sigma_r(j, i) = (THETA*(T(j)-T_M));
       elseif T(j)>A_s_m && T(j)<=A_f_m
              % sigma_r(j, i) = fsolve(@(s) ce(s, T(j)), (sigma_As_r+sigma_Af_r)/2);
       else
              sigma_r(j, i) = (THETA*(T(j)-A_f_m) + sigma_Af_r);
       end  
    end
end

% text('Position', [-10, 6], 'String', ['A_s = ' num2str(A_s) ', A_f = ' num2str(A_f)]);
% xlim([2 3.7]);
% ylim([-2.5 0.5]);
% xlabel('Temperature (Deg. C)','FontName','Times New Roman','FontSize',10);
% ylabel('Recovery Strain','FontName','Times New Roman','FontSize',10);
% title('Figure 19. Recovery strain vs. temperature of free recovery', 'position', [-15, -4],'FontName','Times New Roman','FontSize',10);
% set(gcf,'unit','centimeters','position',[1,2,8.8,6])
% set(gca,'XTick',[]);
% set(gca,'YTick',[]);
