function [U1rmws U2rmws U1dc U2dc Pact Prea Papp PF] = swm1(U1,U2,rootedWW,WW)

% Input-arguments:
%------------------------
% U1:  Signal-array for U1 (Voltage-channel)
% U2:  Signal array for U2 (Current-channel)
% rootedWW: Array containing the windowing-function

% Function Outout
%------------------------
% U1rmws: RMS of U1 (Voltage-channel)
% U2rmws: RMS of U2 (Current-channel)
% U1dc: RMS of U1 (Voltage-channel)
% U2dc: RMS of U2 (Current-channel)
% Pact:   Calculated Active Power 
% Papp:   Calculated Aparent Power
% Prea:   Calculated Reactive Power   
% PF:     Power factor

% ===== Time-domain SWM =================================
WinSVM1 = U1.*rootedWW;
WinSVM2 = U2.*rootedWW;

U1dc = mean(U1.*WW);%.*(rootedWW)^2); 
U2dc = mean(U2.*WW);%.*rootedWW.*rootedWW);
U1rmws  = sqrt(mean( WinSVM1.*WinSVM1 )); % U1 RMS (rmws)
U2rmws  = sqrt(mean( WinSVM2.*WinSVM2 )); % U2 RMS (rmws)
%[mean( WinSVM1.*WinSVM2)]
Pact    = mean( WinSVM1.*WinSVM2 );       % P active
Papp    = U1rmws * U2rmws;                % P apperent
Prea    = sqrt(max(Papp^2,Pact^2)-Pact^2);% P reactive
PF      = Pact/Papp;                      % Power Factor
% ===== Time-domain  SWM ================================

