function PlotResult(LineWidth,State,mem,sec,DOF,varargin)

set(0, 'DefaultLineLineWidth',LineWidth);

if nargin==5
    OSAnalysis = false;
    legendStr = {"MATLAB"};
elseif nargin>5
    OSAnalysis = true;
    DispTop = varargin{1};
    ForceCol = varargin{2};
    DefoColSec = varargin{3};
    ForceColSec = varargin{4};
    legendStr = {"MATLAB","OS"};
end

% GLOBAL RESPONSE
Qk = squeeze(State.Qk(DOF,mem,:));
qk = squeeze(State.qk(DOF,mem,:));
figure(1); hold on;
plot(qk,Qk);
if OSAnalysis
    plot(DispTop(:,1),ForceCol(:,4));
end
xlabel('q (m)'); ylabel('Q (N)'); title(['Q-q plot Member = ',num2str(mem),' DOF = ',num2str(DOF)]);
grid on;
legend(legendStr)

% SECTION RESPONSE
Ssk1 = squeeze(State.Ssk{mem}(1,sec,:));
vsk1 = squeeze(State.vsk{mem}(1,sec,:));
figure(2); hold on;
plot(vsk1,Ssk1);
if OSAnalysis
    plot(DefoColSec(:,1),ForceColSec(:,1));
end
xlabel('{\epsilon}'); ylabel('P (N)');
title(['P-{\epsilon} plot Member = ',num2str(mem),' Section = ',num2str(sec)]);
grid on;
legend(legendStr);

Ssk2 = squeeze(State.Ssk{mem}(2,sec,:));
vsk2 = squeeze(State.vsk{mem}(2,sec,:));
figure(3); hold on;
plot(vsk2,Ssk2);
if OSAnalysis
    plot(DefoColSec(:,2),ForceColSec(:,2));
end
xlabel('{\psi} (1/m)'); ylabel('M (N-m)');
title(['M-{\psi} plot Member = ',num2str(mem),' Section = ',num2str(sec)]);
grid on;
legend(legendStr);