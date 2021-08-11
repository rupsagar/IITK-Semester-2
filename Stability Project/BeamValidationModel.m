clear; format shortg;

% MODEL DATA
NumNodeDOF = 3;
ElementType = 'Frame';
STol = 1e-16;
IterMax = 10;
LoadNodeID = 11;
LoadPattern = [LoadNodeID,0,-1,0];

X = [zeros(1,9),12,24,40,56,72,88,104,120];
Y = [(0:15:120),120*ones(1,8)];
NumNode = numel(X);
XY = [X;Y];
XY = XY(:);

MemCon = [1:(NumNode-1);2:NumNode]';
NumMem = size(MemCon,1);
E = 7060.8*ones(NumMem,1);
A = 6*ones(NumMem,1);
I = 2*ones(NumMem,1);

StrBC = false(3,NumNode);
StrBC(:,[1,NumNode]) = [1,1;1,1;0,0];


% ANALYSIS
% DispHis = 0:-2:-61;
% CtrlNodeID = LoadNodeID;
% CtrlDOFID = 2;
% NumStep = numel(DispHis);
% [XYk,Re] = DCM(NumNodeDOF,ElementType,STol,IterMax,XY,MemCon,StrBC,LoadPattern,NumStep,E,A,I,...
%     LoadNodeID,CtrlDOFID,DispHis);

% NumStep = 250;
% psi = 0;
% ds = 2;
% [XYk,Re] = ALCM(NumNodeDOF,ElementType,STol,IterMax,XY,MemCon,StrBC,LoadPattern,NumStep,E,A,I,psi,ds);

NumStep = 250;
dLambdabar = 1;
[XYk,Re] = GDCM(NumNodeDOF,ElementType,STol,IterMax,XY,MemCon,StrBC,LoadPattern,NumStep,E,A,I,dLambdabar);


% PLOTTING RESULTS
set(0, 'DefaultLineLineWidth',1.0);

delta = abs(XYk(22,:)-120);
P = -Re(30,:);
figure(1); hold on; plot(delta,P);
ValidationData = load('BeamValidation.csv');
hold on; plot(ValidationData(:,1),ValidationData(:,2));
grid on;
xlabel('\delta (cm)'); ylabel('P (kN)');
title('Load displacement curve');
legend('MATLAB','Validation');


% PLOTTING DEFORMED GEOMETRY
StepID = NumStep;
figure(2); hold on;
for mem = 1:NumMem
    MemNode = MemCon(mem,:);
    
    Node1XY1 = XYk([2*MemNode(1)-1,2*MemNode(1)],1);
    Node2XY1 = XYk([2*MemNode(2)-1,2*MemNode(2)],1);
    plot([Node1XY1(1),Node2XY1(1)],[Node1XY1(2),Node2XY1(2)],'b-')
    
    Node1XYk = XYk([2*MemNode(1)-1,2*MemNode(1)],StepID);
    Node2XYk = XYk([2*MemNode(2)-1,2*MemNode(2)],StepID);
    plot([Node1XYk(1),Node2XYk(1)],[Node1XYk(2),Node2XYk(2)],'r--')
end
xlabel('X (cm)'); ylabel('Y (cm)');
grid on;