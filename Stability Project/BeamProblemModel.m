clear; format shortg;

% MODEL DATA
NumNodeDOF = 3;
ElementType = 'Frame';
STol = 1e-16;
IterMax = 10;

NumNode = 21;
LoadNodeID = (NumNode+1)/2;
LoadPattern = [LoadNodeID,0,-1,0];

MidHeight = 9e-3;
Length = 300e-3;
Depth = 5e-3;
Breadth = 120e-3;
Func = @(x)MidHeight*(1-4*x.^2/Length^2);
X = linspace(-Length/2,Length/2,NumNode);
Y = Func(X);
XY = [X;Y];
XY = XY(:);

MemCon = [1:(NumNode-1);2:NumNode]';
NumMem = size(MemCon,1);
E = 2e11*ones(NumMem,1);
A = Breadth*Depth*ones(NumMem,1);
I = Breadth*Depth^3/12*ones(NumMem,1);

StrBC = false(3,NumNode);
StrBC(:,[1,NumNode]) = ones(3,2);


% ANALYSIS
DispHis = 0:-0.04e-3:-16e-3;
CtrlNodeID = LoadNodeID;
CtrlDOFID = 2;
NumStep = numel(DispHis);
[XYk,Re] = DCM(NumNodeDOF,ElementType,STol,IterMax,XY,MemCon,StrBC,LoadPattern,NumStep,E,A,I,...
    LoadNodeID,CtrlDOFID,DispHis);

% NumStep = 25;
% psi = 0;
% ds = 2e-2;
% [XYk,Re] = ALCM(NumNodeDOF,ElementType,STol,IterMax,XY,MemCon,StrBC,LoadPattern,NumStep,E,A,I,psi,ds);

% NumStep = 100;
% dLambdabar = 1e3;
% [XYk,Re] = GDCM(NumNodeDOF,ElementType,STol,IterMax,XY,MemCon,StrBC,LoadPattern,NumStep,E,A,I,dLambdabar);


% PLOTTING RESULTS
fs = 15;
set(0, 'DefaultLineLineWidth',1.2);

delta = abs(XYk((NumNode+1),:)-MidHeight);
P = abs(Re((NumNode-3)*3/2+2,:));
figure(1); hold on; plot(delta,P);
PaperData = load('TrussProblem.csv');
hold on; plot(PaperData(:,1),PaperData(:,2));
grid on;
xlabel('\delta (m)'); ylabel('P (N)');
title('Load displacement curve');
legend('MATLAB','Xenidis et al. 2013');


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
xlabel('X (m)'); ylabel('Y (m)');
xlim([-Length,Length]); ylim([-Length/4,Length/4]);
grid on;