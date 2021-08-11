clear; format shortg;

% MODEL DATA
NumNodeDOF = 3;
ElementType = 'Truss';
STol = 1e-16;
IterMax = 10;
NumNodeChord = 7;
LoadNodeID = (3*NumNodeChord+1)/2;
LoadPattern = [LoadNodeID,0,-1,0];

MidHeight = 9e-3;
Length = 300e-3;
Depth = 5e-3;
AreaChord = 100e-6;
AreaDiagonal = 200e-6;
AreaStrut = 30e-6;
Func = @(x)MidHeight*(1-4*x.^2/Length^2);
X = linspace(-Length/2,Length/2,NumNodeChord);
YBottom = Func(X)-Depth/2;
YTop = Func(X)+Depth/2;
XY = [[X;YBottom],[X;YTop]];
XY = XY(:);
NumNode = numel(XY)/2;

MemCon = [(1:(NumNodeChord-1))',(2:NumNodeChord)';
          ((NumNodeChord+1):(2*NumNodeChord-1))',((NumNodeChord+2):2*NumNodeChord)';
          (1:(NumNodeChord-1))',((NumNodeChord+2):2*NumNodeChord)';
          (2:NumNodeChord)',(NumNodeChord+1:(2*NumNodeChord-1))';
          (1:NumNodeChord)',(NumNodeChord+1:2*NumNodeChord)'];
A = [AreaChord*ones(NumNodeChord-1,1);
    AreaChord*ones(NumNodeChord-1,1);
    AreaDiagonal*ones(NumNodeChord-1,1);
    AreaDiagonal*ones(NumNodeChord-1,1);
    AreaStrut*ones(NumNodeChord,1)];
NumMem = size(MemCon,1);
E = 2e11*ones(NumMem,1);
I = ones(NumMem,1);

StrBC = [false(2,NumNode);true(1,NumNode)];
StrBC(1:2,[1,NumNodeChord,NumNodeChord+1,2*NumNodeChord]) = ones(2,4);


% ANALYSIS
% DispHis = 0:-0.04e-3:-16e-3;
% CtrlNodeID = LoadNodeID;
% CtrlDOFID = 2;
% NumStep = numel(DispHis);
% [XYk,Re] = DCM(NumNodeDOF,ElementType,STol,IterMax,XY,MemCon,StrBC,LoadPattern,NumStep,E,A,I,...
%     LoadNodeID,CtrlDOFID,DispHis);

% NumStep = 200;
% psi = 0;
% ds = 2e-4;
% [XYk,Re] = ALCM(NumNodeDOF,ElementType,STol,IterMax,XY,MemCon,StrBC,LoadPattern,NumStep,E,A,I,psi,ds);

% NumStep = 250;
% dLambdabar = 5e2;
% [XYk,Re] = GDCM(NumNodeDOF,ElementType,STol,IterMax,XY,MemCon,StrBC,LoadPattern,NumStep,E,A,I,dLambdabar);



% PLOTTING RESULTS
fs = 15;
set(0, 'DefaultLineLineWidth',1.2);

delta_m = abs((XYk(NumNodeChord+1,:)+XYk(3*NumNodeChord+1,:))/2-MidHeight);
P = abs(Re((NumNodeChord-2)*2+NumNodeChord-1,:));
figure(1); hold on; plot(delta_m,P);
PaperData = load('TrussProblem.csv');
hold on; plot(PaperData(:,1),PaperData(:,2));
grid on;
xlabel('\delta_m (m)'); ylabel('P (N)');
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
    plot([Node1XYk(1),Node2XYk(1)],[Node1XYk(2),Node2XYk(2)],'r-')
end
xlabel('X (m)'); ylabel('Y (m)');
xlim([-Length,Length]); ylim([-Length/4,Length/4]);
grid on;