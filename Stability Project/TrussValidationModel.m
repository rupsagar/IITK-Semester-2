clear; format shortg;

% MODEL DATA
NumNodeDOF = 3;
ElementType = 'Truss';
STol = 1e-16;
IterMax = 10;
L_0 = 10;
T0 = pi()/3;
LoadNodeID = 2;
LoadPattern = [LoadNodeID,0,-1,0];

XY = [0,L_0*cos(T0),2*L_0*cos(T0);0,L_0*sin(T0),0];
XY = XY(:);
NumNode = numel(XY)/2;

MemCon = [1,2;2,3];
A = [10;10];
I = [1;1];
NumMem = size(MemCon,1);
E = 2e11*ones(NumMem,1);

StrBC = [false(2,NumNode);true(1,NumNode)];
StrBC(1:2,[1,NumNode]) = ones(2,2);


% ANALYSIS
% DispHis = 0:-L_0/20:-2.5*L_0;
% CtrlNodeID = LoadNodeID;
% CtrlDOFID = 2;
% NumStep = numel(DispHis);
% [XYk,Re] = DCM(NumNodeDOF,ElementType,STol,IterMax,XY,MemCon,StrBC,LoadPattern,NumStep,E,A,I,...
%     LoadNodeID,CtrlDOFID,DispHis);

% NumStep = 125;
% psi = 0;
% ds = 2e-1;
% [XYk,Re] = ALCM(NumNodeDOF,ElementType,STol,IterMax,XY,MemCon,StrBC,LoadPattern,NumStep,E,A,I,psi,ds);

NumStep = 150;
dLambdabar = 5e10;
[XYk,Re] = GDCM(NumNodeDOF,ElementType,STol,IterMax,XY,MemCon,StrBC,LoadPattern,NumStep,E,A,I,dLambdabar);


% PLOTTING RESULTS
set(0, 'DefaultLineLineWidth',1.0);

delta = abs(XYk(4,:)-L_0*sin(pi()/3))/L_0;
P = -Re(2,:)/(2*E(1)*A(1));
figure(1); hold on; plot(delta,P,'-');
DispHis = 0:-L_0/20:-2.5*L_0;
a = abs(DispHis)/L_0;
Lambda = (1./sqrt(1-2*a*sin(T0)+a.^2)-1).*(sin(T0)-a);
hold on; plot(a,Lambda);
grid on;
xlabel('Non dimensional displacement'); ylabel('Non dimensional load');
title('Non dimensional load displacement curve');
legend('MATLAB','Analytical');


% PLOTTING DEFORMED GEOMETRY
% StepID = NumStep;
% figure(2); hold on;
% for mem = 1:NumMem
%     MemNode = MemCon(mem,:);
%     
%     Node1XY1 = XYk([2*MemNode(1)-1,2*MemNode(1)],1);
%     Node2XY1 = XYk([2*MemNode(2)-1,2*MemNode(2)],1);
%     plot([Node1XY1(1),Node2XY1(1)],[Node1XY1(2),Node2XY1(2)],'b-')
%     
%     Node1XYk = XYk([2*MemNode(1)-1,2*MemNode(1)],StepID);
%     Node2XYk = XYk([2*MemNode(2)-1,2*MemNode(2)],StepID);
%     plot([Node1XYk(1),Node2XYk(1)],[Node1XYk(2),Node2XYk(2)],'r--')
% end
% xlabel('X'); ylabel('Y');
% grid on;