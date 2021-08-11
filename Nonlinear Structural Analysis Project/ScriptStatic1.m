% CANTILEVER MODEL WITH LOAD-CONTROLLED ANALYSIS
clear;
format shortg;

% UNITS: N,m
% INPUT PARAMETERS
g=9.81;
NumNodeDOF = 3;
NumSecDOF = 2;
NumBasDOF = 3;

Nodes = [0,0;0,3];
StrBC = [1,1,1;0,0,0];
MemCon = [1,2];
SecDim = [0.3,0.5];
BasSys = "Cantilever";
IntType = "Lobatto";
NumSec = 5;
% MatModel = {"Bilinear",415e6,210e9,.003};
MatModel = {"GMP",415e6,210e9,.003,20,0.925,0.15};
NumFib = 20;
ConvCritStr = "Energy Increment";
StrMaxIter = 10;
ConvCritEle = "Energy Increment";
EleMaxIter = 10;
Stol = 1e-16;
Etol = 1e-16;


% CALCULATED PARAMETERS
[NumTotalDOF,StrBCMap,StrKBCMap,NumStrDOF,NumMem,NumMemDOF,MemDOF,MemL,IzSec,apq] = ...
    ParamGeom(Nodes,MemCon,SecDim,StrBC,NumNodeDOF);
[SecComp,FibA,ASec,NumFibParam,FibState,EtanSig,fsec] = ...
    ParamFiber(SecDim,NumMem,NumSec,NumSecDOF,NumFib,MatModel);
[wt,x] = ParamIntgr(IntType,MemL,NumMem,NumSec);
[bx,Kv0,avq] = ParamBasSys(BasSys,NumBasDOF,NumMem,MemL,ASec,IzSec,MatModel,apq);
Ks0 = AssembleStiff(avq,Kv0,NumStrDOF,StrKBCMap,NumTotalDOF,NumMemDOF,NumMem,MemDOF);
Param = ParamCompress(NumNodeDOF,NumSec,NumFib,NumFibParam,MatModel,ConvCritStr,ConvCritEle,...
    StrMaxIter,EleMaxIter,Stol,Etol,NumTotalDOF,StrBCMap,StrKBCMap,NumStrDOF,NumMem,NumMemDOF,...
    MemDOF,MemL,NumBasDOF,bx,Kv0,Ks0,avq,wt,x,NumSecDOF,SecComp,FibA,EtanSig);
State = ParamState0(FibState,fsec);


% LOAD CONTROLLED ANALYSIS - GRAVITY LOAD
LoadPattern = [2,0,-1,0;3,0,-1,0];
PCol = 10000e3;
PIncr = 50e3;
LoadCtrlData = [0,PIncr,PCol];
AnalysisData = {"Static","LoadCtrl",LoadCtrlData,[]};
State = Analyze(Param,LoadPattern,State,AnalysisData);

% LOAD CONTROLLED ANALYSIS - TOP HORIZONTAL LOAD
LoadPattern = [2,1,0,0];
ConstLoad = State.Rk(:,end);
HIncr = 50e3;
LoadCtrlDataSeries = [0,HIncr,4000e3];
for i = 1:size(LoadCtrlDataSeries,1)
    LoadCtrlData = LoadCtrlDataSeries(i,:);
    AnalysisData = {"Static","LoadCtrl",LoadCtrlData,ConstLoad};
    State = Analyze(Param,LoadPattern,State,AnalysisData);
end


% OPENSEES RESULTS
system('"OpenSees.exe" ScriptStatic1.tcl');
DispTop = load('OSResults/DispTop.txt');
ForceCol = load('OSResults/ForceCol.txt');
DefoColSec = load('OSResults/DefoColSec1.txt');
ForceColSec = load('OSResults/ForceColSec1.txt');


% PLOTTING RESULTS
LineWidth = 1.0;
mem = 1;
sec = 1;
DOF = 4;
PlotResult(LineWidth,State,mem,sec,DOF,DispTop,ForceCol,DefoColSec,ForceColSec);