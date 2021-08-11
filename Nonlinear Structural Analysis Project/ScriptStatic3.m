% PORTAL FRAME MODEL WITH LOAD-CONTROLLED STATIC ANALYSIS WITH A GIVEN LOADING HISTORY
clear;
format shortg;

% UNITS: N,m
% INPUT PARAMETERS
g=9.81;
NumNodeDOF = 3;
NumSecDOF = 2;
NumBasDOF = 3;

Nodes = [0,0;0,3;5,3;5,0];
StrBC = [1,1,1;0,0,0;0,0,0;1,1,1];
MemCon = [1,2;2,3;3,4];
SecDim = [0.3,0.5;0.3,0.5;0.3,0.5];
BasSys = "Cantilever";
IntType = ["Lobatto","Lobatto","Lobatto"];
NumSec = [5,5,5];
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
PCol = 20000e3;
PIncr = 50e3;
LoadPattern = [2,0,-1,0;3,0,-1,0];
LoadCtrlData = [0,PIncr,PCol];
AnalysisData = {"Static","LoadCtrl",LoadCtrlData,[]};
State = Analyze(Param,LoadPattern,State,AnalysisData);

% LOAD CONTROLLED ANALYSIS - TOP HORIZONTAL LOAD
LoadPattern = [2,1,0,0];
ConstLoad = State.Rk(:,end);
HIncr = 50e3;
LoadCtrlDataSeries = [0,11750,-10850,12350,-11650,12850]*1e3;
for i = 2:numel(LoadCtrlDataSeries)
    LoadCtrlData = [LoadCtrlDataSeries(i-1),HIncr,LoadCtrlDataSeries(i)];
    AnalysisData = {"Static","LoadCtrl",LoadCtrlData,ConstLoad};
    State = Analyze(Param,LoadPattern,State,AnalysisData);
end


% PLOTTING RESULTS
LineWidth = 1.0;
mem = 1;
sec = 1;
DOF = 4;
PlotResult(LineWidth,State,mem,sec,DOF);