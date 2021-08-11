function State = Analyze(Param,LoadPattern,State,AnalysisData)

% EXTRACT PARAMETERS
[NumNodeDOF,NumSec,NumFib,NumFibParam,MatModel,ConvCritStr,ConvCritEle,StrMaxIter,EleMaxIter,...
    Stol,Etol,NumTotalDOF,StrBCMap,StrKBCMap,NumStrDOF,NumMem,NumMemDOF,MemDOF,MemL,NumBasDOF,bx,...
    Kv0,Ks0,avq,wt,x,NumSecDOF,SecComp,FibA,EtanSig] = ParamExtract(Param);

% PARAMETERS FOR STATIC AND DYNAMIC ANALYSIS (DYNAMIC ANALYSIS NOT YET IMPLEMENTED)
ProblemType = AnalysisData{1};
if strcmp(ProblemType,"Static")
    Algorithm = AnalysisData{2};
    CtrlData = AnalysisData{3};
    ConstLoad = AnalysisData{4}; % CONSTANT LOAD PARAMETERS
    if numel(ConstLoad)==0
        Rc = zeros(NumStrDOF,1);
    else
        Rc = ConstLoad;
    end
    if strcmp(Algorithm,"LoadCtrl") % LOAD CONTROL PARAMETERS
        [NumDataPts,LoadHis] = ParamLoad(LoadPattern,CtrlData,NumTotalDOF,NumNodeDOF,StrBCMap);   
    elseif strcmp(Algorithm,"DispCtrl") % DISPLACEMENT CONTROL PARAMETERS
        [NumDataPts,DispHis,LoadFac,RefLoad,CtrlDOFMap] = ParamDisp(LoadPattern,CtrlData,...
            NumTotalDOF,NumNodeDOF,StrBCMap);
    end
end

% INITIALIZATION
[NumTotalDataPts,kIni,e,rx,Ssi,vsi,Ssk,vsk,Kv,Smk,vmk,Qk,qk,Ks,Rk,rk,i] = ...
    Initialize(NumSecDOF,NumSec,NumFib,NumBasDOF,NumStrDOF,NumMem,NumMemDOF,Kv0,Ks0,State,NumDataPts);

% ANALYSIS STARTS
if strcmp(Algorithm,"LoadCtrl")
    State = AlgoLoadCtrl(NumSec,NumFib,NumFibParam,MatModel,ConvCritStr,ConvCritEle,...
    StrMaxIter,EleMaxIter,Stol,Etol,NumTotalDOF,StrBCMap,StrKBCMap,NumStrDOF,NumMem,NumMemDOF,...
    MemDOF,MemL,NumBasDOF,bx,Ks0,avq,wt,x,NumSecDOF,SecComp,FibA,EtanSig,Rc,State,NumTotalDataPts,kIni,...
    e,rx,Ssi,vsi,Ssk,vsk,Kv,Smk,vmk,Qk,qk,Ks,Rk,rk,i,LoadHis,AnalysisData);

elseif strcmp(Algorithm,"DispCtrl")
    State = AlgoDispCtrl(NumSec,NumFib,NumFibParam,MatModel,ConvCritStr,ConvCritEle,...
    StrMaxIter,EleMaxIter,Stol,Etol,NumTotalDOF,StrBCMap,StrKBCMap,NumStrDOF,NumMem,NumMemDOF,...
    MemDOF,MemL,NumBasDOF,bx,avq,wt,x,NumSecDOF,SecComp,FibA,EtanSig,Rc,State,NumTotalDataPts,kIni,...
    e,rx,Ssi,vsi,Ssk,vsk,Kv,Smk,vmk,Qk,qk,Ks,Rk,rk,i,DispHis,LoadFac,RefLoad,CtrlDOFMap);
end

end