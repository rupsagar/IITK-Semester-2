function Param = ...
    ParamCompress(NumNodeDOF,NumSec,NumFib,NumFibParam,MatModel,ConvCritStr,ConvCritEle,StrMaxIter,EleMaxIter,...
    Stol,Etol,NumTotalDOF,StrBCMap,StrKBCMap,NumStrDOF,NumMem,NumMemDOF,MemDOF,MemL,NumBasDOF,bx,Kv0,Ks0,avq,wt,x,...
    NumSecDOF,SecComp,FibA,EtanSig)

Param = cell(28,1);

Param{1} = NumNodeDOF;
Param{2} = NumSecDOF;
Param{3} = MatModel;
Param{4} = ConvCritStr;
Param{5} = ConvCritEle;
Param{6} = StrMaxIter;
Param{7} = EleMaxIter;
Param{8} = Stol;
Param{9} = Etol;
Param{10} = NumFib;
Param{11} = NumBasDOF;
Param{12} = NumMem;
Param{13} = NumMemDOF;
Param{14} = NumTotalDOF;
Param{15} = NumSec;
Param{16} = NumFibParam;
Param{17} = MemDOF;
Param{18} = MemL;
Param{19} = StrBCMap;
Param{20} = StrKBCMap;
Param{21} = bx;
Param{22} = Kv0;
Param{23} = Ks0;
Param{24} = avq;
Param{25} = wt;
Param{26} = x;
Param{27} = SecComp;
Param{28} = FibA;
Param{29} = EtanSig;
Param{30} = NumStrDOF;

end