function [NumNodeDOF,NumSec,NumFib,NumFibParam,MatModel,ConvCritStr,ConvCritEle,StrMaxIter,EleMaxIter,...
    Stol,Etol,NumTotalDOF,StrBCMap,StrKBCMap,NumStrDOF,NumMem,NumMemDOF,MemDOF,MemL,NumBasDOF,bx,Kv0,Ks0,avq,wt,x,...
    NumSecDOF,SecComp,FibA,EtanSig] = ParamExtract(Param)

NumNodeDOF = Param{1};
NumSecDOF = Param{2};
MatModel = Param{3};
ConvCritStr = Param{4};
ConvCritEle = Param{5};
StrMaxIter = Param{6};
EleMaxIter = Param{7};
Stol = Param{8};
Etol = Param{9};
NumFib = Param{10};
NumBasDOF = Param{11};
NumMem = Param{12};
NumMemDOF = Param{13};
NumTotalDOF = Param{14};
NumSec = Param{15};
NumFibParam = Param{16};
MemDOF = Param{17};
MemL = Param{18};
StrBCMap = Param{19};
StrKBCMap = Param{20};
bx = Param{21};
Kv0 = Param{22};
Ks0 = Param{23};
avq = Param{24};
wt = Param{25};
x = Param{26};
SecComp = Param{27};
FibA = Param{28};
EtanSig = Param{29};
NumStrDOF = Param{30};

end