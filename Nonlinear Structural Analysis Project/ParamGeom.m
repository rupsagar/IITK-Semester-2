function [NumTotalDOF,StrBCMap,StrKBCMap,NumStrDOF,NumMem,NumMemDOF,MemDOF,MemL,IzSec,apq] = ...
    ParamGeom(Nodes,MemCon,SecDim,StrBC,NumNodeDOF)

NumNodes = size(Nodes,1);
NumTotalDOF = NumNodes*NumNodeDOF;

StrBC = StrBC';
StrBC = StrBC(:);
StrBCMap = logical(~StrBC);
StrKBCMap = logical(StrBCMap'.*StrBCMap);

NumStrDOF = sum(StrBCMap);
NumMem = size(MemCon,1);
NumNodePerMem = size(MemCon,2);
NumMemDOF = NumNodeDOF*NumNodePerMem;
MemDOF = zeros(NumMem,NumMemDOF);
MemL = zeros(NumMem,1);
IzSec = zeros(NumMem,1);
apq = zeros(6,6,NumMem);
DOFIndex = zeros(1,NumNodeDOF);

for i = 1:NumNodeDOF
    DOFIndex(i) = NumNodeDOF-i;
end

for i = 1:NumMem
    MemNodes = MemCon(i,:);
    MemDOF(i,:) = [NumNodeDOF*MemNodes(1)-DOFIndex,NumNodeDOF*MemNodes(2)-DOFIndex];
    
    ri = diff(Nodes(MemNodes,:));
    MemL(i) = norm(ri);
    
    BSec = SecDim(i,1);
    DSec = SecDim(i,2);
    IzSec(i) = BSec*DSec^3/12;
    
    CS = ri/MemL(i);
    t = [CS(1),CS(2),0;
        -CS(2),CS(1),0;
        0,0,1];
    apq(:,:,i) = [t,zeros(3);zeros(3),t];
end

end