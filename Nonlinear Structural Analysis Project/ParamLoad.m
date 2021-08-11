function [NumLoadPts,LoadHis] = ParamLoad(LoadPattern,LoadCtrlData,NumTotalDOF,NumNodeDOF,StrBCMap)

ReIni = LoadCtrlData(1);
ReChng = LoadCtrlData(2);
ReFin = LoadCtrlData(3);

if ((ReFin-ReIni)>0&&ReChng<0)||((ReFin-ReIni)<0&&ReChng>0)
    ReChng = -ReChng;
end

LoadFac = ReIni:ReChng:ReFin;
NumLoadPts = numel(LoadFac);

NumLoadPattern = size(LoadPattern,1);

RefLoadAllDOF = zeros(NumTotalDOF,1);

for i = 1:NumLoadPattern
    CtrlNodeID = LoadPattern(i,1);
    RefLoadAllDOF(NumNodeDOF*CtrlNodeID-2:NumNodeDOF*CtrlNodeID) = LoadPattern(i,2:4)';
end
RefLoad = RefLoadAllDOF(StrBCMap);

LoadHis = LoadFac.*RefLoad;

end