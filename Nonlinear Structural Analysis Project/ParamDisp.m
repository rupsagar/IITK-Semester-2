function [NumDispPts,DispHis,LoadFac,RefLoad,CtrlDOFMap] = ...
    ParamDisp(LoadPattern,DispCtrlData,NumTotalDOF,NumNodeDOF,StrBCMap)

DispNode = DispCtrlData(1);
DispDOF = DispCtrlData(2);
DispIni = DispCtrlData(3);
DispChng = DispCtrlData(4);
DispFin = DispCtrlData(5);

if ((DispFin-DispIni)>0&&DispChng<0)||((DispFin-DispIni)<0&&DispChng>0)
    DispChng = -DispChng;
end

DispPts = DispIni:DispChng:DispFin;
NumDispPts = numel(DispPts);

NumLoadPattern = size(LoadPattern,1);

LoadFac = zeros(NumDispPts,1);
RefLoadAllDOF = zeros(NumTotalDOF,1);
CtrlDOFMap = false(NumTotalDOF,1);
DispHis = zeros(NumTotalDOF,NumDispPts);

for i = 1:NumLoadPattern
    DOF = LoadPattern(i,1);
    RefLoadAllDOF(3*DOF-2:3*DOF,:) = LoadPattern(i,2:4)';
end
RefLoad = RefLoadAllDOF(StrBCMap);

CtrlDOFMap(NumNodeDOF*(DispNode-1)+DispDOF) = 1;
DispHis(CtrlDOFMap,:) = DispPts;

CtrlDOFMap = CtrlDOFMap(StrBCMap);
DispHis = DispHis(StrBCMap,:);

end