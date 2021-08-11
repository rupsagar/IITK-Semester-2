function [XYk,Re] = DCM(NumNodeDOF,ElementType,STol,IterMax,XY,MemCon,StrBC,LoadPattern,NumStep,E,A,I,...
    CtrlNodeID,CtrlDOFID,DispHis)

% INITIALIZATION
Init = Initialize(NumNodeDOF,ElementType,XY,MemCon,StrBC,LoadPattern,NumStep,E,A,I);

% ITERATION STARTS
for kk = 2:NumStep
    Init.XYk(:,kk) = Init.XYk(:,kk-1);
    Init.Thetak(:,kk) = Init.Thetak(:,kk-1);
    dLk = 0;
    drk = 0;
    Kt = Init.Ks(:,:,kk-1);
    
    StrNotConv = true;
    while StrNotConv && Init.ii(kk)<IterMax
        drI = Kt\Init.RefLoad;
        drII = Kt\Init.Ru;
        
        % START OF DISPLACEMENT CONTROL METHOD (DCM) ALGORITHM
        Init.drSys(Init.BCMap) = drI;
        drICtrl = Init.drSys(NumNodeDOF*(CtrlNodeID-1)+CtrlDOFID);
        Init.drSys(Init.BCMap) = drII;
        drIICtrl = Init.drSys(NumNodeDOF*(CtrlNodeID-1)+CtrlDOFID);
        if Init.ii(kk)==0
            dLoadFaci = (DispHis(kk)-DispHis(kk-1))/drICtrl-drIICtrl/drICtrl;
        else
            dLoadFaci = -drIICtrl/drICtrl;
        end
        % END OF ALGORITHM
        
        dLk = dLk+dLoadFaci;
        Init.LoadFac(kk) = Init.LoadFac(kk-1)+dLk;
        dri = dLoadFaci*drI+drII;
        Init.drSys(Init.BCMap) = dri;
        drk = drk+dri;
        
        Init.XYk(:,kk) = Init.XYk(:,kk)+Init.drSys(Init.XYThetaMap);
        Init.Thetak(:,kk) = Init.Thetak(:,kk)+Init.drSys(~Init.XYThetaMap);
        
        [Init.Ks(:,:,kk),Ri] = ...
            CorAssemble(ElementType,MemCon,Init.XYk(:,kk),Init.Thetak(:,kk),Init.L0,Init.Beta0,...
            Init.NumMem,Init.NumTotalDOF,Init.NumStrDOF,Init.NumMemDOF,Init.MemDOF,Init.KBCMap,Init.BCMap,E,A,I);
        Init.Ru = Init.LoadFac(kk)*Init.RefLoad-Ri;
        
        StrChk = abs(Init.Ru'*dri);
        if StrChk<STol
            StrNotConv = false;
        else
            Kt = Init.Ks(:,:,kk);
        end
        Init.ii(kk) = Init.ii(kk)+1;
    end
    Init.Re(:,kk) = Ri;
end

XYk = Init.XYk;
Re = Init.Re;