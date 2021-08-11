function [XYk,Re] = GDCM(NumNodeDOF,ElementType,STol,IterMax,XY,MemCon,StrBC,LoadPattern,NumStep,...
    E,A,I,dLambdabar)

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
        
        % START OF GENERALIZED DISPLACEMENT CONTROL METHOD (GDCM) ALGORITHM
        if Init.ii(kk)==0
            if kk==2
                drI1 = Init.Ks(:,:,1)\Init.RefLoad;
                drIprev = drI1;
                sgnprev = 1;
            else
                drIprev = drIprevtemp;
            end
            drII = 0;
            drIprevtemp = drI;
            GSP = (drI1'*drI1)/(drIprev'*drI);
            sgnprev = sgnprev*sign(GSP);
            dLi = sgnprev*dLambdabar*abs(GSP)^0.5;
            
        else
            drII = Kt\Init.Ru;
            dLi = -(drIprev'*drII)/(drIprev'*drI);
        end
        % END OF ALGORITHM
        
        dLk = dLk+dLi;
        Init.LoadFac(kk) = Init.LoadFac(kk-1)+dLk;
        dri = dLi*drI+drII;
        Init.drSys(Init.BCMap) = dri;
        drk = drk+dri;
        
        Init.XYk(:,kk) = Init.XYk(:,kk)+Init.drSys(Init.XYThetaMap);
        Init.Thetak(:,kk) = Init.Thetak(:,kk)+Init.drSys(~Init.XYThetaMap);
        
        [Init.Ks(:,:,kk),Ri] = ...
            CorAssemble(ElementType,MemCon,Init.XYk(:,kk),Init.Thetak(:,kk),Init.L0,Init.Beta0,Init.NumMem,...
            Init.NumTotalDOF,Init.NumStrDOF,Init.NumMemDOF,Init.MemDOF,Init.KBCMap,Init.BCMap,E,A,I);
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