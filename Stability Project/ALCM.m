function [XYk,Re] = ALCM(NumNodeDOF,ElementType,STol,IterMax,XY,MemCon,StrBC,LoadPattern,NumStep,E,A,I,...
    psi,ds)

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
        
        % START OF ARC LENGTH CONTROL METHOD (ALCM) ALGORITHM
        a1 = drI'*drI+psi^2*(Init.RefLoad'*Init.RefLoad);
        a2 = 2*((drk+drII)'*drI+psi^2*dLk*(Init.RefLoad'*Init.RefLoad));
        a3 = (drk+drII)'*(drk+drII)+psi^2*dLk^2*(Init.RefLoad'*Init.RefLoad)-ds^2;
        dLki1 = (-a2+(a2^2-4*a1*a3)^0.5)/2/a1;
        dLki2 = (-a2-(a2^2-4*a1*a3)^0.5)/2/a1;
        if Init.ii(kk)==0
            sgn = sign(det(Kt));
            if sign(dLki1)==sgn
                dLoadFaci = dLki1;
            elseif sign(dLki2)==sgn
                dLoadFaci = dLki2;
            end
        else
            dri1 = dLki1*drI+drII;
            dri2 = dLki2*drI+drII;
            DOT1 = (drk+dri1)'*drk+psi^2*dLk*(dLk+dLki1)*(Init.RefLoad'*Init.RefLoad);
            DOT2 = (drk+dri2)'*drk+psi^2*dLk*(dLk+dLki2)*(Init.RefLoad'*Init.RefLoad);
            if DOT1>DOT2
                dLoadFaci = dLki1;
            elseif DOT1<DOT2
                dLoadFaci = dLki2; 
            end  
        end
        % END OF ALGORITHM
        
        dLk = dLk+dLoadFaci;
        Init.LoadFac(kk) = Init.LoadFac(kk-1)+dLk;
        dri = dLoadFaci*drI+drII;
        Init.drSys(Init.BCMap) = dri;
        drk = drk+dri;
        
        Init.XYk(:,kk) = Init.XYk(:,kk)+Init.drSys(Init.XYThetaMap);
        Init.Thetak(:,kk) = Init.Thetak(:,kk)+Init.drSys(~Init.XYThetaMap);
        
        [Init.Ks(:,:,kk),Ri] = CorAssemble(ElementType,MemCon,Init.XYk(:,kk),Init.Thetak(:,kk),Init.L0,...
            Init.Beta0,Init.NumMem,Init.NumTotalDOF,Init.NumStrDOF,Init.NumMemDOF,Init.MemDOF,Init.KBCMap,...
            Init.BCMap,E,A,I);
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