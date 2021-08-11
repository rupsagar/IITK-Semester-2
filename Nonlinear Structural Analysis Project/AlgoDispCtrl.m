function State = AlgoDispCtrl(NumSec,NumFib,NumFibParam,MatModel,ConvCritStr,ConvCritEle,...
    StrMaxIter,EleMaxIter,Stol,Etol,NumTotalDOF,StrBCMap,StrKBCMap,NumStrDOF,NumMem,NumMemDOF,...
    MemDOF,MemL,NumBasDOF,bx,avq,wt,x,NumSecDOF,SecComp,FibA,EtanSig,Rc,State,NumTotalDataPts,kIni,...
    e,rx,Ssi,vsi,Ssk,vsk,Kv,Smk,vmk,Qk,qk,Ks,Rk,rk,i,DispHis,LoadFac,RefLoad,CtrlDOFMap)

FibState = State.FibState;
fsec = State.fsec;

StrChk0 = 0;
kk = 2;

for k = kIni:NumTotalDataPts
    drk = 0;
    dLoadFack = 0;
    
    drI = Ks(:,:,k-1)\RefLoad;    
    drICtrl = drI(CtrlDOFMap);
    
    dLoadFaci = (DispHis(CtrlDOFMap,kk)-DispHis(CtrlDOFMap,kk-1))/drICtrl;
    dri = dLoadFaci*drI;
    
    StrNotConv = true;
    while StrNotConv && i(k)<=StrMaxIter
        drk = drk+dri;
        
        [EtanSig,e,rx,Ssi,vsi,fsec,Kv,Smk,vmk] = ...
            FBE(k,i,dri,NumSecDOF,NumSec,NumFib,MatModel,ConvCritEle,EleMaxIter,Etol,NumTotalDOF,StrBCMap,...
            NumMem,MemDOF,MemL,NumBasDOF,bx,avq,wt,x,SecComp,FibA,EtanSig,e,rx,Ssi,vsi,FibState,fsec,...
            Ssk,vsk,Kv,Smk,vmk);
        
        Ks(:,:,k) = AssembleStiff(avq,Kv(:,:,:,k),NumStrDOF,StrKBCMap,NumTotalDOF,NumMemDOF,NumMem,MemDOF);
        Ri = AssembleForce(avq,Smk(:,:,k),NumTotalDOF,NumMem,MemDOF,StrBCMap);
                
        dLoadFack = dLoadFack+dLoadFaci;
        LoadFac(kk) = LoadFac(kk-1)+dLoadFack;       
        Ru = LoadFac(kk)*RefLoad+Rc-Ri;
        
        drI = Ks(:,:,k)\RefLoad;
        drII = Ks(:,:,k)\Ru;
        
        drICtrl = drI(CtrlDOFMap);
        drIICtrl = drII(CtrlDOFMap);
        
        dLoadFaci = -drIICtrl/drICtrl;
        dri = dLoadFaci*drI+drII;
        
        [StrNotConv,StrChk0] = ConvTest(ConvCritStr,Stol,StrChk0,StrNotConv,i(k),dri,Ru);
        i(k) = i(k)+1;
    end

    kk = kk+1;
    
    % STATE UPDATE
    [FibState,Ssk,vsk,Qk,qk,Rk,rk] = FibStrStateUpdt(k,FibState,e,Ssi,vsi,Ssk,vsk,Smk,Qk,qk,Rk,rk,...
    drk,EtanSig,Ri,NumTotalDOF,NumMem,NumSec,NumFib,NumFibParam,StrBCMap,MemDOF,avq);

end

State.NumTotalDataPts = NumTotalDataPts;
State.FibState = FibState;
State.fsec = fsec;
State.Ssk = Ssk;
State.vsk = vsk;
State.Kv = Kv;
State.Smk = Smk;
State.vmk = vmk;
State.Qk = Qk;
State.qk = qk;
State.Ks = Ks;
State.Rk = Rk;
State.rk = rk;

end