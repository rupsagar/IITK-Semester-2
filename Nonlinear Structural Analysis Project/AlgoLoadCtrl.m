function State = AlgoLoadCtrl(NumSec,NumFib,NumFibParam,MatModel,ConvCritStr,ConvCritEle,...
    StrMaxIter,EleMaxIter,Stol,Etol,NumTotalDOF,StrBCMap,StrKBCMap,NumStrDOF,NumMem,NumMemDOF,...
    MemDOF,MemL,NumBasDOF,bx,Ks0,avq,wt,x,NumSecDOF,SecComp,FibA,EtanSig,Rc,State,NumTotalDataPts,kIni,...
    e,rx,Ssi,vsi,Ssk,vsk,Kv,Smk,vmk,Qk,qk,Ks,Rk,rk,i,LoadHis,AnalysisData)

FibState = State.FibState;
fsec = State.fsec;

if strcmp(AnalysisData{1},"Static")
    Kbar =  0;
elseif strcmp(AnalysisData{1},"Dynamic") % MODIFICATION IN PROGRESS
    Kbar = AnalysisData{2};
end

StrChk0 = 0;
kk = 2;

for k = kIni:NumTotalDataPts
    drk = 0;
    
    Ru = LoadHis(:,kk)-LoadHis(:,kk-1);
    
    dW = rk(:,k-1)'*Ru;
    if dW>=0
        Kt = Ks(:,:,k-1);
    else
        Kt = Ks0;
    end
    
    Khat =  Kt+Kbar;
    dri = Khat\Ru;
    
    StrNotConv = true;
    while StrNotConv && i(k)<=StrMaxIter
        drk = drk+dri;
        
        [EtanSig,e,rx,Ssi,vsi,fsec,Kv,Smk,vmk] = ...
            FBE(k,i,dri,NumSecDOF,NumSec,NumFib,MatModel,ConvCritEle,EleMaxIter,Etol,NumTotalDOF,StrBCMap,...
            NumMem,MemDOF,MemL,NumBasDOF,bx,avq,wt,x,SecComp,FibA,EtanSig,e,rx,Ssi,vsi,FibState,fsec,...
            Ssk,vsk,Kv,Smk,vmk);
        
        Ks(:,:,k) = AssembleStiff(avq,Kv(:,:,:,k),NumStrDOF,StrKBCMap,NumTotalDOF,NumMemDOF,NumMem,MemDOF);
        Ri = AssembleForce(avq,Smk(:,:,k),NumTotalDOF,NumMem,MemDOF,StrBCMap);
                
        [StrNotConv,StrChk0] = ConvTest(ConvCritStr,Stol,StrChk0,StrNotConv,i(k),dri,Ru);
        
        if strcmp(AnalysisData{1},"Static")
            Ru = LoadHis(:,kk)+Rc-Ri;
        elseif strcmp(AnalysisData{1},"Dynamic") % MODIFICATION IN PROGRESS
%             drdk = gamma/beta/dt*drk-gamma/beta*rd(:,k-1)+dt*(1-gamma/2/beta)*rdd(:,k-1);
%             drddk = drk/beta/dt^2-rd(:,k-1)/beta/dt-rdd(:,k-1)/2/beta;
%             Ru = LoadHis(:,kk)+Rc-(M*(rdd(:,k-1)+drddk)+C*(rd(:,k-1)+drdk)+Ri);
        end      
        
        dri = Ks(:,:,k)\Ru;
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