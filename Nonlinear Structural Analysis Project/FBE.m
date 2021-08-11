function [EtanSig,e,rx,Ssi,vsi,fsec,Kv,Smk,vmk] = ...
            FBE(k,i,dri,NumSecDOF,NumSec,NumFib,MatModel,ConvCritEle,EleMaxIter,Etol,NumTotalDOF,StrBCMap,...
            NumMem,MemDOF,MemL,NumBasDOF,bx,avq,wt,x,SecComp,FibA,EtanSig,e,rx,Ssi,vsi,FibState,fsec,...
            Ssk,vsk,Kv,Smk,vmk)

for mem = 1:NumMem
    dvm = BasicDeform(mem,avq,dri,NumTotalDOF,MemDOF,StrBCMap);
    if i(k)==0
        Ssi{mem} = Ssk{mem}(:,:,k-1);
        vsi{mem} = vsk{mem}(:,:,k-1);
    end
    
    EleChk0 = 0;
    EleNotConv = true;
    j = 0;
    while EleNotConv && j<=EleMaxIter
        if j==0 && i(k)==0
            dSm = Kv(:,:,mem,k-1)*dvm;
            Smk(:,mem,k) = Smk(:,mem,k-1)+dSm;
            vmk(:,mem,k) = vmk(:,mem,k-1)+dvm;
        else
            dSm = Kv(:,:,mem,k)*dvm;
            Smk(:,mem,k) = Smk(:,mem,k)+dSm;
            vmk(:,mem,k) = vmk(:,mem,k)+dvm;
        end

        Fv = zeros(NumBasDOF);
        s = zeros(NumBasDOF,1);
        for sec = 1:NumSec(mem)
            dSs = bx{mem}(x{mem}(sec))*dSm;
            Ssi{mem}(:,sec) = Ssi{mem}(:,sec)+dSs;
            if j==0
                dvs = fsec{mem}(:,:,sec)*dSs;
            else
                dvs = rx{mem}(:,sec)+fsec{mem}(:,:,sec)*dSs;
            end
            vsi{mem}(:,sec) = vsi{mem}(:,sec)+dvs;
            e{mem}(:,sec) = SecComp{mem}*vsi{mem}(:,sec);

            % FIBER STATE DETERMINATION
            EtanSig = FibTrialState(NumFib,MatModel,mem,sec,e,FibState,EtanSig);

            ksec = SecComp{mem}'*EtanSig{mem,1}(:,:,sec)*FibA{mem}*SecComp{mem};
            fsec{mem}(:,:,sec) = ksec\eye(NumSecDOF);
            Ssr = SecComp{mem}'*FibA{mem}*EtanSig{mem,2}(:,1,sec);
            Ssu = Ssi{mem}(:,sec)-Ssr;
            rx{mem}(:,sec) = fsec{mem}(:,:,sec)*Ssu;

            Fv = Fv + MemL(mem)/2*wt{mem}(sec)*bx{mem}(x{mem}(sec))'*fsec{mem}(:,:,sec)*bx{mem}(x{mem}(sec));
            s = s + MemL(mem)/2*wt{mem}(sec)*bx{mem}(x{mem}(sec))'*rx{mem}(:,sec);
        end
        
        [EleNotConv,EleChk0] = ConvTest(ConvCritEle,Etol,EleChk0,EleNotConv,j,dvm,dSm);

        Kv(:,:,mem,k) = Fv\eye(NumBasDOF);
        dvm = -s;

        j = j+1;           
    end

end

end