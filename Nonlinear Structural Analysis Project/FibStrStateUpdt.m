function [FibState,Ssk,vsk,Qk,qk,Rk,rk] = FibStrStateUpdt(k,FibState,e,Ssi,vsi,Ssk,vsk,Smk,Qk,qk,Rk,rk,...
    drk,EtanSig,Ri,NumTotalDOF,NumMem,NumSec,NumFib,NumFibParam,StrBCMap,MemDOF,avq)

rg = zeros(NumTotalDOF,1);

Rk(:,k) = Ri;
rk(:,k) = rk(:,k-1)+drk;
rg(StrBCMap) = rk(:,k);

for mem = 1:NumMem
    for sec = 1:NumSec(mem)
        for fib = 1:NumFib
            FibState{mem}(fib,1,sec) = EtanSig{mem,1}(fib,fib,sec); % Et
            FibState{mem}(fib,2,sec) = e{mem}(fib,sec); % eps
            for fibparam = 3:NumFibParam
                FibState{mem}(fib,fibparam,sec) = EtanSig{mem,2}(fib,fibparam-2,sec);
            end
        end
    end
    Ssk{mem}(:,:,k) = Ssi{mem};
    vsk{mem}(:,:,k) = vsi{mem};
    Qk(:,mem,k) = avq(:,:,mem)'*Smk(:,mem,k);
    qk(:,mem,k) = rg(MemDOF(mem,:));
end

end