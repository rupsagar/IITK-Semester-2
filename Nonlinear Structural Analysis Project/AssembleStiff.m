function Ks = AssembleStiff(avq,Kv,NumStrDOF,StrKBCMap,NumTotalDOF,NumMemDOF,NumMem,MemDOF)

Kg = zeros(NumTotalDOF);
Ks = zeros(NumStrDOF);

for mem = 1:NumMem
    Kq = avq(:,:,mem)'*Kv(:,:,mem)*avq(:,:,mem);
    for i = 1:NumMemDOF
        for j = 1:NumMemDOF
            Kg(MemDOF(mem,j),MemDOF(mem,i)) = Kg(MemDOF(mem,j),MemDOF(mem,i))+Kq(j,i);
        end
    end
end

Ks(:) = Kg(StrKBCMap);
end