function Ri = AssembleForce(avq,Sm,NumTotalDOF,NumMem,MemDOF,StrBCMap)

Rg = zeros(NumTotalDOF,1);

for mem = 1:NumMem
    QMem = avq(:,:,mem)'*Sm(:,mem);
    Rg(MemDOF(mem,:)) = Rg(MemDOF(mem,:))+QMem;
end

Ri = Rg(StrBCMap);

end