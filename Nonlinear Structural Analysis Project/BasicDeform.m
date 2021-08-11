function dvm = BasicDeform(mem,avq,dri,NumTotalDOF,MemDOF,StrBCMap)

drGlo = zeros(NumTotalDOF,1);
drGlo(StrBCMap) = dri;

dqm = drGlo(MemDOF(mem,:));
dvm = avq(:,:,mem)*dqm;

end