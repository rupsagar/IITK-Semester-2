function EtanSig = FibTrialState(NumFib,MatModel,mem,sec,e,FibState,EtanSig)

if strcmp(MatModel{1},"Bilinear")
    EtanSig = MatBilinear(NumFib,MatModel,mem,sec,e,FibState,EtanSig);
elseif strcmp(MatModel{1},"GMP")
    EtanSig = MatGMP(NumFib,MatModel,mem,sec,e,FibState,EtanSig);
elseif strcmp(MatModel{1},"Elastic")
    EtanSig = MatElastic(NumFib,MatModel,mem,sec,e,FibState,EtanSig);
end

end