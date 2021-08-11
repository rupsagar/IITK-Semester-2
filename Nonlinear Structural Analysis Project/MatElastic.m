function EtanSig = MatElastic(NumFib,MatModel,mem,sec,e,FibState,EtanSig)

E = ParamMat(MatModel);

for fib = 1:NumFib
    eps = [FibState{mem}(fib,2,sec),e{mem}(fib,sec)];
    NumStrainPts = numel(eps);
    
    sigma = zeros(1,NumStrainPts);
    sigma(1) = FibState{mem}(fib,3,sec);
    
    deps = eps(2)-eps(1);
    
    sigma(2) = sigma(1)+E*deps;
    
    EtanSig{mem,1}(fib,fib,sec) = E; % storing Et
    EtanSig{mem,2}(fib,1,sec) = sigma(end); % storing sigma
end

end