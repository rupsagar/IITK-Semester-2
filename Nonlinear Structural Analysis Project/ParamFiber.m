function [SecComp,FibA,ASec,NumFibParam,FibState,EtanSig,fsec0] = ...
    ParamFiber(SecDim,NumMem,NumSec,NumSecDOF,NumFib,MatModel)

if strcmp(MatModel{1},"Bilinear")
    NumFibParam = 5; % E,eps,sigma,sigma_lin,sigma_epp
    E = MatModel{3};
elseif strcmp(MatModel{1},"GMP")
    NumFibParam = 13; % E,eps,sigma,depsi,epsr,sigmar,eps0,sigma0,epsxi(1),epsxi(2),R,epsmax(1),epsmax(2)
    Fy = MatModel{2};
    E = MatModel{3};
    R0 = MatModel{5};
elseif strcmp(MatModel{1},"Elastic")
    NumFibParam = 3; % E,eps,sigma
    E = MatModel{2};
end

SecComp = cell(NumMem,1);
FibA = cell(NumMem,1);
ASec = cell(NumMem,1);
FibState = cell(NumMem,1);
EtanSig = cell(NumMem,2);
fsec0 = cell(NumMem,1);

for mem = 1:NumMem
    BSec = SecDim(mem,1);
    DSec = SecDim(mem,2);
    FibZCoord = linspace(DSec/2,-DSec/2,NumFib+1);
    SecComp{mem} = zeros(NumFib,2);
    FibA{mem} = zeros(NumFib);
    ASec{mem} = zeros(NumSec(mem),1);
    FibState{mem} = zeros(NumFib,NumFibParam,NumSec(mem));
    EtanSig{mem,1} = zeros(NumFib,NumFib,NumSec(mem)); % tangent modulus
    EtanSig{mem,2} = zeros(NumFib,NumFibParam-2,NumSec(mem));
    fsec0{mem} = zeros(NumSecDOF,NumSecDOF,NumSec(mem));
    
    for fib = 1:NumFib
        SecComp{mem}(fib,:) = [1,-sum(FibZCoord(fib:fib+1))/2];
        FibA{mem}(fib,fib) = abs(FibZCoord(fib+1)-FibZCoord(fib))*BSec;
        if strcmp(MatModel{1},"Bilinear")||strcmp(MatModel{1},"Elastic")
            FibState{mem}(fib,1,:) = E;
        elseif strcmp(MatModel{1},"GMP")
            FibState{mem}(fib,1,:) = E;
            FibState{mem}(fib,7,:) = Fy/E;
            FibState{mem}(fib,8,:) = Fy;
            FibState{mem}(fib,9,:) = -1;
            FibState{mem}(fib,10,:) = 1;
            FibState{mem}(fib,11,:) = R0;
            FibState{mem}(fib,12,:) = -Fy/E;
            FibState{mem}(fib,13,:) = Fy/E;
        end
        for sec = 1:NumSec(mem)
            EtanSig{mem,1}(fib,fib,sec) = E;
        end
    end

    for sec = 1:NumSec(mem)
        ASec{mem}(sec) = sum(FibA{mem}(:));
        ksec = SecComp{mem}'*EtanSig{mem,1}(:,:,sec)*FibA{mem}*SecComp{mem};
        fsec0{mem}(:,:,sec) = ksec\eye(NumSecDOF);
    end
end

end