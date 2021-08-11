function EtanSig = MatBilinear(NumFib,MatModel,mem,sec,e,FibState,EtanSig)

[Elin,Eepp,Fyepp,NumStrainPts] = ParamMat(MatModel);

for fib = 1:NumFib
    eps = linspace(FibState{mem}(fib,2,sec),e{mem}(fib,sec),NumStrainPts);
    
    sigma = zeros(1,NumStrainPts);
    sigma_lin = zeros(1,NumStrainPts);
    sigma_epp = zeros(1,NumStrainPts);

    sigma(1) = FibState{mem}(fib,3,sec);
    sigma_lin(1) = FibState{mem}(fib,4,sec);
    sigma_epp(1) = FibState{mem}(fib,5,sec);

    for i = 2:NumStrainPts
            deps = eps(i)-eps(i-1); % CHANGE IN ELEMENT STRAIN   
            sigma_lin(i) = sigma_lin(i-1)+Elin*deps;
            sigmaepp_trial = sigma_epp(i-1)+Eepp*deps;
            if abs(sigmaepp_trial)>Fyepp
                sigma_epp(i) = Fyepp*sign(sigmaepp_trial);
                E_ins = Elin;
            else
                sigma_epp(i) = sigmaepp_trial;
                E_ins = Elin+Eepp;
            end
            sigma(i) = sigma_lin(i)+sigma_epp(i); % TOTAL STRESS FOR PARALLEL SPRINGS
    end

    EtanSig{mem,1}(fib,fib,sec) = E_ins; % storing Et
    EtanSig{mem,2}(fib,1,sec) = sigma(end); % storing sigma
    EtanSig{mem,2}(fib,2,sec) = sigma_lin(end); % storing sigma_lin
    EtanSig{mem,2}(fib,3,sec) = sigma_epp(end); % storing sigma_epp
end
end