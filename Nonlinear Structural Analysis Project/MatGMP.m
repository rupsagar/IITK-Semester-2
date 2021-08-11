function EtanSig = MatGMP(NumFib,MatModel,mem,sec,e,FibState,EtanSig)

[Fy,E,b,R0,cR1,cR2,a,NumStrainPts] = ParamMat(MatModel);
eps_y = Fy/E;
a1 = cR1*R0;

for fib = 1:NumFib
    eps = linspace(FibState{mem}(fib,2,sec),e{mem}(fib,sec),NumStrainPts);
    
    sigma = zeros(1,NumStrainPts);
    d_eps = zeros(1,NumStrainPts);
    eps_xi = zeros(1,2);
    eps_max = zeros(1,2);
    
    sigma(1) = FibState{mem}(fib,3,sec);
    d_eps(1) = FibState{mem}(fib,4,sec);
    eps_r = FibState{mem}(fib,5,sec);
    sigma_r = FibState{mem}(fib,6,sec);
    eps_0 = FibState{mem}(fib,7,sec);
    sigma_0 = FibState{mem}(fib,8,sec);
    eps_xi(1) = FibState{mem}(fib,9,sec);
    eps_xi(2) = FibState{mem}(fib,10,sec);
    R = FibState{mem}(fib,11,sec);
    eps_max(1) = FibState{mem}(fib,12,sec);
    eps_max(2) = FibState{mem}(fib,13,sec);
    
    for i = 2:NumStrainPts
        d_eps(i) = eps(i)-eps(i-1);

        % CHECK FOR STRAIN REVERSAL AND IF CURRENT STRAIN INCREMENT IS NON-ZERO
        if sign(d_eps(i-1)*d_eps(i))~=1 && d_eps(i)~=0
            eps_r = eps(i-1);
            sigma_r = sigma(i-1);

            %%% MODIFICATION FOR ISOTROPIC HARDENING
            % CALCULATE MAXIMUM STRAIN
            if sigma(i-1)~=0 || eps(i-1)~=0
                if d_eps(i)<0
                    ct = 1; % COMPRESSION VALUES OF a TO BE TAKEN
                    if eps_r>eps_max(2)
                        eps_max(2) = eps_r;
                    end
                elseif d_eps(i)>0 && (sigma(i-1)~=0 && eps(i-1)~=0)
                    ct = 2; % TENSION VALUES OF a TO BE TAKEN
                    if eps_r<eps_max(1)
                        eps_max(1) = eps_r;
                    end
                end
                eps_diff = (eps_max(2)-eps_max(1))/(2*a(2*ct)*eps_y);
                shift = 1+a(2*ct-1)*eps_diff^0.8;    
            elseif sigma(i-1)==0 && eps(i-1)==0 % NO SHIFT FOR INITIAL POINT 
                shift = 1;
            end
            
            %%% WITHOUT ISOTROPIC HARDENING (USEFUL FOR DEBUGGING - COMMENT ABOVE PART)
%             shift = 1;
             
            % UPDATE YIELD STRESS BASED ON MAXIMUM STRAIN
            Fy_new = Fy*shift;
            eps_y_new = eps_y*shift;
            
            % INTERSECTION OF ASYMPTOTES
            Sgn_d_eps = sign(d_eps(i));
            eps_0 = ((eps_r-b*eps_y_new*Sgn_d_eps)-(sigma_r-Fy_new*Sgn_d_eps)/E)/(1-b);
            sigma_0 = (b*E*(eps_r-eps_y_new*Sgn_d_eps)-(b*sigma_r-Fy_new*Sgn_d_eps))/(1-b);
            
            % UPDATE MAXIMUM NORMALIZED STRAIN REACHED IN SAME SIGN AS EPS_R AND XI
            Sgn_eps_r = sign(eps_r)==sign(eps_xi);
            if abs(eps_r/eps_y)>abs(eps_xi(Sgn_eps_r))
                eps_xi(Sgn_eps_r) = eps_r/eps_y;
            end
            Sgn_sigma_0 = Sgn_d_eps==sign(eps_xi);
            xi = abs(eps_xi(Sgn_sigma_0)-eps_0/eps_y);
            
            R = R0-a1*xi/(cR2+xi);
        end
        
        eps_star = (eps(i)-eps_r)/(eps_0-eps_r);
        sigma_star = b*eps_star+(1-b)*eps_star/((1+eps_star^R)^(1/R));        
        sigma(i) = sigma_star*(sigma_0-sigma_r)+sigma_r;
        E_ins = b*E+(1-b)*E/(1+eps_star^R)^(1+1/R);
    end
    
    EtanSig{mem,1}(fib,fib,sec) = E_ins; % storing Et
    EtanSig{mem,2}(fib,1,sec) = sigma(end); % storing sigma
    EtanSig{mem,2}(fib,2,sec) = d_eps(end); % storing d_eps
    EtanSig{mem,2}(fib,3,sec) = eps_r; % storing eps_r
    EtanSig{mem,2}(fib,4,sec) = sigma_r; % storing sigma_r
    EtanSig{mem,2}(fib,5,sec) = eps_0; % storing eps_0
    EtanSig{mem,2}(fib,6,sec) = sigma_0; % storing sigma_0
    EtanSig{mem,2}(fib,7,sec) = eps_xi(1); % storing eps_xi(1)
    EtanSig{mem,2}(fib,8,sec) = eps_xi(2); % storing eps_xi(2)
    EtanSig{mem,2}(fib,9,sec) = R; % storing R
    EtanSig{mem,2}(fib,10,sec) = eps_max(1); % storing eps_max(1)
    EtanSig{mem,2}(fib,11,sec) = eps_max(2); % storing eps_max(2)
    
end

end