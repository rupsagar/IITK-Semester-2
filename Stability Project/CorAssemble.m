function [Ks,Ri] = CorAssemble(ElementType,MemCon,NodeXYk,Thetak,L0,Beta0,NumMem,NumTotalDOF,...
    NumStrDOF,NumMemDOF,MemDOF,KBCMap,BCMap,E,Area,MOI)

L = zeros(NumMem,1);
Q = zeros(NumMemDOF,NumMem);
Rg = zeros(NumTotalDOF,1);
Kg = zeros(NumTotalDOF);
Ks = zeros(NumStrDOF);

if strcmp(ElementType,'Truss') % FOR TRUSS ELEMENT
    P = zeros(NumMem,1); 
    T = zeros(NumMemDOF,NumMemDOF,NumMem);
    for mem = 1:NumMem
        MemNodes = MemCon(mem,:);
        Node1 = NodeXYk([2*MemNodes(1)-1,2*MemNodes(1)]);
        Node2 = NodeXYk([2*MemNodes(2)-1,2*MemNodes(2)]);
        ri = Node2-Node1;
        L(mem) = norm(ri);
        CS = ri/L(mem);
        t = [CS(1),CS(2),0;
            -CS(2),CS(1),0;
            0,0,1];
        T(:,:,mem) = [t,zeros(3);zeros(3),t];
        P(mem) = E(mem)*Area(mem)*(L(mem)-L0(mem))/L0(mem);
        
        
        KMLocal = E(mem)*Area(mem)/L0(mem)*[1,0,0,-1,0,0;
                                            0,0,0,0,0,0;
                                            0,0,0,0,0,0;
                                            -1,0,0,1,0,0;
                                            0,0,0,0,0,0;
                                            0,0,0,0,0,0];
        KGLocal = P(mem)/L(mem)*[1,0,0,-1,0,0;
                                 0,1,0,0,-1,0;
                                 0,0,0,0,0,0;
                                 -1,0,0,1,0,0;
                                 0,-1,0,0,1,0;
                                 0,0,0,0,0,0];
        
        Kt = T(:,:,mem)'*(KMLocal+KGLocal)*T(:,:,mem);
        Q(:,mem) = T(:,:,mem)'*P(mem)*[-1,0,0,1,0,0]';
        
        Rg(MemDOF(mem,:)) = Rg(MemDOF(mem,:))+Q(:,mem);
        for p = 1:NumMemDOF
            for q = 1:NumMemDOF
                Kg(MemDOF(mem,q),MemDOF(mem,p)) = Kg(MemDOF(mem,q),MemDOF(mem,p))+Kt(q,p);
            end
        end
    end
elseif strcmp(ElementType,'Frame') % FOR FRAME ELEMENT
    betai = zeros(NumMem,1);
    for mem = 1:NumMem
        MemNode = MemCon(mem,:);
        Node1XY = NodeXYk([2*MemNode(1)-1,2*MemNode(1)]);
        Node2XY = NodeXYk([2*MemNode(2)-1,2*MemNode(2)]);
        ri = Node2XY-Node1XY;
        L(mem) = norm(ri);
        uL = L(mem)-L0(mem);

        betai(mem) = atan(ri(2)/ri(1));
        Theta = Thetak(MemNode);
        beta = Theta+Beta0(mem);
        num = @(j)cos(betai(mem))*sin(beta(j))-sin(betai(mem))*cos(beta(j));
        den = @(j)cos(betai(mem))*cos(beta(j))+sin(betai(mem))*sin(beta(j));
        ThetaL = @(j)atan(num(j)/den(j));

        dv = [uL;ThetaL(1);ThetaL(2)];

        Kv = [E(mem)*Area(mem)/L0(mem),0,0;
              0,4*E(mem)*MOI(mem)/L0(mem),2*E(mem)*MOI(mem)/L0(mem);
              0,2*E(mem)*MOI(mem)/L0(mem),4*E(mem)*MOI(mem)/L0(mem)];
        S = Kv*dv;

        CS = ri/L(mem);
        z = [CS(2),-CS(1),0,-CS(2),CS(1),0]';
        r = [-CS(1),-CS(2),0,CS(1),CS(2),0]';
        A = ([0,0,1,0,0,0;0,0,0,0,0,1]-1/L(mem)*[z';z'])';
        B = [r';A'];  

        Kt = B'*Kv*B+S(1)/L(mem)*(z*z')+(S(2)+S(3))/L(mem)^2*(r*z'+z*r');
        Q(:,mem) = B'*S;

        Rg(MemDOF(mem,:)) = Rg(MemDOF(mem,:))+Q(:,mem);
        for p = 1:NumMemDOF
            for q = 1:NumMemDOF
                Kg(MemDOF(mem,q),MemDOF(mem,p)) = Kg(MemDOF(mem,q),MemDOF(mem,p))+Kt(q,p);
            end
        end   
    end
end

Ks(:) = Kg(KBCMap);
Ri = Rg(BCMap);
end