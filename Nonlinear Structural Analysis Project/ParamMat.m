function [varargout] = ParamMat(MatModel)

NumArg = numel(MatModel);

if strcmp(MatModel{1},"Bilinear Plastic")||strcmp(MatModel{1},"Bilinear")
    Fy = MatModel{2};
    E = MatModel{3};
    b = MatModel{4};
    varargout{1} = b*E; % Elin
    varargout{2} = (1-b)*E; % Eepp
    varargout{3} = (1-b)*Fy; % Fyepp
    if NumArg==4
        varargout{4} = 2; % NumStrainPts
    elseif NumArg==5
        varargout{4} = MatModel{5}; % NumStrainPts
    end
    
elseif strcmp(MatModel{1},"GMP")
    varargout{1} = MatModel{2}; % Fy
    varargout{2} = MatModel{3}; % E
    varargout{3} = MatModel{4}; % b
    varargout{4} = MatModel{5}; % R0
    varargout{5} = MatModel{6}; % cR1
    varargout{6} = MatModel{7}; % cR2
    if NumArg==7
        varargout{7} = [0,1,0,1]; % a1,a2,a3,a4
        varargout{8} = 2; % NumStrainPts
    elseif NumArg==8
        varargout{7} = [0,1,0,1]; % a1,a2,a3,a4
        varargout{8} = MatModel{8}; % NumStrainPts
    elseif NumArg==11
        varargout{7} = [MatModel{8},MatModel{9},MatModel{10},MatModel{11}]; % a1,a2,a3,a4
        varargout{8} = 2; % NumStrainPts
    elseif NumArg==12
        varargout{7} = [MatModel{8},MatModel{9},MatModel{10},MatModel{11}]; % a1,a2,a3,a4
        varargout{8} = MatModel{12}; % NumStrainPts
    end
    
elseif strcmp(MatModel{1},"Elastic")
    varargout{1} = MatModel{2}; % E
end

end