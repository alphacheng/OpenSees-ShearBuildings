function [natFreqs,modeShapes] = modalAnalysis(bldg,spring)

m = bldg.storyMass;
k = spring.designStiffness;

M = diag(m);

K = zeros(bldg.nStories);
for i = 1:bldg.nStories
    for j = 1:bldg.nStories
        if     i == j && i ~= bldg.nStories
            K(i,j) = k(i) + k(i+1);
            K(i,j+1) = -k(i+1);
        elseif i == j && i == bldg.nStories
            K(i,j) = k(i);
        end
    end
end
K = K + K' - diag(diag(K)); % Mirror K along the diagonal

sys.mass = M;
sys.stiffness = K;
sys.nDOF = bldg.nStories;

syms wsym
natFreqs = solve(det(-wsym^2*sys.mass + sys.stiffness));
natFreqs = sort(real(double(vpa(natFreqs(vpa(natFreqs)>0)))));

modeShapes = zeros(sys.nDOF);
for i = 1:sys.nDOF
    eigmat = -natFreqs(i)^2*sys.mass + sys.stiffness;  % Set up
    usablematrix = eigmat(2:end,:);     % Remove first equation
    A = usablematrix(:,2:end);          % Extract non-constants
    B = -usablematrix(:,1);             % Move constants
    C = A^-1 * B;                       % Calculate mode shapes
    D = [1 ; C];                        % Include assumed mode shape
    D = D ./ max(abs(D));               % Scale to 1
    modeShapes(:,i) = D;                % Stick into full matrix
end

end
