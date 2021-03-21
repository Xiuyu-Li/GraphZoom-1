% graph reduction

%input:
% mtx file:  input graph Laplacian matrix

%output:
% Gs:        reduced graph Laplacian matrix
% setup:     data structure for storing all information of reduction
% 		     setup.level{i}.R: mapping operator between level i and i-1

% graph reduction

%input:
% mtx file:  input graph Laplacian matrix

%output:
% Gs:        reduced graph Laplacian matrix
% setup:     data structure for storing all information of reduction
% 		     setup.level{i}.R: mapping operator between level i and i-1

% function LamgSetup()
function LamgSetup(GraphPath, yPath, ReductionRatio, tvNum, Fusion, useLabel, SavePath)
    fprintf('Loading Graph to be Reduced......\n');
    ReductionRatio = str2num(ReductionRatio);
    lda = 0.1;                                         % self_loop
    kpower = 2;                                        % power of graph filter
    tvNum = str2num(tvNum);
    
    % tvNum = 4;
    % ReductionRatio = 2;
    % GraphPath = 'cora.mtx';
    % yPath = 'cora_y.mat';
    % Fusion = 'n';
    % useLabel = 1;
    % SavePath = 'save';
    
    y_mat = load(yPath);
    y0 = double(y_mat.data);                           % groundtruth labels
    y0 = kron(y0, ones(1, tvNum));
    useLabel = str2num(useLabel);

    fp = fopen(GraphPath, 'r');
    B = textscan(fp, '%d %d %f', 'headerlines', 3);    % skip first two rows in mtx file
    row = cell2mat(B(1));
    col = cell2mat(B(2));
    val = cell2mat(B(3));
    fclose(fp);

    A = sparse(double(row), double(col), double(val));
    if(~issymmetric(A))
        A = A'+A-diag(diag(A));
    end
    n = length(A);	
    fprintf('###### Running LamgSetup ######\n');
    t = cputime;
    
    if (useLabel~=1)
        lamg  = Solvers.newSolver('lamg', 'randomSeed', 1,  'maxDirectSolverSize', floor(n/ReductionRatio), 'lda', lda, 'kpower', kpower,...
            'tvNum', tvNum);
    else
        y_size = size(y0);
        tvNum = y_size(2);
        tvMax = y_size(2)
        lamg  = Solvers.newSolver('lamg', 'randomSeed', 1,  'maxDirectSolverSize', floor(n/ReductionRatio), 'lda', lda, 'kpower', kpower,...
            'y0', y0, 'tvNum', tvNum, 'tvMax', tvMax, 'useLabel', true);
    end

    %tStart = tic;
    setup = lamg.setup('laplacian', A);
    %tSetup = toc(tStart);
    disp(setup)

    %setRandomSeed(now);
    
    lv = length(setup.level);
    assert(lv>=2, 'ERROR: Reduction ratio is too small, plese try ReductionRatio > 2 !!!!!!\n');
    X = setup.level{2}.R; % R is m-by-n
    if Fusion~='f'
        writeMtx(X, nnz(X), strcat(SavePath,'/Projection_1.mtx'));
        A = setup.level{2}.R*A*setup.level{2}.R';
    end
    
    i = 3;
    while(lv > 2 & i <= lv)
        X = setup.level{i}.R * X;
        if Fusion~='f'
            A = setup.level{i}.R*A*setup.level{i}.R';
            writeMtx(setup.level{i}.R, nnz(setup.level{i}.R), strcat(SavePath,'/Projection_',num2str(i-1),'.mtx'));
        end
        i = i+1;
    end
    cpu_time = cputime - t;
    writeMtx(X, nnz(X), strcat(SavePath,'/Mapping.mtx'));
    if Fusion~='f'
        dlmwrite(sprintf(strcat(SavePath,'/NumLevels.txt')), lv);
        writematrix(A, nnz(A), strcat(SavePath,'/Gs.mtx'));
    end
    dlmwrite(sprintf(strcat(SavePath,'/CPUtime.txt')), cpu_time);
end
