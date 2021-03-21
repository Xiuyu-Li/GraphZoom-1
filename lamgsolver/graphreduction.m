% graph reduction

% Gs: reduced graph
% H: m-by-n matrix

function [Gs, H, setup] = graphreduction(A, ratio)
	n = length(A);
	A = A-spdiags(sum(A, 2), 0, n, n);
	lamg    = Solvers.newSolver('lamg', 'randomSeed', 1,  'maxDirectSolverSize', floor(n/ratio));

    tStart = tic;
    setup = lamg.setup('laplacian', A);
	tSetup = toc(tStart);
	disp(setup);
    
	setRandomSeed(now);

	i = 3;
	X = setup.level{2}.R';
	lv = length(setup.level);
	while(lv > 2 & i <= lv)
		X = X*setup.level{i}.R';
		i = i+1;
	end

	Gs = setup.level{lv}.A;
	H = X';
 end
