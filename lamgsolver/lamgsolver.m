% solve equation Ax=b using LAMG

function x=lamgsolver(A, b)
	lamg    = Solvers.newSolver('lamg', 'randomSeed', 1);
    tStart = tic;

	sumA = abs(sum(A));
	kk = find(sumA > 1e-8);
	if(length(kk)>0)
    	setup = lamg.setup('sdd', A);
	else
    	setup = lamg.setup('laplacian', A);
	end
	    	setup = lamg.setup('laplacian', A);
	tSetup = toc(tStart);
    setRandomSeed(now);
	% Turn on debugging printouts during the run
	core.logging.Logger.setLevel('lin.api.AcfComputer', core.logging.LogLevel.DEBUG);
	tStart = tic;
	[x, ~, ~, details] = lamg.solve(setup, b, 'errorReductionTol', 1e-12);
	tSolve = toc(tStart);
	%disp(setup.setupLaplacian);
	tMvm    = mvmTime(A, 5);
	nnz     = numel(nonzeros(A));
	fprintf('|A*x-b|/|b|: %.2e\n', norm(A*x-b)/norm(b));
 end
