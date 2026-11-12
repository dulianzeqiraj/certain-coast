function ok = run_tests()
%RUN_TESTS Run every test of the CERTAIN-COAST toolbox, in MATLAB or GNU Octave.
%   run_tests runs the test files tests/test_*.m, prints one line per test and raises an error if any test fails
%   (so that continuous integration stops). OK = run_tests returns true when all tests pass, without the error.

%   Copyright 2026 Dulian Zeqiraj. Licensed under the Apache License, Version 2.0.

root = fileparts(mfilename('fullpath'));
addpath(root);
addpath(fullfile(root, 'tests'));
files = dir(fullfile(root, 'tests', 'test_*.m'));
nPass = 0;
failed = {};
fprintf('CERTAIN-COAST %s\n', certaincoast.version());
for i = 1:numel(files)
    [~, name] = fileparts(files(i).name);
    tests = feval(name);
    for j = 1:numel(tests)
        label = sprintf('%s/%s', name, func2str(tests{j}));
        t0 = tic;
        try
            tests{j}();
            nPass = nPass + 1;
            fprintf('  pass  %-70s %6.1f s\n', label, toc(t0));
        catch err
            failed{end + 1} = label; %#ok<AGROW>
            fprintf('  FAIL  %-70s %s\n', label, err.message);
        end
    end
end
fprintf('%d passed, %d failed\n', nPass, numel(failed));
ok = isempty(failed);
if ~ok && nargout == 0
    error('certaincoast:tests', '%d test(s) failed: %s', numel(failed), strjoin(failed, ', '));
end
end
