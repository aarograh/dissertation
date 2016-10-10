close all; clear all; clc;

addpath tests
fid = fopen('testList.txt');
tests = 0;
failed = 0;
EOF = 0;
verbose = false;

nextline = fgets(fid);
ntests = strread(nextline);
for i=1:ntests
    nextline = fgets(fid);
    if ischar(nextline)
        testname = strtrim(strrep(nextline,sprintf('\n'),''));
        fh = str2func(testname);
        tests = tests + 1;
        display(sprintf('Starting test %i/%i: %15s',tests,ntests,testname))
        result = fh(verbose);
        if ~result
            failed = failed + 1;
        end
    else
        EOF = 1;
    end
end

display(sprintf('\n%i of %i tests failed!',failed,tests));
rmpath tests