% inspired by  https://rse.shef.ac.uk/blog/2022-05-05-concise-guide-to-reproducible-matlab/

% ├── data/             # directory to contain your data
% │   └── raw/          # all your raw, unprocessed data
% │   └── processed/    # any data that you have altered
% ├── figures/          # keep the figures that you produce from your code here
% ├── reports/          # a clear folder for your papers or reports for the project
% ├── src/              # source code
% │   └── @MyClass/     # a class directory
% │   └── functions/    # a directory containing authored functions or classes with some commonality
% │   └── utils/        # a directory containing utility functions
% │   └── reqs/         # a directory containing MATHWORKS etc user submitted functions
% ├── docs/             # documentation
% ├── tests/            # software tests for the project
% ├── README            # readme file! (essential)
% └── my.prj            # MATLAB project file (read on for more details)

mkdir(['data' filesep 'raw'])
mkdir(['data' filesep 'processed'])
mkdir figures
mkdir reports      
mkdir(['src' filesep 'utils'])
mkdir(['src' filesep 'functions'])
mkdir(['src' filesep 'req'])
mkdir(['src' filesep 'utils'])      
mkdir docs            
mkdir tests 

addpath(genpath('data'))
addpath(genpath('src'))
