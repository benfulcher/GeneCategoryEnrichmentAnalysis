function ermineJResults = RunErmineJ(inputFileName,numIterations)
% Runs ermineJ on an input file:

% ---Info on running ermineJ in the commandline---:
% cf. http://erminej.chibi.ubc.ca/help/tutorials/erminej-cli/
%
% -a: sets the annotation file
% -b: sets bigger is better for gene scores
% -c: sets the gene set (class file): i.e., GO XML
% -d sets the data directory
% -e sets the column for the scores in your gene score file
% -g sets how to deal with scores for replicate genes
% -i number of iterations
% -j include gene symbols for all gene sets
% -M multiple test correction method
% -m method for computing raw class statistics for GSR
% -n method for computing gene set significance
% -o output filename
% -s gene score file
% -x maximum class size (100)
% -y minimum class size

if nargin < 2
    numIterations = 20000;
end


% Set the ermineJ home directory:
% homeDir = '/home/benfulcher/Monash076/Ben/MouseConnectome/ermineJ/ermineJ-3.0.2';
% [status,cmdOut] = system('ERMINEJ_HOME=/Users/benfulcher/Downloads/ermineJ-3.0.2');

% Get settings for writing out to file:
shellScriptPath = '/Users/benfulcher/Downloads/ermineJ-3.0.2/bin/ermineJ.sh';
propertiesFilePath = which('ermineJBP.properties');
fprintf(1,'Using ermineJ settings in %s.\n',propertiesFilePath);
inputFilePath = which(inputFileName);
if isempty(inputFilePath)
    error('%s not found',inputFileName);
end
outputFile = fullfile(pwd,'tmp_ermineJ.out');

% Set the JAVA_HOME variable:
setenv('JAVA_HOME','/Library/Java/JavaVirtualMachines/jdk1.8.0_73.jdk/Contents/Home');

% Delete the temp file if it already exists:
if exist(outputFile,'file')==2
    delete(outputFile);
end

% Construct the command:
command = sprintf('%s -C %s -i %u -b -j -s %s -o %s',shellScriptPath,...
                    propertiesFilePath,numIterations,inputFilePath,outputFile);

%-------------------------------------------------------------------------------
% Execute the command:
%-------------------------------------------------------------------------------
fprintf(1,'Running ermineJ on %s for biological processes using %u iterations...',...
                                            inputFileName,numIterations);
[status,cmdOut] = system(command);
fprintf(1,' Done.\n');

%-------------------------------------------------------------------------------
% Read in the data:
%-------------------------------------------------------------------------------
if ~exist(outputFile,'file')
    warning('ermineJ didn''t write output??')
    keyboard
end
ermineJResults = ReadInErmineJ(outputFile);

% Write out to screen:
fisSig_005 = find(ermineJResults.pVal_corr < 0.05);
fisSig_01 = find(ermineJResults.pVal_corr>=0.05 & ermineJResults.pVal_corr < 0.1);
numSig_005 = length(fisSig_005);
numSig_01 = length(fisSig_01);
if numSig_01==0 && numSig_005==0
    fprintf(1,'No significant GO enrichment for %s at q < 0.1\n',inputFileName);
else
    fprintf(1,'---%s---\n',inputFileName);
    % The FDR 0.05 club:
    for i = 1:numSig_005
        fprintf(1,'%s (%s) [p=%.2g]\n',ermineJResults.GOName{fisSig_005(i)},...
                    ermineJResults.GOID{fisSig_005(i)},ermineJResults.pVal_corr(fisSig_005(i)));
    end
    % What about outside contenders:
    if numSig_01 > 0
        fprintf(1,'--------\n');
        for i = 1:numSig_01
            fprintf(1,'%s (%s) [p=%.2g]\n',ermineJResults.GOName{fisSig_01(i)},...
                        ermineJResults.GOID{fisSig_01(i)},ermineJResults.pVal_corr(fisSig_01(i)));
        end
    end
end

% $ERMINEJ_HOME/bin/ermineJ.sh -C $ERMINEJ_HOME/ermineJBP.properties -b -j -s
% $ERMINEJ_HOME/ermineJ.data/ermineJInputFile_GGblock_17349genes_k44_norm_energy_reciprocalConnected_distance_alldiv-expFitAll_tstat.txt
% -o ermineJ_reciprocalConnectedBP.out
