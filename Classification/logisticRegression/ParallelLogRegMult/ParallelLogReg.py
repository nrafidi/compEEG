import os, QSubHelpers

def writeDataSizes(
  srcRoot, 
  jobDir,
  outputDir,
  jobNameSuffixFmt,
  dataLoaderFmt,
  matlabDeclarations,
  parameters):
  
  if (matlabDeclarations is None):
    matlabDeclarations = list();
    
  logRegRoot = os.path.join(srcRoot, 'shared', 'logisticRegression');
  addpathFmt = 'addpath(genpath(\'{0}\'));'
  
  matlabDeclarations.append(addpathFmt.format(logRegRoot));
  
  jobNameFmt = 'plllogreg_getsz_' + jobNameSuffixFmt;
  
  writeSizeFuncFmt = 'WriteSize( {0}, \'{1}\' );'.format(
    dataLoaderFmt, 
    os.path.join(outputDir, 'output_' + jobNameFmt + '.txt'));
  
  writeDataSizeJobDict = QSubHelpers.createGridJobs(
    jobDir,
    outputDir,
    jobNameFmt,
    writeSizeFuncFmt,
    matlabDeclarations,
    parameters);
    
  QSubHelpers.runScripts(
    jobDir,
    outputDir,
    writeDataSizeJobDict,
    100,
    None,
    None);
    
  return os.path.join(outputDir, 'output_' + jobNameFmt + '.txt');

def readDataSize( outputPath ):
  with open(outputPath, 'r') as sizeFile:
    for line in sizeFile:
      sizeStrings = line.split('\t');
      sizes = list();
      for sizeStr in sizeStrings:
        sizes.append(long(sizeStr));
      return sizes;

def createLeaveOneOutJobs(
  srcRoot,
  jobDir,
  jobNameFmt,
  finalJobName,
  outputDir,
  numObservations,
  dataLoader,
  lambdas,
  shouldRegBias,
  matlabDeclarations):
  
  jobNamePart = jobNameFmt.format('NA', 'NA');
  
  foldIndices = range(1, numObservations + 1);
  
  if (matlabDeclarations is None):
    matlabDeclarations = list();
    
  logRegRoot = os.path.join(srcRoot, 'shared', 'logisticRegression');
  addpathFmt = 'addpath(genpath(\'{0}\'));'
  
  matlabDeclarations.append(addpathFmt.format(logRegRoot));
  
  if (shouldRegBias):
    shouldRegBiasStr = 'true';
  else:
    shouldRegBiasStr = 'false';
  
  cvFuncFmt = 'LogRegSingleLambdaSingleFoldCVFromLoader({0}, {1}, {2}, {3}, {4})'.format(
    dataLoader,
    '{0}',
    numObservations,
    '{1}',
    shouldRegBiasStr);
    
  looJobDictionary = QSubHelpers.createGridJobs(
    jobDir,
    outputDir,
    jobNameFmt,
    cvFuncFmt,
    matlabDeclarations,
    [ lambdas, foldIndices ],
    50);
  
  mergeJobDictionary = QSubHelpers.createMergeAllJobOutputsJob(
    jobDir,
    outputDir,
    outputDir,
    jobNameFmt,
    'loo_mrg_' + jobNamePart,
    [ lambdas, foldIndices ]);
    
  finalizeFunc = 'LogRegChooseLambdaAndTrainAll({0}, {1}, \'{2}\')'.format(
    dataLoader,
    shouldRegBiasStr,
    os.path.join(outputDir, 'output_loo_mrg_' + jobNamePart + '.mat'));
    
  finalizeJobDictionary = QSubHelpers.createGridJobs(
    jobDir,
    outputDir,
    finalJobName,
    finalizeFunc,
    matlabDeclarations,
    ());
    
  finalJobManifestName = 'looLogRegFinalJobManifest_' + jobNamePart + '.txt';
  return [looJobDictionary, mergeJobDictionary, finalizeJobDictionary];
  
def runLeaveOneOutJobs(
  jobDir,
  outputDir,
  jobDictionaries,
  maxSimulataneousJobs,
  resourceStr,
  pool):
  
  for jobDictionary in jobDictionaries:
    QSubHelpers.runScripts(
      jobDir,
      outputDir,
      jobDictionary,
      maxSimulataneousJobs,
      resourceStr,
      pool);
       
def createLeaveOneOutJobsForEachParameterCombination(
  srcRoot,
  jobDir,
  outputDir,
  jobNameFmt,
  dataLoaderFmt,
  parameters,
  lambdas,
  shouldRegBias,
  matlabDeclarations):

  combinedJobDictionaries = list();
  
  formatArgsForFinalJob = list();
  for (indexParameterList, parameterList) in enumerate(parameters):
    formatArgsForFinalJob.append('{' + '{0}'.format(indexParameterList) + '}');
  formatArgsForFinalJob.append('NA');
  formatArgsForFinalJob.append('NA');
  
  jobNameSuffixFmt = jobNameFmt.format(*tuple(formatArgsForFinalJob));
  finalJobNameFormat = 'fnl_loo_' + jobNameSuffixFmt;
  
  dataSizeOutputFmt = writeDataSizes(
    srcRoot,
    jobDir,
    outputDir,
    jobNameSuffixFmt,
    dataLoaderFmt,
    matlabDeclarations,
    parameters);
  
  combinedDictionaries = None;
  
  for parameterCombination in QSubHelpers.parameterCombinationIterator(parameters):
    dataLoader = dataLoaderFmt.format(*parameterCombination);
    
    parameterCombinationWithPlaceholders = list();
    for parameter in parameterCombination:
      parameterCombinationWithPlaceholders.append(parameter);
    parameterCombinationWithPlaceholders.append('{0}');
    parameterCombinationWithPlaceholders.append('{1}');
    
    jobNameFmtPartial = jobNameFmt.format(*tuple(parameterCombinationWithPlaceholders));
    
    print ('Creating jobs for parameter combination: ' + str(parameterCombination));
    
    dataSize = readDataSize( dataSizeOutputFmt.format(*parameterCombination) );
    numObservations = dataSize[0];
    
    print ('num observations: {0}'.format(numObservations));
    
    jobDictionaries = createLeaveOneOutJobs(
      srcRoot,
      jobDir,
      jobNameFmtPartial,
      finalJobNameFormat.format(*parameterCombination),
      outputDir,
      numObservations,
      dataLoader,
      lambdas,
      shouldRegBias,
      matlabDeclarations);
      
    if (not combinedDictionaries):
      combinedDictionaries = jobDictionaries;
    else:
      for (indexDictionary, jobDictionary) in enumerate(jobDictionaries):
        for (jobName, outputs) in jobDictionary.items():
          combinedDictionaries[indexDictionary][jobName] = outputs;
  
  return (combinedDictionaries, finalJobNameFormat);