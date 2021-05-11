% Script to follow PARAFAC modeling (w/ notes!)
% NOTE: This is ALL samples (both reservoirs; all sites) for samples with a
% 1 OR 2 dilution
% 11 May 2021, A Hounshell

% Load in samples w/ pfileloader
% DO NOT: divide by fluorescence - will do this later one!
pfileloader3
% Ex = 240:5:450
% Em = 300:2:600
% DO NOT DIVIDE BY FLUORESCENCE

% Visualize EEMs
eemview(MyData,'X',[4 4],[],[],'Labels')
% So much scatter! Hard to visualize samples; remove scatter with smootheem
% below and re-visualize

% Remove Rayleigh and Raman scatter
% Remove 1st order Rayleigh scattering: 15 nm on each side
% Remove 1st order Raman scattering: 10 nm on each side
% Remove 2nd order Rayleigh scattering: 20 nm above; 18 nm below
% Remove 2nd order Raman scattering: 15 nm on each side
% No visualization (visualize below!)
Xs = smootheem(MyData,[15 15],[10 10],[20 18],[15 15],0)

% Re-visualize EEMs with scattering removed
eemview(Xs,'X',[4 4],[],[],'Labels')

% Normalize EEMs (this is why we don't divide by fluorescence above!)
Xnorm = normeem(Xs)

% Visualize normalized EEMs (as a check!)
eemview(Xnorm,'X',[4 4],[],[],'Labels')

% Apply a PARAFAC model!
% Use all data, apply components 2-7, constrain to non-negativity 
Test1 = outliertest(Xnorm,[],2:7,'nonnegativity',[])
% Then look through the pop-up window!
% Sample leverages look pretty good (<0.30 for all models); looks like
% Model 4 might be the best (based on %explained and Core Consistency)

% Look at other metrics
% Sum of squared error - 4 components looks good; 2 & 3 components too low;
% 4-7 all look similar
specsse(Test1,2:7)
specsse(Test1,4:7)

% Compare the spectra for each
comparespectra(Test1,3:6)

% Look at the components
fingerprint(Test1,4)

% 4 components looks good!! (Also the same # of components as for the other models)
% Let's look at module residuals
eemview({Test1 Test1}, {'X','Model4','error_residuals'},[3 3],1)
% Missing some fluorescence in the M-like region (for a few samples), but
% is not ubiquitous (and probably why it wasn't picked up)
% Going to move forward with the 4-component model

% Randomly initiate a 4-component model to verify model results
[LSModel4, Convg4, DSit4] = randinitanal(Xnorm,4,10,'nonnegativity')

% Check model output to verify it meets the previous model
spectralloadings(LSModel4,4)
comparespectra(LSModel4,4)
fingerprint(LSModel4,4)

% Conduct split half-analysis for model validation
S1 = splitds(Xnorm,[],4,'alternating',{[1 2],[3 4], [1 3], [2 4], [1 4], [2 3]})

% Generate 4 component PARAFAC model on each split (this will take
% awhile!!)
A1 = splitanalysis(S1,4,100,'nonnegativity')

% Check the split-half validation - including the 'final' LSModel4 from
% above
splitvalidation(A1,4,[1 2;3 4;5 6],{'AB','CD','AC','BD','AD','BC'},LSModel4)
% All models validated!!

% Reverse normalization
LSModel4r = normeem(LSModel4,'reverse',4)

% Convert raw scores to Fmax
[F4,Scores4] = scores2fmax(LSModel4r,4)

% Export model!
modelout(LSModel4r,4,'20210325_all_1to2Dil.xls')

% Also export out sample names
writecell(MyData.Labels,'20210325_Names.xls')

% Save workspace as: '20210325_all_1to2Dil.mat'
