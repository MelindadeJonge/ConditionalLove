function[] = run_static_model(id,repN,niter,thinning,chains)
addpath(fullfile('MATLAB','HMSC_Class'))
addpath('data')
addpath('output')

%% Import & prepare data
all_occurrences = importdata(fullfile('data','Occurrences.csv'));
all_covariates = importdata(fullfile('data','Covariates.csv'));
covariate_names = all_covariates.colheaders(2:end);
species_names = all_occurrences.colheaders(2:end);
all_occurrences = all_occurrences.data;
all_covariates = all_covariates.data;

% Select only plot per environmental grid cell
plot_combinations = importdata(fullfile('data','UniquePlotCombinations.csv'));
plot_combinations = plot_combinations.data(:,id);
selection = ismember(all_occurrences(:,1),plot_combinations);
occurrences_selected = all_occurrences(selection,2:end);
covariates_selected = all_covariates(selection,2:end);

% Scaling of the covariates and adding intercept and quadratic terms
covariates_mean = mean(covariates_selected);
covariates_centered = covariates_selected - covariates_mean; 
covariates_sd = std(covariates_centered);
covariates_scaled = covariates_centered./covariates_sd;
covariates_squared = covariates_scaled.^2;
[n, ~] = size(covariates_scaled);
covariates = [ones(n, 1), covariates_scaled, covariates_squared];
covariate_names2 = strcat(covariate_names,'2');
covariate_names =  ['intercept',covariate_names,covariate_names2];

% Preparing random effect for CWD dependent latent loading
[n, nc] = size(covariates);
piCell = num2cell(1:n)';
piCell = cellfun(@num2str,piCell,'Uniformoutput',false);
model_folder = fullfile('output',strcat('StaticModel',int2str(id)));
mkdir(model_folder);
save(fullfile(model_folder,'DataUsed.mat'),'thinning','repN','niter','chains','species_names','covariate_names','occurrences_selected','covariates','piCell')

%% Run model, parallel chains
parfor i = 1:chains
   folder = fullfile(model_folder,strcat('Chain',int2str(i)));
   mkdir(folder);
    m = Hmsc(folder, false, false, false, 0, 0, false, [false], [false]);
    m.setData(occurences_selected,'probit',covariates,piCell,[], [], [], [], [], []);
    m.setCovNames(covariate_names);
    m.setSpeciesNames(species_names);
    covScaleFlag = [2,zeros(1,(nc-1))];
    m.setCovScaling(covScaleFlag); %Scaling 0: no scaling, 2: intercept, 1: normal scaling.
    m.setPriorsDefault();
    
    m.setMCMCOptions(niter, thinning); %(draws,thinning) draws * thinning iterations per run. for each run X draws are sampled from the posterior with thinning of X.
    m.setMCMCAdapt([10,0], false, [0,0], false); %This influences the runs (10 in this case) in which the number of latent variables still adapted.
    m.setMCMCSaveOptions(false, true); %save to (RAM,file)
    m.sampleMCMC(repN, false, 0, [], 2); %sampleMCMC(m, nRuns, append, startPar,lastRun, verbose)
    m.postRamClear();
end
end
