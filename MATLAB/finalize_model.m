function[] = finalize_model(model,id,burnin,repN,chains,niter,thinning,extrathinning)
addpath(fullfile('MATLAB','HMSC_Class'))
addpath('output')
model_folder = fullfile('output',strcat(model,'Model',int2str(id)));
load(fullfile(model_folder,'DataUsed.mat'));
[~, nc] = size(covariates);
for i=1:chains
    folder = fullfile(model_folder,strcat('Chain',int2str(i)));
    if strcmp(model,'Conditional')
        m = Hmsc(folder, false, false, false, 0, 0, false, [false], [true]);
        m.setData(occurrences_selected,'probit',covariates,piCell,[], [], [], XrCell, [], []);
        m.setFactorCovScaling({[2,0]})
    elseif strcmp(model,'Static')
        m = Hmsc(folder, false, false, false, 0, 0, false, [false], [false]);
        m.setData(occurrences_selected,'probit',covariates,piCell,[], [], [], [], [], []);
    end        
    m.setCovNames(covariate_names);
    m.setSpeciesNames(species_names);
    covScaleFlag = [2,zeros(1,(nc-1))];
    m.setCovScaling(covScaleFlag); %Scaling 0: no scaling, 2: intercept, 1: normal scaling.
    m.setPriorsDefault();
    
    m.setMCMCOptions(niter, thinning); %(draws,thinning) draws * thinning iterations per run. for each run X draws are sampled from the posterior with thinning of X.
    m.setMCMCAdapt([10,0], false, [0,0], false); %This influences the runs (10 in this case) in which the number of latent variables still adapted.
    for j=(burnin+1):repN
        m.postFileToRam(j);
    end
    m.setPostThinning((burnin+1):repN,extrathinning);
    m.postRamClear();
    save(fullfile(folder,strcat(model,int2str(id),'Chain',int2str(i),'.mat')),'-v7.3')
    mkdir(fullfile(folder,'Correlations'));
    
    % Calculate and save species associations to csv file(s)
    if strcmp(model,'Conditional')
        Xr=m.Xr{1};
        Xr=Xr(:,2);
        XrGrad = min(Xr):0.1:max(Xr);
        XrGrad=[ones(61,1),XrGrad'];
        
        % Residual correlations along the gradient
        for grad = 1:length(XrGrad)
            [correlations, combination] = m.computeCorrelationsFullP(1,XrGrad(grad,:));                      
            csvwrite(fullfile(folder,'Correlations',strcat('Correlations',num2str(grad),'.csv')),correlations);
        end
        
        % Resdiual correlations at the 5th, 25th, 50th 75th and 95th
        % percentiles of the gradient
        CWDquants = quantile(m.Xr{1},[0.05,0.25,0.5,0.75,0.95]);
        for grad = 1:length(CWDquants)
            [correlations, combination] = m.computeCorrelationsFullP(1,CWDquants(grad,:));                      
            csvwrite(fullfile(folder,'Correlations',strcat('QuantileCorrelations',num2str(grad),'.csv')),correlations);
            % Extract parameters for effective sample size estimation
            Omega = m.getPostOmega(1,CWDquants(grad,:));
            Omega = Omega{1};
            csvwrite(fullfile(folder,strcat('Omega',num2str(grad),'.csv')),Omega);

        end

        csvwrite(fullfile(folder,'Correlations','XrGrad.csv'),XrGrad);
        csvwrite(fullfile(folder,'Correlations','combinations.csv'),combination);
        
        % Extract parameters for effective sample size estimation
        Beta = m.getPostBeta();
        Beta = Beta{1};
        csvwrite(fullfile(folder,'Beta.csv'),Beta);

    elseif strcmp(model,'Static')
        [correlations,combinations] = m.computeCorrelationsFullP(1,[]);
        csvwrite(fullfile(folder,'Correlations','Correlations.csv'),correlations);
        csvwrite(fullfile(folder,'Correlations','SpeciesCombinations.csv'),combinations)
        
        % Extract parameters for effective sample size estimation
        Omega = m.getPostOmega(1,[]);
        Omega = Omega{1};
        csvwrite(fullfile(folder,strcat('Omega',num2str(grad),'.csv')),Omega);
        Beta = m.getPostBeta();
        Beta = Beta{1};
        csvwrite(fullfile(folder,'Beta.csv'),Beta);
    end
end
end