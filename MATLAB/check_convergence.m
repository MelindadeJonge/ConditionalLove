function[] = check_convergence(model,id,start_run)
addpath(fullfile('MATLAB','HMSC_Class'))
addpath('output')
model_folder = fullfile('output',strcat(model,'Model',int2str(id)));
figures_folder = fullfile('figures',strcat(model,'Model',int2str(id)),'convergence');
mkdir(figures_folder)
load(fullfile(model_folder,'DataUsed.mat'));
covariate_names = regexprep(covariate_names,'"','');
species_names = regexprep(species_names,'/','');

%Load posteriors
[beta_post, beta_names] = extract_beta(model_folder,repN,start_run,niter,chains,covariate_names,species_names);

%Calculate and plot the sampling variance corrected Gelman-Rubin statistics (Rc)
itertotal = (repN-start_run)*niter;
Rc = GelmanRubin(beta_post,niter,itertotal);
figure1=figure;
xaxis = 1:niter:itertotal;
xaxis = xaxis + start_run*niter;
plot(xaxis,Rc);
hold on;
axis([start_run*niter itertotal 0.9 max(max(Rc))]);
hline=refline([0,1.05]);
hline.Color='k';
hline.LineWidth = 0.5;
xlabel('last iteration in chain');
ylabel('Potential scale reduction factor');
filename = fullfile(figures_folder,'ShrinkageAllVars.jpg');
saveas(figure1,filename,'jpg')
hold off

%Make traceplots for variables that have a GR statistic below 1.05
not_converged = Rc(:,(repN-start_run))>2;
nplots = sum(not_converged);
for i=1:chains
    beta_plot{i} = beta_post{i}(:,not_converged);
end
beta_names_plot = beta_names(not_converged);
for i=1:nplots
	figure1=figure;
	plot(1:itertotal,beta_plot{1}(:,i),1:itertotal,beta_plot{2}(:,i),1:itertotal,beta_plot{3}(:,i));
	xlabel('MCMC iteration');
	ylabel(strcat(beta_names_plot(i)));
	legend('Chain1','Chain2','Chain3')
	filename = fullfile(figures_folder,char(beta_names_plot(i)));
	saveas(figure1,filename,'jpg')
end
end

%% -------------------------------------------------------------------------------------------------------------
% Gelman-Rubin statistics for MCMC chain convergence
%	This function calculates the Gelman-Rubin statistics for a set of MCMC chains that was given by the ChainsToArray function
%	Input A: cell containing M arrays (M = number of chains)
%	arrays should be structured as n rows and nc columns with nc the number of variables and n the number of iterations in a chain.
%
% Derived from:
%	Brooks, S. P., & Gelman, A. (1998). General methods for monitoring convergence of iterative simulations.
%	Journal of computational and graphical statistics, 7(4), 434-455.
%	Gelman, A., & Rubin, D. (1992). Inference from Iterative Simulation Using Multiple Sequences.
%	Statistical Science, 7(4), 457-472. Retrieved from http://www.jstor.org/stable/2246093
function [Rc] = GelmanRubin(A,sectionsize,itertotal)
[~,nchain] = size(A);
[nvar, ~] = size(A{1});
nsections = itertotal/sectionsize;
for section = 1:nsections
	lastit = section*sectionsize;
	for var = 1:nvar
		X = [];
		for chain = 1:nchain
			X = [X, A{chain}(1:lastit,var)];
        end
        [n, m] = size(X);
        overallmean = mean(mean(X));
        chainmean = mean(X);
        chainvariance = std(X).^2;
        B = n/(m-1)*(sum((chainmean-overallmean).^2));
        W = 1/m * sum(chainvariance);
        Var=(n-1)/n*W + B/n;
        Vpooled = Var + B/(m*n);
        %Correcting for sampling variance
        varssquare = std(chainvariance).^2 ;
        covsxsquare = cov(chainvariance,chainmean.^2);
        covsx = cov(chainvariance,chainmean);
        covsx(1,2);
        varv = ((n-1)/n)^2*1/m*varssquare + ((m+1)/(m*n))^2*2/(m-1)*B^2 + 2*((m+1)*(n-1))/(m*n^2) * n/m * (covsxsquare(1,2) - 2*overallmean*covsx(1,2));
        df = (2*Vpooled^2)/varv;
        Rc(var,section) = sqrt((df+3)/(df+1)*Vpooled/W);
	end
end
end

%% ---------------------------------------------------------------------------------------------------
% Combine beta posteriors from multiple chains into array
function[beta_post, beta_names] = extract_beta(folder,repN,start_run,niter,chains,covariate_names,species_names)
beta_post = cell(1,chains);
runs_wanted = repN - start_run;
[~, ns] = size(species_names);
[~, nc] = size(covariate_names);
a=1; nf = 1;
for chain = 1:chains
	for i=1:runs_wanted
		load(fullfile(folder,strcat('Chain',num2str(chain)),'posteriors',strcat('postPar_run_',num2str(i+start_run),'.mat')))
		for j=1:niter
			try
				beta_post{chain}((j+(i-1)*niter),:) = parVec{j}.beta(:)';
				a=a+1;
			catch
			end
		end
	end
end

for Species=1:ns
	for Cov = 1:nc
		coeffnames(Species,Cov) = string(strcat(covariate_names{Cov},species_names{Species}));
	end
	coeffnames = regexprep(coeffnames,'"','');
end
beta_names = reshape(coeffnames.',[],1);
end

%% ---------------------------------------------------------------------------------------------------
%Function to generate traceplots, NOTE! only for 3 chains! But can be easily updated to three.
function [] = plotTraceplots(X,names,burnin,Nitertotal,Niter,Sample,folder)

end

