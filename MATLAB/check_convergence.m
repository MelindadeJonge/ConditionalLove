function[] = check_convergence(model,id,startN,repN,niter,chains,q)
addpath(fullfile('MATLAB','HMSC_Class'))
addpath('output')
model_folder = fullfile('output',strcat(model,'Model',int2str(id)));
figures_folder = fullfile('figures',strcat(model,'Model',int2str(id)),'convergence');
mkdir(figures_folder)
load(fullfile(model_folder,'DataUsed.mat'));
covariate_names = regexprep(covariate_names,'"','');
species_names = regexprep(species_names,'/','');

species_exclude = [];
species_include = setdiff(1:length(species_names),species_exclude);

%Load posteriors for beta and omega
if strcmp(model,'Static')
    [beta_post, beta_names, omega_post, omega_names] = extract_beta(model_folder,repN,startN,niter,chains,covariate_names,species_names,0,species_include);
else
   	x = cell2mat(XrCell{1,1}(:,3));
   	Env = quantile(x,q)
   	[beta_post, beta_names, omega_post, omega_names] = extract_beta(model_folder,repN,startN,niter,chains,covariate_names,species_names,Env,species_include);
end


%Calculate and plot the sampling variance corrected Gelman-Rubin statistics (Rc) for the beta parameters and plot them
itertotal = (repN-startN)*niter;
Rc = GelmanRubin(beta_post,niter,itertotal);
figure1=figure;
xaxis = 1:niter:itertotal;
plot(xaxis,Rc);
hold on;
hline=refline([0,1.1]);
hline.Color='k';
hline.LineWidth = 0.5;
xlabel('Model iteration');
ylabel('Potential scale reduction factor');
filename = fullfile(figures_folder,'ShrinkageAllVars.jpg');
saveas(figure1,filename,'jpg')
hold off
%Make traceplots for beta parameters that have a GR statistic below 1.1
not_converged = Rc(:,(repN-startN))>1.1;
nplots = sum(not_converged);
for i=1:chains
    beta_plot{i} = beta_post{i}(:,not_converged);
end
beta_names_plot = beta_names(not_converged);
for i=1:nplots
	figure1=figure;
  xaxis = ((startN)*niter+1):repN*niter;
	plot(xaxis,beta_plot{1}(:,i),xaxis,beta_plot{2}(:,i),xaxis,beta_plot{3}(:,i));
	xlabel('MCMC iteration');
	ylabel(strcat(beta_names_plot(i)));
	legend('Chain1','Chain2','Chain3')
	filename = fullfile(figures_folder,char(beta_names_plot(i)));
	saveas(figure1,filename,'jpg')
end

if strcmp(model,'Static')
 	RcO = GelmanRubin(omega_post,niter,itertotal);
 	xaxis = 1:niter:itertotal;
	plot(xaxis,RcO);
	hold on;
	hline=refline([0,1.1]);
	hline.Color='k';
	hline.LineWidth = 0.5;
	xlabel('last iteration in chain');
	ylabel('Potential scale reduction factor');
	filename = fullfile(figures_folder,'ShrinkageOmega.jpg');
	saveas(figure1,filename,'jpg')
	hold off
	
	not_convergedO = RcO(:,(repN-startN))>1.1;
	omega_names_plot = omega_names(not_convergedO);
	for j=1:chains
		omega_plot{j} = omega_post{j}(:,not_convergedO);
	end
	for j=1:sum(not_convergedO)
		figure1=figure;
    xaxis = ((startN)*niter+1):repN*niter;
		plot(xaxis,omega_plot{1}(:,j),xaxis,omega_plot{2}(:,j),xaxis,omega_plot{3}(:,j));
		xlabel('MCMC iteration');
		ylabel(strcat(omega_names_plot(j)));
		legend('Chain1','Chain2','Chain3')
		filename = fullfile(figures_folder,char(strcat(omega_names_plot(j))));
		saveas(figure1,filename,'jpg')
	end
else
 	RcO = GelmanRubin(omega_post,niter,itertotal);
 	xaxis = 1:niter:itertotal;
	plot(xaxis,RcO);
	hold on;
	hline=refline([0,1.1]);
	hline.Color='k';
	hline.LineWidth = 0.5;
	xlabel('Model iteration');
	ylabel('Potential scale reduction factor');
	filename = fullfile(figures_folder,'ShrinkageOmega.jpg');
	saveas(figure1,filename,'jpg')
	hold off
	
	not_convergedO = RcO(:,(repN-startN))>1.1;
	omega_names_plot = omega_names(not_convergedO);
	for j=1:chains
		omega_plot{j} = omega_post{j}(:,not_convergedO);
	end
	for j=1:sum(not_convergedO)
		figure1=figure;
    		xaxis = ((startN)*niter+1):repN*niter;
		plot(xaxis,omega_plot{1}(:,j),xaxis,omega_plot{2}(:,j),xaxis,omega_plot{3}(:,j));
		xlabel('MCMC iteration');
		ylabel(strcat(omega_names_plot(j)));
		legend('Chain1','Chain2','Chain3')
			filename = fullfile(figures_folder,char(strcat(omega_names_plot(j),'CWD',int2str(q))));
		saveas(figure1,filename,'jpg')
	end

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
[~, nvar] = size(A{1});
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
%% ---------------------------------------------------------------------------------------------------
% Combine beta & omega posteriors from multiple chains into array
function[beta_post, beta_names, omega_post, omega_names] = extract_beta(folder,repN,start,niter,chains,covariate_names,species_names,Env,species_include)
species_names = species_names(species_include);
[~, ns] = size(species_names);
[~, nc] = size(covariate_names);
beta_post = {zeros((repN-start)*niter,nc*ns),zeros((repN-start)*niter,nc*ns),zeros((repN-start)*niter,nc*ns)};
omega_post = {zeros((repN-start)*niter,(ns*ns-ns)/2+ns),zeros((repN-start)*niter,(ns*ns-ns)/2+ns),zeros((repN-start)*niter,(ns*ns-ns)/2+ns)};  
for chain = 1:chains
	for i=1:(repN - start)
		p = load(fullfile(folder,strcat('Chain',num2str(chain)),'posteriors',strcat('postPar_run_',num2str(i+start),'.mat')));
		parVec = p.parVec;
		for j=1:niter
				bp = parVec{j}.beta(:,species_include);
				beta_post{chain}((j+(i-1)*niter),:) = bp(:)';
				if Env(1) == 0
					lambda = parVec{j}.lambda{1};
					lam = lambda(:,species_include);
					RR = tril(lam'*lam);
					RR = RR(:);
					omega_post{chain}((j+(i-1)*niter),:) = RR(RR~=0);;
				else
					lambda = parVec{j}.lambda{1};
					lam = lambda(:,species_include,1) + lambda(:,species_include,2)*Env;
					RR = tril(lam'*lam);
					RR = RR(:);
					omega_post{chain}((j+(i-1)*niter),:) = RR(RR~=0);;
				end
		end
		disp(strcat('Run ', int2str(i), 'of chain ',int2str(chain), 'loaded'))
	end
end

for Species=1:ns
	for Cov = 1:nc
		coeffnames(Species,Cov) = string(strcat(covariate_names{Cov},species_names{Species}));
	end
	for Species2=1:ns
		omega_names(Species,Species2) = string(strcat(species_names{Species},'_vs_',species_names{Species2}));
	end
	coeffnames = regexprep(coeffnames,'"','');
end
beta_names = reshape(coeffnames.',[],1);
ind = tril(ones(ns));
ind = ind(:);
omega_names = omega_names(:);
omega_names = omega_names(ind~=0);
end

