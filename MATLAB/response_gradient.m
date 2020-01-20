% Array job for predictions
function [] = response_gradient(model,id,chains)
%% Prep environment
addpath(fullfile('MATLAB','HMSC_Class'))
addpath('output')
model_folder = fullfile('output',strcat(model,'Model',int2str(id)));

for i=1:chains
   load(fullfile(model_folder,strcat('Chain',int2str(i)),strcat(model,int2str(id),'Chain',int2str(i),'.mat')));
   nvar = (m.nc-1)/2;
   CWDGrad = min(m.X(:,3)):0.1:max(m.X(:,3));
   [~, nNew] = size(CWDGrad);
   Env = [m.X;repmat(median(m.X),nNew,1)];
   Env((m.ny+1):(m.ny+nNew),3) = CWDGrad;
   Env((m.ny+1):(m.ny+nNew),(3+nvar)) = CWDGrad.^2;
   piNew = ones(1,nNew)*(m.ny+1);
   piCell = num2cell(([1:(m.ny),piNew])');
   piCell = cellfun(@num2str,piCell,'Uniformoutput',false);
   Xr = [[ones(m.ny,1), m.X(:,3)];  [ones(nNew,1), repmat(median(m.X(:,3)),nNew,1)]];
   XrCell = {[piCell(:,1),num2cell(Xr)]};
   PoO = zeros(nNew,m.ns);

   if strcmp(model,'Conditional')
       predList = m.predict(500,Env,[],[],piCell,[],XrCell,true);
	   PoO_all = mean(cat(3, predList{:}), 3);
       PoO = PoO + PoO_all((m.ny+1):end,:);
   elseif strcmp(model,'Static')
       predList = m.predict(500,Env,[],[],piCell,[],[],true);
       PoO_all = mean(cat(3, predList{:}), 3);
       PoO = PoO + PoO_all((m.ny+1):end,:);   
   end
end
PoO = PoO./3;
csvwrite(fullfile(model_folder,'CWDResponses.csv'),PoO);

end
