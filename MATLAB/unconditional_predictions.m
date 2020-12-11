%  job for predictions
function [] = unconditional_predictions(n,id,nchain)
addpath('HMSC_Class');
model_folder = fullfile('output',strcat('ConditionalModel',int2str(id)));
mkdir(fullfile(model_folder,'Predictions'));

for chain =1:nchain
    Data =load(fullfile(model_folder,'DataUsed.mat'));
    m = load(fullfile(model_folder,strcat('Chain',num2str(chain),'/Conditional',num2str(id),'Chain',num2str(chain),'.mat')));
    m=m.m;

    %Start making the predictions
    piNew = num2cell(1:(m.ny+m.ny))';
    piNew = cellfun(@num2str,piNew,'Uniformoutput',false);
    XNew = [m.X;m.X];
    Xr = [ones(m.ny*2,1),XNew(:,3)];
    XrNew = {[piNew(:,1),num2cell(Xr)]};
    predListU = m.predict(n,XNew,false,false,piNew,[],XrNew,false);
    predList{chain} = cat(3,predListU{:});
end

%Restructure predictions and save for each species seperately
for sp=1:m.ns
    pchain1 = predList{1}(:,sp,:);
    pchain2 = predList{2}(:,sp,:);
    pchain3 = predList{3}(:,sp,:);
    p = [squeeze(pchain1),squeeze(pchain2),squeeze(pchain3)];
    filename = strcat(model_folder,'/Predictions','/Species',int2str(sp),'.csv');
    csvwrite(filename,p);
end

end