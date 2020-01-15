function [correlations,combinations] = computeCorrelationsFullP(m,level,Xr)


AOmega = m.getPostOmega(level, Xr);
npost=size(AOmega{1},1)
if npost == 0
	AOmega=m.getPostOmega(level);
	npost = size(AOmega{1},1);
end
mR=zeros(m.ns,m.ns);
Rsig=zeros(m.ns,m.ns);

combinations = nchoosek(1:length(m.spNames),2);

for i=1:npost
	Omega=reshape(AOmega{1}(i,:),m.ns,m.ns);
	%[~,R]=cov2corr(Omega);
	R = zeros(m.ns,m.ns);
	for k=1:size(Omega,1)
		for l =1:size(Omega,1)
			R(k,l) = Omega(k,l)/sqrt(Omega(k,k)*Omega(l,l));
		end
	end
	for Combi = 1:length(combinations)
		correlations(Combi,i) = R(combinations(Combi,1),combinations(Combi,2));
	end
end



