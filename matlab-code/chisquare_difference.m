%% The goodness-of-fit test
% The Simon task (Congruency effect x SOA)
clear
NumP=6;
% Load the result of each dataset (in .mat) 
matData={};
data={};
for k=1:NumP
	% OPTION 1: Create a mat filename, and load it into a structure called matData.
	matFileName = sprintf('Simon_CE_SOA_SIMPLEX_%02d.mat', k);
	if isfile(matFileName)
		matData{1, k} = load(matFileName);
        data{1, k} = matData{1, k}.fitted_Param;
        data{2, k} = matData{1, k}.fval;
	else
		fprintf('File %s does not exist.\n', matFileName);
    end
end

df(1, 1)=abs(length(data{1, 1}(~isnan(data{1, 1})))-length(data{1, 2}(~isnan(data{1, 2}))));
chidiff(1, 1)=abs(min(data{2, 1})-min(data{2, 2}));
chi2(1, 1)=1-chi2cdf(chidiff(1, 1), df(1, 1));

df(1, 2)=abs(length(data{1, 2}(~isnan(data{1, 2})))-length(data{1, 3}(~isnan(data{1, 3}))));
chidiff(1, 2)=abs(min(data{2, 2})-min(data{2, 3}));
chi2(1, 2)=1-chi2cdf(chidiff(1, 2), df(1, 2));

df(1, 3)=abs(length(data{1, 2}(~isnan(data{1, 2})))-length(data{1, 4}(~isnan(data{1, 4}))));
chidiff(1, 3)=abs(min(data{2, 2})-min(data{2, 4}));
chi2(1, 3)=1-chi2cdf(chidiff(1, 3), df(1, 3));

df(1, 4)=abs(length(data{1, 2}(~isnan(data{1, 2})))-length(data{1, 5}(~isnan(data{1, 5}))));
chidiff(1, 4)=abs(min(data{2, 2})-min(data{2, 5}));
chi2(1, 4)=1-chi2cdf(chidiff(1, 4), df(1, 4));

df(1, 5)=abs(length(data{1, 2}(~isnan(data{1, 2})))-length(data{1, 6}(~isnan(data{1, 6}))));
chidiff(1, 5)=abs(min(data{2, 2})-min(data{2, 6}));
chi2(1, 5)=1-chi2cdf(chidiff(1, 5), df(1, 5));

%% The flanker task (Congruency effect x SOA)
clear
NumP=4;
% Load the result of each dataset (in .mat) 
matData={};
data={};
for k=1:NumP
	% OPTION 1: Create a mat filename, and load it into a structure called matData.
	matFileName = sprintf('flanker_CE_SOA_SIMPLEX_%02d.mat', k);
	if isfile(matFileName)
		matData{1, k} = load(matFileName);
        data{1, k} = matData{1, k}.fitted_Param;
        data{2, k} = matData{1, k}.fval;
	else
		fprintf('File %s does not exist.\n', matFileName);
    end
end

df(1, 1)=abs(length(data{1, 1}(~isnan(data{1, 1})))-length(data{1, 2}(~isnan(data{1, 2}))));
chidiff(1, 1)=abs(min(data{2, 1})-min(data{2, 2}));
chi2(1, 1)=1-chi2cdf(chidiff(1, 1), df(1, 1));

df(1, 2)=abs(length(data{1, 1}(~isnan(data{1, 1})))-length(data{1, 2}(~isnan(data{1, 2}))));
chidiff(1, 2)=abs(min(data{2, 1})-min(data{2, 2}));
chi2(1, 2)=1-chi2cdf(chidiff(1, 2), df(1, 2));

df(1, 3)=abs(length(data{1, 1}(~isnan(data{1, 1})))-length(data{1, 4}(~isnan(data{1, 4}))));
chidiff(1, 3)=abs(min(data{2, 1})-min(data{2, 4}));
chi2(1, 3)=1-chi2cdf(chidiff(1, 3), df(1, 3));

%% The Simon task (CSE x SOA)
clear
NumP=7;
% Load the result of each dataset (in .mat) 
matData={};
data={};
for k=1:NumP
	% OPTION 1: Create a mat filename, and load it into a structure called matData.
	matFileName = sprintf('Simon_CSE_SOA_SIMPLEX_%02d.mat', k);
	if isfile(matFileName)
		matData{1, k} = load(matFileName);
        data{1, k} = matData{1, k}.fitted_Param;
        data{2, k} = matData{1, k}.fval;
	else
		fprintf('File %s does not exist.\n', matFileName);
    end
end

df(1, 1)=abs(length(data{1, 1}(~isnan(data{1, 1})))-length(data{1, 4}(~isnan(data{1, 4}))));
chidiff(1, 1)=abs(min(data{2, 1})-min(data{2, 4}));
chi2(1, 1)=1-chi2cdf(chidiff(1, 1), df(1, 1));

df(1, 2)=abs(length(data{1, 2}(~isnan(data{1, 2})))-length(data{1, 4}(~isnan(data{1, 4}))));
chidiff(1, 2)=abs(min(data{2, 2})-min(data{2, 4}));
chi2(1, 2)=1-chi2cdf(chidiff(1, 1), df(1, 2));

df(1, 3)=abs(length(data{1, 3}(~isnan(data{1, 3})))-length(data{1, 4}(~isnan(data{1, 4}))));
chidiff(1, 3)=abs(min(data{2, 3})-min(data{2, 4}));
chi2(1, 3)=1-chi2cdf(chidiff(1, 3), df(1, 3));

df(1, 4)=abs(length(data{1, 5}(~isnan(data{1, 5})))-length(data{1, 4}(~isnan(data{1, 4}))));
chidiff(1, 4)=abs(min(data{2, 5})-min(data{2, 4}));
chi2(1, 4)=1-chi2cdf(chidiff(1, 4), df(1, 4));

df(1, 5)=abs(length(data{1, 6}(~isnan(data{1, 6})))-length(data{1, 4}(~isnan(data{1, 4}))));
chidiff(1, 5)=abs(min(data{2, 6})-min(data{2, 4}));
chi2(1, 5)=1-chi2cdf(chidiff(1, 5), df(1, 5));

%% The flanker task (CSE x SOA)
clear
NumP=2;
% Load the result of each dataset (in .mat) 
matData={};
data={};
for k=1:NumP
	% OPTION 1: Create a mat filename, and load it into a structure called matData.
	matFileName = sprintf('flanker_CSE_SOA_SIMPLEX_%02d.mat', k);
	if isfile(matFileName)
		matData{1, k} = load(matFileName);
        data{1, k} = matData{1, k}.fitted_Param;
        data{2, k} = matData{1, k}.fval;
	else
		fprintf('File %s does not exist.\n', matFileName);
    end
end

df(1, 1)=abs(length(data{1, 1}(~isnan(data{1, 1})))-length(data{1, 2}(~isnan(data{1, 2}))));
chidiff(1, 1)=abs(min(data{2, 1})-min(data{2, 2}));
chi2(1, 1)=1-chi2cdf(chidiff(1, 1), df(1, 1));

