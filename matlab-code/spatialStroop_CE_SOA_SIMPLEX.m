function spatialStroop_CE_SOA_SIMPLEX(m)
tic;
options = optimset('MaxIter', 100, 'MaxFunEvals', 3500);
load('spatialStroop_CE_SOA_new_st')

f = @(w)wrapper_spatialStroop_CE_SOA(m, w);
[fitted_Param,fval] = fminsearch(f, new_st(m, :), options);

objFctHandle = @spatialStroop_CE_SOA_SIMPLEX;
objFctSettings = m;
if ~isempty(objFctHandle)
  objFctName = strrep(func2str(objFctHandle), '/', '_');
else
  objFctName = 'error';
end

resultFileName = sprintf('%s_%02d.mat', objFctName, objFctSettings);
% save data
sem = setfilesemaphore(resultFileName, 'set always');
save(resultFileName, 'fitted_Param', 'fval');
removefilesemaphore(sem);

toc;
end