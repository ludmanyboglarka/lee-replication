% This document serves as storing the empirical spatialStroop data in a m file (congruency sequence x SOA)
clear
warning('off');
% Empirical spatialStroop data
data=readtable('spatialStroopSOA.xlsx');
data=table2array(data);

NumP= 40;

Ppt = cell(NumP, 1);
for k=1:NumP
    Ppt{k} = data(data(:,1)==k,:);
end

cC = cell(NumP, 1);
cI = cell(NumP, 1);
iC = cell(NumP, 1);
iI = cell(NumP, 1);
for n=1:NumP
    for k=1:size(Ppt{n},1)-1
        if Ppt{n}(k,7)==1 && Ppt{n}(k+1, 7)==1 && Ppt{n}(k+1,2)-Ppt{n}(k,2)==1
            cC{n} = [cC{n} ; Ppt{n}(k+1,1) Ppt{n}(k+1,2) Ppt{n}(k+1,5) Ppt{n}(k+1,6), Ppt{n}(k+1,8)];
        elseif Ppt{n}(k,7)==0 && Ppt{n}(k+1, 7)==1 && Ppt{n}(k+1,2)-Ppt{n}(k,2)==1
            iC{n} = [iC{n} ; Ppt{n}(k+1,1) Ppt{n}(k+1,2) Ppt{n}(k+1,5) Ppt{n}(k+1,6), Ppt{n}(k+1,8)];
        elseif Ppt{n}(k,7)==1 && Ppt{n}(k+1, 7)==0 && Ppt{n}(k+1,2)-Ppt{n}(k,2)==1
            cI{n} = [cI{n} ; Ppt{n}(k+1,1) Ppt{n}(k+1,2) Ppt{n}(k+1,5) Ppt{n}(k+1,6), Ppt{n}(k+1,8)];
        elseif Ppt{n}(k,7)==0 && Ppt{n}(k+1, 7)==0 && Ppt{n}(k+1,2)-Ppt{n}(k,2)==1
            iI{n} = [iI{n} ; Ppt{n}(k+1,1) Ppt{n}(k+1,2) Ppt{n}(k+1,5) Ppt{n}(k+1,6), Ppt{n}(k+1,8)];
        end
    end
end

cC = cell2mat(cC);
de3=cC(:,3)<200;
cC(de3,:) = [];
de4=cC(:,3)>1500;
cC(de4,:) = [];

% delete the first trial of each block
de1 = cC(:, 2) ==[1 121 241 361 481 601 721 841];
de1=[find(de1(:, 1)==1);find(de1(:, 2)==1);find(de1(:, 3)==1);find(de1(:, 4)==1);find(de1(:, 5)==1);find(de1(:, 6)==1);find(de1(:, 7)==1);find(de1(:, 8)==1)];
cC(de1,:) = [];

cI = cell2mat(cI);
de7=cI(:,3)<200;
cI(de7,:) = [];
de8=cI(:,3)>1500;
cI(de8,:) = [];
de2 = cI(:, 2) ==[1 121 241 361 481 601 721 841];
de2=[find(de2(:, 1)==1);find(de2(:, 2)==1);find(de2(:, 3)==1);find(de2(:, 4)==1);find(de2(:, 5)==1);find(de2(:, 6)==1);find(de2(:, 7)==1);find(de2(:, 8)==1)];

cI(de2,:) = [];

iC = cell2mat(iC);
de11=iC(:,3)<200;
iC(de11,:) = [];
de12=iC(:,3)>1500;
iC(de12,:) = [];
de3 = iC(:, 2) ==[1 121 241 361 481 601 721 841];
de3=[find(de3(:, 1)==1);find(de3(:, 2)==1);find(de3(:, 3)==1);find(de3(:, 4)==1);find(de3(:, 5)==1);find(de3(:, 6)==1);find(de3(:, 7)==1);find(de3(:, 8)==1)];

iC(de3,:) = [];

iI = cell2mat(iI);
de15=iI(:,3)<200;
iI(de15,:) = [];
de16=iI(:,3)>1500;
iI(de16,:) = [];
de4 = iI(:, 2) ==[1 121 241 361 481 601 721 841];
de4=[find(de4(:, 1)==1);find(de4(:, 2)==1);find(de4(:, 3)==1);find(de4(:, 4)==1);find(de4(:, 5)==1);find(de4(:, 6)==1);find(de4(:, 7)==1);find(de4(:, 8)==1)];

iI(de4,:) = [];

cC1=cell(NumP, 1);
cC2=cell(NumP, 1);
for n=1:NumP
    k=cC(:, 1)==n;
    cC1{n} = cC(k, :);
    cC2{n} = [cC1{n}(:, 3), cC1{n}(:, 4), cC1{n}(:, 5)];
end
cC = cC2;

cI1=cell(NumP, 1);
cI2=cell(NumP, 1);
for n=1:NumP
    l=cI(:, 1)==n;
    cI1{n} = cI(l, :);
    cI2{n} = [cI1{n}(:, 3), cI1{n}(:, 4), cI1{n}(:, 5)];
end
cI=cI2;

iC1=cell(NumP, 1);
iC2=cell(NumP, 1);
for n=1:NumP
    k=iC(:, 1)==n;
    iC1{n} = iC(k, :);
    iC2{n} = [iC1{n}(:, 3), iC1{n}(:, 4), iC1{n}(:, 5)];
end
iC=iC2;

iI1=cell(NumP, 1);
iI2=cell(NumP, 1);
for n=1:NumP
    k=iI(:, 1)==n;
    iI1{n} = iI(k, :);
    iI2{n} = [iI1{n}(:, 3), iI1{n}(:, 4), iI1{n}(:, 5)];
end
iI=iI2;

% We have different SOAs in the spatialStroop task, -100s (target shown before
% distractors), 0, 50, 100s
% Five quantiles for each condition

RTcc = cell(n, 4);
for n=1:NumP
    RTcc{n, 1} = cC{n}(cC{n}(:, 3)==100,1:2);
    RTcc{n, 2} = cC{n}(cC{n}(:, 3)==0,1:2);
    RTcc{n, 3} = cC{n}(cC{n}(:, 3)==-50,1:2);
    RTcc{n, 4} = cC{n}(cC{n}(:, 3)==-100,1:2);
end

RTci = cell(n, 4);
for n=1:NumP
    RTci{n, 1} = cI{n}(cI{n}(:, 3)==100,1:2);
    RTci{n, 2} = cI{n}(cI{n}(:, 3)==0,1:2);
    RTci{n, 3} = cI{n}(cI{n}(:, 3)==-50,1:2); 
    RTci{n, 4} = cI{n}(cI{n}(:, 3)==-100,1:2);
end

RTic = cell(n, 4);
for n=1:NumP
    RTic{n, 1} = iC{n}(iC{n}(:, 3)==100,1:2);
    RTic{n, 2} = iC{n}(iC{n}(:, 3)==0,1:2);
    RTic{n, 3} = iC{n}(iC{n}(:, 3)==-50,1:2);
    RTic{n, 4} = iC{n}(iC{n}(:, 3)==-100,1:2);
end

RTii = cell(n, 4);
for n=1:NumP
    RTii{n, 1} = iI{n}(iI{n}(:, 3)==100,1:2);
    RTii{n, 2} = iI{n}(iI{n}(:, 3)==0,1:2);
    RTii{n, 3} = iI{n}(iI{n}(:, 3)==-50,1:2);
    RTii{n, 4} = iI{n}(iI{n}(:, 3)==-100,1:2);
end

cC_acc=cell(NumP, 4);
cI_acc=cell(NumP, 4);
for n=1:NumP
    for m=1:4
        cC_acc{n, m}=sum(RTcc{n, m}(:, 2))/numel(RTcc{n, m}(:, 2));
        cI_acc{n, m}=sum(RTci{n, m}(:, 2))/numel(RTci{n, m}(:, 2));
    end
end

iC_acc=cell(NumP, 4);
iI_acc=cell(NumP, 4);
for n=1:NumP
    for m=1:4
        iC_acc{n, m}=sum(RTic{n, m}(:, 2))/numel(RTic{n, m}(:, 2));
        iI_acc{n, m}=sum(RTii{n, m}(:, 2))/numel(RTii{n, m}(:, 2));
    end
end

cC_corr=cell(NumP, 4);
cC_err=cell(NumP, 4);
for n=1:NumP
    for m=1:4
        cC_corr{n, m} = RTcc{n, m}(RTcc{n, m}(:, 2)==1, :);
        cC_err{n, m} = RTcc{n, m}(RTcc{n, m}(:, 2)==0, :);
    end
end

cI_corr=cell(NumP, 4);
cI_err=cell(NumP, 4);
for n=1:NumP
    for m=1:4
        cI_corr{n, m} = RTci{n, m}(RTci{n, m}(:, 2)==1, :);
        cI_err{n, m} = RTci{n, m}(RTci{n, m}(:, 2)==0, :);
    end
end

iC_corr=cell(NumP, 4);
iC_err=cell(NumP, 4);
for n=1:NumP
    for m=1:4
        iC_corr{n, m} = RTic{n, m}(RTic{n, m}(:, 2)==1, :);
        iC_err{n, m} = RTic{n, m}(RTic{n, m}(:, 2)==0, :);
    end
end

iI_corr=cell(NumP, 4);
iI_err=cell(NumP, 4);
for n=1:NumP
    for m=1:4
        iI_corr{n, m} = RTii{n, m}(RTii{n, m}(:, 2)==1, :);
        iI_err{n, m} = RTii{n, m}(RTii{n, m}(:, 2)==0, :);
    end
end

bin_cc=cell(NumP, 4);
bin_cc_err=cell(NumP, 4);
for n=1:NumP
    for m=1:4
        bin_cc{n, m}=[prctile(cC_corr{n, m}(:, 1), 10), prctile(cC_corr{n, m}(:, 1),30), prctile(cC_corr{n, m}(:, 1),50), prctile(cC_corr{n, m}(:, 1),70), ...
            prctile(cC_corr{n, m}(:, 1),90), prctile(cC_corr{n, m}(:, 1),100)];
        bin_cc_err{n, m}=cC_err{n, m}(:, 1);
    end
end

bin_ci=cell(NumP, 4);
bin_ci_err=cell(NumP, 4);
for n=1:NumP
    for m=1:4
        bin_ci{n, m}=[prctile(cI_corr{n, m}(:, 1), 10), prctile(cI_corr{n, m}(:, 1),30), prctile(cI_corr{n, m}(:, 1),50), prctile(cI_corr{n, m}(:, 1),70), ...
            prctile(cI_corr{n, m}(:, 1),90), prctile(cI_corr{n, m}(:, 1),100)];
        bin_ci_err{n, m}=cI_err{n, m}(:, 1);
    end
end

for m=1:4
    mean_bin_cc(m, :)=[mean(cellfun(@(c) c(1,1), bin_cc(:, m))), mean(cellfun(@(c) c(1,2), bin_cc(:, m))), mean(cellfun(@(c) c(1,3), bin_cc(:, m))), ...
        mean(cellfun(@(c) c(1,4), bin_cc(:, m))), mean(cellfun(@(c) c(1,5), bin_cc(:, m))), mean(cellfun(@(c) c(1,6), bin_cc(:, m)))]; 
    mean_bin_ci(m, :)=[mean(cellfun(@(c) c(1,1), bin_ci(:, m))), mean(cellfun(@(c) c(1,2), bin_ci(:, m))), mean(cellfun(@(c) c(1,3), bin_ci(:, m))), ...
        mean(cellfun(@(c) c(1,4), bin_ci(:, m))), mean(cellfun(@(c) c(1,5), bin_ci(:, m))), mean(cellfun(@(c) c(1,6), bin_ci(:, m)))]; 
end

idx_cc=cell(NumP, m);
idx_ci=cell(NumP, m);
for n=1:NumP
    for m=1:4
        idx_cc{n, m} = [cC_corr{n, m}(:, 1)<=mean_bin_cc(m, 1), cC_corr{n, m}(:, 1)>mean_bin_cc(m, 1)&cC_corr{n, m}(:, 1)<=mean_bin_cc(m, 2), ...
            cC_corr{n, m}(:, 1)>mean_bin_cc(m, 2)&cC_corr{n, m}(:, 1)<=mean_bin_cc(m, 3), cC_corr{n, m}(:, 1)>mean_bin_cc(m, 3)&cC_corr{n, m}(:, 1)<=mean_bin_cc(m, 4), ...
            cC_corr{n, m}(:, 1)>mean_bin_cc(m, 4)&cC_corr{n, m}(:, 1)<=mean_bin_cc(m, 5), cC_corr{n, m}(:, 1)>mean_bin_cc(m, 5)&cC_corr{n, m}(:, 1)<=mean_bin_cc(m, 6)];
        idx_ci{n, m} = [cI_corr{n, m}(:, 1)<=mean_bin_ci(m, 1), cI_corr{n, m}(:, 1)>mean_bin_ci(m, 1)&cI_corr{n, m}(:, 1)<=mean_bin_ci(m, 2), cI_corr{n, m}(:, 1)>mean_bin_ci(m, 2)&...
            cI_corr{n, m}(:, 1)<=mean_bin_ci(m, 3), cI_corr{n, m}(:, 1)>mean_bin_ci(m, 3)&cI_corr{n, m}(:, 1)<=mean_bin_ci(m, 4), cI_corr{n, m}(:, 1)>mean_bin_ci(m, 4)&...
            cI_corr{n, m}(:, 1)<=mean_bin_ci(m, 5), cI_corr{n, m}(:, 1)>mean_bin_ci(m, 5)&cI_corr{n, m}(:, 1)<=mean_bin_ci(m, 6)];
    end
end

for n=1:NumP
    for m=1:4
        RT_cc{n,m}={cC_corr{n, m}(idx_cc{n, m}(:, 1),:), cC_corr{n, m}(idx_cc{n, m}(:, 2),:), cC_corr{n, m}(idx_cc{n, m}(:, 3),:), cC_corr{n, m}(idx_cc{n, m}(:, 4),:), ...
            cC_corr{n, m}(idx_cc{n, m}(:, 5),:), cC_corr{n, m}(idx_cc{n, m}(:, 6),:)};
        RT_ci{n,m}={cI_corr{n, m}(idx_ci{n, m}(:, 1),:), cI_corr{n, m}(idx_ci{n, m}(:, 2),:), cI_corr{n, m}(idx_ci{n, m}(:, 3),:), cI_corr{n, m}(idx_ci{n, m}(:, 4),:), ...
            cI_corr{n, m}(idx_ci{n, m}(:, 5),:), cI_corr{n, m}(idx_ci{n, m}(:, 6),:)};
    end
end

bin_ic=cell(NumP, 4);
bin_ic_err=cell(NumP, 4);
for n=1:NumP
    for m=1:4
        bin_ic{n, m}=[prctile(iC_corr{n, m}(:, 1), 10), prctile(iC_corr{n, m}(:, 1),30), prctile(iC_corr{n, m}(:, 1),50), prctile(iC_corr{n, m}(:, 1),70), ...
            prctile(iC_corr{n, m}(:, 1),90), prctile(iC_corr{n, m}(:, 1),100)];
        bin_ic_err{n, m}=iC_err{n, m}(:, 1);
    end
end

bin_ii=cell(NumP, 4);
bin_ii_err=cell(NumP, 4);
for n=1:NumP
    for m=1:4
        bin_ii{n, m}=[prctile(iI_corr{n, m}(:, 1), 10), prctile(iI_corr{n, m}(:, 1),30), prctile(iI_corr{n, m}(:, 1),50), prctile(iI_corr{n, m}(:, 1),70), ...
            prctile(iI_corr{n, m}(:, 1),90), prctile(iI_corr{n, m}(:, 1),100)];
        bin_ii_err{n, m}=iI_err{n, m}(:, 1);
    end
end

for m=1:4
    mean_bin_ic(m, :)=[mean(cellfun(@(c) c(1,1), bin_ic(:, m))), mean(cellfun(@(c) c(1,2), bin_ic(:, m))), mean(cellfun(@(c) c(1,3), bin_ic(:, m))), ...
        mean(cellfun(@(c) c(1,4), bin_ic(:, m))), mean(cellfun(@(c) c(1,5), bin_ic(:, m))), mean(cellfun(@(c) c(1,6), bin_ic(:, m)))]; 
    mean_bin_ii(m, :)=[mean(cellfun(@(c) c(1,1), bin_ii(:, m))), mean(cellfun(@(c) c(1,2), bin_ii(:, m))), mean(cellfun(@(c) c(1,3), bin_ii(:, m))), ...
        mean(cellfun(@(c) c(1,4), bin_ii(:, m))), mean(cellfun(@(c) c(1,5), bin_ii(:, m))), mean(cellfun(@(c) c(1,6), bin_ii(:, m)))]; 
end

idx_ic=cell(NumP, m);
idx_ii=cell(NumP, m);
for n=1:NumP
    for m=1:4
        idx_ic{n, m} = [iC_corr{n, m}(:, 1)<=mean_bin_ic(m, 1), iC_corr{n, m}(:, 1)>mean_bin_ic(m, 1)&iC_corr{n, m}(:, 1)<=mean_bin_ic(m, 2), ...
            iC_corr{n, m}(:, 1)>mean_bin_ic(m, 2)&iC_corr{n, m}(:, 1)<=mean_bin_ic(m, 3), iC_corr{n, m}(:, 1)>mean_bin_ic(m, 3)&iC_corr{n, m}(:, 1)<=mean_bin_ic(m, 4), ...
            iC_corr{n, m}(:, 1)>mean_bin_ic(m, 4)&iC_corr{n, m}(:, 1)<=mean_bin_ic(m, 5), iC_corr{n, m}(:, 1)>mean_bin_ic(m, 5)&iC_corr{n, m}(:, 1)<=mean_bin_ic(m, 6)];
        idx_ii{n, m} = [iI_corr{n, m}(:, 1)<=mean_bin_ii(m, 1), iI_corr{n, m}(:, 1)>mean_bin_ii(m, 1)&iI_corr{n, m}(:, 1)<=mean_bin_ii(m, 2), iI_corr{n, m}(:, 1)>mean_bin_ii(m, 2)&...
            iI_corr{n, m}(:, 1)<=mean_bin_ii(m, 3), iI_corr{n, m}(:, 1)>mean_bin_ii(m, 3)&iI_corr{n, m}(:, 1)<=mean_bin_ii(m, 4), iI_corr{n, m}(:, 1)>mean_bin_ii(m, 4)&...
            iI_corr{n, m}(:, 1)<=mean_bin_ii(m, 5), iI_corr{n, m}(:, 1)>mean_bin_ii(m, 5)&iI_corr{n, m}(:, 1)<=mean_bin_ii(m, 6)];
    end
end

RT_ic=cell(NumP, 4);
RT_ii=cell(NumP, 4);
for n=1:NumP
    for m=1:4
        RT_ic{n,m}={iC_corr{n, m}(idx_ic{n, m}(:, 1),:), iC_corr{n, m}(idx_ic{n, m}(:, 2),:), iC_corr{n, m}(idx_ic{n, m}(:, 3),:), iC_corr{n, m}(idx_ic{n, m}(:, 4),:), ...
            iC_corr{n, m}(idx_ic{n, m}(:, 5),:), iC_corr{n, m}(idx_ic{n, m}(:, 6),:)};
        RT_ii{n,m}={iI_corr{n, m}(idx_ii{n, m}(:, 1),:), iI_corr{n, m}(idx_ii{n, m}(:, 2),:), iI_corr{n, m}(idx_ii{n, m}(:, 3),:), iI_corr{n, m}(idx_ii{n, m}(:, 4),:), ...
            iI_corr{n, m}(idx_ii{n, m}(:, 5),:), iI_corr{n, m}(idx_ii{n, m}(:, 6),:)};
    end
end

Num_cc=cell(NumP, 4);
for n=1:NumP
    for m=1:4
        Num_cc{n, m}=cellfun(@(x) numel(x), RT_cc{n,m})/2;
    end
end

frq_cc=[];
for m=1:4
    frq_cc = [frq_cc ; mean(cell2mat(Num_cc(:, m))).*mean(cell2mat(cC_acc(:, m)))];
end

frq_cc_err=[];
for n=1:NumP
    for m=1:4
        frq_cc_err(n,m) = mean(numel(cell2mat(bin_cc_err(n, m)))*(1-mean(cell2mat(cC_acc(:, m)))));
    end
end

frq_cc_err=mean(frq_cc_err);

Num_ci=cell(NumP, 4);
for n=1:NumP
    for m=1:4
        Num_ci{n, m}=cellfun(@(x) numel(x), RT_ci{n,m})/2;
    end
end

frq_ci=[];
for m=1:4
    frq_ci = [frq_ci ; mean(cell2mat(Num_ci(:, m))).*mean(cell2mat(cI_acc(:, m)))];
end

frq_ci_err=[];
for n=1:NumP
    for m=1:4
        frq_ci_err(n,m) = mean(numel(cell2mat(bin_ci_err(n, m)))*(1-mean(cell2mat(cI_acc(:, m)))));
    end
end

frq_ci_err=mean(frq_ci_err);

Num_ic=cell(NumP, 4);
for n=1:NumP
    for m=1:4
        Num_ic{n, m}=cellfun(@(x) numel(x), RT_ic{n,m})/2;
    end
end

frq_ic=[];
for m=1:4
    frq_ic = [frq_ic ; mean(cell2mat(Num_ic(:, m))).*mean(cell2mat(iC_acc(:, m)))];
end

frq_ic_err=[];
for n=1:NumP
    for m=1:4
        frq_ic_err(n,m) = mean(numel(cell2mat(bin_ic_err(n, m)))*(1-mean(cell2mat(iC_acc(:, m)))));
    end
end

frq_ic_err=mean(frq_ic_err);

Num_ii=cell(NumP, 4);
for n=1:NumP
    for m=1:4
        Num_ii{n, m}=cellfun(@(x) numel(x), RT_ii{n,m})/2;
    end
end

frq_ii=[];
for m=1:4
    frq_ii = [frq_ii ; mean(cell2mat(Num_ii(:, m))).*mean(cell2mat(iI_acc(:, m)))];
end

frq_ii_err=[];
for n=1:NumP
    for m=1:4
        frq_ii_err(n,m) = mean(numel(cell2mat(bin_ii_err(n, m)))*(1-mean(cell2mat(iI_acc(:, m)))));
    end
end

frq_ii_err=mean(frq_ii_err);

save('spatialStroop_SOA_CSE_ReadData.mat')