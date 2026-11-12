function deviance=spatialStroop_CSE_SOA(A_01_pc, A_02_pc, A_03_pc, A_04_pc, A_01_pi, A_02_pi, A_03_pi, A_04_pi, ...
    kc_pc, kc_pi, ki_pc, ki_pi, da1_pc, da2_pc, da3_pc, da4_pc, da1_pi, da2_pi, da3_pi, da4_pi, ...
    dc, a_pc, a_pi, Ter1, Ter2, Ter3, Ter4, std, c)

load('spatialStroop_SOA_CSE_ReadData.mat');

Number_of_simulation = 50000;

upperbd1 = Ter1+0.5*std;
lowerbd1 = Ter1-0.5*std;

upperbd2 = Ter2+0.5*std;
lowerbd2 = Ter2-0.5*std;

upperbd3 = Ter3+0.5*std;
lowerbd3 = Ter3-0.5*std;

upperbd4 = Ter4+0.5*std;
lowerbd4 = Ter4-0.5*std;

R=[];
R(1, :)=randraw('uniform', [lowerbd1, upperbd1], 1, Number_of_simulation);
R(2, :)=randraw('uniform', [lowerbd2, upperbd2], 1, Number_of_simulation);
R(3, :)=randraw('uniform', [lowerbd3, upperbd3], 1, Number_of_simulation);
R(4, :)=randraw('uniform', [lowerbd4, upperbd4], 1, Number_of_simulation);

cc_pred1 = [];
cc_pred2 = [];
cc_pred3 = [];
cc_pred4 = [];
ci_pred1 = [];
ci_pred2 = [];
ci_pred3 = [];
ci_pred4 = [];
ic_pred1 = [];
ic_pred2 = [];
ic_pred3 = [];
ic_pred4 = [];
ii_pred1 = [];
ii_pred2 = [];
ii_pred3 = [];
ii_pred4 = [];

parfor x=1:Number_of_simulation
    cc_pred1(x, :) = RDMC_SOA(A_01_pc, kc_pc, da1_pc, dc, c);
    cc_pred2(x, :) = RDMC_SOA(A_02_pc, kc_pc, da2_pc, dc, c);
    cc_pred3(x, :) = RDMC_SOA(A_03_pc, kc_pc, da3_pc, dc, c);
    cc_pred4(x, :) = RDMC_SOA(A_04_pc, kc_pc, da4_pc, dc, c);
    ci_pred1(x, :) = RDMC_SOA(A_01_pc, ki_pc, -da1_pc, dc, c);
    ci_pred2(x, :) = RDMC_SOA(A_02_pc, ki_pc, -da2_pc, dc, c);
    ci_pred3(x, :) = RDMC_SOA(A_03_pc, ki_pc, -da3_pc, dc, c);
    ci_pred4(x, :) = RDMC_SOA(A_04_pc, ki_pc, -da4_pc, dc, c);
    ic_pred1(x, :) = RDMC_SOA(A_01_pi, kc_pi, da1_pi, dc, c);
    ic_pred2(x, :) = RDMC_SOA(A_02_pi, kc_pi, da2_pi, dc, c);
    ic_pred3(x, :) = RDMC_SOA(A_03_pi, kc_pi, da3_pi, dc, c);
    ic_pred4(x, :) = RDMC_SOA(A_04_pi, kc_pi, da4_pi, dc, c);
    ii_pred1(x, :) = RDMC_SOA(A_01_pi, ki_pi, -da1_pi, dc, c);
    ii_pred2(x, :) = RDMC_SOA(A_02_pi, ki_pi, -da2_pi, dc, c);
    ii_pred3(x, :) = RDMC_SOA(A_03_pi, ki_pi, -da3_pi, dc, c);
    ii_pred4(x, :) = RDMC_SOA(A_04_pi, ki_pi, -da4_pi, dc, c);
end

cc_pred={cc_pred1, cc_pred2, cc_pred3, cc_pred4};
ci_pred={ci_pred1, ci_pred2, ci_pred3, ci_pred4};
ic_pred={ic_pred1, ic_pred2, ic_pred3, ic_pred4};
ii_pred={ii_pred1, ii_pred2, ii_pred3, ii_pred4};

cC_predicted={};
for n=1:4
    for x=1:Number_of_simulation
        if isempty(find(cc_pred{1, n}(x,:)<=-a_pc,1))==0
            cC_predicted{1, n}(x, 1) = find(cc_pred{1, n}(x,:)<=-a_pc,1)+R(n, x);
            cC_predicted{1, n}(x, 2) = 0;
        elseif isempty(find(cc_pred{1, n}(x,:)<=-a_pc,1))==1 && isempty(find(cc_pred{1, n}(x,:)>=a_pc,1))==0
            cC_predicted{1, n}(x, 1) = find(cc_pred{1, n}(x,:)>=a_pc,1)+R(n, x);
            cC_predicted{1, n}(x, 2) = 1;
        else
            cC_predicted{1, n}(x, 1) = 1500;
            cC_predicted{1, n}(x, 2) = 1;
        end
    end
end

for n=1:4
    cC_predicted{1, n}(cC_predicted{1, n}(:, 1)>=1500, :)=[];
end

cI_predicted={};
for n=1:4
    for x=1:Number_of_simulation
        if isempty(find(ci_pred{1, n}(x,:)<=-a_pc,1))==0
            cI_predicted{1, n}(x, 1) = find(ci_pred{1, n}(x,:)<=-a_pc,1)+R(n, x);
            cI_predicted{1, n}(x, 2) = 0;
        elseif isempty(find(ci_pred{1, n}(x,:)<=-a_pc,1))==1 && isempty(find(ci_pred{1, n}(x,:)>=a_pc,1))==0
            cI_predicted{1, n}(x, 1) = find(ci_pred{1, n}(x,:)>=a_pc,1)+R(n, x);
            cI_predicted{1, n}(x, 2) = 1;
        else
            cI_predicted{1, n}(x, 1) = 1500;
            cI_predicted{1, n}(x, 2) = 1;
        end
    end
end

for n=1:4
    cI_predicted{1, n}(cI_predicted{1, n}(:, 1)>=1500, :)=[];
end

cC_acc_pred=[];
cI_acc_pred=[];
cC_corr_pred={};
cC_err_pred={};
cI_corr_pred={};
cI_err_pred={};
for n=1:4
    cC_acc_pred(1, n)=sum(cC_predicted{1, n}(:, 2))/numel(cC_predicted{1, n}(:, 2));
    cI_acc_pred(1, n)=sum(cI_predicted{1, n}(:, 2))/numel(cI_predicted{1, n}(:, 2));
end

for n=1:4
    cC_corr_pred{1, n} = cC_predicted{1, n}(cC_predicted{1, n}(:, 2)==1, :);
    cI_corr_pred{1, n} = cI_predicted{1, n}(cI_predicted{1, n}(:, 2)==1, :);
    cC_err_pred{1, n} = cC_predicted{1, n}(cC_predicted{1, n}(:, 2)==0, :);
    cI_err_pred{1, n} = cI_predicted{1, n}(cI_predicted{1, n}(:, 2)==0, :);
end

bincc=[];
binci=[];
for n=1:4
    bincc(n, :) = mean(cell2mat(bin_cc(:, n)));
    binci(n, :) = mean(cell2mat(bin_ci(:, n)));
end

% Plot the congruent and incongruent CDF
bin_cc_err_pred={};
for n=1:4
    bin_cc_err_pred{1, n}=cC_err_pred{1, n}(:, 1);
end

idx_cc_pred={};
RT_cc_pred={};
for n=1:4
    idx_cc_pred{1, n} = [cC_corr_pred{1, n}(:, 1)<=bincc(n, 1), cC_corr_pred{1, n}(:, 1)>bincc(n, 1)&cC_corr_pred{1, n}(:, 1)<=bincc(n, 2), ...
        cC_corr_pred{1, n}(:, 1)>bincc(n, 2)&cC_corr_pred{1, n}(:, 1)<=bincc(n, 3), ...
        cC_corr_pred{1, n}(:, 1)>bincc(n, 3)&cC_corr_pred{1, n}(:, 1)<=bincc(n, 4), ...
        cC_corr_pred{1, n}(:, 1)>bincc(n, 4)&cC_corr_pred{1, n}(:, 1)<=bincc(n, 5), ...
        cC_corr_pred{1, n}(:, 1)>bincc(n, 5)&cC_corr_pred{1, n}(:, 1)<=bincc(n, 6)];
    RT_cc_pred{1, n}={cC_corr_pred{1, n}(idx_cc_pred{1, n}(:, 1), :), cC_corr_pred{1, n}(idx_cc_pred{1, n}(:, 2), :), ...
        cC_corr_pred{1, n}(idx_cc_pred{1, n}(:, 3), :), cC_corr_pred{1, n}(idx_cc_pred{1, n}(:, 4), :), ...
        cC_corr_pred{1, n}(idx_cc_pred{1, n}(:, 5), :), cC_corr_pred{1, n}(idx_cc_pred{1, n}(:, 6), :)}; 
    
end 

bin_ci_err_pred={};
for n=1:4
    bin_ci_err_pred{1, n}=cI_err_pred{1, n}(:, 1);
end

idx_ci_pred={};
RT_ci_pred={};
for n=1:4
    idx_ci_pred{1, n} = [cI_corr_pred{1, n}(:, 1)<=binci(n, 1), cI_corr_pred{1, n}(:, 1)>binci(n, 1)&cI_corr_pred{1, n}(:, 1)<=binci(n, 2), ...
        cI_corr_pred{1, n}(:, 1)>binci(n, 2)&cI_corr_pred{1, n}(:, 1)<=binci(n, 3), ...
        cI_corr_pred{1, n}(:, 1)>binci(n, 3)&cI_corr_pred{1, n}(:, 1)<=binci(n, 4), ...
        cI_corr_pred{1, n}(:, 1)>binci(n, 4)&cI_corr_pred{1, n}(:, 1)<=binci(n, 5), ...
        cI_corr_pred{1, n}(:, 1)>binci(n, 5)&cI_corr_pred{1, n}(:, 1)<=binci(n, 6)];
    RT_ci_pred{1, n}={cI_corr_pred{1, n}(idx_ci_pred{1, n}(:, 1), :), cI_corr_pred{1, n}(idx_ci_pred{1, n}(:, 2), :), ...
        cI_corr_pred{1, n}(idx_ci_pred{1, n}(:, 3), :), cI_corr_pred{1, n}(idx_ci_pred{1, n}(:, 4), :), ...
        cI_corr_pred{1, n}(idx_ci_pred{1, n}(:, 5), :), cI_corr_pred{1, n}(idx_ci_pred{1, n}(:, 6), :)}; 
    
end 

iC_predicted={};
for n=1:4
    for x=1:Number_of_simulation
        if isempty(find(ic_pred{1, n}(x,:)<=-a_pi,1))==0
            iC_predicted{1, n}(x, 1) = find(ic_pred{1, n}(x,:)<=-a_pi,1)+R(n, x);
            iC_predicted{1, n}(x, 2) = 0;
        elseif isempty(find(ic_pred{1, n}(x,:)<=-a_pi,1))==1 && isempty(find(ic_pred{1, n}(x,:)>=a_pi,1))==0
            iC_predicted{1, n}(x, 1) = find(ic_pred{1, n}(x,:)>=a_pi,1)+R(n, x);
            iC_predicted{1, n}(x, 2) = 1;
        else
            iC_predicted{1, n}(x, 1) = 1500;
            iC_predicted{1, n}(x, 2) = 1;
        end
    end
end

for n=1:4
    iC_predicted{1, n}(iC_predicted{1, n}(:, 1)>=1500, :)=[];
end

iI_predicted={};
for n=1:4
    for x=1:Number_of_simulation
        if isempty(find(ii_pred{1, n}(x,:)<=-a_pi,1))==0
            iI_predicted{1, n}(x, 1) = find(ii_pred{1, n}(x,:)<=-a_pi,1)+R(n, x);
            iI_predicted{1, n}(x, 2) = 0;
        elseif isempty(find(ii_pred{1, n}(x,:)<=-a_pi,1))==1 && isempty(find(ii_pred{1, n}(x,:)>=a_pi,1))==0
            iI_predicted{1, n}(x, 1) = find(ii_pred{1, n}(x,:)>=a_pi,1)+R(n, x);
            iI_predicted{1, n}(x, 2) = 1;
        else
            iI_predicted{1, n}(x, 1) = 1500;
            iI_predicted{1, n}(x, 2) = 1;
        end
    end
end

for n=1:4
    iI_predicted{1, n}(iI_predicted{1, n}(:, 1)>=1500, :)=[];
end

iC_acc_pred=[];
iI_acc_pred=[];
iC_corr_pred={};
iC_err_pred={};
iI_corr_pred={};
iI_err_pred={};
for n=1:4
    iC_acc_pred(1, n)=sum(iC_predicted{1, n}(:, 2))/numel(iC_predicted{1, n}(:, 2));
    iI_acc_pred(1, n)=sum(iI_predicted{1, n}(:, 2))/numel(iI_predicted{1, n}(:, 2));
end

for n=1:4
    iC_corr_pred{1, n} = iC_predicted{1, n}(iC_predicted{1, n}(:, 2)==1, :);
    iI_corr_pred{1, n} = iI_predicted{1, n}(iI_predicted{1, n}(:, 2)==1, :);
    iC_err_pred{1, n} = iC_predicted{1, n}(iC_predicted{1, n}(:, 2)==0, :);
    iI_err_pred{1, n} = iI_predicted{1, n}(iI_predicted{1, n}(:, 2)==0, :);
end

binic=[];
binii=[];
for n=1:4
    binic(n, :) = mean(cell2mat(bin_ic(:, n)));
    binii(n, :) = mean(cell2mat(bin_ii(:, n)));
end

% Plot the congruent and incongruent CDF
bin_ic_err_pred={};
for n=1:4
    bin_ic_err_pred{1, n}=iC_err_pred{1, n}(:, 1);
end

idx_ic_pred={};
RT_ic_pred={};
for n=1:4
    idx_ic_pred{1, n} = [iC_corr_pred{1, n}(:, 1)<=binic(n, 1), iC_corr_pred{1, n}(:, 1)>binic(n, 1)&iC_corr_pred{1, n}(:, 1)<=binic(n, 2), ...
        iC_corr_pred{1, n}(:, 1)>binic(n, 2)&iC_corr_pred{1, n}(:, 1)<=binic(n, 3), ...
        iC_corr_pred{1, n}(:, 1)>binic(n, 3)&iC_corr_pred{1, n}(:, 1)<=binic(n, 4), ...
        iC_corr_pred{1, n}(:, 1)>binic(n, 4)&iC_corr_pred{1, n}(:, 1)<=binic(n, 5), ...
        iC_corr_pred{1, n}(:, 1)>binic(n, 5)&iC_corr_pred{1, n}(:, 1)<=binic(n, 6)];
    RT_ic_pred{1, n}={iC_corr_pred{1, n}(idx_ic_pred{1, n}(:, 1), :), iC_corr_pred{1, n}(idx_ic_pred{1, n}(:, 2), :), ...
        iC_corr_pred{1, n}(idx_ic_pred{1, n}(:, 3), :), iC_corr_pred{1, n}(idx_ic_pred{1, n}(:, 4), :), ...
        iC_corr_pred{1, n}(idx_ic_pred{1, n}(:, 5), :), iC_corr_pred{1, n}(idx_ic_pred{1, n}(:, 6), :)}; 
    
end 

bin_ii_err_pred={};
for n=1:4
    bin_ii_err_pred{1, n}=iI_err_pred{1, n}(:, 1);
end

idx_ii_pred={};
RT_ii_pred={};
for n=1:4
    idx_ii_pred{1, n} = [iI_corr_pred{1, n}(:, 1)<=binii(n, 1), iI_corr_pred{1, n}(:, 1)>binii(n, 1)&iI_corr_pred{1, n}(:, 1)<=binii(n, 2), ...
        iI_corr_pred{1, n}(:, 1)>binii(n, 2)&iI_corr_pred{1, n}(:, 1)<=binii(n, 3), ...
        iI_corr_pred{1, n}(:, 1)>binii(n, 3)&iI_corr_pred{1, n}(:, 1)<=binii(n, 4), ...
        iI_corr_pred{1, n}(:, 1)>binii(n, 4)&iI_corr_pred{1, n}(:, 1)<=binii(n, 5), ...
        iI_corr_pred{1, n}(:, 1)>binii(n, 5)&iI_corr_pred{1, n}(:, 1)<=binii(n, 6)];
    RT_ii_pred{1, n}={iI_corr_pred{1, n}(idx_ii_pred{1, n}(:, 1), :), iI_corr_pred{1, n}(idx_ii_pred{1, n}(:, 2), :), ...
        iI_corr_pred{1, n}(idx_ii_pred{1, n}(:, 3), :), iI_corr_pred{1, n}(idx_ii_pred{1, n}(:, 4), :), ...
        iI_corr_pred{1, n}(idx_ii_pred{1, n}(:, 5), :), iI_corr_pred{1, n}(idx_ii_pred{1, n}(:, 6), :)}; 
    
end 

cC_acc=cell2mat(cC_acc);
Num_cc={};
frq_cc={};
for n=1:4
    Num_cc{1, n}= cellfun(@(x) numel(x), vertcat(RT_cc{:, n}))/2;
    frq_cc{1, n} = Num_cc{1, n}.*cC_acc(:, n);
end

frq_cc_err=[];
for n=1:4
    for j=1:40
        frq_cc_err(j, n) = numel(bin_cc_err{j, n}(:, 1))*(1-cC_acc(j, n));
    end
end
frq_cc_err=mean(frq_cc_err);

frq_cc_tl=[];
frq_cc_all=[];
for n=1:4
    frq_cc{1, n} = mean(frq_cc{1, n});
    frq_cc_tl(1, n)=sum(frq_cc{1, n});
    frq_cc_all(1, n) = sum(frq_cc_tl(1, n))+frq_cc_err(1, n);
end

cI_acc=cell2mat(cI_acc);
Num_ci={};
frq_ci={};
for n=1:4
    Num_ci{1, n}= cellfun(@(x) numel(x), vertcat(RT_ci{:, n}))/2;
    frq_ci{1, n} = Num_ci{1, n}.*cI_acc(:, n);
end

frq_ci_err=[];
for n=1:4
    for j=1:40
        frq_ci_err(j, n) = numel(bin_ci_err{j, n}(:, 1))*(1-cI_acc(j, n));
    end
end
frq_ci_err=mean(frq_ci_err);

frq_ci_tl=[];
frq_ci_all=[];
for n=1:4
    frq_ci{1, n} = mean(frq_ci{1, n});
    frq_ci_tl(1, n)=sum(frq_ci{1, n});
    frq_ci_all(1, n) = sum(frq_ci_tl(1, n))+frq_ci_err(1, n);
end

iC_acc=cell2mat(iC_acc);
Num_ic={};
frq_ic={};
for n=1:4
    Num_ic{1, n}= cellfun(@(x) numel(x), vertcat(RT_ic{:, n}))/2;
    frq_ic{1, n} = Num_ic{1, n}.*iC_acc(:, n);
end

frq_ic_err=[];
for n=1:4
    for j=1:40
        frq_ic_err(j, n) = numel(bin_ic_err{j, n}(:, 1))*(1-iC_acc(j, n));
    end
end
frq_ic_err=mean(frq_ic_err);

frq_ic_tl=[];
frq_ic_all=[];
for n=1:4
    frq_ic{1, n} = mean(frq_ic{1, n});
    frq_ic_tl(1, n)=sum(frq_ic{1, n});
    frq_ic_all(1, n) = sum(frq_ic_tl(1, n))+frq_ic_err(1, n);
end

iI_acc=cell2mat(iI_acc);
Num_ii={};
frq_ii={};
for n=1:4
    Num_ii{1, n}= cellfun(@(x) numel(x), vertcat(RT_ii{:, n}))/2;
    frq_ii{1, n} = Num_ii{1, n}.*iI_acc(:, n);
end

frq_ii_err=[];
for n=1:4
    for j=1:40
        frq_ii_err(j, n) = numel(bin_ii_err{j, n}(:, 1))*(1-iI_acc(j, n));
    end
end
frq_ii_err=mean(frq_ii_err);

frq_ii_tl=[];
frq_ii_all=[];
for n=1:4
    frq_ii{1, n} = mean(frq_ii{1, n});
    frq_ii_tl(1, n)=sum(frq_ii{1, n});
    frq_ii_all(1, n) = sum(frq_ii_tl(1, n))+frq_ii_err(1, n);
end

Num_cc_pred={};
Q_cc_pred={};
Q_cc_pred_err=[];
for n=1:4
    Num_cc_pred{1, n}= cellfun(@(x) numel(x), vertcat(RT_cc_pred{:, n}))/2;
    Q_cc_pred{1, n} = Num_cc_pred{1, n}.*cC_acc_pred(:, n);
    Q_cc_pred_err(:, n) = numel(bin_cc_err_pred(:, n))*(1-cC_acc_pred(:, n));
end

Q_cc_pred_all=[];
for n=1:4
    Q_cc_pred_all(1, n) = sum(Q_cc_pred{1, n})+Q_cc_pred_err(1, n); 
end

prp_cc_pred=[];
prp_cc_err_pred=[];

for n=1:4
    prp_cc_pred(n, :) = Q_cc_pred{1, n}/Q_cc_pred_all(1, n);
    prp_cc_err_pred(1, n) = Q_cc_pred_err(1, n)/Q_cc_pred_all(1, n);  
end

frq_cc_pred={};
frq_cc_pred_err=[];
for n=1:4
    frq_cc_pred{1, n} = frq_cc_all(1, n).*prp_cc_pred(n, :);
    frq_cc_pred_err(1, n) = frq_cc_all(1, n).*prp_cc_err_pred(1, n);
end

chi_cc={};
for n=1:4
    chi_cc{1, n}=frq_cc{1, n}.*log(frq_cc{1, n}./frq_cc_pred{1, n});
end

Num_ci_pred={};
Q_ci_pred={};
Q_ci_pred_err=[];
for n=1:4
    Num_ci_pred{1, n}= cellfun(@(x) numel(x), vertcat(RT_ci_pred{:, n}))/2;
    Q_ci_pred{1, n} = Num_ci_pred{1, n}.*cI_acc_pred(:, n);
    Q_ci_pred_err(:, n) = numel(bin_ci_err_pred(:, n))*(1-cI_acc_pred(:, n));
end

Q_ci_pred_all=[];
for n=1:4
    Q_ci_pred_all(1, n) = sum(Q_ci_pred{1, n})+Q_ci_pred_err(1, n); 
end

prp_ci_pred=[];
prp_ci_err_pred=[];
for n=1:4
    prp_ci_pred(n, :) = Q_ci_pred{1, n}/Q_ci_pred_all(1, n);
    prp_ci_err_pred(1, n) = Q_ci_pred_err(1, n)/Q_ci_pred_all(1, n); 
end

frq_ci_pred={};
frq_ci_pred_err=[];
for n=1:4
    frq_ci_pred{1, n} = frq_ci_all(1, n).*prp_ci_pred(n, :);
    frq_ci_pred_err(1, n) = frq_ci_all(1, n).*prp_ci_err_pred(1, n);
end

chi_ci={};
for n=1:4
    chi_ci{1, n}=frq_ci{1, n}.*log(frq_ci{1, n}./frq_ci_pred{1, n});
end

chi_cc_err=[];
for n=1:4
    if frq_cc_err(1, n)==0 || frq_cc_pred_err(1, n)==0
        chi_cc_err(1, n)=0;
    else 
        chi_cc_err(1, n)=frq_cc_err(1, n).*log(frq_cc_err(1, n)./frq_cc_pred_err(1, n));
    end
end

chi_ci_err=[];
for n=1:4
    if frq_ci_err(1, n)==0 || frq_ci_pred_err(1, n)==0
        chi_ci_err(1, n)=0;
    else 
        chi_ci_err(1, n)=frq_ci_err(1, n).*log(frq_ci_err(1, n)./frq_ci_pred_err(1, n));
    end
end

Num_ic_pred={};
Q_ic_pred={};
Q_ic_pred_err=[];
for n=1:4
    Num_ic_pred{1, n}= cellfun(@(x) numel(x), vertcat(RT_ic_pred{:, n}))/2;
    Q_ic_pred{1, n} = Num_ic_pred{1, n}.*cC_acc_pred(:, n);
    Q_ic_pred_err(:, n) = numel(bin_ic_err_pred(:, n))*(1-cC_acc_pred(:, n));
end

Q_ic_pred_all=[];
for n=1:4
    Q_ic_pred_all(1, n) = sum(Q_ic_pred{1, n})+Q_ic_pred_err(1, n); 
end

prp_ic_pred=[];
prp_ic_err_pred=[];

for n=1:4
    prp_ic_pred(n, :) = Q_ic_pred{1, n}/Q_ic_pred_all(1, n);
    prp_ic_err_pred(1, n) = Q_ic_pred_err(1, n)/Q_ic_pred_all(1, n);  
end

frq_ic_pred={};
frq_ic_pred_err=[];
for n=1:4
    frq_ic_pred{1, n} = frq_ic_all(1, n).*prp_ic_pred(n, :);
    frq_ic_pred_err(1, n) = frq_ic_all(1, n).*prp_ic_err_pred(1, n);
end

chi_ic={};
for n=1:4
    chi_ic{1, n}=frq_ic{1, n}.*log(frq_ic{1, n}./frq_ic_pred{1, n});
end

Num_ii_pred={};
Q_ii_pred={};
Q_ii_pred_err=[];
for n=1:4
    Num_ii_pred{1, n}= cellfun(@(x) numel(x), vertcat(RT_ii_pred{:, n}))/2;
    Q_ii_pred{1, n} = Num_ii_pred{1, n}.*cI_acc_pred(:, n);
    Q_ii_pred_err(:, n) = numel(bin_ii_err_pred(:, n))*(1-cI_acc_pred(:, n));
end

Q_ii_pred_all=[];
for n=1:4
    Q_ii_pred_all(1, n) = sum(Q_ii_pred{1, n})+Q_ii_pred_err(1, n); 
end

prp_ii_pred=[];
prp_ii_err_pred=[];
for n=1:4
    prp_ii_pred(n, :) = Q_ii_pred{1, n}/Q_ii_pred_all(1, n);
    prp_ii_err_pred(1, n) = Q_ii_pred_err(1, n)/Q_ii_pred_all(1, n); 
end

frq_ii_pred={};
frq_ii_pred_err=[];
for n=1:4
    frq_ii_pred{1, n} = frq_ii_all(1, n).*prp_ii_pred(n, :);
    frq_ii_pred_err(1, n) = frq_ii_all(1, n).*prp_ii_err_pred(1, n);
end

chi_ii={};
for n=1:4
    chi_ii{1, n}=frq_ii{1, n}.*log(frq_ii{1, n}./frq_ii_pred{1, n});
end

chi_ic_err=[];
for n=1:4
    if frq_ic_err(1, n)==0 || frq_ic_pred_err(1, n)==0
        chi_ic_err(1, n)=0;
    else 
        chi_ic_err(1, n)=frq_ic_err(1, n).*log(frq_ic_err(1, n)./frq_ic_pred_err(1, n));
    end
end

chi_ii_err=[];
for n=1:4
    if frq_ii_err(1, n)==0 || frq_ii_pred_err(1, n)==0
        chi_ii_err(1, n)=0;
    else 
        chi_ii_err(1, n)=frq_ii_err(1, n).*log(frq_ii_err(1, n)./frq_ii_pred_err(1, n));
    end
end

chi_p=[];
for n=1:4
    chi_p(n, :)=chi_cc{1, n}+chi_ci{1, n}+chi_ic{1, n}+chi_ii{1, n};
end

deviance=2*(sum(sum(chi_p))+sum(chi_cc_err)+sum(chi_ci_err)+sum(chi_ic_err)+sum(chi_ii_err));