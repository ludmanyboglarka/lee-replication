function deviance=spatialStroop_CE_SOA(A_01, A_02, A_03, A_04, kc1, kc2, kc3, kc4, ki1, ki2, ki3, ki4, da1, da2, da3, da4, dc1, dc2, dc3, dc4, ...
     a1, a2, a3, a4, Ter1, Ter2, Ter3, Ter4, std, c1, c2, c3, c4)

load('spatialStroop_SOA_ReadData.mat');

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

a=[a1, a2, a3, a4];

c_pred1 = [];
c_pred2 = [];
c_pred3 = [];
c_pred4 = [];
i_pred1 = [];
i_pred2 = [];
i_pred3 = [];
i_pred4 = [];

parfor x=1:Number_of_simulation
    c_pred1(x, :) = RDMC_SOA(A_01, kc1, da1, dc1, c1);
    c_pred2(x, :) = RDMC_SOA(A_02, kc2, da2, dc2, c2);
    c_pred3(x, :) = RDMC_SOA(A_03, kc3, da3, dc3, c3);
    c_pred4(x, :) = RDMC_SOA(A_04, kc4, da4, dc4, c4);
    i_pred1(x, :) = RDMC_SOA(A_01, ki1, -da1, dc1, c1);
    i_pred2(x, :) = RDMC_SOA(A_02, ki2, -da2, dc2, c2);
    i_pred3(x, :) = RDMC_SOA(A_03, ki3, -da3, dc3, c3);
    i_pred4(x, :) = RDMC_SOA(A_04, ki4, -da4, dc4, c4);
end

c_pred={c_pred1, c_pred2, c_pred3, c_pred4};
i_pred={i_pred1, i_pred2, i_pred3, i_pred4};

C_predicted={};
for n=1:4
    for x=1:Number_of_simulation
        if isempty(find(c_pred{1, n}(x,:)<=-a(1, n),1))==0
            C_predicted{1, n}(x, 1) = find(c_pred{1, n}(x,:)<=-a(1, n),1)+R(n, x);
            C_predicted{1, n}(x, 2) = 0;
        elseif isempty(find(c_pred{1, n}(x,:)<=-a(1, n),1))==1 && isempty(find(c_pred{1, n}(x,:)>=a(1, n),1))==0
            C_predicted{1, n}(x, 1) = find(c_pred{1, n}(x,:)>=a(1, n),1)+R(n, x);
            C_predicted{1, n}(x, 2) = 1;
        else
            C_predicted{1, n}(x, 1) = 1500;
            C_predicted{1, n}(x, 2) = 1;
        end
    end
end

for n=1:4
    C_predicted{1, n}(C_predicted{1, n}(:, 1)>=1500, :)=[];
end

I_predicted={};
for n=1:4
    for x=1:Number_of_simulation
        if isempty(find(i_pred{1, n}(x,:)<=-a(1, n),1))==0
            I_predicted{1, n}(x, 1) = find(i_pred{1, n}(x,:)<=-a(1, n),1)+R(n, x);
            I_predicted{1, n}(x, 2) = 0;
        elseif isempty(find(i_pred{1, n}(x,:)<=-a(1, n),1))==1 && isempty(find(i_pred{1, n}(x,:)>=a(1, n),1))==0
            I_predicted{1, n}(x, 1) = find(i_pred{1, n}(x,:)>=a(1, n),1)+R(n, x);
            I_predicted{1, n}(x, 2) = 1;
        else
            I_predicted{1, n}(x, 1) = 1500;
            I_predicted{1, n}(x, 2) = 1;
        end
    end
end

for n=1:4
    I_predicted{1, n}(I_predicted{1, n}(:, 1)>=1500, :)=[];
end

C_acc_pred=[];
I_acc_pred=[];
C_corr_pred={};
C_err_pred={};
I_corr_pred={};
I_err_pred={};
for n=1:4
    C_acc_pred(1, n)=sum(C_predicted{1, n}(:, 2))/numel(C_predicted{1, n}(:, 2));
    I_acc_pred(1, n)=sum(I_predicted{1, n}(:, 2))/numel(I_predicted{1, n}(:, 2));
end

for n=1:4
    C_corr_pred{1, n} = C_predicted{1, n}(C_predicted{1, n}(:, 2)==1, :);
    I_corr_pred{1, n} = I_predicted{1, n}(I_predicted{1, n}(:, 2)==1, :);
    C_err_pred{1, n} = C_predicted{1, n}(C_predicted{1, n}(:, 2)==0, :);
    I_err_pred{1, n} = I_predicted{1, n}(I_predicted{1, n}(:, 2)==0, :);
end

binc=[];
bini=[];
for n=1:4
    binc(n, :) = mean(cell2mat(bin_c(:, n)));
    bini(n, :) = mean(cell2mat(bin_i(:, n)));
end

% Plot the congruent and incongruent CDF
bin_c_err_pred={};
for n=1:4
    bin_c_err_pred{1, n}=C_err_pred{1, n}(:, 1);
end

idx_c_pred={};
RT_c_pred={};
for n=1:4
    idx_c_pred{1, n} = [C_corr_pred{1, n}(:, 1)<=binc(n, 1), C_corr_pred{1, n}(:, 1)>binc(n, 1)&C_corr_pred{1, n}(:, 1)<=binc(n, 2), ...
        C_corr_pred{1, n}(:, 1)>binc(n, 2)&C_corr_pred{1, n}(:, 1)<=binc(n, 3), ...
        C_corr_pred{1, n}(:, 1)>binc(n, 3)&C_corr_pred{1, n}(:, 1)<=binc(n, 4), ...
        C_corr_pred{1, n}(:, 1)>binc(n, 4)&C_corr_pred{1, n}(:, 1)<=binc(n, 5), ...
        C_corr_pred{1, n}(:, 1)>binc(n, 5)&C_corr_pred{1, n}(:, 1)<=binc(n, 6)];
    RT_c_pred{1, n}={C_corr_pred{1, n}(idx_c_pred{1, n}(:, 1), :), C_corr_pred{1, n}(idx_c_pred{1, n}(:, 2), :), ...
        C_corr_pred{1, n}(idx_c_pred{1, n}(:, 3), :), C_corr_pred{1, n}(idx_c_pred{1, n}(:, 4), :), ...
        C_corr_pred{1, n}(idx_c_pred{1, n}(:, 5), :), C_corr_pred{1, n}(idx_c_pred{1, n}(:, 6), :)}; 
    
end 

bin_i_err_pred={};
for n=1:4
    bin_i_err_pred{1, n}=I_err_pred{1, n}(:, 1);
end

idx_i_pred={};
RT_i_pred={};
for n=1:4
    idx_i_pred{1, n} = [I_corr_pred{1, n}(:, 1)<=bini(n, 1), I_corr_pred{1, n}(:, 1)>bini(n, 1)&I_corr_pred{1, n}(:, 1)<=bini(n, 2), ...
        I_corr_pred{1, n}(:, 1)>bini(n, 2)&I_corr_pred{1, n}(:, 1)<=bini(n, 3), ...
        I_corr_pred{1, n}(:, 1)>bini(n, 3)&I_corr_pred{1, n}(:, 1)<=bini(n, 4), ...
        I_corr_pred{1, n}(:, 1)>bini(n, 4)&I_corr_pred{1, n}(:, 1)<=bini(n, 5), ...
        I_corr_pred{1, n}(:, 1)>bini(n, 5)&I_corr_pred{1, n}(:, 1)<=bini(n, 6)];
    RT_i_pred{1, n}={I_corr_pred{1, n}(idx_i_pred{1, n}(:, 1), :), I_corr_pred{1, n}(idx_i_pred{1, n}(:, 2), :), ...
        I_corr_pred{1, n}(idx_i_pred{1, n}(:, 3), :), I_corr_pred{1, n}(idx_i_pred{1, n}(:, 4), :), ...
        I_corr_pred{1, n}(idx_i_pred{1, n}(:, 5), :), I_corr_pred{1, n}(idx_i_pred{1, n}(:, 6), :)}; 
    
end 

C_acc=cell2mat(C_acc);
Num_c={};
frq_c={};
for n=1:4
    Num_c{1, n}= cellfun(@(x) numel(x), vertcat(RT_c{:, n}))/2;
    frq_c{1, n} = Num_c{1, n}.*C_acc(:, n);
end

frq_c_err=[];
for n=1:4
    for j=1:40
        frq_c_err(j, n) = numel(bin_c_err{j, n}(:, 1))*(1-C_acc(j, n));
    end
end
frq_c_err=mean(frq_c_err);

frq_c_tl=[];
frq_c_all=[];
for n=1:4
    frq_c{1, n} = mean(frq_c{1, n});
    frq_c_tl(1, n)=sum(frq_c{1, n});
    frq_c_all(1, n) = sum(frq_c_tl(1, n))+frq_c_err(1, n);
end

I_acc=cell2mat(I_acc);
Num_i={};
frq_i={};
for n=1:4
    Num_i{1, n}= cellfun(@(x) numel(x), vertcat(RT_i{:, n}))/2;
    frq_i{1, n} = Num_i{1, n}.*I_acc(:, n);
end

frq_i_err=[];
for n=1:4
    for j=1:40
        frq_i_err(j, n) = numel(bin_i_err{j, n}(:, 1))*(1-I_acc(j, n));
    end
end
frq_i_err=mean(frq_i_err);

frq_i_tl=[];
frq_i_all=[];
for n=1:4
    frq_i{1, n} = mean(frq_i{1, n});
    frq_i_tl(1, n)=sum(frq_i{1, n});
    frq_i_all(1, n) = sum(frq_i_tl(1, n))+frq_i_err(1, n);
end

Num_c_pred={};
Q_c_pred={};
Q_c_pred_err=[];
for n=1:4
    Num_c_pred{1, n}= cellfun(@(x) numel(x), vertcat(RT_c_pred{:, n}))/2;
    Q_c_pred{1, n} = Num_c_pred{1, n}.*C_acc_pred(:, n);
    Q_c_pred_err(:, n) = numel(bin_c_err_pred(:, n))*(1-C_acc_pred(:, n));
end

Q_c_pred_all=[];
for n=1:4
    Q_c_pred_all(1, n) = sum(Q_c_pred{1, n})+Q_c_pred_err(1, n); 
end

prp_c_pred=[];
prp_c_err_pred=[];

for n=1:4
    prp_c_pred(n, :) = Q_c_pred{1, n}/Q_c_pred_all(1, n);
    prp_c_err_pred(1, n) = Q_c_pred_err(1, n)/Q_c_pred_all(1, n);  
end

frq_c_pred={};
frq_c_pred_err=[];
for n=1:4
    frq_c_pred{1, n} = frq_c_all(1, n).*prp_c_pred(n, :);
    frq_c_pred_err(1, n) = frq_c_all(1, n).*prp_c_err_pred(1, n);
end

chi_c={};
for n=1:4
    chi_c{1, n}=frq_c{1, n}.*log(frq_c{1, n}./frq_c_pred{1, n});
end

Num_i_pred={};
Q_i_pred={};
Q_i_pred_err=[];
for n=1:4
    Num_i_pred{1, n}= cellfun(@(x) numel(x), vertcat(RT_i_pred{:, n}))/2;
    Q_i_pred{1, n} = Num_i_pred{1, n}.*I_acc_pred(:, n);
    Q_i_pred_err(:, n) = numel(bin_i_err_pred(:, n))*(1-I_acc_pred(:, n));
end

Q_i_pred_all=[];
for n=1:4
    Q_i_pred_all(1, n) = sum(Q_i_pred{1, n})+Q_i_pred_err(1, n); 
end

prp_i_pred=[];
prp_i_err_pred=[];
for n=1:4
    prp_i_pred(n, :) = Q_i_pred{1, n}/Q_i_pred_all(1, n);
    prp_i_err_pred(1, n) = Q_i_pred_err(1, n)/Q_i_pred_all(1, n); 
end

frq_i_pred={};
frq_i_pred_err=[];
for n=1:4
    frq_i_pred{1, n} = frq_i_all(1, n).*prp_i_pred(n, :);
    frq_i_pred_err(1, n) = frq_i_all(1, n).*prp_i_err_pred(1, n);
end

chi_i={};
for n=1:4
    chi_i{1, n}=frq_i{1, n}.*log(frq_i{1, n}./frq_i_pred{1, n});
end

chi_c_err=[];
for n=1:4
    if frq_c_err(1, n)==0 || frq_c_pred_err(1, n)==0
        chi_c_err(1, n)=0;
    else 
        chi_c_err(1, n)=frq_c_err(1, n).*log(frq_c_err(1, n)./frq_c_pred_err(1, n));
    end
end

chi_i_err=[];
for n=1:4
    if frq_i_err(1, n)==0 || frq_i_pred_err(1, n)==0
        chi_i_err(1, n)=0;
    else 
        chi_i_err(1, n)=frq_i_err(1, n).*log(frq_i_err(1, n)./frq_i_pred_err(1, n));
    end
end

chi_p=[];
for n=1:4
    chi_p(n, :)=chi_c{1, n}+chi_i{1, n};
end

deviance=2*(sum(sum(chi_p))+sum(chi_c_err)+sum(chi_i_err));