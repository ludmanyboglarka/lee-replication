%% The CDF, delta plot, and CAF of the best performing model (fitting congruency sequence x SOA)

% spatialStroop CDF (SOA manipulation)
clear
% Empirical spatialStroop data
warning off
spatialStroop=readtable('spatialStroopSOA.xlsx');
spatialStroop=table2array(spatialStroop);

NumP= 40;

Ppt = cell(NumP, 1);
for k=1:NumP
    Ppt{k} = spatialStroop(spatialStroop(:,1)==k,:);
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
    cC2{n} = [cC1{n}(:, 3), cC1{n}(:, 5)];
end
cC = cC2;

cI1=cell(NumP, 1);
cI2=cell(NumP, 1);
for n=1:NumP
    l=cI(:, 1)==n;
    cI1{n} = cI(l, :);
    cI2{n} = [cI1{n}(:, 3), cI1{n}(:, 5)];
end
cI=cI2;

iC1=cell(NumP, 1);
iC2=cell(NumP, 1);
for n=1:NumP
    k=iC(:, 1)==n;
    iC1{n} = iC(k, :);
    iC2{n} = [iC1{n}(:, 3), iC1{n}(:, 5)];
end
iC=iC2;

iI1=cell(NumP, 1);
iI2=cell(NumP, 1);
for n=1:NumP
    k=iI(:, 1)==n;
    iI1{n} = iI(k, :);
    iI2{n} = [iI1{n}(:, 3), iI1{n}(:, 5)];
end
iI=iI2;

% We have different SOAs in the spatialStroop task, 100s (target shown before
% distractors), 0, -50, -100s
% Five quantiles for each condition

RTcc = cell(n, 4);
for n=1:NumP
    RTcc{n, 1} = cC{n}(cC{n}(:, 2)==100,1);
    RTcc{n, 2} = cC{n}(cC{n}(:, 2)==0,1);
    RTcc{n, 3} = cC{n}(cC{n}(:, 2)==-50,1);
    RTcc{n, 4} = cC{n}(cC{n}(:, 2)==-100,1);
end

RTci = cell(n, 4);
for n=1:NumP
    RTci{n, 1} = cI{n}(cI{n}(:, 2)==100,1);
    RTci{n, 2} = cI{n}(cI{n}(:, 2)==0,1);
    RTci{n, 3} = cI{n}(cI{n}(:, 2)==-50,1); 
    RTci{n, 4} = cI{n}(cI{n}(:, 2)==-100,1);
end

RTic = cell(n, 4);
for n=1:NumP
    RTic{n, 1} = iC{n}(iC{n}(:, 2)==100,1);
    RTic{n, 2} = iC{n}(iC{n}(:, 2)==0,1);
    RTic{n, 3} = iC{n}(iC{n}(:, 2)==-50,1);
    RTic{n, 4} = iC{n}(iC{n}(:, 2)==-100,1);
end

RTii = cell(n, 4);
for n=1:NumP
    RTii{n, 1} = iI{n}(iI{n}(:, 2)==100,1);
    RTii{n, 2} = iI{n}(iI{n}(:, 2)==0,1);
    RTii{n, 3} = iI{n}(iI{n}(:, 2)==-50,1);
    RTii{n, 4} = iI{n}(iI{n}(:, 2)==-100,1);
end

cC_SOA_1 = cell2mat(RTcc(:, 1));
cC_SOA_2 = cell2mat(RTcc(:, 2));
cC_SOA_3 = cell2mat(RTcc(:, 3));
cC_SOA_4 = cell2mat(RTcc(:, 4));

cI_SOA_1 = cell2mat(RTci(:, 1));
cI_SOA_2 = cell2mat(RTci(:, 2));
cI_SOA_3 = cell2mat(RTci(:, 3));
cI_SOA_4 = cell2mat(RTci(:, 4));

iC_SOA_1 = cell2mat(RTic(:, 1));
iC_SOA_2 = cell2mat(RTic(:, 2));
iC_SOA_3 = cell2mat(RTic(:, 3));
iC_SOA_4 = cell2mat(RTic(:, 4));

iI_SOA_1 = cell2mat(RTii(:, 1));
iI_SOA_2 = cell2mat(RTii(:, 2));
iI_SOA_3 = cell2mat(RTii(:, 3));
iI_SOA_4 = cell2mat(RTii(:, 4));


% Simulated data
% Load the estimated parameters from the best performing model

% Model Prediction
% cC
Number_of_simulation = 50000;

load('spatialStroop_CSE_SOA_SIMPLEX_04.mat')

A_01_pc=fitted_Param(1);
A_01_pi=fitted_Param(1);
kc_pc=fitted_Param(3);
ki_pc=fitted_Param(4);
da1_pc=fitted_Param(5);
da1_pi=fitted_Param(6);
dc=fitted_Param(8);
a_pc=fitted_Param(9);
Ter1=fitted_Param(10);
Ter2=fitted_Param(11);
Ter3=fitted_Param(12);
Ter4=fitted_Param(13);
std=fitted_Param(14);
A_02_pc=fitted_Param(2);
A_03_pc=fitted_Param(2);
A_04_pc=fitted_Param(2);
A_02_pi=fitted_Param(2);
A_03_pi=fitted_Param(2);
A_04_pi=fitted_Param(2);
kc_pi=fitted_Param(3);
ki_pi=fitted_Param(4);
da2_pc=fitted_Param(5);
da3_pc=fitted_Param(5);
da4_pc=fitted_Param(5);
da2_pi=fitted_Param(6);
da3_pi=fitted_Param(7);
da4_pi=fitted_Param(7);
a_pi=fitted_Param(9);
c=0;

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

cC_SOA_1 = cell2mat(RTcc(:, 1));
cC_SOA_2 = cell2mat(RTcc(:, 2));
cC_SOA_3 = cell2mat(RTcc(:, 3));
cC_SOA_4 = cell2mat(RTcc(:, 4));

cI_SOA_1 = cell2mat(RTci(:, 1));
cI_SOA_2 = cell2mat(RTci(:, 2));
cI_SOA_3 = cell2mat(RTci(:, 3));
cI_SOA_4 = cell2mat(RTci(:, 4));

iC_SOA_1 = cell2mat(RTic(:, 1));
iC_SOA_2 = cell2mat(RTic(:, 2));
iC_SOA_3 = cell2mat(RTic(:, 3));
iC_SOA_4 = cell2mat(RTic(:, 4));

iI_SOA_1 = cell2mat(RTii(:, 1));
iI_SOA_2 = cell2mat(RTii(:, 2));
iI_SOA_3 = cell2mat(RTii(:, 3));
iI_SOA_4 = cell2mat(RTii(:, 4));

% CDFs for previous congruent condition
figure('Position', [1, 1, 800, 400])
figure(1)
t=tiledlayout(2, 2);
title(t, 'Previous Congruent');
nexttile
h1 = cdfplot(cC_SOA_1);
hold on
h2 = cdfplot(cI_SOA_1);
hold on 
h3 = cdfplot(cC_predicted{1, 1}(:, 1));
hold on 
h4 = cdfplot(cI_predicted{1, 1}(:, 1));
xlim([200 1000]);
xlabel('RT (ms)');
ylabel('CDF');
set(h1, 'LineStyle', ':', 'Color', [0 0.5 0], 'LineWidth', 1.5);
set(h2, 'LineStyle', ':', 'Color', [.5 0 .5], 'LineWidth', 1.5);
set(h3, 'LineStyle', '-', 'Color', [0 0.5 0], 'LineWidth', 1.5);
set(h4, 'LineStyle', '-', 'Color', [.5 0 .5], 'LineWidth', 1.5);
legend('Empirical cC', 'Empirical cI', 'Predicted cC', 'Predicted cI', 'Location', 'best')
title('SOA = +100 ms')
hold off

nexttile
h3 = cdfplot(cC_SOA_2);
hold on 
h4 = cdfplot(cI_SOA_2);
hold on 
h5 = cdfplot(cC_predicted{1, 2}(:, 1));
hold on 
h6 = cdfplot(cI_predicted{1, 2}(:, 1));
xlim([200 1000]);
xlabel('RT (ms)');
ylabel('CDF');
set(h3, 'LineStyle', ':', 'Color', [0 0.5 0], 'LineWidth', 1.5);
set(h4, 'LineStyle', ':', 'Color', [.5 0 .5], 'LineWidth', 1.5);
set(h5, 'LineStyle', '-', 'Color', [0 0.5 0], 'LineWidth', 1.5);
set(h6, 'LineStyle', '-', 'Color', [.5 0 .5], 'LineWidth', 1.5);
legend('Empirical cC', 'Empirical cI', 'Predicted cC', 'Predicted cI', 'Location', 'best')
title('SOA = 0 ms')
hold off

nexttile
h3 = cdfplot(cC_SOA_3);
hold on 
h4 = cdfplot(cI_SOA_3);
hold on 
h5 = cdfplot(cC_predicted{1, 3}(:, 1));
hold on 
h6 = cdfplot(cI_predicted{1, 3}(:, 1));
xlim([200 1000]);
xlabel('RT (ms)');
ylabel('CDF');
set(h3, 'LineStyle', ':', 'Color', [0 0.5 0], 'LineWidth', 1.5);
set(h4, 'LineStyle', ':', 'Color', [.5 0 .5], 'LineWidth', 1.5);
set(h5, 'LineStyle', '-', 'Color', [0 0.5 0], 'LineWidth', 1.5);
set(h6, 'LineStyle', '-', 'Color', [.5 0 .5], 'LineWidth', 1.5);
legend('Empirical cC', 'Empirical cI', 'Predicted cC', 'Predicted cI', 'Location', 'best')
title('SOA = -50 ms')
hold off

nexttile
h3 = cdfplot(cC_SOA_4);
hold on 
h4 = cdfplot(cI_SOA_4);
hold on 
h5 = cdfplot(cC_predicted{1, 4}(:, 1));
hold on 
h6 = cdfplot(cI_predicted{1, 4}(:, 1));
xlim([200 1000]);
xlabel('RT (ms)');
ylabel('CDF');
set(h3, 'LineStyle', ':', 'Color', [0 0.5 0], 'LineWidth', 1.5);
set(h4, 'LineStyle', ':', 'Color', [.5 0 .5], 'LineWidth', 1.5);
set(h5, 'LineStyle', '-', 'Color', [0 0.5 0], 'LineWidth', 1.5);
set(h6, 'LineStyle', '-', 'Color', [.5 0 .5], 'LineWidth', 1.5);
legend('Empirical cC', 'Empirical cI', 'Predicted cC', 'Predicted cI', 'Location', 'best')
title('SOA = -100 ms')
hold off

% CDFs for previous incongruent condition
figure('Position', [1, 1, 800, 400])
figure(2)
h=tiledlayout(2, 2);
title(h, 'Previous Incongruent');
nexttile
h1 = cdfplot(iC_SOA_1);
hold on
h2 = cdfplot(iI_SOA_1);
hold on 
h3 = cdfplot(iC_predicted{1, 1}(:, 1));
hold on 
h4 = cdfplot(iI_predicted{1, 1}(:, 1));
xlim([200 1000]);
xlabel('RT (ms)');
ylabel('CDF');
set(h1, 'LineStyle', ':', 'Color', [0 0.5 0], 'LineWidth', 1.5);
set(h2, 'LineStyle', ':', 'Color', [.5 0 .5], 'LineWidth', 1.5);
set(h3, 'LineStyle', '-', 'Color', [0 0.5 0], 'LineWidth', 1.5);
set(h4, 'LineStyle', '-', 'Color', [.5 0 .5], 'LineWidth', 1.5);
legend('Empirical iC', 'Empirical iI', 'Predicted iC', 'Predicted iI', 'Location', 'best')
title('SOA = +100 ms')
hold off

nexttile
h3 = cdfplot(iC_SOA_2);
hold on 
h4 = cdfplot(iI_SOA_2);
hold on 
h5 = cdfplot(iC_predicted{1, 2}(:, 1));
hold on 
h6 = cdfplot(iI_predicted{1, 2}(:, 1));
xlim([200 1000]);
xlabel('RT (ms)');
ylabel('CDF');
set(h3, 'LineStyle', ':', 'Color', [0 0.5 0], 'LineWidth', 1.5);
set(h4, 'LineStyle', ':', 'Color', [.5 0 .5], 'LineWidth', 1.5);
set(h5, 'LineStyle', '-', 'Color', [0 0.5 0], 'LineWidth', 1.5);
set(h6, 'LineStyle', '-', 'Color', [.5 0 .5], 'LineWidth', 1.5);
legend('Empirical iC', 'Empirical iI', 'Predicted iC', 'Predicted iI', 'Location', 'best')
title('SOA = 0 ms')
hold off

nexttile
h3 = cdfplot(iC_SOA_3);
hold on 
h4 = cdfplot(iI_SOA_3);
hold on 
h5 = cdfplot(iC_predicted{1, 3}(:, 1));
hold on 
h6 = cdfplot(iI_predicted{1, 3}(:, 1));
xlim([200 1000]);
xlabel('RT (ms)');
ylabel('CDF');
set(h3, 'LineStyle', ':', 'Color', [0 0.5 0], 'LineWidth', 1.5);
set(h4, 'LineStyle', ':', 'Color', [.5 0 .5], 'LineWidth', 1.5);
set(h5, 'LineStyle', '-', 'Color', [0 0.5 0], 'LineWidth', 1.5);
set(h6, 'LineStyle', '-', 'Color', [.5 0 .5], 'LineWidth', 1.5);
legend('Empirical iC', 'Empirical iI', 'Predicted iC', 'Predicted iI', 'Location', 'best')
title('SOA = -50 ms')
hold off

nexttile
h3 = cdfplot(iC_SOA_4);
hold on 
h4 = cdfplot(iI_SOA_4);
hold on 
h5 = cdfplot(iC_predicted{1, 4}(:, 1));
hold on 
h6 = cdfplot(iI_predicted{1, 4}(:, 1));
xlim([200 1000]);
xlabel('RT (ms)');
ylabel('CDF');
set(h3, 'LineStyle', ':', 'Color', [0 0.5 0], 'LineWidth', 1.5);
set(h4, 'LineStyle', ':', 'Color', [.5 0 .5], 'LineWidth', 1.5);
set(h5, 'LineStyle', '-', 'Color', [0 0.5 0], 'LineWidth', 1.5);
set(h6, 'LineStyle', '-', 'Color', [.5 0 .5], 'LineWidth', 1.5);
legend('Empirical iC', 'Empirical iI', 'Predicted iC', 'Predicted iI', 'Location', 'best')
title('SOA = -100 ms')
hold off

exportgraphics(t,'SpatialStroop_CSE_SOA_CDF.pdf')
exportgraphics(h,'SpatialStroop_CSE_SOA_CDF.pdf','Append',true)

%% spatialStroop DP (CSE x SOA manipulation)

clear
% Empirical spatialStroop data
warning off
spatialStroop=readtable('spatialStroopSOA.xlsx');
spatialStroop=table2array(spatialStroop);

NumP= 40;

Ppt = cell(NumP, 1);
for k=1:NumP
    Ppt{k} = spatialStroop(spatialStroop(:,1)==k,:);
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
de=cC(:,4)==0;
cC(de,:) = [];
de3=cC(:,3)<200;
cC(de3,:) = [];
de4=cC(:,3)>1500;
cC(de4,:) = [];

% delete the first trial of each block
de1 = cC(:, 2) ==[1 121 241 361 481 601 721 841];
de1=[find(de1(:, 1)==1);find(de1(:, 2)==1);find(de1(:, 3)==1);find(de1(:, 4)==1);find(de1(:, 5)==1);find(de1(:, 6)==1);find(de1(:, 7)==1);find(de1(:, 8)==1)];
cC(de1,:) = [];

cI = cell2mat(cI);
de=cI(:,4)==0;
cI(de,:) = [];
de7=cI(:,3)<200;
cI(de7,:) = [];
de8=cI(:,3)>1500;
cI(de8,:) = [];
de2 = cI(:, 2) ==[1 121 241 361 481 601 721 841];
de2=[find(de2(:, 1)==1);find(de2(:, 2)==1);find(de2(:, 3)==1);find(de2(:, 4)==1);find(de2(:, 5)==1);find(de2(:, 6)==1);find(de2(:, 7)==1);find(de2(:, 8)==1)];
cI(de2,:) = [];

iC = cell2mat(iC);
de=iC(:,4)==0;
iC(de,:) = [];
de11=iC(:,3)<200;
iC(de11,:) = [];
de12=iC(:,3)>1500;
iC(de12,:) = [];
de3 = iC(:, 2) ==[1 121 241 361 481 601 721 841];
de3=[find(de3(:, 1)==1);find(de3(:, 2)==1);find(de3(:, 3)==1);find(de3(:, 4)==1);find(de3(:, 5)==1);find(de3(:, 6)==1);find(de3(:, 7)==1);find(de3(:, 8)==1)];
iC(de3,:) = [];

iI = cell2mat(iI);
de=iI(:,4)==0;
iI(de,:) = [];
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
    cC2{n} = [cC1{n}(:, 3), cC1{n}(:, 5)];
end
cC = cC2;

cI1=cell(NumP, 1);
cI2=cell(NumP, 1);
for n=1:NumP
    l=cI(:, 1)==n;
    cI1{n} = cI(l, :);
    cI2{n} = [cI1{n}(:, 3), cI1{n}(:, 5)];
end
cI=cI2;

iC1=cell(NumP, 1);
iC2=cell(NumP, 1);
for n=1:NumP
    k=iC(:, 1)==n;
    iC1{n} = iC(k, :);
    iC2{n} = [iC1{n}(:, 3), iC1{n}(:, 5)];
end
iC=iC2;

iI1=cell(NumP, 1);
iI2=cell(NumP, 1);
for n=1:NumP
    k=iI(:, 1)==n;
    iI1{n} = iI(k, :);
    iI2{n} = [iI1{n}(:, 3), iI1{n}(:, 5)];
end
iI=iI2;

% We have different SOAs in the spatialStroop task, -100s (target shown before
% distractors), 0, 50, 100s
% Five quantiles for each condition

RTcc = cell(n, 4);
for n=1:NumP
    RTcc{n, 1} = cC{n}(cC{n}(:, 2)==100,1);
    RTcc{n, 1} = [prctile(RTcc{n, 1}(:,1),10), prctile(RTcc{n, 1}(:,1),30), prctile(RTcc{n, 1}(:,1),50), prctile(RTcc{n, 1}(:,1),70), ...
        prctile(RTcc{n, 1}(:,1),90)];
    RTcc{n, 2} = cC{n}(cC{n}(:, 2)==0,1);
    RTcc{n, 2} = [prctile(RTcc{n, 2}(:,1),10), prctile(RTcc{n, 2}(:,1),30), prctile(RTcc{n, 2}(:,1),50), prctile(RTcc{n, 2}(:,1),70), ...
        prctile(RTcc{n, 2}(:,1),90)];
    RTcc{n, 3} = cC{n}(cC{n}(:, 2)==-50,1);
    RTcc{n, 3} = [prctile(RTcc{n, 3}(:,1),10), prctile(RTcc{n, 3}(:,1),30), prctile(RTcc{n, 3}(:,1),50), prctile(RTcc{n, 3}(:,1),70), ...
        prctile(RTcc{n, 3}(:,1),90)];
    RTcc{n, 4} = cC{n}(cC{n}(:, 2)==-100,1);
    RTcc{n, 4} = [prctile(RTcc{n, 4}(:,1),10), prctile(RTcc{n, 4}(:,1),30), prctile(RTcc{n, 4}(:,1),50), prctile(RTcc{n, 4}(:,1),70), ...
        prctile(RTcc{n, 4}(:,1),90)];
end

RTci = cell(n, 4);
for n=1:NumP
    RTci{n, 1} = cI{n}(cI{n}(:, 2)==100,1);
    RTci{n, 1} = [prctile(RTci{n, 1}(:,1),10), prctile(RTci{n, 1}(:,1),30), prctile(RTci{n, 1}(:,1),50), prctile(RTci{n, 1}(:,1),70), ...
        prctile(RTci{n, 1}(:,1),90)];
    RTci{n, 2} = cI{n}(cI{n}(:, 2)==0,1);
    RTci{n, 2} = [prctile(RTci{n, 2}(:,1),10), prctile(RTci{n, 2}(:,1),30), prctile(RTci{n, 2}(:,1),50), prctile(RTci{n, 2}(:,1),70), ...
        prctile(RTci{n, 2}(:,1),90)];
    RTci{n, 3} = cI{n}(cI{n}(:, 2)==-50,1);
    RTci{n, 3} = [prctile(RTci{n, 3}(:,1),10), prctile(RTci{n, 3}(:,1),30), prctile(RTci{n, 3}(:,1),50), prctile(RTci{n, 3}(:,1),70), ...
        prctile(RTci{n, 3}(:,1),90)];
    RTci{n, 4} = cI{n}(cI{n}(:, 2)==-100,1);
    RTci{n, 4} = [prctile(RTci{n, 4}(:,1),10), prctile(RTci{n, 4}(:,1),30), prctile(RTci{n, 4}(:,1),50), prctile(RTci{n, 4}(:,1),70), ...
        prctile(RTci{n, 4}(:,1),90)];
end

RTic = cell(n, 4);
for n=1:NumP
    RTic{n, 1} = iC{n}(iC{n}(:, 2)==100,1);
    RTic{n, 1} = [prctile(RTic{n, 1}(:,1),10), prctile(RTic{n, 1}(:,1),30), prctile(RTic{n, 1}(:,1),50), prctile(RTic{n, 1}(:,1),70), ...
        prctile(RTic{n, 1}(:,1),90)];
    RTic{n, 2} = iC{n}(iC{n}(:, 2)==0,1);
    RTic{n, 2} = [prctile(RTic{n, 2}(:,1),10), prctile(RTic{n, 2}(:,1),30), prctile(RTic{n, 2}(:,1),50), prctile(RTic{n, 2}(:,1),70), ...
        prctile(RTic{n, 2}(:,1),90)];
    RTic{n, 3} = iC{n}(iC{n}(:, 2)==-50,1);
    RTic{n, 3} = [prctile(RTic{n, 3}(:,1),10), prctile(RTic{n, 3}(:,1),30), prctile(RTic{n, 3}(:,1),50), prctile(RTic{n, 3}(:,1),70), ...
        prctile(RTic{n, 3}(:,1),90)];
    RTic{n, 4} = iC{n}(iC{n}(:, 2)==-100,1);
    RTic{n, 4} = [prctile(RTic{n, 4}(:,1),10), prctile(RTic{n, 4}(:,1),30), prctile(RTic{n, 4}(:,1),50), prctile(RTic{n, 4}(:,1),70), ...
        prctile(RTic{n, 4}(:,1),90)];
end

RTii = cell(n, 4);
for n=1:NumP
    RTii{n, 1} = iI{n}(iI{n}(:, 2)==100,1);
    RTii{n, 1} = [prctile(RTii{n, 1}(:,1),10), prctile(RTii{n, 1}(:,1),30), prctile(RTii{n, 1}(:,1),50), prctile(RTii{n, 1}(:,1),70), ...
        prctile(RTii{n, 1}(:,1),90)];
    RTii{n, 2} = iI{n}(iI{n}(:, 2)==0,1);
    RTii{n, 2} = [prctile(RTii{n, 2}(:,1),10), prctile(RTii{n, 2}(:,1),30), prctile(RTii{n, 2}(:,1),50), prctile(RTii{n, 2}(:,1),70), ...
        prctile(RTii{n, 2}(:,1),90)];
    RTii{n, 3} = iI{n}(iI{n}(:, 2)==-50,1);
    RTii{n, 3} = [prctile(RTii{n, 3}(:,1),10), prctile(RTii{n, 3}(:,1),30), prctile(RTii{n, 3}(:,1),50), prctile(RTii{n, 3}(:,1),70), ...
        prctile(RTii{n, 3}(:,1),90)];
    RTii{n, 4} = iI{n}(iI{n}(:, 2)==-100,1);
    RTii{n, 4} = [prctile(RTii{n, 4}(:,1),10), prctile(RTii{n, 4}(:,1),30), prctile(RTii{n, 4}(:,1),50), prctile(RTii{n, 4}(:,1),70), ...
        prctile(RTii{n, 4}(:,1),90)];
end

CE_c=cell(NumP, 4);
for n=1:NumP
    CE_c{n, 1} = [RTci{n, 1}(1) - RTcc{n, 1}(1), RTci{n, 1}(2) - RTcc{n, 1}(2), RTci{n, 1}(3) - RTcc{n, 1}(3), RTci{n, 1}(4) - RTcc{n, 1}(4), ...
        RTci{n, 1}(5) - RTcc{n, 1}(5)];
    CE_c{n, 2} = [RTci{n, 2}(1) - RTcc{n, 2}(1), RTci{n, 2}(2) - RTcc{n, 2}(2), RTci{n, 2}(3) - RTcc{n, 2}(3), RTci{n, 2}(4) - RTcc{n, 2}(4), ...
        RTci{n, 2}(5) - RTcc{n, 2}(5)];
    CE_c{n, 3} = [RTci{n, 3}(1) - RTcc{n, 3}(1), RTci{n, 3}(2) - RTcc{n, 3}(2), RTci{n, 3}(3) - RTcc{n, 3}(3), RTci{n, 3}(4) - RTcc{n, 3}(4), ...
        RTci{n, 3}(5) - RTcc{n, 3}(5)];
    CE_c{n, 4} = [RTci{n, 4}(1) - RTcc{n, 4}(1), RTci{n, 4}(2) - RTcc{n, 4}(2), RTci{n, 4}(3) - RTcc{n, 4}(3), RTci{n, 4}(4) - RTcc{n, 4}(4), ...
        RTci{n, 4}(5) - RTcc{n, 4}(5)];
end

RT_cc_1=[];
RT_cc_2=[];
RT_cc_3=[];
RT_cc_4=[];
for n=1:NumP
    RT_cc_1=[RT_cc_1 ; RTcc{n, 1}];
    RT_cc_2=[RT_cc_2 ; RTcc{n, 2}];
    RT_cc_3=[RT_cc_3 ; RTcc{n, 3}];
    RT_cc_4=[RT_cc_4 ; RTcc{n, 4}];
end
meanRTcc=[mean(RT_cc_1) ; mean(RT_cc_2); mean(RT_cc_3); mean(RT_cc_4)];

RT_ci_1=[];
RT_ci_2=[];
RT_ci_3=[];
RT_ci_4=[];
for n=1:NumP
    RT_ci_1=[RT_ci_1 ; RTci{n, 1}];
    RT_ci_2=[RT_ci_2 ; RTci{n, 2}];
    RT_ci_3=[RT_ci_3 ; RTci{n, 3}];
    RT_ci_4=[RT_ci_4 ; RTci{n, 4}];
end
meanRTci=[mean(RT_ci_1) ; mean(RT_ci_2); mean(RT_ci_3); mean(RT_ci_4)];

RT_c=[];
for n=1:4
    RT_c = [RT_c ; (meanRTcc(n, 1)+meanRTci(n, 1))/2, (meanRTcc(n, 2)+meanRTci(n, 2))/2, (meanRTcc(n, 3)+meanRTci(n, 3))/2, (meanRTcc(n, 4)+meanRTci(n, 4))/2, ...
        (meanRTcc(n, 5)+meanRTci(n, 5))/2];
end

CE_i=cell(NumP, 4);
for n=1:NumP
    CE_i{n, 1} = [RTii{n, 1}(1) - RTic{n, 1}(1), RTii{n, 1}(2) - RTic{n, 1}(2), RTii{n, 1}(3) - RTic{n, 1}(3), RTii{n, 1}(4) - RTic{n, 1}(4), ...
        RTii{n, 1}(5) - RTic{n, 1}(5)];
    CE_i{n, 2} = [RTii{n, 2}(1) - RTic{n, 2}(1), RTii{n, 2}(2) - RTic{n, 2}(2), RTii{n, 2}(3) - RTic{n, 2}(3), RTii{n, 2}(4) - RTic{n, 2}(4), ...
        RTii{n, 2}(5) - RTic{n, 2}(5)];
    CE_i{n, 3} = [RTii{n, 3}(1) - RTic{n, 3}(1), RTii{n, 3}(2) - RTic{n, 3}(2), RTii{n, 3}(3) - RTic{n, 3}(3), RTii{n, 3}(4) - RTic{n, 3}(4), ...
        RTii{n, 3}(5) - RTic{n, 3}(5)];
    CE_i{n, 4} = [RTii{n, 4}(1) - RTic{n, 4}(1), RTii{n, 4}(2) - RTic{n, 4}(2), RTii{n, 4}(3) - RTic{n, 4}(3), RTii{n, 4}(4) - RTic{n, 4}(4), ...
        RTii{n, 4}(5) - RTic{n, 4}(5)];
end

RT_ic_1=[];
RT_ic_2=[];
RT_ic_3=[];
RT_ic_4=[];
for n=1:NumP
    RT_ic_1=[RT_ic_1 ; RTic{n, 1}];
    RT_ic_2=[RT_ic_2 ; RTic{n, 2}];
    RT_ic_3=[RT_ic_3 ; RTic{n, 3}];
    RT_ic_4=[RT_ic_4 ; RTic{n, 4}];
end
meanRTic=[mean(RT_ic_1) ; mean(RT_ic_2); mean(RT_ic_3); mean(RT_ic_4)];

RT_ii_1=[];
RT_ii_2=[];
RT_ii_3=[];
RT_ii_4=[];
for n=1:NumP
    RT_ii_1=[RT_ii_1 ; RTii{n, 1}];
    RT_ii_2=[RT_ii_2 ; RTii{n, 2}];
    RT_ii_3=[RT_ii_3 ; RTii{n, 3}];
    RT_ii_4=[RT_ii_4 ; RTii{n, 4}];
end
meanRTii=[mean(RT_ii_1) ; mean(RT_ii_2); mean(RT_ii_3); mean(RT_ii_4)];

RT_i=[];
for n=1:4
    RT_i = [RT_i ; (meanRTic(n, 1)+meanRTii(n, 1))/2, (meanRTic(n, 2)+meanRTii(n, 2))/2, (meanRTic(n, 3)+meanRTii(n, 3))/2, (meanRTic(n, 4)+meanRTii(n, 4))/2, ...
        (meanRTic(n, 5)+meanRTii(n, 5))/2];
end

CE_1_c=[];
CE_2_c=[];
CE_3_c=[];
CE_4_c=[];
for n=1:NumP
    CE_1_c=[CE_1_c ; CE_c{n, 1}];
    CE_2_c=[CE_2_c ; CE_c{n, 2}];
    CE_3_c=[CE_3_c ; CE_c{n, 3}];
    CE_4_c=[CE_4_c ; CE_c{n, 4}];
end

meanCE_c=[mean(CE_1_c) ; mean(CE_2_c); mean(CE_3_c); mean(CE_4_c)];

CE_1_i=[];
CE_2_i=[];
CE_3_i=[];
CE_4_i=[];
for n=1:NumP
    CE_1_i=[CE_1_i ; CE_i{n, 1}];
    CE_2_i=[CE_2_i ; CE_i{n, 2}];
    CE_3_i=[CE_3_i ; CE_i{n, 3}];
    CE_4_i=[CE_4_i ; CE_i{n, 4}];
end

meanCE_i=[mean(CE_1_i) ; mean(CE_2_i); mean(CE_3_i); mean(CE_4_i)];

SEM_c = [std(CE_1_c)/sqrt(length(CE_1_c)); std(CE_2_c)/sqrt(length(CE_2_c)); std(CE_3_c)/sqrt(length(CE_3_c)); std(CE_4_c)/sqrt(length(CE_4_c))];
ts_c = tinv([0.025  0.975],length(CE_1_c)-1);

inv_c={};
for r=1:4
    inv_c{1, r} = ts_c.*SEM_c(r, :)';
end

CI_c={};
for r=1:4
    for j=1:5
        CI_c{1, r}(j, :) = meanCE_c(r, j) + inv_c{1, r}(j, :);
    end
end

SEM_i = [std(CE_1_i)/sqrt(length(CE_1_i)); std(CE_2_i)/sqrt(length(CE_2_i)); std(CE_3_i)/sqrt(length(CE_3_i)); std(CE_4_i)/sqrt(length(CE_4_i))];
ts_i = tinv([0.025  0.975],length(CE_1_i)-1);

inv_i={};
for r=1:4
    inv_i{1, r} = ts_i.*SEM_i(r, :)';
end

CI_i={};
for r=1:4
    for j=1:5
        CI_i{1, r}(j, :) = meanCE_i(r, j) + inv_i{1, r}(j, :);
    end
end

% Simulated data
% Load the estimated parameters from the best performing model

% Model Prediction
% cC
Number_of_simulation = 50000;

load('spatialStroop_CSE_SOA_SIMPLEX_04.mat')

A_01_pc=fitted_Param(1);
A_01_pi=fitted_Param(1);
kc_pc=fitted_Param(3);
ki_pc=fitted_Param(4);
da1_pc=fitted_Param(5);
da1_pi=fitted_Param(6);
dc=fitted_Param(8);
a_pc=fitted_Param(9);
Ter1=fitted_Param(10);
Ter2=fitted_Param(11);
Ter3=fitted_Param(12);
Ter4=fitted_Param(13);
std=fitted_Param(14);
A_02_pc=fitted_Param(2);
A_03_pc=fitted_Param(2);
A_04_pc=fitted_Param(2);
A_02_pi=fitted_Param(2);
A_03_pi=fitted_Param(2);
A_04_pi=fitted_Param(2);
kc_pi=fitted_Param(3);
ki_pi=fitted_Param(4);
da2_pc=fitted_Param(5);
da3_pc=fitted_Param(5);
da4_pc=fitted_Param(5);
da2_pi=fitted_Param(6);
da3_pi=fitted_Param(7);
da4_pi=fitted_Param(7);
a_pi=fitted_Param(9);
c=0;

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
    cC_predicted{1, n}(cC_predicted{1, n}(:, 2)==0, :)=[];
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
    cI_predicted{1, n}(cI_predicted{1, n}(:, 2)==0, :)=[];
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
    iC_predicted{1, n}(iC_predicted{1, n}(:, 2)==0, :)=[];
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
    iI_predicted{1, n}(iI_predicted{1, n}(:, 2)==0, :)=[];
end

RTcc_predicted={};
for n=1:4
    RTcc_predicted{1, n} = [prctile(cC_predicted{1, n}(:,1),10), prctile(cC_predicted{1, n}(:,1),30), prctile(cC_predicted{1, n}(:,1),50),...
            prctile(cC_predicted{1, n}(:,1),70), prctile(cC_predicted{1, n}(:,1),90)];
end

RTci_predicted={};
for n=1:4
    RTci_predicted{1, n} = [prctile(cI_predicted{1, n}(:,1),10), prctile(cI_predicted{1, n}(:,1),30), prctile(cI_predicted{1, n}(:,1),50),...
            prctile(cI_predicted{1, n}(:,1),70), prctile(cI_predicted{1, n}(:,1),90)];
end

CE_c_predicted=[];
for n=1:4
    CE_c_predicted(n, :) = [RTci_predicted{1, n}(1) - RTcc_predicted{1, n}(1), RTci_predicted{1, n}(2) - RTcc_predicted{1, n}(2), ...
        RTci_predicted{1, n}(3) - RTcc_predicted{1, n}(3), RTci_predicted{1, n}(4) - RTcc_predicted{1, n}(4), RTci_predicted{1, n}(5) - RTcc_predicted{1, n}(5)];
end

RT_c_predicted=[];
for n=1:4
    RT_c_predicted(n, :) = (RTcc_predicted{1, n}+RTci_predicted{1, n})/2;
end

RTic_predicted={};
for n=1:4
    RTic_predicted{1, n} = [prctile(iC_predicted{1, n}(:,1),10), prctile(iC_predicted{1, n}(:,1),30), prctile(iC_predicted{1, n}(:,1),50),...
            prctile(iC_predicted{1, n}(:,1),70), prctile(iC_predicted{1, n}(:,1),90)];
end

RTii_predicted={};
for n=1:4
    RTii_predicted{1, n} = [prctile(iI_predicted{1, n}(:,1),10), prctile(iI_predicted{1, n}(:,1),30), prctile(iI_predicted{1, n}(:,1),50),...
            prctile(iI_predicted{1, n}(:,1),70), prctile(iI_predicted{1, n}(:,1),90)];
end

CE_i_predicted=[];
for n=1:4
    CE_i_predicted(n, :) = [RTii_predicted{1, n}(1) - RTic_predicted{1, n}(1), RTii_predicted{1, n}(2) - RTic_predicted{1, n}(2), ...
        RTii_predicted{1, n}(3) - RTic_predicted{1, n}(3), RTii_predicted{1, n}(4) - RTic_predicted{1, n}(4), RTii_predicted{1, n}(5) - RTic_predicted{1, n}(5)];
end

RT_i_predicted=[];
for n=1:4
    RT_i_predicted(n, :) = (RTic_predicted{1, n}+RTii_predicted{1, n})/2;
end

% Plot the delta plot of each SOA 
figure('Position', [300 300 1000 700])
t= tiledlayout(2,2);
nexttile
g1 = plot(RT_c(1, :), meanCE_c(1, :));
hold on 
h1 = plot(RT_i(1, :), meanCE_i(1, :));
hold on
g2 = plot(RT_c_predicted(1, :), CE_c_predicted(1, :));
hold on 
h2 = plot(RT_i_predicted(1, :), CE_i_predicted(1, :));
set(g1, 'LineStyle', ':', 'Marker', 'o', 'Color', [0 0.5 0], 'LineWidth', 1.5);
set(g2, 'LineStyle', '-', 'Marker', 'o', 'Color', [0 0.5 0], 'LineWidth', 1.5);
set(h1, 'LineStyle', ':', 'Marker', 'o', 'Color', [.5 0 .5], 'LineWidth', 1.5);
set(h2, 'LineStyle', '-', 'Marker', 'o', 'Color', [.5 0 .5], 'LineWidth', 1.5);
legend([g1 g2 h1 h2], 'Empirical PC', 'Predicted PC', 'Empirical PI', 'Predicted PI', 'Location', 'best')
xlim([350 650]);
ylim([0 120]);
xlabel('RT (ms)')
ylabel('Congruency effects (ms)')
title('SOA = +100 ms');
hold off

nexttile
g1 = plot(RT_c(2, :), meanCE_c(2, :));
hold on 
h1 = plot(RT_i(2, :), meanCE_i(2, :));
hold on
g2 = plot(RT_c_predicted(2, :), CE_c_predicted(2, :));
hold on 
h2 = plot(RT_i_predicted(2, :), CE_i_predicted(2, :));
set(g1, 'LineStyle', ':', 'Marker', 'o', 'Color', [0 0.5 0], 'LineWidth', 1.5);
set(g2, 'LineStyle', '-', 'Marker', 'o', 'Color', [0 0.5 0], 'LineWidth', 1.5);
set(h1, 'LineStyle', ':', 'Marker', 'o', 'Color', [.5 0 .5], 'LineWidth', 1.5);
set(h2, 'LineStyle', '-', 'Marker', 'o', 'Color', [.5 0 .5], 'LineWidth', 1.5);
legend([g1 g2 h1 h2], 'Empirical PC', 'Predicted PC', 'Empirical PI', 'Predicted PI', 'Location', 'best')
xlim([350 650]);
ylim([0 120]);
xlabel('RT (ms)')
ylabel('Congruency effects (ms)')
title('SOA = 0 ms');
hold off

nexttile
g1 = plot(RT_c(3, :), meanCE_c(3, :));
hold on 
h1 = plot(RT_i(3, :), meanCE_i(3, :));
hold on
g2 = plot(RT_c_predicted(3, :), CE_c_predicted(3, :));
hold on 
h2 = plot(RT_i_predicted(3, :), CE_i_predicted(3, :));
set(g1, 'LineStyle', ':', 'Marker', 'o', 'Color', [0 0.5 0], 'LineWidth', 1.5);
set(g2, 'LineStyle', '-', 'Marker', 'o', 'Color', [0 0.5 0], 'LineWidth', 1.5);
set(h1, 'LineStyle', ':', 'Marker', 'o', 'Color', [.5 0 .5], 'LineWidth', 1.5);
set(h2, 'LineStyle', '-', 'Marker', 'o', 'Color', [.5 0 .5], 'LineWidth', 1.5);
legend([g1 g2 h1 h2], 'Empirical PC', 'Predicted PC', 'Empirical PI', 'Predicted PI', 'Location', 'best')
xlim([350 650]);
ylim([0 120]);
xlabel('RT (ms)')
ylabel('Congruency effects (ms)')
title('SOA = -50 ms');
hold off

nexttile
g1 = plot(RT_c(4, :), meanCE_c(4, :));
hold on 
h1 = plot(RT_i(4, :), meanCE_i(4, :));
hold on
g2 = plot(RT_c_predicted(4, :), CE_c_predicted(4, :));
hold on 
h2 = plot(RT_i_predicted(4, :), CE_i_predicted(4, :));
set(g1, 'LineStyle', ':', 'Marker', 'o', 'Color', [0 0.5 0], 'LineWidth', 1.5);
set(g2, 'LineStyle', '-', 'Marker', 'o', 'Color', [0 0.5 0], 'LineWidth', 1.5);
set(h1, 'LineStyle', ':', 'Marker', 'o', 'Color', [.5 0 .5], 'LineWidth', 1.5);
set(h2, 'LineStyle', '-', 'Marker', 'o', 'Color', [.5 0 .5], 'LineWidth', 1.5);
legend([g1 g2 h1 h2], 'Empirical PC', 'Predicted PC', 'Empirical PI', 'Predicted PI', 'Location', 'best')
xlim([350 650]);
ylim([0 120]);
xlabel('RT (ms)')
ylabel('Congruency effects (ms)')
title('SOA = -100 ms');
hold off

exportgraphics(t,'SpatialStroop_CSE_SOA_DP.pdf')

%% spatialStroop CAF (CSE x SOA)
clear
% Empirical spatialStroop data
warning off
spatialStroop=readtable('spatialStroopSOA.xlsx');
spatialStroop=table2array(spatialStroop);

NumP= 40;

Ppt = cell(NumP, 1);
for k=1:NumP
    Ppt{k} = spatialStroop(spatialStroop(:,1)==k,:);
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

% We have different SOAs in the spatialStroop task, 100s (target shown before
% distractors), 0, -50, -100s
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

bin_cc = cell(NumP, 4);
for n=1:NumP
    bin_cc{n, 1}=[prctile(RTcc{n, 1}(:,1),20), prctile(RTcc{n, 1}(:,1),40), prctile(RTcc{n, 1}(:,1),60), prctile(RTcc{n, 1}(:,1),80), prctile(RTcc{n, 1}(:,1),100)];
    bin_cc{n, 2}=[prctile(RTcc{n, 2}(:,1),20), prctile(RTcc{n, 2}(:,1),40), prctile(RTcc{n, 2}(:,1),60), prctile(RTcc{n, 2}(:,1),80), prctile(RTcc{n, 2}(:,1),100)];
    bin_cc{n, 3}=[prctile(RTcc{n, 3}(:,1),20), prctile(RTcc{n, 3}(:,1),40), prctile(RTcc{n, 3}(:,1),60), prctile(RTcc{n, 3}(:,1),80), prctile(RTcc{n, 3}(:,1),100)];
    bin_cc{n, 4}=[prctile(RTcc{n, 4}(:,1),20), prctile(RTcc{n, 4}(:,1),40), prctile(RTcc{n, 4}(:,1),60), prctile(RTcc{n, 4}(:,1),80), prctile(RTcc{n, 4}(:,1),100)];
end

corrcc = cell(NumP, 4);
incorrcc = cell(NumP, 4);

for n=1:NumP 
    for m=1:4
        for x=1:size(RTcc{n, m},1)
            corrcc{n, m} = [corrcc{n, m} ; RTcc{n, m}(x,1)<=bin_cc{n, m}(:,1) && RTcc{n, m}(x,2)==1, ...
               RTcc{n, m}(x,1)>bin_cc{n, m}(:,1) & RTcc{n, m}(x,1)<=bin_cc{n, m}(:,2) & RTcc{n, m}(x,2)==1, ...
               RTcc{n, m}(x,1)>bin_cc{n, m}(:,2) & RTcc{n, m}(x,1)<=bin_cc{n, m}(:,3) & RTcc{n, m}(x,2)==1, ...
               RTcc{n, m}(x,1)>bin_cc{n, m}(:,3) & RTcc{n, m}(x,1)<=bin_cc{n, m}(:,4) & RTcc{n, m}(x,2)==1, ...
               RTcc{n, m}(x,1)>bin_cc{n, m}(:,4) & RTcc{n, m}(x,1)<=bin_cc{n, m}(:,5) & RTcc{n, m}(x,2)==1];
            incorrcc{n, m} = [incorrcc{n, m} ; RTcc{n, m}(x,1)<=bin_cc{n, m}(:,1) && RTcc{n, m}(x,2)~=1, ...
               RTcc{n, m}(x,1)>bin_cc{n, m}(:,1) & RTcc{n, m}(x,1)<=bin_cc{n, m}(:,2) & RTcc{n, m}(x,2)~=1, ...
               RTcc{n, m}(x,1)>bin_cc{n, m}(:,2) & RTcc{n, m}(x,1)<=bin_cc{n, m}(:,3) & RTcc{n, m}(x,2)~=1, ...
               RTcc{n, m}(x,1)>bin_cc{n, m}(:,3) & RTcc{n, m}(x,1)<=bin_cc{n, m}(:,4) & RTcc{n, m}(x,2)~=1, ...
               RTcc{n, m}(x,1)>bin_cc{n, m}(:,4) & RTcc{n, m}(x,1)<=bin_cc{n, m}(:,5) & RTcc{n, m}(x,2)~=1];  
        end
    end
end

bin_ci = cell(NumP, 4);
for n=1:NumP
    bin_ci{n, 1}=[prctile(RTci{n, 1}(:,1),20), prctile(RTci{n, 1}(:,1),40), prctile(RTci{n, 1}(:,1),60), prctile(RTci{n, 1}(:,1),80), prctile(RTci{n, 1}(:,1),100)];
    bin_ci{n, 2}=[prctile(RTci{n, 2}(:,1),20), prctile(RTci{n, 2}(:,1),40), prctile(RTci{n, 2}(:,1),60), prctile(RTci{n, 2}(:,1),80), prctile(RTci{n, 2}(:,1),100)];
    bin_ci{n, 3}=[prctile(RTci{n, 3}(:,1),20), prctile(RTci{n, 3}(:,1),40), prctile(RTci{n, 3}(:,1),60), prctile(RTci{n, 3}(:,1),80), prctile(RTci{n, 3}(:,1),100)];
    bin_ci{n, 4}=[prctile(RTci{n, 4}(:,1),20), prctile(RTci{n, 4}(:,1),40), prctile(RTci{n, 4}(:,1),60), prctile(RTci{n, 4}(:,1),80), prctile(RTci{n, 4}(:,1),100)];
end

corrci = cell(NumP, 4);
incorrci = cell(NumP, 4);

for n=1:NumP 
    for m=1:4
        for x=1:size(RTci{n, m},1)
            corrci{n, m} = [corrci{n, m} ; RTci{n, m}(x,1)<=bin_ci{n, m}(:,1) && RTci{n, m}(x,2)==1, ...
               RTci{n, m}(x,1)>bin_ci{n, m}(:,1) & RTci{n, m}(x,1)<=bin_ci{n, m}(:,2) & RTci{n, m}(x,2)==1, ...
               RTci{n, m}(x,1)>bin_ci{n, m}(:,2) & RTci{n, m}(x,1)<=bin_ci{n, m}(:,3) & RTci{n, m}(x,2)==1, ...
               RTci{n, m}(x,1)>bin_ci{n, m}(:,3) & RTci{n, m}(x,1)<=bin_ci{n, m}(:,4) & RTci{n, m}(x,2)==1, ...
               RTci{n, m}(x,1)>bin_ci{n, m}(:,4) & RTci{n, m}(x,1)<=bin_ci{n, m}(:,5) & RTci{n, m}(x,2)==1];
            incorrci{n, m} = [incorrci{n, m} ; RTci{n, m}(x,1)<=bin_ci{n, m}(:,1) && RTci{n, m}(x,2)~=1, ...
               RTci{n, m}(x,1)>bin_ci{n, m}(:,1) & RTci{n, m}(x,1)<=bin_ci{n, m}(:,2) & RTci{n, m}(x,2)~=1, ...
               RTci{n, m}(x,1)>bin_ci{n, m}(:,2) & RTci{n, m}(x,1)<=bin_ci{n, m}(:,3) & RTci{n, m}(x,2)~=1, ...
               RTci{n, m}(x,1)>bin_ci{n, m}(:,3) & RTci{n, m}(x,1)<=bin_ci{n, m}(:,4) & RTci{n, m}(x,2)~=1, ...
               RTci{n, m}(x,1)>bin_ci{n, m}(:,4) & RTci{n, m}(x,1)<=bin_ci{n, m}(:,5) & RTci{n, m}(x,2)~=1];
            
        end
    end
end

p_cc=[];
p_cc(1, :) = sum(cell2mat(corrcc(:, 1)))./(sum(cell2mat(corrcc(:, 1)))+sum(cell2mat(incorrcc(:, 1))));
p_cc(2, :) = sum(cell2mat(corrcc(:, 2)))./(sum(cell2mat(corrcc(:, 2)))+sum(cell2mat(incorrcc(:, 2))));
p_cc(3, :) = sum(cell2mat(corrcc(:, 3)))./(sum(cell2mat(corrcc(:, 3)))+sum(cell2mat(incorrcc(:, 3))));
p_cc(4, :) = sum(cell2mat(corrcc(:, 4)))./(sum(cell2mat(corrcc(:, 4)))+sum(cell2mat(incorrcc(:, 4))));

p_ci=[];
p_ci(1, :) = sum(cell2mat(corrci(:, 1)))./(sum(cell2mat(corrci(:, 1)))+sum(cell2mat(incorrci(:, 1))));
p_ci(2, :) = sum(cell2mat(corrci(:, 2)))./(sum(cell2mat(corrci(:, 2)))+sum(cell2mat(incorrci(:, 2))));
p_ci(3, :) = sum(cell2mat(corrci(:, 3)))./(sum(cell2mat(corrci(:, 3)))+sum(cell2mat(incorrci(:, 3))));
p_ci(4, :) = sum(cell2mat(corrci(:, 4)))./(sum(cell2mat(corrci(:, 4)))+sum(cell2mat(incorrci(:, 4))));

bin_ic = cell(NumP, 4);
for n=1:NumP
    bin_ic{n, 1}=[prctile(RTic{n, 1}(:,1),20), prctile(RTic{n, 1}(:,1),40), prctile(RTic{n, 1}(:,1),60), prctile(RTic{n, 1}(:,1),80), prctile(RTic{n, 1}(:,1),100)];
    bin_ic{n, 2}=[prctile(RTic{n, 2}(:,1),20), prctile(RTic{n, 2}(:,1),40), prctile(RTic{n, 2}(:,1),60), prctile(RTic{n, 2}(:,1),80), prctile(RTic{n, 2}(:,1),100)];
    bin_ic{n, 3}=[prctile(RTic{n, 3}(:,1),20), prctile(RTic{n, 3}(:,1),40), prctile(RTic{n, 3}(:,1),60), prctile(RTic{n, 3}(:,1),80), prctile(RTic{n, 3}(:,1),100)];
    bin_ic{n, 4}=[prctile(RTic{n, 4}(:,1),20), prctile(RTic{n, 4}(:,1),40), prctile(RTic{n, 4}(:,1),60), prctile(RTic{n, 4}(:,1),80), prctile(RTic{n, 4}(:,1),100)];
end

corric = cell(NumP, 4);
incorric = cell(NumP, 4);

for n=1:NumP 
    for m=1:4
        for x=1:size(RTic{n, m},1)
            corric{n, m} = [corric{n, m} ; RTic{n, m}(x,1)<=bin_ic{n, m}(:,1) && RTic{n, m}(x,2)==1, ...
               RTic{n, m}(x,1)>bin_ic{n, m}(:,1) & RTic{n, m}(x,1)<=bin_ic{n, m}(:,2) & RTic{n, m}(x,2)==1, ...
               RTic{n, m}(x,1)>bin_ic{n, m}(:,2) & RTic{n, m}(x,1)<=bin_ic{n, m}(:,3) & RTic{n, m}(x,2)==1, ...
               RTic{n, m}(x,1)>bin_ic{n, m}(:,3) & RTic{n, m}(x,1)<=bin_ic{n, m}(:,4) & RTic{n, m}(x,2)==1, ...
               RTic{n, m}(x,1)>bin_ic{n, m}(:,4) & RTic{n, m}(x,1)<=bin_ic{n, m}(:,5) & RTic{n, m}(x,2)==1];
            incorric{n, m} = [incorric{n, m} ; RTic{n, m}(x,1)<=bin_ic{n, m}(:,1) && RTic{n, m}(x,2)~=1, ...
               RTic{n, m}(x,1)>bin_ic{n, m}(:,1) & RTic{n, m}(x,1)<=bin_ic{n, m}(:,2) & RTic{n, m}(x,2)~=1, ...
               RTic{n, m}(x,1)>bin_ic{n, m}(:,2) & RTic{n, m}(x,1)<=bin_ic{n, m}(:,3) & RTic{n, m}(x,2)~=1, ...
               RTic{n, m}(x,1)>bin_ic{n, m}(:,3) & RTic{n, m}(x,1)<=bin_ic{n, m}(:,4) & RTic{n, m}(x,2)~=1, ...
               RTic{n, m}(x,1)>bin_ic{n, m}(:,4) & RTic{n, m}(x,1)<=bin_ic{n, m}(:,5) & RTic{n, m}(x,2)~=1];  
        end
    end
end

bin_ii = cell(NumP, 4);
for n=1:NumP
    bin_ii{n, 1}=[prctile(RTii{n, 1}(:,1),20), prctile(RTii{n, 1}(:,1),40), prctile(RTii{n, 1}(:,1),60), prctile(RTii{n, 1}(:,1),80), prctile(RTii{n, 1}(:,1),100)];
    bin_ii{n, 2}=[prctile(RTii{n, 2}(:,1),20), prctile(RTii{n, 2}(:,1),40), prctile(RTii{n, 2}(:,1),60), prctile(RTii{n, 2}(:,1),80), prctile(RTii{n, 2}(:,1),100)];
    bin_ii{n, 3}=[prctile(RTii{n, 3}(:,1),20), prctile(RTii{n, 3}(:,1),40), prctile(RTii{n, 3}(:,1),60), prctile(RTii{n, 3}(:,1),80), prctile(RTii{n, 3}(:,1),100)];
    bin_ii{n, 4}=[prctile(RTii{n, 4}(:,1),20), prctile(RTii{n, 4}(:,1),40), prctile(RTii{n, 4}(:,1),60), prctile(RTii{n, 4}(:,1),80), prctile(RTii{n, 4}(:,1),100)];
end

corrii = cell(NumP, 4);
incorrii = cell(NumP, 4);

for n=1:NumP 
    for m=1:4
        for x=1:size(RTii{n, m},1)
            corrii{n, m} = [corrii{n, m} ; RTii{n, m}(x,1)<=bin_ii{n, m}(:,1) && RTii{n, m}(x,2)==1, ...
               RTii{n, m}(x,1)>bin_ii{n, m}(:,1) & RTii{n, m}(x,1)<=bin_ii{n, m}(:,2) & RTii{n, m}(x,2)==1, ...
               RTii{n, m}(x,1)>bin_ii{n, m}(:,2) & RTii{n, m}(x,1)<=bin_ii{n, m}(:,3) & RTii{n, m}(x,2)==1, ...
               RTii{n, m}(x,1)>bin_ii{n, m}(:,3) & RTii{n, m}(x,1)<=bin_ii{n, m}(:,4) & RTii{n, m}(x,2)==1, ...
               RTii{n, m}(x,1)>bin_ii{n, m}(:,4) & RTii{n, m}(x,1)<=bin_ii{n, m}(:,5) & RTii{n, m}(x,2)==1];
            incorrii{n, m} = [incorrii{n, m} ; RTii{n, m}(x,1)<=bin_ii{n, m}(:,1) && RTii{n, m}(x,2)~=1, ...
               RTii{n, m}(x,1)>bin_ii{n, m}(:,1) & RTii{n, m}(x,1)<=bin_ii{n, m}(:,2) & RTii{n, m}(x,2)~=1, ...
               RTii{n, m}(x,1)>bin_ii{n, m}(:,2) & RTii{n, m}(x,1)<=bin_ii{n, m}(:,3) & RTii{n, m}(x,2)~=1, ...
               RTii{n, m}(x,1)>bin_ii{n, m}(:,3) & RTii{n, m}(x,1)<=bin_ii{n, m}(:,4) & RTii{n, m}(x,2)~=1, ...
               RTii{n, m}(x,1)>bin_ii{n, m}(:,4) & RTii{n, m}(x,1)<=bin_ii{n, m}(:,5) & RTii{n, m}(x,2)~=1];
            
        end
    end
end

p_ic=[];
p_ic(1, :) = sum(cell2mat(corric(:, 1)))./(sum(cell2mat(corric(:, 1)))+sum(cell2mat(incorric(:, 1))));
p_ic(2, :) = sum(cell2mat(corric(:, 2)))./(sum(cell2mat(corric(:, 2)))+sum(cell2mat(incorric(:, 2))));
p_ic(3, :) = sum(cell2mat(corric(:, 3)))./(sum(cell2mat(corric(:, 3)))+sum(cell2mat(incorric(:, 3))));
p_ic(4, :) = sum(cell2mat(corric(:, 4)))./(sum(cell2mat(corric(:, 4)))+sum(cell2mat(incorric(:, 4))));

p_ii=[];
p_ii(1, :) = sum(cell2mat(corrii(:, 1)))./(sum(cell2mat(corrii(:, 1)))+sum(cell2mat(incorrii(:, 1))));
p_ii(2, :) = sum(cell2mat(corrii(:, 2)))./(sum(cell2mat(corrii(:, 2)))+sum(cell2mat(incorrii(:, 2))));
p_ii(3, :) = sum(cell2mat(corrii(:, 3)))./(sum(cell2mat(corrii(:, 3)))+sum(cell2mat(incorrii(:, 3))));
p_ii(4, :) = sum(cell2mat(corrii(:, 4)))./(sum(cell2mat(corrii(:, 4)))+sum(cell2mat(incorrii(:, 4))));

% Simulated data
% Load the estimated parameters from the best performing model

% Model Prediction
% cC
Number_of_simulation = 50000;

load('spatialStroop_CSE_SOA_SIMPLEX_04.mat')

A_01_pc=fitted_Param(1);
A_01_pi=fitted_Param(1);
kc_pc=fitted_Param(3);
ki_pc=fitted_Param(4);
da1_pc=fitted_Param(5);
da1_pi=fitted_Param(6);
dc=fitted_Param(8);
a_pc=fitted_Param(9);
Ter1=fitted_Param(10);
Ter2=fitted_Param(11);
Ter3=fitted_Param(12);
Ter4=fitted_Param(13);
std=fitted_Param(14);
A_02_pc=fitted_Param(2);
A_03_pc=fitted_Param(2);
A_04_pc=fitted_Param(2);
A_02_pi=fitted_Param(2);
A_03_pi=fitted_Param(2);
A_04_pi=fitted_Param(2);
kc_pi=fitted_Param(3);
ki_pi=fitted_Param(4);
da2_pc=fitted_Param(5);
da3_pc=fitted_Param(5);
da4_pc=fitted_Param(5);
da2_pi=fitted_Param(6);
da3_pi=fitted_Param(7);
da4_pi=fitted_Param(7);
a_pi=fitted_Param(9);
c=0;

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

bin_cc_pred = {};
for m=1:4
    bin_cc_pred{1, m}=[prctile(cC_predicted{1, m}(:,1),20), prctile(cC_predicted{1, m}(:,1),40), prctile(cC_predicted{1, m}(:,1),60), ...
        prctile(cC_predicted{1, m}(:,1),80), prctile(cC_predicted{1, m}(:,1),100)];
end

corrcc_pred = {};
incorrcc_pred = {};

for m=1:4
    corrcc_pred{1, m} = [cC_predicted{1, m}(:,1)<=bin_cc_pred{1, m}(:,1) & cC_predicted{1, m}(:,2)==1, ...
        cC_predicted{1, m}(:,1)>bin_cc_pred{1, m}(:,1) & cC_predicted{1, m}(:,1)<=bin_cc_pred{1, m}(:,2) & cC_predicted{1, m}(:,2)==1, ...
        cC_predicted{1, m}(:,1)>bin_cc_pred{1, m}(:,2) & cC_predicted{1, m}(:,1)<=bin_cc_pred{1, m}(:,3) & cC_predicted{1, m}(:,2)==1, ...
        cC_predicted{1, m}(:,1)>bin_cc_pred{1, m}(:,3) & cC_predicted{1, m}(:,1)<=bin_cc_pred{1, m}(:,4) & cC_predicted{1, m}(:,2)==1, ...
        cC_predicted{1, m}(:,1)>bin_cc_pred{1, m}(:,4) & cC_predicted{1, m}(:,1)<=bin_cc_pred{1, m}(:,5) & cC_predicted{1, m}(:,2)==1];
    incorrcc_pred{1, m} = [cC_predicted{1, m}(:,1)<=bin_cc_pred{1, m}(:,1) & cC_predicted{1, m}(:,2)==0, ...
        cC_predicted{1, m}(:,1)>bin_cc_pred{1, m}(:,1) & cC_predicted{1, m}(:,1)<=bin_cc_pred{1, m}(:,2) & cC_predicted{1, m}(:,2)==0, ...
        cC_predicted{1, m}(:,1)>bin_cc_pred{1, m}(:,2) & cC_predicted{1, m}(:,1)<=bin_cc_pred{1, m}(:,3) & cC_predicted{1, m}(:,2)==0, ...
        cC_predicted{1, m}(:,1)>bin_cc_pred{1, m}(:,3) & cC_predicted{1, m}(:,1)<=bin_cc_pred{1, m}(:,4) & cC_predicted{1, m}(:,2)==0, ...
        cC_predicted{1, m}(:,1)>bin_cc_pred{1, m}(:,4) & cC_predicted{1, m}(:,1)<=bin_cc_pred{1, m}(:,5) & cC_predicted{1, m}(:,2)==0];   
end


bin_ci_pred = {};
for m=1:4
    bin_ci_pred{1, m}=[prctile(cI_predicted{1, m}(:,1),20), prctile(cI_predicted{1, m}(:,1),40), prctile(cI_predicted{1, m}(:,1),60), ...
        prctile(cI_predicted{1, m}(:,1),80), prctile(cI_predicted{1, m}(:,1),100)];
end

corrci_pred = {};
incorrci_pred = {};

for m=1:4
    corrci_pred{1, m} = [cI_predicted{1, m}(:,1)<=bin_ci_pred{1, m}(:,1) & cI_predicted{1, m}(:,2)==1, ...
        cI_predicted{1, m}(:,1)>bin_ci_pred{1, m}(:,1) & cI_predicted{1, m}(:,1)<=bin_ci_pred{1, m}(:,2) & cI_predicted{1, m}(:,2)==1, ...
        cI_predicted{1, m}(:,1)>bin_ci_pred{1, m}(:,2) & cI_predicted{1, m}(:,1)<=bin_ci_pred{1, m}(:,3) & cI_predicted{1, m}(:,2)==1, ...
        cI_predicted{1, m}(:,1)>bin_ci_pred{1, m}(:,3) & cI_predicted{1, m}(:,1)<=bin_ci_pred{1, m}(:,4) & cI_predicted{1, m}(:,2)==1, ...
        cI_predicted{1, m}(:,1)>bin_ci_pred{1, m}(:,4) & cI_predicted{1, m}(:,1)<=bin_ci_pred{1, m}(:,5) & cI_predicted{1, m}(:,2)==1];
    incorrci_pred{1, m} = [cI_predicted{1, m}(:,1)<=bin_ci_pred{1, m}(:,1) & cI_predicted{1, m}(:,2)==0, ...
        cI_predicted{1, m}(:,1)>bin_ci_pred{1, m}(:,1) & cI_predicted{1, m}(:,1)<=bin_ci_pred{1, m}(:,2) & cI_predicted{1, m}(:,2)==0, ...
        cI_predicted{1, m}(:,1)>bin_ci_pred{1, m}(:,2) & cI_predicted{1, m}(:,1)<=bin_ci_pred{1, m}(:,3) & cI_predicted{1, m}(:,2)==0, ...
        cI_predicted{1, m}(:,1)>bin_ci_pred{1, m}(:,3) & cI_predicted{1, m}(:,1)<=bin_ci_pred{1, m}(:,4) & cI_predicted{1, m}(:,2)==0, ...
        cI_predicted{1, m}(:,1)>bin_ci_pred{1, m}(:,4) & cI_predicted{1, m}(:,1)<=bin_ci_pred{1, m}(:,5) & cI_predicted{1, m}(:,2)==0];   
end

p_cc_pred=[];
for m=1:4
    p_cc_pred(m, :) = sum(cell2mat(corrcc_pred(:, m)))./(sum(cell2mat(corrcc_pred(:, m)))+sum(cell2mat(incorrcc_pred(:, m))));
end

p_ci_pred=[];
for m=1:4
    p_ci_pred(m, :) = sum(cell2mat(corrci_pred(:, m)))./(sum(cell2mat(corrci_pred(:, m)))+sum(cell2mat(incorrci_pred(:, m))));
end

bin_ic_pred = {};
for m=1:4
    bin_ic_pred{1, m}=[prctile(iC_predicted{1, m}(:,1),20), prctile(iC_predicted{1, m}(:,1),40), prctile(iC_predicted{1, m}(:,1),60), ...
        prctile(iC_predicted{1, m}(:,1),80), prctile(iC_predicted{1, m}(:,1),100)];
end

corric_pred = {};
incorric_pred = {};

for m=1:4
    corric_pred{1, m} = [iC_predicted{1, m}(:,1)<=bin_ic_pred{1, m}(:,1) & iC_predicted{1, m}(:,2)==1, ...
        iC_predicted{1, m}(:,1)>bin_ic_pred{1, m}(:,1) & iC_predicted{1, m}(:,1)<=bin_ic_pred{1, m}(:,2) & iC_predicted{1, m}(:,2)==1, ...
        iC_predicted{1, m}(:,1)>bin_ic_pred{1, m}(:,2) & iC_predicted{1, m}(:,1)<=bin_ic_pred{1, m}(:,3) & iC_predicted{1, m}(:,2)==1, ...
        iC_predicted{1, m}(:,1)>bin_ic_pred{1, m}(:,3) & iC_predicted{1, m}(:,1)<=bin_ic_pred{1, m}(:,4) & iC_predicted{1, m}(:,2)==1, ...
        iC_predicted{1, m}(:,1)>bin_ic_pred{1, m}(:,4) & iC_predicted{1, m}(:,1)<=bin_ic_pred{1, m}(:,5) & iC_predicted{1, m}(:,2)==1];
    incorric_pred{1, m} = [iC_predicted{1, m}(:,1)<=bin_ic_pred{1, m}(:,1) & iC_predicted{1, m}(:,2)==0, ...
        iC_predicted{1, m}(:,1)>bin_ic_pred{1, m}(:,1) & iC_predicted{1, m}(:,1)<=bin_ic_pred{1, m}(:,2) & iC_predicted{1, m}(:,2)==0, ...
        iC_predicted{1, m}(:,1)>bin_ic_pred{1, m}(:,2) & iC_predicted{1, m}(:,1)<=bin_ic_pred{1, m}(:,3) & iC_predicted{1, m}(:,2)==0, ...
        iC_predicted{1, m}(:,1)>bin_ic_pred{1, m}(:,3) & iC_predicted{1, m}(:,1)<=bin_ic_pred{1, m}(:,4) & iC_predicted{1, m}(:,2)==0, ...
        iC_predicted{1, m}(:,1)>bin_ic_pred{1, m}(:,4) & iC_predicted{1, m}(:,1)<=bin_ic_pred{1, m}(:,5) & iC_predicted{1, m}(:,2)==0];   
end


bin_ii_pred = {};
for m=1:4
    bin_ii_pred{1, m}=[prctile(iI_predicted{1, m}(:,1),20), prctile(iI_predicted{1, m}(:,1),40), prctile(iI_predicted{1, m}(:,1),60), ...
        prctile(iI_predicted{1, m}(:,1),80), prctile(iI_predicted{1, m}(:,1),100)];
end

corrii_pred = {};
incorrii_pred = {};

for m=1:4
    corrii_pred{1, m} = [iI_predicted{1, m}(:,1)<=bin_ii_pred{1, m}(:,1) & iI_predicted{1, m}(:,2)==1, ...
        iI_predicted{1, m}(:,1)>bin_ii_pred{1, m}(:,1) & iI_predicted{1, m}(:,1)<=bin_ii_pred{1, m}(:,2) & iI_predicted{1, m}(:,2)==1, ...
        iI_predicted{1, m}(:,1)>bin_ii_pred{1, m}(:,2) & iI_predicted{1, m}(:,1)<=bin_ii_pred{1, m}(:,3) & iI_predicted{1, m}(:,2)==1, ...
        iI_predicted{1, m}(:,1)>bin_ii_pred{1, m}(:,3) & iI_predicted{1, m}(:,1)<=bin_ii_pred{1, m}(:,4) & iI_predicted{1, m}(:,2)==1, ...
        iI_predicted{1, m}(:,1)>bin_ii_pred{1, m}(:,4) & iI_predicted{1, m}(:,1)<=bin_ii_pred{1, m}(:,5) & iI_predicted{1, m}(:,2)==1];
    incorrii_pred{1, m} = [iI_predicted{1, m}(:,1)<=bin_ii_pred{1, m}(:,1) & iI_predicted{1, m}(:,2)==0, ...
        iI_predicted{1, m}(:,1)>bin_ii_pred{1, m}(:,1) & iI_predicted{1, m}(:,1)<=bin_ii_pred{1, m}(:,2) & iI_predicted{1, m}(:,2)==0, ...
        iI_predicted{1, m}(:,1)>bin_ii_pred{1, m}(:,2) & iI_predicted{1, m}(:,1)<=bin_ii_pred{1, m}(:,3) & iI_predicted{1, m}(:,2)==0, ...
        iI_predicted{1, m}(:,1)>bin_ii_pred{1, m}(:,3) & iI_predicted{1, m}(:,1)<=bin_ii_pred{1, m}(:,4) & iI_predicted{1, m}(:,2)==0, ...
        iI_predicted{1, m}(:,1)>bin_ii_pred{1, m}(:,4) & iI_predicted{1, m}(:,1)<=bin_ii_pred{1, m}(:,5) & iI_predicted{1, m}(:,2)==0];   
end

p_ic_pred=[];
for m=1:4
    p_ic_pred(m, :) = sum(cell2mat(corric_pred(:, m)))./(sum(cell2mat(corric_pred(:, m)))+sum(cell2mat(incorric_pred(:, m))));
end

p_ii_pred=[];
for m=1:4
    p_ii_pred(m, :) = sum(cell2mat(corrii_pred(:, m)))./(sum(cell2mat(corrii_pred(:, m)))+sum(cell2mat(incorrii_pred(:, m))));
end

% Plot CAF of each SOA
figure(1)
figure('Position', [300 300 1000 700])
t= tiledlayout(2,2);
title(t, 'Previous Congruent');
nexttile
RT20 = [p_cc(1, 1), p_cc_pred(1, 1), p_ci(1, 1), p_ci_pred(1, 1)];
RT40 = [p_cc(1, 2), p_cc_pred(1, 2), p_ci(1, 2), p_ci_pred(1, 2)];
RT60 = [p_cc(1, 3), p_cc_pred(1, 3), p_ci(1, 3), p_ci_pred(1, 3)];
RT80 = [p_cc(1, 4), p_cc_pred(1, 4), p_ci(1, 4), p_ci_pred(1, 4)];
RT100 = [p_cc(1, 5), p_cc_pred(1, 5), p_ci(1, 5), p_ci_pred(1, 5)];
plot([0.8, 1.1, 1.4, 1.7, 2.0], [RT20(1), RT40(1), RT60(1), RT80(1), RT100(1)], 'LineStyle', ...
    ':', 'Color', [0 0.5 0], 'LineWidth', 1.5);
hold on;
plot([0.8, 1.1, 1.4, 1.7, 2.0], [RT20(2), RT40(2), RT60(2), RT80(2), RT100(2)], 'LineStyle', ...
    '-', 'Color', [0 0.5 0], 'Marker', 'o', 'LineWidth', 1.5);
hold on 
plot([0.8, 1.1, 1.4, 1.7, 2.0], [RT20(3), RT40(3), RT60(3), RT80(3), RT100(3)], 'LineStyle', ...
    ':', 'Color', [0.5 0 0.5], 'LineWidth', 1.5);
hold on 
plot([0.8, 1.1, 1.4, 1.7, 2.0], [RT20(4), RT40(4), RT60(4), RT80(4), RT100(4)], 'LineStyle', ...
    '-', 'Color', [0.5 0 0.5], 'Marker', 'o', 'LineWidth', 1.5);
grid on;
xlim([0.5, 2.3]);
ylim([0, 1]);
xlabel('RT bins');
ylabel('Accuracy');
ax = gca;
ax.XTick = [0.8, 1.1, 1.4, 1.7, 2.0];
ax.XTickLabels = {'RT20','RT40', 'RT60', 'RT80', 'RT100'};
grid on;
legend('Empirical cC', 'Predicted cC', 'Empirical cI', 'Predicted cI', 'Location', 'best')
title('SOA = +100 ms');
hold off

nexttile
RT20 = [p_cc(2, 1), p_cc_pred(2, 1), p_ci(2, 1), p_ci_pred(2, 1)];
RT40 = [p_cc(2, 2), p_cc_pred(2, 2), p_ci(2, 2), p_ci_pred(2, 2)];
RT60 = [p_cc(2, 3), p_cc_pred(2, 3), p_ci(2, 3), p_ci_pred(2, 3)];
RT80 = [p_cc(2, 4), p_cc_pred(2, 4), p_ci(2, 4), p_ci_pred(2, 4)];
RT100 = [p_cc(2, 5), p_cc_pred(2, 5), p_ci(2, 5), p_ci_pred(2, 5)];
plot([0.8, 1.1, 1.4, 1.7, 2.0], [RT20(1), RT40(1), RT60(1), RT80(1), RT100(1)], 'LineStyle', ...
    ':', 'Color', [0 0.5 0], 'LineWidth', 1.5);
hold on;
plot([0.8, 1.1, 1.4, 1.7, 2.0], [RT20(2), RT40(2), RT60(2), RT80(2), RT100(2)], 'LineStyle', ...
    '-', 'Color', [0 0.5 0], 'Marker', 'o', 'LineWidth', 1.5);
hold on 
plot([0.8, 1.1, 1.4, 1.7, 2.0], [RT20(3), RT40(3), RT60(3), RT80(3), RT100(3)], 'LineStyle', ...
    ':', 'Color', [0.5 0 0.5], 'LineWidth', 1.5);
hold on 
plot([0.8, 1.1, 1.4, 1.7, 2.0], [RT20(4), RT40(4), RT60(4), RT80(4), RT100(4)], 'LineStyle', ...
    '-', 'Color', [0.5 0 0.5], 'Marker', 'o', 'LineWidth', 1.5);
grid on;
xlim([0.5, 2.3]);
ylim([0, 1]);
xlabel('RT bins');
ylabel('Accuracy');
ax = gca;
ax.XTick = [0.8, 1.1, 1.4, 1.7, 2.0];
ax.XTickLabels = {'RT20','RT40', 'RT60', 'RT80', 'RT100'};
grid on;
legend('Empirical cC', 'Predicted cC', 'Empirical cI', 'Predicted cI', 'Location', 'best')
title('SOA = 0 ms');
hold off

nexttile
RT20 = [p_cc(3, 1), p_cc_pred(3, 1), p_ci(3, 1), p_ci_pred(3, 1)];
RT40 = [p_cc(3, 2), p_cc_pred(3, 2), p_ci(3, 2), p_ci_pred(3, 2)];
RT60 = [p_cc(3, 3), p_cc_pred(3, 3), p_ci(3, 3), p_ci_pred(3, 3)];
RT80 = [p_cc(3, 4), p_cc_pred(3, 4), p_ci(3, 4), p_ci_pred(3, 4)];
RT100 = [p_cc(3, 5), p_cc_pred(3, 5), p_ci(3, 5), p_ci_pred(3, 5)];
plot([0.8, 1.1, 1.4, 1.7, 2.0], [RT20(1), RT40(1), RT60(1), RT80(1), RT100(1)], 'LineStyle', ...
    ':', 'Color', [0 0.5 0], 'LineWidth', 1.5);
hold on;
plot([0.8, 1.1, 1.4, 1.7, 2.0], [RT20(2), RT40(2), RT60(2), RT80(2), RT100(2)], 'LineStyle', ...
    '-', 'Color', [0 0.5 0], 'Marker', 'o', 'LineWidth', 1.5);
hold on 
plot([0.8, 1.1, 1.4, 1.7, 2.0], [RT20(3), RT40(3), RT60(3), RT80(3), RT100(3)], 'LineStyle', ...
    ':', 'Color', [0.5 0 0.5], 'LineWidth', 1.5);
hold on 
plot([0.8, 1.1, 1.4, 1.7, 2.0], [RT20(4), RT40(4), RT60(4), RT80(4), RT100(4)], 'LineStyle', ...
    '-', 'Color', [0.5 0 0.5], 'Marker', 'o', 'LineWidth', 1.5);
grid on;
xlim([0.5, 2.3]);
ylim([0, 1]);
xlabel('RT bins');
ylabel('Accuracy');
ax = gca;
ax.XTick = [0.8, 1.1, 1.4, 1.7, 2.0];
ax.XTickLabels = {'RT20','RT40', 'RT60', 'RT80', 'RT100'};
grid on;
legend('Empirical cC', 'Predicted cC', 'Empirical cI', 'Predicted cI', 'Location', 'best')
title('SOA = -50 ms');
hold off

nexttile
RT20 = [p_cc(4, 1), p_cc_pred(4, 1), p_ci(4, 1), p_ci_pred(4, 1)];
RT40 = [p_cc(4, 2), p_cc_pred(4, 2), p_ci(4, 2), p_ci_pred(4, 2)];
RT60 = [p_cc(4, 3), p_cc_pred(4, 3), p_ci(4, 3), p_ci_pred(4, 3)];
RT80 = [p_cc(4, 4), p_cc_pred(4, 4), p_ci(4, 4), p_ci_pred(4, 4)];
RT100 = [p_cc(4, 5), p_cc_pred(4, 5), p_ci(4, 5), p_ci_pred(4, 5)];
plot([0.8, 1.1, 1.4, 1.7, 2.0], [RT20(1), RT40(1), RT60(1), RT80(1), RT100(1)], 'LineStyle', ...
    ':', 'Color', [0 0.5 0], 'LineWidth', 1.5);
hold on;
plot([0.8, 1.1, 1.4, 1.7, 2.0], [RT20(2), RT40(2), RT60(2), RT80(2), RT100(2)], 'LineStyle', ...
    '-', 'Color', [0 0.5 0], 'Marker', 'o', 'LineWidth', 1.5);
hold on 
plot([0.8, 1.1, 1.4, 1.7, 2.0], [RT20(3), RT40(3), RT60(3), RT80(3), RT100(3)], 'LineStyle', ...
    ':', 'Color', [0.5 0 0.5], 'LineWidth', 1.5);
hold on 
plot([0.8, 1.1, 1.4, 1.7, 2.0], [RT20(4), RT40(4), RT60(4), RT80(4), RT100(4)], 'LineStyle', ...
    '-', 'Color', [0.5 0 0.5], 'Marker', 'o', 'LineWidth', 1.5);
grid on;
xlim([0.5, 2.3]);
ylim([0, 1]);
xlabel('RT bins');
ylabel('Accuracy');
ax = gca;
ax.XTick = [0.8, 1.1, 1.4, 1.7, 2.0];
ax.XTickLabels = {'RT20','RT40', 'RT60', 'RT80', 'RT100'};
grid on;
legend('Empirical cC', 'Predicted cC', 'Empirical cI', 'Predicted cI', 'Location', 'best')
title('SOA = -100 ms');
hold off

figure(2)
figure('Position', [300 300 1000 700])
h= tiledlayout(2,2);
title(h, 'Previous Incongruent');
nexttile
RT20 = [p_ic(1, 1), p_ic_pred(1, 1), p_ii(1, 1), p_ii_pred(1, 1)];
RT40 = [p_ic(1, 2), p_ic_pred(1, 2), p_ii(1, 2), p_ii_pred(1, 2)];
RT60 = [p_ic(1, 3), p_ic_pred(1, 3), p_ii(1, 3), p_ii_pred(1, 3)];
RT80 = [p_ic(1, 4), p_ic_pred(1, 4), p_ii(1, 4), p_ii_pred(1, 4)];
RT100 = [p_ic(1, 5), p_ic_pred(1, 5), p_ii(1, 5), p_ii_pred(1, 5)];
plot([0.8, 1.1, 1.4, 1.7, 2.0], [RT20(1), RT40(1), RT60(1), RT80(1), RT100(1)], 'LineStyle', ...
    ':', 'Color', [0 0.5 0], 'LineWidth', 1.5);
hold on;
plot([0.8, 1.1, 1.4, 1.7, 2.0], [RT20(2), RT40(2), RT60(2), RT80(2), RT100(2)], 'LineStyle', ...
    '-', 'Color', [0 0.5 0], 'Marker', 'o', 'LineWidth', 1.5);
hold on 
plot([0.8, 1.1, 1.4, 1.7, 2.0], [RT20(3), RT40(3), RT60(3), RT80(3), RT100(3)], 'LineStyle', ...
    ':', 'Color', [0.5 0 0.5], 'LineWidth', 1.5);
hold on 
plot([0.8, 1.1, 1.4, 1.7, 2.0], [RT20(4), RT40(4), RT60(4), RT80(4), RT100(4)], 'LineStyle', ...
    '-', 'Color', [0.5 0 0.5], 'Marker', 'o', 'LineWidth', 1.5);
grid on;
xlim([0.5, 2.3]);
ylim([0, 1]);
xlabel('RT bins');
ylabel('Aicuracy');
ax = gca;
ax.XTick = [0.8, 1.1, 1.4, 1.7, 2.0];
ax.XTickLabels = {'RT20','RT40', 'RT60', 'RT80', 'RT100'};
grid on;
legend('Empirical iC', 'Predicted iC', 'Empirical iI', 'Predicted iI', 'Location', 'best')
title('SOA = +100 ms');
hold off

nexttile
RT20 = [p_ic(2, 1), p_ic_pred(2, 1), p_ii(2, 1), p_ii_pred(2, 1)];
RT40 = [p_ic(2, 2), p_ic_pred(2, 2), p_ii(2, 2), p_ii_pred(2, 2)];
RT60 = [p_ic(2, 3), p_ic_pred(2, 3), p_ii(2, 3), p_ii_pred(2, 3)];
RT80 = [p_ic(2, 4), p_ic_pred(2, 4), p_ii(2, 4), p_ii_pred(2, 4)];
RT100 = [p_ic(2, 5), p_ic_pred(2, 5), p_ii(2, 5), p_ii_pred(2, 5)];
plot([0.8, 1.1, 1.4, 1.7, 2.0], [RT20(1), RT40(1), RT60(1), RT80(1), RT100(1)], 'LineStyle', ...
    ':', 'Color', [0 0.5 0], 'LineWidth', 1.5);
hold on;
plot([0.8, 1.1, 1.4, 1.7, 2.0], [RT20(2), RT40(2), RT60(2), RT80(2), RT100(2)], 'LineStyle', ...
    '-', 'Color', [0 0.5 0], 'Marker', 'o', 'LineWidth', 1.5);
hold on 
plot([0.8, 1.1, 1.4, 1.7, 2.0], [RT20(3), RT40(3), RT60(3), RT80(3), RT100(3)], 'LineStyle', ...
    ':', 'Color', [0.5 0 0.5], 'LineWidth', 1.5);
hold on 
plot([0.8, 1.1, 1.4, 1.7, 2.0], [RT20(4), RT40(4), RT60(4), RT80(4), RT100(4)], 'LineStyle', ...
    '-', 'Color', [0.5 0 0.5], 'Marker', 'o', 'LineWidth', 1.5);
grid on;
xlim([0.5, 2.3]);
ylim([0, 1]);
xlabel('RT bins');
ylabel('Aicuracy');
ax = gca;
ax.XTick = [0.8, 1.1, 1.4, 1.7, 2.0];
ax.XTickLabels = {'RT20','RT40', 'RT60', 'RT80', 'RT100'};
grid on;
legend('Empirical iC', 'Predicted iC', 'Empirical iI', 'Predicted iI', 'Location', 'best')
title('SOA = 0 ms');
hold off

nexttile
RT20 = [p_ic(3, 1), p_ic_pred(3, 1), p_ii(3, 1), p_ii_pred(3, 1)];
RT40 = [p_ic(3, 2), p_ic_pred(3, 2), p_ii(3, 2), p_ii_pred(3, 2)];
RT60 = [p_ic(3, 3), p_ic_pred(3, 3), p_ii(3, 3), p_ii_pred(3, 3)];
RT80 = [p_ic(3, 4), p_ic_pred(3, 4), p_ii(3, 4), p_ii_pred(3, 4)];
RT100 = [p_ic(3, 5), p_ic_pred(3, 5), p_ii(3, 5), p_ii_pred(3, 5)];
plot([0.8, 1.1, 1.4, 1.7, 2.0], [RT20(1), RT40(1), RT60(1), RT80(1), RT100(1)], 'LineStyle', ...
    ':', 'Color', [0 0.5 0], 'LineWidth', 1.5);
hold on;
plot([0.8, 1.1, 1.4, 1.7, 2.0], [RT20(2), RT40(2), RT60(2), RT80(2), RT100(2)], 'LineStyle', ...
    '-', 'Color', [0 0.5 0], 'Marker', 'o', 'LineWidth', 1.5);
hold on 
plot([0.8, 1.1, 1.4, 1.7, 2.0], [RT20(3), RT40(3), RT60(3), RT80(3), RT100(3)], 'LineStyle', ...
    ':', 'Color', [0.5 0 0.5], 'LineWidth', 1.5);
hold on 
plot([0.8, 1.1, 1.4, 1.7, 2.0], [RT20(4), RT40(4), RT60(4), RT80(4), RT100(4)], 'LineStyle', ...
    '-', 'Color', [0.5 0 0.5], 'Marker', 'o', 'LineWidth', 1.5);
grid on;
xlim([0.5, 2.3]);
ylim([0, 1]);
xlabel('RT bins');
ylabel('Aicuracy');
ax = gca;
ax.XTick = [0.8, 1.1, 1.4, 1.7, 2.0];
ax.XTickLabels = {'RT20','RT40', 'RT60', 'RT80', 'RT100'};
grid on;
legend('Empirical iC', 'Predicted iC', 'Empirical iI', 'Predicted iI', 'Location', 'best')
title('SOA = -50 ms');
hold off

nexttile
RT20 = [p_ic(4, 1), p_ic_pred(4, 1), p_ii(4, 1), p_ii_pred(4, 1)];
RT40 = [p_ic(4, 2), p_ic_pred(4, 2), p_ii(4, 2), p_ii_pred(4, 2)];
RT60 = [p_ic(4, 3), p_ic_pred(4, 3), p_ii(4, 3), p_ii_pred(4, 3)];
RT80 = [p_ic(4, 4), p_ic_pred(4, 4), p_ii(4, 4), p_ii_pred(4, 4)];
RT100 = [p_ic(4, 5), p_ic_pred(4, 5), p_ii(4, 5), p_ii_pred(4, 5)];
plot([0.8, 1.1, 1.4, 1.7, 2.0], [RT20(1), RT40(1), RT60(1), RT80(1), RT100(1)], 'LineStyle', ...
    ':', 'Color', [0 0.5 0], 'LineWidth', 1.5);
hold on;
plot([0.8, 1.1, 1.4, 1.7, 2.0], [RT20(2), RT40(2), RT60(2), RT80(2), RT100(2)], 'LineStyle', ...
    '-', 'Color', [0 0.5 0], 'Marker', 'o', 'LineWidth', 1.5);
hold on 
plot([0.8, 1.1, 1.4, 1.7, 2.0], [RT20(3), RT40(3), RT60(3), RT80(3), RT100(3)], 'LineStyle', ...
    ':', 'Color', [0.5 0 0.5], 'LineWidth', 1.5);
hold on 
plot([0.8, 1.1, 1.4, 1.7, 2.0], [RT20(4), RT40(4), RT60(4), RT80(4), RT100(4)], 'LineStyle', ...
    '-', 'Color', [0.5 0 0.5], 'Marker', 'o', 'LineWidth', 1.5);
grid on;
xlim([0.5, 2.3]);
ylim([0, 1]);
xlabel('RT bins');
ylabel('Aicuracy');
ax = gca;
ax.XTick = [0.8, 1.1, 1.4, 1.7, 2.0];
ax.XTickLabels = {'RT20','RT40', 'RT60', 'RT80', 'RT100'};
grid on;
legend('Empirical iC', 'Predicted iC', 'Empirical iI', 'Predicted iI', 'Location', 'best')
title('SOA = -100 ms');
hold off

exportgraphics(t,'SpatialStroop_CSE_SOA_CAF.pdf')
exportgraphics(h,'SpatialStroop_CSE_SOA_CAF.pdf','Append',true)
