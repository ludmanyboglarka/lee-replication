%% Model fitting of the best rank from the model set
% CDF (CE x SOA)
clear
% Load empirical spatialStroop/flanker data
warning off
spatialStroop=readtable('spatialStroopSOA.xlsx');
spatialStroop=table2array(spatialStroop);

NumP= 40;

Ppt = cell(NumP, 1);
for k=1:NumP
    Ppt{k} = spatialStroop(spatialStroop(:,1)==k,:);
end

C = cell(NumP, 1);
I = cell(NumP, 1);
for n=1:NumP
    for k=1:size(Ppt{n},1)
        if Ppt{n}(k,7)==1
            C{n} = [C{n} ; Ppt{n}(k,1) Ppt{n}(k,2) Ppt{n}(k,5) Ppt{n}(k,6), Ppt{n}(k,8)]; 
        elseif Ppt{n}(k,7)==0 
            I{n} = [I{n} ; Ppt{n}(k,1) Ppt{n}(k,2) Ppt{n}(k,5) Ppt{n}(k,6), Ppt{n}(k,8)];
        end
    end
end

C = cell2mat(C);
de3=C(:,3)<200;
C(de3,:) = [];
de4=C(:,3)>1500;
C(de4,:) = [];

I = cell2mat(I);
de7=I(:,3)<200;
I(de7,:) = [];
de8=I(:,3)>1500;
I(de8,:) = [];

C1=cell(NumP, 1);
C2=cell(NumP, 1);
for n=1:NumP
    k=C(:, 1)==n;
    C1{n} = C(k, :);
    C2{n} = [C1{n}(:, 3), C1{n}(:, 5)];
end
C = C2;

I1=cell(NumP, 1);
I2=cell(NumP, 1);
for n=1:NumP
    l=I(:, 1)==n;
    I1{n} = I(l, :);
    I2{n} = [I1{n}(:, 3), I1{n}(:, 5)];
end
I=I2;

% We have different SOAs in the spatialStroop task, -100s (target shown before
% distractors), 0, 50, 100s
% Five RT quantiles for each condition

RTc = cell(n, 4);
for n=1:NumP
    RTc{n, 1} = C{n}(C{n}(:, 2)==100,1);
    RTc{n, 2} = C{n}(C{n}(:, 2)==0,1);
    RTc{n, 3} = C{n}(C{n}(:, 2)==-50,1);
    RTc{n, 4} = C{n}(C{n}(:, 2)==-100,1);
end

RTi = cell(n, 4);
for n=1:NumP
    RTi{n, 1} = I{n}(I{n}(:, 2)==100,1);
    RTi{n, 2} = I{n}(I{n}(:, 2)==0,1);
    RTi{n, 3} = I{n}(I{n}(:, 2)==-50,1); 
    RTi{n, 4} = I{n}(I{n}(:, 2)==-100,1);
end

% Simulated data
% Load the estimated parameters from the best performing model

% Model Prediction
% cC
Number_of_simulation = 50000;

load('spatialStroop_CE_SOA_SIMPLEX_02.mat')

A_01=fitted_Param(1);
kc1=fitted_Param(2);
ki1=fitted_Param(3);
da1=fitted_Param(4);
dc1=fitted_Param(5);
a1=fitted_Param(6);
Ter1=fitted_Param(7);
Ter2=fitted_Param(8);
Ter3=fitted_Param(9);
Ter4=fitted_Param(10);
std=fitted_Param(11);
A_02=A_01;
A_03=A_01;
A_04=A_01;
kc2=kc1;
kc3=kc1;
kc4=kc1;
ki2=ki1;
ki3=ki1;
ki4=ki1;
da2=da1;
da3=da1;
da4=da1;
dc2=dc1;
dc3=dc1;
dc4=dc1;
a2=a1;
a3=a1;
a4=a1;
c1=0;
c2=c1;
c3=c1;
c4=c1;

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

C_SOA_1 = cell2mat(RTc(:, 1));
C_SOA_2 = cell2mat(RTc(:, 2));
C_SOA_3 = cell2mat(RTc(:, 3));
C_SOA_4 = cell2mat(RTc(:, 4));

I_SOA_1 = cell2mat(RTi(:, 1));
I_SOA_2 = cell2mat(RTi(:, 2));
I_SOA_3 = cell2mat(RTi(:, 3));
I_SOA_4 = cell2mat(RTi(:, 4));

% Plot CDFs
figure('Position', [1, 1, 800, 400])
t=tiledlayout(2, 2);
nexttile
h1 = cdfplot(C_SOA_1);
hold on
h2 = cdfplot(I_SOA_1);
hold on 
h3 = cdfplot(C_predicted{1, 1}(:, 1));
hold on 
h4 = cdfplot(I_predicted{1, 1}(:, 1));
xlim([200 1000]);
xlabel('RT (ms)');
ylabel('CDF');
set(h1, 'LineStyle', ':', 'Color', [0 0.5 0], 'LineWidth', 1.5);
set(h2, 'LineStyle', ':', 'Color', [.5 0 .5], 'LineWidth', 1.5);
set(h3, 'LineStyle', '-', 'Color', [0 0.5 0], 'LineWidth', 1.5);
set(h4, 'LineStyle', '-', 'Color', [.5 0 .5], 'LineWidth', 1.5);
legend('Empirical C', 'Empirical I', 'Predicted C', 'Predicted I', 'Location', 'best')
title('SOA = +100 ms')
hold off

nexttile
h3 = cdfplot(C_SOA_2);
hold on 
h4 = cdfplot(I_SOA_2);
hold on 
h5 = cdfplot(C_predicted{1, 2}(:, 1));
hold on 
h6 = cdfplot(I_predicted{1, 2}(:, 1));
xlim([200 1000]);
xlabel('RT (ms)');
ylabel('CDF');
set(h3, 'LineStyle', ':', 'Color', [0 0.5 0], 'LineWidth', 1.5);
set(h4, 'LineStyle', ':', 'Color', [.5 0 .5], 'LineWidth', 1.5);
set(h5, 'LineStyle', '-', 'Color', [0 0.5 0], 'LineWidth', 1.5);
set(h6, 'LineStyle', '-', 'Color', [.5 0 .5], 'LineWidth', 1.5);
legend('Empirical C', 'Empirical I', 'Predicted C', 'Predicted I', 'Location', 'best')
title('SOA = 0 ms')
hold off

nexttile
h3 = cdfplot(C_SOA_3);
hold on 
h4 = cdfplot(I_SOA_3);
hold on 
h5 = cdfplot(C_predicted{1, 3}(:, 1));
hold on 
h6 = cdfplot(I_predicted{1, 3}(:, 1));
xlim([200 1000]);
xlabel('RT (ms)');
ylabel('CDF');
set(h3, 'LineStyle', ':', 'Color', [0 0.5 0], 'LineWidth', 1.5);
set(h4, 'LineStyle', ':', 'Color', [.5 0 .5], 'LineWidth', 1.5);
set(h5, 'LineStyle', '-', 'Color', [0 0.5 0], 'LineWidth', 1.5);
set(h6, 'LineStyle', '-', 'Color', [.5 0 .5], 'LineWidth', 1.5);
legend('Empirical C', 'Empirical I', 'Predicted C', 'Predicted I', 'Location', 'best')
title('SOA = -50 ms')
hold off

nexttile
h3 = cdfplot(C_SOA_4);
hold on 
h4 = cdfplot(I_SOA_4);
hold on 
h5 = cdfplot(C_predicted{1, 4}(:, 1));
hold on 
h6 = cdfplot(I_predicted{1, 4}(:, 1));
xlim([200 1000]);
xlabel('RT (ms)');
ylabel('CDF');
set(h3, 'LineStyle', ':', 'Color', [0 0.5 0], 'LineWidth', 1.5);
set(h4, 'LineStyle', ':', 'Color', [.5 0 .5], 'LineWidth', 1.5);
set(h5, 'LineStyle', '-', 'Color', [0 0.5 0], 'LineWidth', 1.5);
set(h6, 'LineStyle', '-', 'Color', [.5 0 .5], 'LineWidth', 1.5);
legend('Empirical C', 'Empirical I', 'Predicted C', 'Predicted I', 'Location', 'best')
title('SOA = -100 ms')
hold off

exportgraphics(t,'SpatialStroop_CE_SOA_CDF.pdf')


%% Delta plots (CE x SOA)

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

C = cell(NumP, 1);
I = cell(NumP, 1);
for n=1:NumP
    for k=1:size(Ppt{n},1)
        if Ppt{n}(k,7)==1
            C{n} = [C{n} ; Ppt{n}(k,1) Ppt{n}(k,2) Ppt{n}(k,5) Ppt{n}(k,6), Ppt{n}(k,8)]; 
        elseif Ppt{n}(k,7)==0 
            I{n} = [I{n} ; Ppt{n}(k,1) Ppt{n}(k,2) Ppt{n}(k,5) Ppt{n}(k,6), Ppt{n}(k,8)];
        end
    end
end

C = cell2mat(C);
de3=C(:,3)<200;
C(de3,:) = [];
de4=C(:,3)>1500;
C(de4,:) = [];
de2=C(:,4)==0;
C(de2,:) = [];

I = cell2mat(I);
de7=I(:,3)<200;
I(de7,:) = [];
de8=I(:,3)>1500;
I(de8,:) = [];
de6=I(:,4)==0;
I(de6,:) = [];

C1=cell(NumP, 1);
C2=cell(NumP, 1);
for n=1:NumP
    k=C(:, 1)==n;
    C1{n} = C(k, :);
    C2{n} = [C1{n}(:, 3), C1{n}(:, 5)];
end
C = C2;

I1=cell(NumP, 1);
I2=cell(NumP, 1);
for n=1:NumP
    l=I(:, 1)==n;
    I1{n} = I(l, :);
    I2{n} = [I1{n}(:, 3), I1{n}(:, 5)];
end
I=I2;

% We have different SOAs in the spatialStroop task, -100s (target shown before
% distractors), 0, 50, 100s
% Five quantiles for each condition

RTc = cell(n, 4);
for n=1:NumP
    RTc{n, 1} = C{n}(C{n}(:, 2)==100,1);
    RTc{n, 1} = [prctile(RTc{n, 1}(:,1),10), prctile(RTc{n, 1}(:,1),30), prctile(RTc{n, 1}(:,1),50), prctile(RTc{n, 1}(:,1),70), ...
        prctile(RTc{n, 1}(:,1),90)];
    RTc{n, 2} = C{n}(C{n}(:, 2)==0,1);
    RTc{n, 2} = [prctile(RTc{n, 2}(:,1),10), prctile(RTc{n, 2}(:,1),30), prctile(RTc{n, 2}(:,1),50), prctile(RTc{n, 2}(:,1),70), ...
        prctile(RTc{n, 2}(:,1),90)];
    RTc{n, 3} = C{n}(C{n}(:, 2)==-50,1);
    RTc{n, 3} = [prctile(RTc{n, 3}(:,1),10), prctile(RTc{n, 3}(:,1),30), prctile(RTc{n, 3}(:,1),50), prctile(RTc{n, 3}(:,1),70), ...
        prctile(RTc{n, 3}(:,1),90)];
    RTc{n, 4} = C{n}(C{n}(:, 2)==-100,1);
    RTc{n, 4} = [prctile(RTc{n, 4}(:,1),10), prctile(RTc{n, 4}(:,1),30), prctile(RTc{n, 4}(:,1),50), prctile(RTc{n, 4}(:,1),70), ...
        prctile(RTc{n, 4}(:,1),90)];
end

RTi = cell(n, 4);
for n=1:NumP
    RTi{n, 1} = I{n}(I{n}(:, 2)==100,1);
    RTi{n, 1} = [prctile(RTi{n, 1}(:,1),10), prctile(RTi{n, 1}(:,1),30), prctile(RTi{n, 1}(:,1),50), prctile(RTi{n, 1}(:,1),70), ...
        prctile(RTi{n, 1}(:,1),90)];
    RTi{n, 2} = I{n}(I{n}(:, 2)==0,1);
    RTi{n, 2} = [prctile(RTi{n, 2}(:,1),10), prctile(RTi{n, 2}(:,1),30), prctile(RTi{n, 2}(:,1),50), prctile(RTi{n, 2}(:,1),70), ...
        prctile(RTi{n, 2}(:,1),90)];
    RTi{n, 3} = I{n}(I{n}(:, 2)==-50,1);
    RTi{n, 3} = [prctile(RTi{n, 3}(:,1),10), prctile(RTi{n, 3}(:,1),30), prctile(RTi{n, 3}(:,1),50), prctile(RTi{n, 3}(:,1),70), ...
        prctile(RTi{n, 3}(:,1),90)];
    RTi{n, 4} = I{n}(I{n}(:, 2)==-100,1);
    RTi{n, 4} = [prctile(RTi{n, 4}(:,1),10), prctile(RTi{n, 4}(:,1),30), prctile(RTi{n, 4}(:,1),50), prctile(RTi{n, 4}(:,1),70), ...
        prctile(RTi{n, 4}(:,1),90)];
end

CE=cell(NumP, 4);
for n=1:NumP
    CE{n, 1} = [RTi{n, 1}(1) - RTc{n, 1}(1), RTi{n, 1}(2) - RTc{n, 1}(2), RTi{n, 1}(3) - RTc{n, 1}(3), RTi{n, 1}(4) - RTc{n, 1}(4), ...
        RTi{n, 1}(5) - RTc{n, 1}(5)];
    CE{n, 2} = [RTi{n, 2}(1) - RTc{n, 2}(1), RTi{n, 2}(2) - RTc{n, 2}(2), RTi{n, 2}(3) - RTc{n, 2}(3), RTi{n, 2}(4) - RTc{n, 2}(4), ...
        RTi{n, 2}(5) - RTc{n, 2}(5)];
    CE{n, 3} = [RTi{n, 3}(1) - RTc{n, 3}(1), RTi{n, 3}(2) - RTc{n, 3}(2), RTi{n, 3}(3) - RTc{n, 3}(3), RTi{n, 3}(4) - RTc{n, 3}(4), ...
        RTi{n, 3}(5) - RTc{n, 3}(5)];
    CE{n, 4} = [RTi{n, 4}(1) - RTc{n, 4}(1), RTi{n, 4}(2) - RTc{n, 4}(2), RTi{n, 4}(3) - RTc{n, 4}(3), RTi{n, 4}(4) - RTc{n, 4}(4), ...
        RTi{n, 4}(5) - RTc{n, 4}(5)];
end

RT_c_1=[];
RT_c_2=[];
RT_c_3=[];
RT_c_4=[];
for n=1:NumP
    RT_c_1=[RT_c_1 ; RTc{n, 1}];
    RT_c_2=[RT_c_2 ; RTc{n, 2}];
    RT_c_3=[RT_c_3 ; RTc{n, 3}];
    RT_c_4=[RT_c_4 ; RTc{n, 4}];
end
meanRTc=[mean(RT_c_1) ; mean(RT_c_2); mean(RT_c_3); mean(RT_c_4)];

RT_i_1=[];
RT_i_2=[];
RT_i_3=[];
RT_i_4=[];
for n=1:NumP
    RT_i_1=[RT_i_1 ; RTi{n, 1}];
    RT_i_2=[RT_i_2 ; RTi{n, 2}];
    RT_i_3=[RT_i_3 ; RTi{n, 3}];
    RT_i_4=[RT_i_4 ; RTi{n, 4}];
end
meanRTi=[mean(RT_i_1) ; mean(RT_i_2); mean(RT_i_3); mean(RT_i_4)];

RT=[];
for n=1:4
    RT = [RT ; (meanRTc(n, 1)+meanRTi(n, 1))/2, (meanRTc(n, 2)+meanRTi(n, 2))/2, (meanRTc(n, 3)+meanRTi(n, 3))/2, (meanRTc(n, 4)+meanRTi(n, 4))/2, ...
        (meanRTc(n, 5)+meanRTi(n, 5))/2];
end

CE_1=[];
CE_2=[];
CE_3=[];
CE_4=[];
for n=1:NumP
    CE_1=[CE_1 ; CE{n, 1}];
    CE_2=[CE_2 ; CE{n, 2}];
    CE_3=[CE_3 ; CE{n, 3}];
    CE_4=[CE_4 ; CE{n, 4}];
end

meanCE=[mean(CE_1) ; mean(CE_2); mean(CE_3); mean(CE_4)];

SEM = [std(CE_1)/sqrt(length(CE_1)); std(CE_2)/sqrt(length(CE_2)); std(CE_3)/sqrt(length(CE_3)); std(CE_4)/sqrt(length(CE_4))];
ts = tinv([0.025  0.975],length(CE_1)-1);

inv={};
for r=1:4
    inv{1, r} = ts.*SEM(r, :)';
end

CI={};
for r=1:4
    for j=1:5
        CI{1, r}(j, :) = meanCE(r, j) + inv{1, r}(j, :);
    end
end

% Simulated data
% Load the estimated parameters from the best performing model

% Model Prediction
% cC
Number_of_simulation = 50000;

load('spatialStroop_CE_SOA_SIMPLEX_02.mat')

A_01=fitted_Param(1);
kc1=fitted_Param(2);
ki1=fitted_Param(3);
da1=fitted_Param(4);
dc1=fitted_Param(5);
a1=fitted_Param(6);
Ter1=fitted_Param(7);
Ter2=fitted_Param(8);
Ter3=fitted_Param(9);
Ter4=fitted_Param(10);
std=fitted_Param(11);
A_02=A_01;
A_03=A_01;
A_04=A_01;
kc2=kc1;
kc3=kc1;
kc4=kc1;
ki2=ki1;
ki3=ki1;
ki4=ki1;
da2=da1;
da3=da1;
da4=da1;
dc2=dc1;
dc3=dc1;
dc4=dc1;
a2=a1;
a3=a1;
a4=a1;
c1=0;
c2=c1;
c3=c1;
c4=c1;

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
    C_predicted{1, n}(C_predicted{1, n}(:, 2)==0,:)=[];
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
    I_predicted{1, n}(I_predicted{1, n}(:, 2)==0,:)=[];
end

RTc_predicted={};
for n=1:4
    RTc_predicted{1, n} = [prctile(C_predicted{1, n}(:,1),10), prctile(C_predicted{1, n}(:,1),30), prctile(C_predicted{1, n}(:,1),50),...
            prctile(C_predicted{1, n}(:,1),70), prctile(C_predicted{1, n}(:,1),90)];
end

RTi_predicted={};
for n=1:4
    RTi_predicted{1, n} = [prctile(I_predicted{1, n}(:,1),10), prctile(I_predicted{1, n}(:,1),30), prctile(I_predicted{1, n}(:,1),50),...
            prctile(I_predicted{1, n}(:,1),70), prctile(I_predicted{1, n}(:,1),90)];
end

CE_predicted=[];
for n=1:4
    CE_predicted(n, :) = [RTi_predicted{1, n}(1) - RTc_predicted{1, n}(1), RTi_predicted{1, n}(2) - RTc_predicted{1, n}(2), ...
        RTi_predicted{1, n}(3) - RTc_predicted{1, n}(3), RTi_predicted{1, n}(4) - RTc_predicted{1, n}(4), RTi_predicted{1, n}(5) - RTc_predicted{1, n}(5)];
end

RT_predicted=[];
for n=1:4
    RT_predicted(n, :) = (RTc_predicted{1, n}+RTi_predicted{1, n})/2;
end

% Plot the delta plot of each SOA 
figure('Position', [300 300 1000 700])
t= tiledlayout(2,2);
nexttile
y = meanCE(1, :);
x = RT(1, :);
g1 = plot(RT(1, :), meanCE(1, :));
hold on 
g2 = plot(RT_predicted(1, :), CE_predicted(1, :));
set(g1, 'LineStyle', ':', 'Marker', 'o', 'Color', [0 0 0], 'LineWidth', 1.5);
set(g2, 'LineStyle', '-', 'Marker', 'o', 'Color', [0 0 0], 'LineWidth', 1.5);
xlim([350 650]);
ylim([0 120]);
xlabel('RT (ms)')
ylabel('Congruency effects (ms)')
title('SOA = +100 ms');
hold off

nexttile
y = meanCE(2, :);
x = RT(2, :);
g1 = plot(RT(2, :), meanCE(2, :));
hold on 
g2 = plot(RT_predicted(2, :), CE_predicted(2, :));
set(g1, 'LineStyle', ':', 'Marker', 'o', 'Color', [0 0 0], 'LineWidth', 1.5);
set(g2, 'LineStyle', '-', 'Marker', 'o', 'Color', [0 0 0], 'LineWidth', 1.5);
xlim([350 650]);
ylim([0 120]);
xlabel('RT (ms)')
ylabel('Congruency effects (ms)')
title('SOA = 0 ms');
hold off

nexttile
y = meanCE(3, :);
x = RT(3, :);
g1 = plot(RT(3, :), meanCE(3, :));
hold on 
g2 = plot(RT_predicted(3, :), CE_predicted(3, :));
set(g1, 'LineStyle', ':', 'Marker', 'o', 'Color', [0 0 0], 'LineWidth', 1.5);
set(g2, 'LineStyle', '-', 'Marker', 'o', 'Color', [0 0 0], 'LineWidth', 1.5);
xlim([350 650]);
ylim([0 120]);
xlabel('RT (ms)')
ylabel('Congruency effects (ms)')
title('SOA = -50 ms');
hold off

nexttile
y = meanCE(4, :);
x = RT(4, :);
g1 = plot(RT(4, :), meanCE(4, :));
hold on 
g2 = plot(RT_predicted(4, :), CE_predicted(4, :));
set(g1, 'LineStyle', ':', 'Marker', 'o', 'Color', [0 0 0], 'LineWidth', 1.5);
set(g2, 'LineStyle', '-', 'Marker', 'o', 'Color', [0 0 0], 'LineWidth', 1.5);
xlim([350 650]);
ylim([0 120]);
xlabel('RT (ms)')
ylabel('Congruency effects (ms)')
title('SOA = -100 ms');
hold off

exportgraphics(t,'SpatialStroop_CE_SOA_DP.pdf')
%% CAF (CE x SOA)

clear
% Empirical spatialStroop data
spatialStroop=readtable('spatialStroopSOA.xlsx');
spatialStroop=table2array(spatialStroop);
warning('off')
NumP= 40;

Ppt = cell(NumP, 1);
for k=1:NumP
    Ppt{k} = spatialStroop(spatialStroop(:,1)==k,:);
end

C = cell(NumP, 1);
I = cell(NumP, 1);

for n=1:NumP
    for k=1:size(Ppt{n},1)-1
        if Ppt{n}(k,7)==1 
            C{n} = [C{n} ; Ppt{n}(k,1) Ppt{n}(k,2) Ppt{n}(k,5) Ppt{n}(k,6), Ppt{n}(k,8)];
        elseif Ppt{n}(k,7)==0
            I{n} = [I{n} ; Ppt{n}(k,1) Ppt{n}(k,2) Ppt{n}(k,5) Ppt{n}(k,6), Ppt{n}(k,8)];
        end
    end
end

C = cell2mat(C);
de_3=C(:,3)<200;
C(de_3,:) = [];
de_4=C(:,3)>1500;
C(de_4,:) = [];

I = cell2mat(I);
de7=I(:,3)<200;
I(de7,:) = [];
de8=I(:,3)>1500;
I(de8,:) = [];

C1=cell(NumP, 1);
C2=cell(NumP, 1);
for n=1:NumP
    k=C(:, 1)==n;
    C1{n} = C(k, :);
    C2{n} = [C1{n}(:, 3), C1{n}(:, 4), C1{n}(:, 5)];
end
C = C2;

I1=cell(NumP, 1);
I2=cell(NumP, 1);
for n=1:NumP
    l=I(:, 1)==n;
    I1{n} = I(l, :);
    I2{n} = [I1{n}(:, 3), I1{n}(:, 4), I1{n}(:, 5)];
end
I=I2;

RTc = cell(n, 4);
for n=1:NumP
    RTc{n, 1} = C{n}(C{n}(:, 3)==100,:);
    RTc{n, 2} = C{n}(C{n}(:, 3)==0,:);
    RTc{n, 3} = C{n}(C{n}(:, 3)==-50,:);
    RTc{n, 4} = C{n}(C{n}(:, 3)==-100,:);
end

RTi = cell(n, 4);
for n=1:NumP
    RTi{n, 1} = I{n}(I{n}(:, 3)==100,:);
    RTi{n, 2} = I{n}(I{n}(:, 3)==0,:);
    RTi{n, 3} = I{n}(I{n}(:, 3)==-50,:); 
    RTi{n, 4} = I{n}(I{n}(:, 3)==-100,:);
end

% Four RT bins (0-20%, 20-40%, 40-60%, 60-80%, 80-100%)
bin_c = cell(NumP, 4);
for n=1:NumP
    bin_c{n, 1}=[prctile(RTc{n, 1}(:,1),20), prctile(RTc{n, 1}(:,1),40), prctile(RTc{n, 1}(:,1),60), prctile(RTc{n, 1}(:,1),80), prctile(RTc{n, 1}(:,1),100)];
    bin_c{n, 2}=[prctile(RTc{n, 2}(:,1),20), prctile(RTc{n, 2}(:,1),40), prctile(RTc{n, 2}(:,1),60), prctile(RTc{n, 2}(:,1),80), prctile(RTc{n, 2}(:,1),100)];
    bin_c{n, 3}=[prctile(RTc{n, 3}(:,1),20), prctile(RTc{n, 3}(:,1),40), prctile(RTc{n, 3}(:,1),60), prctile(RTc{n, 3}(:,1),80), prctile(RTc{n, 3}(:,1),100)];
    bin_c{n, 4}=[prctile(RTc{n, 4}(:,1),20), prctile(RTc{n, 4}(:,1),40), prctile(RTc{n, 4}(:,1),60), prctile(RTc{n, 4}(:,1),80), prctile(RTc{n, 4}(:,1),100)];
end

corrc = cell(NumP, 4);
incorrc = cell(NumP, 4);

for n=1:NumP 
    for m=1:4
        for x=1:size(RTc{n, m},1)
            corrc{n, m} = [corrc{n, m} ; RTc{n, m}(x,1)<=bin_c{n, m}(:,1) && RTc{n, m}(x,2)==1, ...
               RTc{n, m}(x,1)>bin_c{n, m}(:,1) & RTc{n, m}(x,1)<=bin_c{n, m}(:,2) & RTc{n, m}(x,2)==1, ...
               RTc{n, m}(x,1)>bin_c{n, m}(:,2) & RTc{n, m}(x,1)<=bin_c{n, m}(:,3) & RTc{n, m}(x,2)==1, ...
               RTc{n, m}(x,1)>bin_c{n, m}(:,3) & RTc{n, m}(x,1)<=bin_c{n, m}(:,4) & RTc{n, m}(x,2)==1, ...
               RTc{n, m}(x,1)>bin_c{n, m}(:,4) & RTc{n, m}(x,1)<=bin_c{n, m}(:,5) & RTc{n, m}(x,2)==1];
            incorrc{n, m} = [incorrc{n, m} ; RTc{n, m}(x,1)<=bin_c{n, m}(:,1) && RTc{n, m}(x,2)~=1, ...
               RTc{n, m}(x,1)>bin_c{n, m}(:,1) & RTc{n, m}(x,1)<=bin_c{n, m}(:,2) & RTc{n, m}(x,2)~=1, ...
               RTc{n, m}(x,1)>bin_c{n, m}(:,2) & RTc{n, m}(x,1)<=bin_c{n, m}(:,3) & RTc{n, m}(x,2)~=1, ...
               RTc{n, m}(x,1)>bin_c{n, m}(:,3) & RTc{n, m}(x,1)<=bin_c{n, m}(:,4) & RTc{n, m}(x,2)~=1, ...
               RTc{n, m}(x,1)>bin_c{n, m}(:,4) & RTc{n, m}(x,1)<=bin_c{n, m}(:,5) & RTc{n, m}(x,2)~=1];
            
        end
    end
end

bin_i = cell(NumP, 4);
for n=1:NumP
    bin_i{n, 1}=[prctile(RTi{n, 1}(:,1),20), prctile(RTi{n, 1}(:,1),40), prctile(RTi{n, 1}(:,1),60), prctile(RTi{n, 1}(:,1),80), prctile(RTi{n, 1}(:,1),100)];
    bin_i{n, 2}=[prctile(RTi{n, 2}(:,1),20), prctile(RTi{n, 2}(:,1),40), prctile(RTi{n, 2}(:,1),60), prctile(RTi{n, 2}(:,1),80), prctile(RTi{n, 2}(:,1),100)];
    bin_i{n, 3}=[prctile(RTi{n, 3}(:,1),20), prctile(RTi{n, 3}(:,1),40), prctile(RTi{n, 3}(:,1),60), prctile(RTi{n, 3}(:,1),80), prctile(RTi{n, 3}(:,1),100)];
    bin_i{n, 4}=[prctile(RTi{n, 4}(:,1),20), prctile(RTi{n, 4}(:,1),40), prctile(RTi{n, 4}(:,1),60), prctile(RTi{n, 4}(:,1),80), prctile(RTi{n, 4}(:,1),100)];
end

corri = cell(NumP, 4);
incorri = cell(NumP, 4);

for n=1:NumP 
    for m=1:4
        for x=1:size(RTi{n, m},1)
            corri{n, m} = [corri{n, m} ; RTi{n, m}(x,1)<=bin_i{n, m}(:,1) && RTi{n, m}(x,2)==1, ...
               RTi{n, m}(x,1)>bin_i{n, m}(:,1) & RTi{n, m}(x,1)<=bin_i{n, m}(:,2) & RTi{n, m}(x,2)==1, ...
               RTi{n, m}(x,1)>bin_i{n, m}(:,2) & RTi{n, m}(x,1)<=bin_i{n, m}(:,3) & RTi{n, m}(x,2)==1, ...
               RTi{n, m}(x,1)>bin_i{n, m}(:,3) & RTi{n, m}(x,1)<=bin_i{n, m}(:,4) & RTi{n, m}(x,2)==1, ...
               RTi{n, m}(x,1)>bin_i{n, m}(:,4) & RTi{n, m}(x,1)<=bin_i{n, m}(:,5) & RTi{n, m}(x,2)==1];
            incorri{n, m} = [incorri{n, m} ; RTi{n, m}(x,1)<=bin_i{n, m}(:,1) && RTi{n, m}(x,2)~=1, ...
               RTi{n, m}(x,1)>bin_i{n, m}(:,1) & RTi{n, m}(x,1)<=bin_i{n, m}(:,2) & RTi{n, m}(x,2)~=1, ...
               RTi{n, m}(x,1)>bin_i{n, m}(:,2) & RTi{n, m}(x,1)<=bin_i{n, m}(:,3) & RTi{n, m}(x,2)~=1, ...
               RTi{n, m}(x,1)>bin_i{n, m}(:,3) & RTi{n, m}(x,1)<=bin_i{n, m}(:,4) & RTi{n, m}(x,2)~=1, ...
               RTi{n, m}(x,1)>bin_i{n, m}(:,4) & RTi{n, m}(x,1)<=bin_i{n, m}(:,5) & RTi{n, m}(x,2)~=1];
            
        end
    end
end

p_c=[];
p_c(1, :) = sum(cell2mat(corrc(:, 1)))./(sum(cell2mat(corrc(:, 1)))+sum(cell2mat(incorrc(:, 1))));
p_c(2, :) = sum(cell2mat(corrc(:, 2)))./(sum(cell2mat(corrc(:, 2)))+sum(cell2mat(incorrc(:, 2))));
p_c(3, :) = sum(cell2mat(corrc(:, 3)))./(sum(cell2mat(corrc(:, 3)))+sum(cell2mat(incorrc(:, 3))));
p_c(4, :) = sum(cell2mat(corrc(:, 4)))./(sum(cell2mat(corrc(:, 4)))+sum(cell2mat(incorrc(:, 4))));

p_i=[];
p_i(1, :) = sum(cell2mat(corri(:, 1)))./(sum(cell2mat(corri(:, 1)))+sum(cell2mat(incorri(:, 1))));
p_i(2, :) = sum(cell2mat(corri(:, 2)))./(sum(cell2mat(corri(:, 2)))+sum(cell2mat(incorri(:, 2))));
p_i(3, :) = sum(cell2mat(corri(:, 3)))./(sum(cell2mat(corri(:, 3)))+sum(cell2mat(incorri(:, 3))));
p_i(4, :) = sum(cell2mat(corri(:, 4)))./(sum(cell2mat(corri(:, 4)))+sum(cell2mat(incorri(:, 4))));

% Simulated data
% Load the estimated parameters from the best performing model

% Model Prediction
% C
Number_of_simulation = 50000;

load('spatialStroop_CE_SOA_SIMPLEX_02.mat')

A_01=fitted_Param(1);
kc1=fitted_Param(2);
ki1=fitted_Param(3);
da1=fitted_Param(4);
dc1=fitted_Param(5);
a1=fitted_Param(6);
Ter1=fitted_Param(7);
Ter2=fitted_Param(8);
Ter3=fitted_Param(9);
Ter4=fitted_Param(10);
std=fitted_Param(11);
A_02=A_01;
A_03=A_01;
A_04=A_01;
kc2=kc1;
kc3=kc1;
kc4=kc1;
ki2=ki1;
ki3=ki1;
ki4=ki1;
da2=da1;
da3=da1;
da4=da1;
dc2=dc1;
dc3=dc1;
dc4=dc1;
a2=a1;
a3=a1;
a4=a1;
c1=0;
c2=c1;
c3=c1;
c4=c1;

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

bin_c_pred = {};
for m=1:4
    bin_c_pred{1, m}=[prctile(C_predicted{1, m}(:,1),20), prctile(C_predicted{1, m}(:,1),40), prctile(C_predicted{1, m}(:,1),60), ...
        prctile(C_predicted{1, m}(:,1),80), prctile(C_predicted{1, m}(:,1),100)];
end

corrc_pred = {};
incorrc_pred = {};

for m=1:4
    corrc_pred{1, m} = [C_predicted{1, m}(:,1)<=bin_c_pred{1, m}(:,1) & C_predicted{1, m}(:,2)==1, ...
        C_predicted{1, m}(:,1)>bin_c_pred{1, m}(:,1) & C_predicted{1, m}(:,1)<=bin_c_pred{1, m}(:,2) & C_predicted{1, m}(:,2)==1, ...
        C_predicted{1, m}(:,1)>bin_c_pred{1, m}(:,2) & C_predicted{1, m}(:,1)<=bin_c_pred{1, m}(:,3) & C_predicted{1, m}(:,2)==1, ...
        C_predicted{1, m}(:,1)>bin_c_pred{1, m}(:,3) & C_predicted{1, m}(:,1)<=bin_c_pred{1, m}(:,4) & C_predicted{1, m}(:,2)==1, ...
        C_predicted{1, m}(:,1)>bin_c_pred{1, m}(:,4) & C_predicted{1, m}(:,1)<=bin_c_pred{1, m}(:,5) & C_predicted{1, m}(:,2)==1];
    incorrc_pred{1, m} = [C_predicted{1, m}(:,1)<=bin_c_pred{1, m}(:,1) & C_predicted{1, m}(:,2)==0, ...
        C_predicted{1, m}(:,1)>bin_c_pred{1, m}(:,1) & C_predicted{1, m}(:,1)<=bin_c_pred{1, m}(:,2) & C_predicted{1, m}(:,2)==0, ...
        C_predicted{1, m}(:,1)>bin_c_pred{1, m}(:,2) & C_predicted{1, m}(:,1)<=bin_c_pred{1, m}(:,3) & C_predicted{1, m}(:,2)==0, ...
        C_predicted{1, m}(:,1)>bin_c_pred{1, m}(:,3) & C_predicted{1, m}(:,1)<=bin_c_pred{1, m}(:,4) & C_predicted{1, m}(:,2)==0, ...
        C_predicted{1, m}(:,1)>bin_c_pred{1, m}(:,4) & C_predicted{1, m}(:,1)<=bin_c_pred{1, m}(:,5) & C_predicted{1, m}(:,2)==0];   
end


bin_i_pred = {};
for m=1:4
    bin_i_pred{1, m}=[prctile(I_predicted{1, m}(:,1),20), prctile(I_predicted{1, m}(:,1),40), prctile(I_predicted{1, m}(:,1),60), ...
        prctile(I_predicted{1, m}(:,1),80), prctile(I_predicted{1, m}(:,1),100)];
end

corri_pred = {};
incorri_pred = {};

for m=1:4
    corri_pred{1, m} = [I_predicted{1, m}(:,1)<=bin_i_pred{1, m}(:,1) & I_predicted{1, m}(:,2)==1, ...
        I_predicted{1, m}(:,1)>bin_i_pred{1, m}(:,1) & I_predicted{1, m}(:,1)<=bin_i_pred{1, m}(:,2) & I_predicted{1, m}(:,2)==1, ...
        I_predicted{1, m}(:,1)>bin_i_pred{1, m}(:,2) & I_predicted{1, m}(:,1)<=bin_i_pred{1, m}(:,3) & I_predicted{1, m}(:,2)==1, ...
        I_predicted{1, m}(:,1)>bin_i_pred{1, m}(:,3) & I_predicted{1, m}(:,1)<=bin_i_pred{1, m}(:,4) & I_predicted{1, m}(:,2)==1, ...
        I_predicted{1, m}(:,1)>bin_i_pred{1, m}(:,4) & I_predicted{1, m}(:,1)<=bin_i_pred{1, m}(:,5) & I_predicted{1, m}(:,2)==1];
    incorri_pred{1, m} = [I_predicted{1, m}(:,1)<=bin_i_pred{1, m}(:,1) & I_predicted{1, m}(:,2)==0, ...
        I_predicted{1, m}(:,1)>bin_i_pred{1, m}(:,1) & I_predicted{1, m}(:,1)<=bin_i_pred{1, m}(:,2) & I_predicted{1, m}(:,2)==0, ...
        I_predicted{1, m}(:,1)>bin_i_pred{1, m}(:,2) & I_predicted{1, m}(:,1)<=bin_i_pred{1, m}(:,3) & I_predicted{1, m}(:,2)==0, ...
        I_predicted{1, m}(:,1)>bin_i_pred{1, m}(:,3) & I_predicted{1, m}(:,1)<=bin_i_pred{1, m}(:,4) & I_predicted{1, m}(:,2)==0, ...
        I_predicted{1, m}(:,1)>bin_i_pred{1, m}(:,4) & I_predicted{1, m}(:,1)<=bin_i_pred{1, m}(:,5) & I_predicted{1, m}(:,2)==0];   
end

p_c_pred=[];
for m=1:4
    p_c_pred(m, :) = sum(cell2mat(corrc_pred(:, m)))./(sum(cell2mat(corrc_pred(:, m)))+sum(cell2mat(incorrc_pred(:, m))));
end

p_i_pred=[];
for m=1:4
    p_i_pred(m, :) = sum(cell2mat(corri_pred(:, m)))./(sum(cell2mat(corri_pred(:, m)))+sum(cell2mat(incorri_pred(:, m))));
end

% Plot CAF of each SOA
figure('Position', [300 300 1000 700])
t= tiledlayout(2,2);
nexttile
RT20 = [p_c(1, 1), p_c_pred(1, 1), p_i(1, 1), p_i_pred(1, 1)];
RT40 = [p_c(1, 2), p_c_pred(1, 2), p_i(1, 2), p_i_pred(1, 2)];
RT60 = [p_c(1, 3), p_c_pred(1, 3), p_i(1, 3), p_i_pred(1, 3)];
RT80 = [p_c(1, 4), p_c_pred(1, 4), p_i(1, 4), p_i_pred(1, 4)];
RT100 = [p_c(1, 5), p_c_pred(1, 5), p_i(1, 5), p_i_pred(1, 5)];
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
legend('Empirical C', 'Predicted C', 'Empirical I', 'Predicted I', 'Location', 'best')
title('SOA = +100 ms');
hold off

nexttile
RT20 = [p_c(2, 1), p_c_pred(2, 1), p_i(2, 1), p_i_pred(2, 1)];
RT40 = [p_c(2, 2), p_c_pred(2, 2), p_i(2, 2), p_i_pred(2, 2)];
RT60 = [p_c(2, 3), p_c_pred(2, 3), p_i(2, 3), p_i_pred(2, 3)];
RT80 = [p_c(2, 4), p_c_pred(2, 4), p_i(2, 4), p_i_pred(2, 4)];
RT100 = [p_c(2, 5), p_c_pred(2, 5), p_i(2, 5), p_i_pred(2, 5)];
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
legend('Empirical C', 'Predicted C', 'Empirical I', 'Predicted I', 'Location', 'best')
title('SOA = 0 ms');
hold off

nexttile
RT20 = [p_c(3, 1), p_c_pred(3, 1), p_i(3, 1), p_i_pred(3, 1)];
RT40 = [p_c(3, 2), p_c_pred(3, 2), p_i(3, 2), p_i_pred(3, 2)];
RT60 = [p_c(3, 3), p_c_pred(3, 3), p_i(3, 3), p_i_pred(3, 3)];
RT80 = [p_c(3, 4), p_c_pred(3, 4), p_i(3, 4), p_i_pred(3, 4)];
RT100 = [p_c(3, 5), p_c_pred(3, 5), p_i(3, 5), p_i_pred(3, 5)];
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
legend('Empirical C', 'Predicted C', 'Empirical I', 'Predicted I', 'Location', 'best')
title('SOA = -50 ms');
hold off

nexttile
RT20 = [p_c(4, 1), p_c_pred(4, 1), p_i(4, 1), p_i_pred(4, 1)];
RT40 = [p_c(4, 2), p_c_pred(4, 2), p_i(4, 2), p_i_pred(4, 2)];
RT60 = [p_c(4, 3), p_c_pred(4, 3), p_i(4, 3), p_i_pred(4, 3)];
RT80 = [p_c(4, 4), p_c_pred(4, 4), p_i(4, 4), p_i_pred(4, 4)];
RT100 = [p_c(4, 5), p_c_pred(4, 5), p_i(4, 5), p_i_pred(4, 5)];
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
legend('Empirical C', 'Predicted C', 'Empirical I', 'Predicted I', 'Location', 'best')
title('SOA = -100 ms');
hold off

exportgraphics(t,'SpatialStroop_CE_SOA_CAF.pdf')
