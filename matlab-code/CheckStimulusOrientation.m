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

C_hor = cell(NumP, 1);
I_hor = cell(NumP, 1);
C_ver = cell(NumP, 1);
I_ver = cell(NumP, 1);
for n=1:NumP
    for k=1:size(Ppt{n},1)
        if Ppt{n}(k,7)==1 && Ppt{n}(k,3) <=2
            C_hor{n} = [C_hor{n} ; Ppt{n}(k,1) Ppt{n}(k,2) Ppt{n}(k,5) Ppt{n}(k,6), Ppt{n}(k,8)]; 
        elseif Ppt{n}(k,7)==0 && Ppt{n}(k,3) <=2
            I_hor{n} = [I_hor{n} ; Ppt{n}(k,1) Ppt{n}(k,2) Ppt{n}(k,5) Ppt{n}(k,6), Ppt{n}(k,8)];
        elseif Ppt{n}(k,7)==1 && Ppt{n}(k,3) >2
            C_ver{n} = [C_ver{n} ; Ppt{n}(k,1) Ppt{n}(k,2) Ppt{n}(k,5) Ppt{n}(k,6), Ppt{n}(k,8)]; 
        elseif Ppt{n}(k,7)==0 && Ppt{n}(k,3) >2
            I_ver{n} = [I_ver{n} ; Ppt{n}(k,1) Ppt{n}(k,2) Ppt{n}(k,5) Ppt{n}(k,6), Ppt{n}(k,8)];
        end
    end
end

C_hor = cell2mat(C_hor);
de3=C_hor(:,3)<200;
C_hor(de3,:) = [];
de4=C_hor(:,3)>1500;
C_hor(de4,:) = [];
de2=C_hor(:,4)==0;
C_hor(de2,:) = [];

I_hor = cell2mat(I_hor);
de7=I_hor(:,3)<200;
I_hor(de7,:) = [];
de8=I_hor(:,3)>1500;
I_hor(de8,:) = [];
de6=I_hor(:,4)==0;
I_hor(de6,:) = [];

C_ver = cell2mat(C_ver);
de3=C_ver(:,3)<200;
C_ver(de3,:) = [];
de4=C_ver(:,3)>1500;
C_ver(de4,:) = [];
de2=C_ver(:,4)==0;
C_ver(de2,:) = [];

I_ver = cell2mat(I_ver);
de7=I_ver(:,3)<200;
I_ver(de7,:) = [];
de8=I_ver(:,3)>1500;
I_ver(de8,:) = [];
de6=I_ver(:,4)==0;
I_ver(de6,:) = [];

C1=cell(NumP, 1);
C2=cell(NumP, 1);
for n=1:NumP
    k=C_hor(:, 1)==n;
    C1{n} = C_hor(k, :);
    C2{n} = [C1{n}(:, 3), C1{n}(:, 5)];
end
C_hor = C2;

I1=cell(NumP, 1);
I2=cell(NumP, 1);
for n=1:NumP
    l=I_hor(:, 1)==n;
    I1{n} = I_hor(l, :);
    I2{n} = [I1{n}(:, 3), I1{n}(:, 5)];
end
I_hor=I2;

C1=cell(NumP, 1);
C2=cell(NumP, 1);
for n=1:NumP
    k=C_ver(:, 1)==n;
    C1{n} = C_ver(k, :);
    C2{n} = [C1{n}(:, 3), C1{n}(:, 5)];
end
C_ver = C2;

I1=cell(NumP, 1);
I2=cell(NumP, 1);
for n=1:NumP
    l=I_ver(:, 1)==n;
    I1{n} = I_ver(l, :);
    I2{n} = [I1{n}(:, 3), I1{n}(:, 5)];
end
I_ver=I2;

% We have different SOAs in the spatialStroop task, 100s (target shown before
% distractors), 0, -50, -100s separated by the orientation
% Five quantiles for each condition

RTc_hor = cell(n, 4);
for n=1:NumP
    RTc_hor{n, 1} = C_hor{n}(C_hor{n}(:, 2)==100,1);
    RTc_hor{n, 1} = [prctile(RTc_hor{n, 1}(:,1),10), prctile(RTc_hor{n, 1}(:,1),30), prctile(RTc_hor{n, 1}(:,1),50), prctile(RTc_hor{n, 1}(:,1),70), ...
        prctile(RTc_hor{n, 1}(:,1),90)];
    RTc_hor{n, 2} = C_hor{n}(C_hor{n}(:, 2)==0,1);
    RTc_hor{n, 2} = [prctile(RTc_hor{n, 2}(:,1),10), prctile(RTc_hor{n, 2}(:,1),30), prctile(RTc_hor{n, 2}(:,1),50), prctile(RTc_hor{n, 2}(:,1),70), ...
        prctile(RTc_hor{n, 2}(:,1),90)];
    RTc_hor{n, 3} = C_hor{n}(C_hor{n}(:, 2)==-50,1);
    RTc_hor{n, 3} = [prctile(RTc_hor{n, 3}(:,1),10), prctile(RTc_hor{n, 3}(:,1),30), prctile(RTc_hor{n, 3}(:,1),50), prctile(RTc_hor{n, 3}(:,1),70), ...
        prctile(RTc_hor{n, 3}(:,1),90)];
    RTc_hor{n, 4} = C_hor{n}(C_hor{n}(:, 2)==-100,1);
    RTc_hor{n, 4} = [prctile(RTc_hor{n, 4}(:,1),10), prctile(RTc_hor{n, 4}(:,1),30), prctile(RTc_hor{n, 4}(:,1),50), prctile(RTc_hor{n, 4}(:,1),70), ...
        prctile(RTc_hor{n, 4}(:,1),90)];
end

RTi_hor = cell(n, 4);
for n=1:NumP
    RTi_hor{n, 1} = I_hor{n}(I_hor{n}(:, 2)==100,1);
    RTi_hor{n, 1} = [prctile(RTi_hor{n, 1}(:,1),10), prctile(RTi_hor{n, 1}(:,1),30), prctile(RTi_hor{n, 1}(:,1),50), prctile(RTi_hor{n, 1}(:,1),70), ...
        prctile(RTi_hor{n, 1}(:,1),90)];
    RTi_hor{n, 2} = I_hor{n}(I_hor{n}(:, 2)==0,1);
    RTi_hor{n, 2} = [prctile(RTi_hor{n, 2}(:,1),10), prctile(RTi_hor{n, 2}(:,1),30), prctile(RTi_hor{n, 2}(:,1),50), prctile(RTi_hor{n, 2}(:,1),70), ...
        prctile(RTi_hor{n, 2}(:,1),90)];
    RTi_hor{n, 3} = I_hor{n}(I_hor{n}(:, 2)==-50,1);
    RTi_hor{n, 3} = [prctile(RTi_hor{n, 3}(:,1),10), prctile(RTi_hor{n, 3}(:,1),30), prctile(RTi_hor{n, 3}(:,1),50), prctile(RTi_hor{n, 3}(:,1),70), ...
        prctile(RTi_hor{n, 3}(:,1),90)];
    RTi_hor{n, 4} = I_hor{n}(I_hor{n}(:, 2)==-100,1);
    RTi_hor{n, 4} = [prctile(RTi_hor{n, 4}(:,1),10), prctile(RTi_hor{n, 4}(:,1),30), prctile(RTi_hor{n, 4}(:,1),50), prctile(RTi_hor{n, 4}(:,1),70), ...
        prctile(RTi_hor{n, 4}(:,1),90)];
end

RTc_ver = cell(n, 4);
for n=1:NumP
    RTc_ver{n, 1} = C_ver{n}(C_ver{n}(:, 2)==100,1);
    RTc_ver{n, 1} = [prctile(RTc_ver{n, 1}(:,1),10), prctile(RTc_ver{n, 1}(:,1),30), prctile(RTc_ver{n, 1}(:,1),50), prctile(RTc_ver{n, 1}(:,1),70), ...
        prctile(RTc_ver{n, 1}(:,1),90)];
    RTc_ver{n, 2} = C_ver{n}(C_ver{n}(:, 2)==0,1);
    RTc_ver{n, 2} = [prctile(RTc_ver{n, 2}(:,1),10), prctile(RTc_ver{n, 2}(:,1),30), prctile(RTc_ver{n, 2}(:,1),50), prctile(RTc_ver{n, 2}(:,1),70), ...
        prctile(RTc_ver{n, 2}(:,1),90)];
    RTc_ver{n, 3} = C_ver{n}(C_ver{n}(:, 2)==-50,1);
    RTc_ver{n, 3} = [prctile(RTc_ver{n, 3}(:,1),10), prctile(RTc_ver{n, 3}(:,1),30), prctile(RTc_ver{n, 3}(:,1),50), prctile(RTc_ver{n, 3}(:,1),70), ...
        prctile(RTc_ver{n, 3}(:,1),90)];
    RTc_ver{n, 4} = C_ver{n}(C_ver{n}(:, 2)==-100,1);
    RTc_ver{n, 4} = [prctile(RTc_ver{n, 4}(:,1),10), prctile(RTc_ver{n, 4}(:,1),30), prctile(RTc_ver{n, 4}(:,1),50), prctile(RTc_ver{n, 4}(:,1),70), ...
        prctile(RTc_ver{n, 4}(:,1),90)];
end

RTi_ver = cell(n, 4);
for n=1:NumP
    RTi_ver{n, 1} = I_ver{n}(I_ver{n}(:, 2)==100,1);
    RTi_ver{n, 1} = [prctile(RTi_ver{n, 1}(:,1),10), prctile(RTi_ver{n, 1}(:,1),30), prctile(RTi_ver{n, 1}(:,1),50), prctile(RTi_ver{n, 1}(:,1),70), ...
        prctile(RTi_ver{n, 1}(:,1),90)];
    RTi_ver{n, 2} = I_ver{n}(I_ver{n}(:, 2)==0,1);
    RTi_ver{n, 2} = [prctile(RTi_ver{n, 2}(:,1),10), prctile(RTi_ver{n, 2}(:,1),30), prctile(RTi_ver{n, 2}(:,1),50), prctile(RTi_ver{n, 2}(:,1),70), ...
        prctile(RTi_ver{n, 2}(:,1),90)];
    RTi_ver{n, 3} = I_ver{n}(I_ver{n}(:, 2)==-50,1);
    RTi_ver{n, 3} = [prctile(RTi_ver{n, 3}(:,1),10), prctile(RTi_ver{n, 3}(:,1),30), prctile(RTi_ver{n, 3}(:,1),50), prctile(RTi_ver{n, 3}(:,1),70), ...
        prctile(RTi_ver{n, 3}(:,1),90)];
    RTi_ver{n, 4} = I_ver{n}(I_ver{n}(:, 2)==-100,1);
    RTi_ver{n, 4} = [prctile(RTi_ver{n, 4}(:,1),10), prctile(RTi_ver{n, 4}(:,1),30), prctile(RTi_ver{n, 4}(:,1),50), prctile(RTi_ver{n, 4}(:,1),70), ...
        prctile(RTi_ver{n, 4}(:,1),90)];
end

CE_hor=cell(NumP, 4);
for n=1:NumP
    CE_hor{n, 1} = [RTi_hor{n, 1}(1) - RTc_hor{n, 1}(1), RTi_hor{n, 1}(2) - RTc_hor{n, 1}(2), RTi_hor{n, 1}(3) - RTc_hor{n, 1}(3), RTi_hor{n, 1}(4) - RTc_hor{n, 1}(4), ...
        RTi_hor{n, 1}(5) - RTc_hor{n, 1}(5)];
    CE_hor{n, 2} = [RTi_hor{n, 2}(1) - RTc_hor{n, 2}(1), RTi_hor{n, 2}(2) - RTc_hor{n, 2}(2), RTi_hor{n, 2}(3) - RTc_hor{n, 2}(3), RTi_hor{n, 2}(4) - RTc_hor{n, 2}(4), ...
        RTi_hor{n, 2}(5) - RTc_hor{n, 2}(5)];
    CE_hor{n, 3} = [RTi_hor{n, 3}(1) - RTc_hor{n, 3}(1), RTi_hor{n, 3}(2) - RTc_hor{n, 3}(2), RTi_hor{n, 3}(3) - RTc_hor{n, 3}(3), RTi_hor{n, 3}(4) - RTc_hor{n, 3}(4), ...
        RTi_hor{n, 3}(5) - RTc_hor{n, 3}(5)];
    CE_hor{n, 4} = [RTi_hor{n, 4}(1) - RTc_hor{n, 4}(1), RTi_hor{n, 4}(2) - RTc_hor{n, 4}(2), RTi_hor{n, 4}(3) - RTc_hor{n, 4}(3), RTi_hor{n, 4}(4) - RTc_hor{n, 4}(4), ...
        RTi_hor{n, 4}(5) - RTc_hor{n, 4}(5)];
end

CE_ver=cell(NumP, 4);
for n=1:NumP
    CE_ver{n, 1} = [RTi_ver{n, 1}(1) - RTc_ver{n, 1}(1), RTi_ver{n, 1}(2) - RTc_ver{n, 1}(2), RTi_ver{n, 1}(3) - RTc_ver{n, 1}(3), RTi_ver{n, 1}(4) - RTc_ver{n, 1}(4), ...
        RTi_ver{n, 1}(5) - RTc_ver{n, 1}(5)];
    CE_ver{n, 2} = [RTi_ver{n, 2}(1) - RTc_ver{n, 2}(1), RTi_ver{n, 2}(2) - RTc_ver{n, 2}(2), RTi_ver{n, 2}(3) - RTc_ver{n, 2}(3), RTi_ver{n, 2}(4) - RTc_ver{n, 2}(4), ...
        RTi_ver{n, 2}(5) - RTc_ver{n, 2}(5)];
    CE_ver{n, 3} = [RTi_ver{n, 3}(1) - RTc_ver{n, 3}(1), RTi_ver{n, 3}(2) - RTc_ver{n, 3}(2), RTi_ver{n, 3}(3) - RTc_ver{n, 3}(3), RTi_ver{n, 3}(4) - RTc_ver{n, 3}(4), ...
        RTi_ver{n, 3}(5) - RTc_ver{n, 3}(5)];
    CE_ver{n, 4} = [RTi_ver{n, 4}(1) - RTc_ver{n, 4}(1), RTi_ver{n, 4}(2) - RTc_ver{n, 4}(2), RTi_ver{n, 4}(3) - RTc_ver{n, 4}(3), RTi_ver{n, 4}(4) - RTc_ver{n, 4}(4), ...
        RTi_ver{n, 4}(5) - RTc_ver{n, 4}(5)];
end

RT_c_hor_1=[];
RT_c_hor_2=[];
RT_c_hor_3=[];
RT_c_hor_4=[];
for n=1:NumP
    RT_c_hor_1=[RT_c_hor_1 ; RTc_hor{n, 1}];
    RT_c_hor_2=[RT_c_hor_2 ; RTc_hor{n, 2}];
    RT_c_hor_3=[RT_c_hor_3 ; RTc_hor{n, 3}];
    RT_c_hor_4=[RT_c_hor_4 ; RTc_hor{n, 4}];
end
meanRTc_hor=[mean(RT_c_hor_1) ; mean(RT_c_hor_2); mean(RT_c_hor_3); mean(RT_c_hor_4)];

RT_i_hor_1=[];
RT_i_hor_2=[];
RT_i_hor_3=[];
RT_i_hor_4=[];
for n=1:NumP
    RT_i_hor_1=[RT_i_hor_1 ; RTi_hor{n, 1}];
    RT_i_hor_2=[RT_i_hor_2 ; RTi_hor{n, 2}];
    RT_i_hor_3=[RT_i_hor_3 ; RTi_hor{n, 3}];
    RT_i_hor_4=[RT_i_hor_4 ; RTi_hor{n, 4}];
end
meanRTi_hor=[mean(RT_i_hor_1) ; mean(RT_i_hor_2); mean(RT_i_hor_3); mean(RT_i_hor_4)];

RT_c_ver_1=[];
RT_c_ver_2=[];
RT_c_ver_3=[];
RT_c_ver_4=[];
for n=1:NumP
    RT_c_ver_1=[RT_c_ver_1 ; RTc_ver{n, 1}];
    RT_c_ver_2=[RT_c_ver_2 ; RTc_ver{n, 2}];
    RT_c_ver_3=[RT_c_ver_3 ; RTc_ver{n, 3}];
    RT_c_ver_4=[RT_c_ver_4 ; RTc_ver{n, 4}];
end
meanRTc_ver=[mean(RT_c_ver_1) ; mean(RT_c_ver_2); mean(RT_c_ver_3); mean(RT_c_ver_4)];

RT_i_ver_1=[];
RT_i_ver_2=[];
RT_i_ver_3=[];
RT_i_ver_4=[];
for n=1:NumP
    RT_i_ver_1=[RT_i_ver_1 ; RTi_ver{n, 1}];
    RT_i_ver_2=[RT_i_ver_2 ; RTi_ver{n, 2}];
    RT_i_ver_3=[RT_i_ver_3 ; RTi_ver{n, 3}];
    RT_i_ver_4=[RT_i_ver_4 ; RTi_ver{n, 4}];
end
meanRTi_ver=[mean(RT_i_ver_1) ; mean(RT_i_ver_2); mean(RT_i_ver_3); mean(RT_i_ver_4)];

RT_hor=[];
for n=1:4
    RT_hor = [RT_hor ; (meanRTc_hor(n, 1)+meanRTi_hor(n, 1))/2, (meanRTc_hor(n, 2)+meanRTi_hor(n, 2))/2, (meanRTc_hor(n, 3)+meanRTi_hor(n, 3))/2, (meanRTc_hor(n, 4)+meanRTi_hor(n, 4))/2, ...
        (meanRTc_hor(n, 5)+meanRTi_hor(n, 5))/2];
end

CE_hor_1=[];
CE_hor_2=[];
CE_hor_3=[];
CE_hor_4=[];
for n=1:NumP
    CE_hor_1=[CE_hor_1 ; CE_hor{n, 1}];
    CE_hor_2=[CE_hor_2 ; CE_hor{n, 2}];
    CE_hor_3=[CE_hor_3 ; CE_hor{n, 3}];
    CE_hor_4=[CE_hor_4 ; CE_hor{n, 4}];
end

meanCE_hor=[mean(CE_hor_1) ; mean(CE_hor_2); mean(CE_hor_3); mean(CE_hor_4)];

RT_ver=[];
for n=1:4
    RT_ver = [RT_ver ; (meanRTc_ver(n, 1)+meanRTi_ver(n, 1))/2, (meanRTc_ver(n, 2)+meanRTi_ver(n, 2))/2, (meanRTc_ver(n, 3)+meanRTi_ver(n, 3))/2, (meanRTc_ver(n, 4)+meanRTi_ver(n, 4))/2, ...
        (meanRTc_ver(n, 5)+meanRTi_ver(n, 5))/2];
end

CE_ver_1=[];
CE_ver_2=[];
CE_ver_3=[];
CE_ver_4=[];
for n=1:NumP
    CE_ver_1=[CE_ver_1 ; CE_ver{n, 1}];
    CE_ver_2=[CE_ver_2 ; CE_ver{n, 2}];
    CE_ver_3=[CE_ver_3 ; CE_ver{n, 3}];
    CE_ver_4=[CE_ver_4 ; CE_ver{n, 4}];
end

meanCE_ver=[mean(CE_ver_1) ; mean(CE_ver_2); mean(CE_ver_3); mean(CE_ver_4)];


% Plot the delta plot of each SOA 
figure('Position', [300 300 1200 700])
t= tiledlayout(1,2);
nexttile
g1 = plot(RT_hor(1, :), meanCE_hor(1, :));
hold on 
g2 = plot(RT_hor(2, :), meanCE_hor(2, :));
hold on 
g3 = plot(RT_hor(3, :), meanCE_hor(3, :));
hold on 
g4 = plot(RT_hor(4, :), meanCE_hor(4, :));
set(g1, 'LineStyle', '-', 'Marker', 'o', 'Color', [0 0 0], 'LineWidth', 1.5);
set(g2, 'LineStyle', '-', 'Marker', 'square', 'Color', [0 0 0], 'LineWidth', 1.5);
set(g3, 'LineStyle', '-', 'Marker', 'diamond', 'Color', [0 0 0], 'LineWidth', 1.5);
set(g4, 'LineStyle', '-', 'Marker', '^', 'Color', [0 0 0], 'LineWidth', 1.5);
xlim([350 650]);
ylim([0 120]);
xlabel('RT (ms)')
ylabel('Congruency effects (ms)')
legend('+100 ms', '0 ms', '-50 ms', '-100 ms')
hold off

nexttile
g1 = plot(RT_ver(1, :), meanCE_ver(1, :));
hold on 
g2 = plot(RT_ver(2, :), meanCE_ver(2, :));
hold on 
g3 = plot(RT_ver(3, :), meanCE_ver(3, :));
hold on 
g4 = plot(RT_ver(4, :), meanCE_ver(4, :));
set(g1, 'LineStyle', '-', 'Marker', 'o', 'Color', [0 0 0], 'LineWidth', 1.5);
set(g2, 'LineStyle', '-', 'Marker', 'square', 'Color', [0 0 0], 'LineWidth', 1.5);
set(g3, 'LineStyle', '-', 'Marker', 'diamond', 'Color', [0 0 0], 'LineWidth', 1.5);
set(g4, 'LineStyle', '-', 'Marker', '^', 'Color', [0 0 0], 'LineWidth', 1.5);
xlim([350 650]);
ylim([0 120]);
xlabel('RT (ms)')
ylabel('Congruency effects (ms)')
legend('+100 ms', '0 ms', '-50 ms', '-100 ms')
hold off

exportgraphics(t,'spatialStroopDP_orientation.pdf')

%%
clear
% Empirical flanker data
warning off
flanker=readtable('FlankerSOA.xlsx');
flanker=table2array(flanker);

NumP= 40;

Ppt = cell(NumP, 1);
for k=1:NumP
    Ppt{k} = flanker(flanker(:,1)==k,:);
end

C_hor = cell(NumP, 1);
I_hor = cell(NumP, 1);
C_ver = cell(NumP, 1);
I_ver = cell(NumP, 1);
for n=1:NumP
    for k=1:size(Ppt{n},1)
        if Ppt{n}(k,7)==1 && Ppt{n}(k,3) <=2
            C_hor{n} = [C_hor{n} ; Ppt{n}(k,1) Ppt{n}(k,2) Ppt{n}(k,5) Ppt{n}(k,6), Ppt{n}(k,8)]; 
        elseif Ppt{n}(k,7)==0 && Ppt{n}(k,3) <=2
            I_hor{n} = [I_hor{n} ; Ppt{n}(k,1) Ppt{n}(k,2) Ppt{n}(k,5) Ppt{n}(k,6), Ppt{n}(k,8)];
        elseif Ppt{n}(k,7)==1 && Ppt{n}(k,3) >2
            C_ver{n} = [C_ver{n} ; Ppt{n}(k,1) Ppt{n}(k,2) Ppt{n}(k,5) Ppt{n}(k,6), Ppt{n}(k,8)]; 
        elseif Ppt{n}(k,7)==0 && Ppt{n}(k,3) >2
            I_ver{n} = [I_ver{n} ; Ppt{n}(k,1) Ppt{n}(k,2) Ppt{n}(k,5) Ppt{n}(k,6), Ppt{n}(k,8)];
        end
    end
end

C_hor = cell2mat(C_hor);
de3=C_hor(:,3)<200;
C_hor(de3,:) = [];
de4=C_hor(:,3)>1500;
C_hor(de4,:) = [];
de2=C_hor(:,4)==0;
C_hor(de2,:) = [];

I_hor = cell2mat(I_hor);
de7=I_hor(:,3)<200;
I_hor(de7,:) = [];
de8=I_hor(:,3)>1500;
I_hor(de8,:) = [];
de6=I_hor(:,4)==0;
I_hor(de6,:) = [];

C_ver = cell2mat(C_ver);
de3=C_ver(:,3)<200;
C_ver(de3,:) = [];
de4=C_ver(:,3)>1500;
C_ver(de4,:) = [];
de2=C_ver(:,4)==0;
C_ver(de2,:) = [];

I_ver = cell2mat(I_ver);
de7=I_ver(:,3)<200;
I_ver(de7,:) = [];
de8=I_ver(:,3)>1500;
I_ver(de8,:) = [];
de6=I_ver(:,4)==0;
I_ver(de6,:) = [];

C1=cell(NumP, 1);
C2=cell(NumP, 1);
for n=1:NumP
    k=C_hor(:, 1)==n;
    C1{n} = C_hor(k, :);
    C2{n} = [C1{n}(:, 3), C1{n}(:, 5)];
end
C_hor = C2;

I1=cell(NumP, 1);
I2=cell(NumP, 1);
for n=1:NumP
    l=I_hor(:, 1)==n;
    I1{n} = I_hor(l, :);
    I2{n} = [I1{n}(:, 3), I1{n}(:, 5)];
end
I_hor=I2;

C1=cell(NumP, 1);
C2=cell(NumP, 1);
for n=1:NumP
    k=C_ver(:, 1)==n;
    C1{n} = C_ver(k, :);
    C2{n} = [C1{n}(:, 3), C1{n}(:, 5)];
end
C_ver = C2;

I1=cell(NumP, 1);
I2=cell(NumP, 1);
for n=1:NumP
    l=I_ver(:, 1)==n;
    I1{n} = I_ver(l, :);
    I2{n} = [I1{n}(:, 3), I1{n}(:, 5)];
end
I_ver=I2;

% We have different SOAs in the spatialStroop task, 100s (target shown before
% distractors), 0, -50, -100s separated by the orientation
% Five quantiles for each condition

RTc_hor = cell(n, 4);
for n=1:NumP
    RTc_hor{n, 1} = C_hor{n}(C_hor{n}(:, 2)==100,1);
    RTc_hor{n, 1} = [prctile(RTc_hor{n, 1}(:,1),10), prctile(RTc_hor{n, 1}(:,1),30), prctile(RTc_hor{n, 1}(:,1),50), prctile(RTc_hor{n, 1}(:,1),70), ...
        prctile(RTc_hor{n, 1}(:,1),90)];
    RTc_hor{n, 2} = C_hor{n}(C_hor{n}(:, 2)==0,1);
    RTc_hor{n, 2} = [prctile(RTc_hor{n, 2}(:,1),10), prctile(RTc_hor{n, 2}(:,1),30), prctile(RTc_hor{n, 2}(:,1),50), prctile(RTc_hor{n, 2}(:,1),70), ...
        prctile(RTc_hor{n, 2}(:,1),90)];
    RTc_hor{n, 3} = C_hor{n}(C_hor{n}(:, 2)==-50,1);
    RTc_hor{n, 3} = [prctile(RTc_hor{n, 3}(:,1),10), prctile(RTc_hor{n, 3}(:,1),30), prctile(RTc_hor{n, 3}(:,1),50), prctile(RTc_hor{n, 3}(:,1),70), ...
        prctile(RTc_hor{n, 3}(:,1),90)];
    RTc_hor{n, 4} = C_hor{n}(C_hor{n}(:, 2)==-100,1);
    RTc_hor{n, 4} = [prctile(RTc_hor{n, 4}(:,1),10), prctile(RTc_hor{n, 4}(:,1),30), prctile(RTc_hor{n, 4}(:,1),50), prctile(RTc_hor{n, 4}(:,1),70), ...
        prctile(RTc_hor{n, 4}(:,1),90)];
end

RTi_hor = cell(n, 4);
for n=1:NumP
    RTi_hor{n, 1} = I_hor{n}(I_hor{n}(:, 2)==100,1);
    RTi_hor{n, 1} = [prctile(RTi_hor{n, 1}(:,1),10), prctile(RTi_hor{n, 1}(:,1),30), prctile(RTi_hor{n, 1}(:,1),50), prctile(RTi_hor{n, 1}(:,1),70), ...
        prctile(RTi_hor{n, 1}(:,1),90)];
    RTi_hor{n, 2} = I_hor{n}(I_hor{n}(:, 2)==0,1);
    RTi_hor{n, 2} = [prctile(RTi_hor{n, 2}(:,1),10), prctile(RTi_hor{n, 2}(:,1),30), prctile(RTi_hor{n, 2}(:,1),50), prctile(RTi_hor{n, 2}(:,1),70), ...
        prctile(RTi_hor{n, 2}(:,1),90)];
    RTi_hor{n, 3} = I_hor{n}(I_hor{n}(:, 2)==-50,1);
    RTi_hor{n, 3} = [prctile(RTi_hor{n, 3}(:,1),10), prctile(RTi_hor{n, 3}(:,1),30), prctile(RTi_hor{n, 3}(:,1),50), prctile(RTi_hor{n, 3}(:,1),70), ...
        prctile(RTi_hor{n, 3}(:,1),90)];
    RTi_hor{n, 4} = I_hor{n}(I_hor{n}(:, 2)==-100,1);
    RTi_hor{n, 4} = [prctile(RTi_hor{n, 4}(:,1),10), prctile(RTi_hor{n, 4}(:,1),30), prctile(RTi_hor{n, 4}(:,1),50), prctile(RTi_hor{n, 4}(:,1),70), ...
        prctile(RTi_hor{n, 4}(:,1),90)];
end

RTc_ver = cell(n, 4);
for n=1:NumP
    RTc_ver{n, 1} = C_ver{n}(C_ver{n}(:, 2)==100,1);
    RTc_ver{n, 1} = [prctile(RTc_ver{n, 1}(:,1),10), prctile(RTc_ver{n, 1}(:,1),30), prctile(RTc_ver{n, 1}(:,1),50), prctile(RTc_ver{n, 1}(:,1),70), ...
        prctile(RTc_ver{n, 1}(:,1),90)];
    RTc_ver{n, 2} = C_ver{n}(C_ver{n}(:, 2)==0,1);
    RTc_ver{n, 2} = [prctile(RTc_ver{n, 2}(:,1),10), prctile(RTc_ver{n, 2}(:,1),30), prctile(RTc_ver{n, 2}(:,1),50), prctile(RTc_ver{n, 2}(:,1),70), ...
        prctile(RTc_ver{n, 2}(:,1),90)];
    RTc_ver{n, 3} = C_ver{n}(C_ver{n}(:, 2)==-50,1);
    RTc_ver{n, 3} = [prctile(RTc_ver{n, 3}(:,1),10), prctile(RTc_ver{n, 3}(:,1),30), prctile(RTc_ver{n, 3}(:,1),50), prctile(RTc_ver{n, 3}(:,1),70), ...
        prctile(RTc_ver{n, 3}(:,1),90)];
    RTc_ver{n, 4} = C_ver{n}(C_ver{n}(:, 2)==-100,1);
    RTc_ver{n, 4} = [prctile(RTc_ver{n, 4}(:,1),10), prctile(RTc_ver{n, 4}(:,1),30), prctile(RTc_ver{n, 4}(:,1),50), prctile(RTc_ver{n, 4}(:,1),70), ...
        prctile(RTc_ver{n, 4}(:,1),90)];
end

RTi_ver = cell(n, 4);
for n=1:NumP
    RTi_ver{n, 1} = I_ver{n}(I_ver{n}(:, 2)==100,1);
    RTi_ver{n, 1} = [prctile(RTi_ver{n, 1}(:,1),10), prctile(RTi_ver{n, 1}(:,1),30), prctile(RTi_ver{n, 1}(:,1),50), prctile(RTi_ver{n, 1}(:,1),70), ...
        prctile(RTi_ver{n, 1}(:,1),90)];
    RTi_ver{n, 2} = I_ver{n}(I_ver{n}(:, 2)==0,1);
    RTi_ver{n, 2} = [prctile(RTi_ver{n, 2}(:,1),10), prctile(RTi_ver{n, 2}(:,1),30), prctile(RTi_ver{n, 2}(:,1),50), prctile(RTi_ver{n, 2}(:,1),70), ...
        prctile(RTi_ver{n, 2}(:,1),90)];
    RTi_ver{n, 3} = I_ver{n}(I_ver{n}(:, 2)==-50,1);
    RTi_ver{n, 3} = [prctile(RTi_ver{n, 3}(:,1),10), prctile(RTi_ver{n, 3}(:,1),30), prctile(RTi_ver{n, 3}(:,1),50), prctile(RTi_ver{n, 3}(:,1),70), ...
        prctile(RTi_ver{n, 3}(:,1),90)];
    RTi_ver{n, 4} = I_ver{n}(I_ver{n}(:, 2)==-100,1);
    RTi_ver{n, 4} = [prctile(RTi_ver{n, 4}(:,1),10), prctile(RTi_ver{n, 4}(:,1),30), prctile(RTi_ver{n, 4}(:,1),50), prctile(RTi_ver{n, 4}(:,1),70), ...
        prctile(RTi_ver{n, 4}(:,1),90)];
end

CE_hor=cell(NumP, 4);
for n=1:NumP
    CE_hor{n, 1} = [RTi_hor{n, 1}(1) - RTc_hor{n, 1}(1), RTi_hor{n, 1}(2) - RTc_hor{n, 1}(2), RTi_hor{n, 1}(3) - RTc_hor{n, 1}(3), RTi_hor{n, 1}(4) - RTc_hor{n, 1}(4), ...
        RTi_hor{n, 1}(5) - RTc_hor{n, 1}(5)];
    CE_hor{n, 2} = [RTi_hor{n, 2}(1) - RTc_hor{n, 2}(1), RTi_hor{n, 2}(2) - RTc_hor{n, 2}(2), RTi_hor{n, 2}(3) - RTc_hor{n, 2}(3), RTi_hor{n, 2}(4) - RTc_hor{n, 2}(4), ...
        RTi_hor{n, 2}(5) - RTc_hor{n, 2}(5)];
    CE_hor{n, 3} = [RTi_hor{n, 3}(1) - RTc_hor{n, 3}(1), RTi_hor{n, 3}(2) - RTc_hor{n, 3}(2), RTi_hor{n, 3}(3) - RTc_hor{n, 3}(3), RTi_hor{n, 3}(4) - RTc_hor{n, 3}(4), ...
        RTi_hor{n, 3}(5) - RTc_hor{n, 3}(5)];
    CE_hor{n, 4} = [RTi_hor{n, 4}(1) - RTc_hor{n, 4}(1), RTi_hor{n, 4}(2) - RTc_hor{n, 4}(2), RTi_hor{n, 4}(3) - RTc_hor{n, 4}(3), RTi_hor{n, 4}(4) - RTc_hor{n, 4}(4), ...
        RTi_hor{n, 4}(5) - RTc_hor{n, 4}(5)];
end

CE_ver=cell(NumP, 4);
for n=1:NumP
    CE_ver{n, 1} = [RTi_ver{n, 1}(1) - RTc_ver{n, 1}(1), RTi_ver{n, 1}(2) - RTc_ver{n, 1}(2), RTi_ver{n, 1}(3) - RTc_ver{n, 1}(3), RTi_ver{n, 1}(4) - RTc_ver{n, 1}(4), ...
        RTi_ver{n, 1}(5) - RTc_ver{n, 1}(5)];
    CE_ver{n, 2} = [RTi_ver{n, 2}(1) - RTc_ver{n, 2}(1), RTi_ver{n, 2}(2) - RTc_ver{n, 2}(2), RTi_ver{n, 2}(3) - RTc_ver{n, 2}(3), RTi_ver{n, 2}(4) - RTc_ver{n, 2}(4), ...
        RTi_ver{n, 2}(5) - RTc_ver{n, 2}(5)];
    CE_ver{n, 3} = [RTi_ver{n, 3}(1) - RTc_ver{n, 3}(1), RTi_ver{n, 3}(2) - RTc_ver{n, 3}(2), RTi_ver{n, 3}(3) - RTc_ver{n, 3}(3), RTi_ver{n, 3}(4) - RTc_ver{n, 3}(4), ...
        RTi_ver{n, 3}(5) - RTc_ver{n, 3}(5)];
    CE_ver{n, 4} = [RTi_ver{n, 4}(1) - RTc_ver{n, 4}(1), RTi_ver{n, 4}(2) - RTc_ver{n, 4}(2), RTi_ver{n, 4}(3) - RTc_ver{n, 4}(3), RTi_ver{n, 4}(4) - RTc_ver{n, 4}(4), ...
        RTi_ver{n, 4}(5) - RTc_ver{n, 4}(5)];
end

RT_c_hor_1=[];
RT_c_hor_2=[];
RT_c_hor_3=[];
RT_c_hor_4=[];
for n=1:NumP
    RT_c_hor_1=[RT_c_hor_1 ; RTc_hor{n, 1}];
    RT_c_hor_2=[RT_c_hor_2 ; RTc_hor{n, 2}];
    RT_c_hor_3=[RT_c_hor_3 ; RTc_hor{n, 3}];
    RT_c_hor_4=[RT_c_hor_4 ; RTc_hor{n, 4}];
end
meanRTc_hor=[mean(RT_c_hor_1) ; mean(RT_c_hor_2); mean(RT_c_hor_3); mean(RT_c_hor_4)];

RT_i_hor_1=[];
RT_i_hor_2=[];
RT_i_hor_3=[];
RT_i_hor_4=[];
for n=1:NumP
    RT_i_hor_1=[RT_i_hor_1 ; RTi_hor{n, 1}];
    RT_i_hor_2=[RT_i_hor_2 ; RTi_hor{n, 2}];
    RT_i_hor_3=[RT_i_hor_3 ; RTi_hor{n, 3}];
    RT_i_hor_4=[RT_i_hor_4 ; RTi_hor{n, 4}];
end
meanRTi_hor=[mean(RT_i_hor_1) ; mean(RT_i_hor_2); mean(RT_i_hor_3); mean(RT_i_hor_4)];

RT_c_ver_1=[];
RT_c_ver_2=[];
RT_c_ver_3=[];
RT_c_ver_4=[];
for n=1:NumP
    RT_c_ver_1=[RT_c_ver_1 ; RTc_ver{n, 1}];
    RT_c_ver_2=[RT_c_ver_2 ; RTc_ver{n, 2}];
    RT_c_ver_3=[RT_c_ver_3 ; RTc_ver{n, 3}];
    RT_c_ver_4=[RT_c_ver_4 ; RTc_ver{n, 4}];
end
meanRTc_ver=[mean(RT_c_ver_1) ; mean(RT_c_ver_2); mean(RT_c_ver_3); mean(RT_c_ver_4)];

RT_i_ver_1=[];
RT_i_ver_2=[];
RT_i_ver_3=[];
RT_i_ver_4=[];
for n=1:NumP
    RT_i_ver_1=[RT_i_ver_1 ; RTi_ver{n, 1}];
    RT_i_ver_2=[RT_i_ver_2 ; RTi_ver{n, 2}];
    RT_i_ver_3=[RT_i_ver_3 ; RTi_ver{n, 3}];
    RT_i_ver_4=[RT_i_ver_4 ; RTi_ver{n, 4}];
end
meanRTi_ver=[mean(RT_i_ver_1) ; mean(RT_i_ver_2); mean(RT_i_ver_3); mean(RT_i_ver_4)];

RT_hor=[];
for n=1:4
    RT_hor = [RT_hor ; (meanRTc_hor(n, 1)+meanRTi_hor(n, 1))/2, (meanRTc_hor(n, 2)+meanRTi_hor(n, 2))/2, (meanRTc_hor(n, 3)+meanRTi_hor(n, 3))/2, (meanRTc_hor(n, 4)+meanRTi_hor(n, 4))/2, ...
        (meanRTc_hor(n, 5)+meanRTi_hor(n, 5))/2];
end

CE_hor_1=[];
CE_hor_2=[];
CE_hor_3=[];
CE_hor_4=[];
for n=1:NumP
    CE_hor_1=[CE_hor_1 ; CE_hor{n, 1}];
    CE_hor_2=[CE_hor_2 ; CE_hor{n, 2}];
    CE_hor_3=[CE_hor_3 ; CE_hor{n, 3}];
    CE_hor_4=[CE_hor_4 ; CE_hor{n, 4}];
end

meanCE_hor=[mean(CE_hor_1) ; mean(CE_hor_2); mean(CE_hor_3); mean(CE_hor_4)];

RT_ver=[];
for n=1:4
    RT_ver = [RT_ver ; (meanRTc_ver(n, 1)+meanRTi_ver(n, 1))/2, (meanRTc_ver(n, 2)+meanRTi_ver(n, 2))/2, (meanRTc_ver(n, 3)+meanRTi_ver(n, 3))/2, (meanRTc_ver(n, 4)+meanRTi_ver(n, 4))/2, ...
        (meanRTc_ver(n, 5)+meanRTi_ver(n, 5))/2];
end

CE_ver_1=[];
CE_ver_2=[];
CE_ver_3=[];
CE_ver_4=[];
for n=1:NumP
    CE_ver_1=[CE_ver_1 ; CE_ver{n, 1}];
    CE_ver_2=[CE_ver_2 ; CE_ver{n, 2}];
    CE_ver_3=[CE_ver_3 ; CE_ver{n, 3}];
    CE_ver_4=[CE_ver_4 ; CE_ver{n, 4}];
end

meanCE_ver=[mean(CE_ver_1) ; mean(CE_ver_2); mean(CE_ver_3); mean(CE_ver_4)];


% Plot the delta plot of each SOA 
figure('Position', [300 300 1200 700])
t= tiledlayout(1,2);
nexttile
g1 = plot(RT_hor(1, :), meanCE_hor(1, :));
hold on 
g2 = plot(RT_hor(2, :), meanCE_hor(2, :));
hold on 
g3 = plot(RT_hor(3, :), meanCE_hor(3, :));
hold on 
g4 = plot(RT_hor(4, :), meanCE_hor(4, :));
set(g1, 'LineStyle', '-', 'Marker', 'o', 'Color', [0 0 0], 'LineWidth', 1.5);
set(g2, 'LineStyle', '-', 'Marker', 'square', 'Color', [0 0 0], 'LineWidth', 1.5);
set(g3, 'LineStyle', '-', 'Marker', 'diamond', 'Color', [0 0 0], 'LineWidth', 1.5);
set(g4, 'LineStyle', '-', 'Marker', '^', 'Color', [0 0 0], 'LineWidth', 1.5);
xlim([350 650]);
ylim([-20 140]);
xlabel('RT (ms)')
ylabel('Congruency effects (ms)')
legend('+100 ms', '0 ms', '-50 ms', '-100 ms')
hold off

nexttile
g1 = plot(RT_ver(1, :), meanCE_ver(1, :));
hold on 
g2 = plot(RT_ver(2, :), meanCE_ver(2, :));
hold on 
g3 = plot(RT_ver(3, :), meanCE_ver(3, :));
hold on 
g4 = plot(RT_ver(4, :), meanCE_ver(4, :));
set(g1, 'LineStyle', '-', 'Marker', 'o', 'Color', [0 0 0], 'LineWidth', 1.5);
set(g2, 'LineStyle', '-', 'Marker', 'square', 'Color', [0 0 0], 'LineWidth', 1.5);
set(g3, 'LineStyle', '-', 'Marker', 'diamond', 'Color', [0 0 0], 'LineWidth', 1.5);
set(g4, 'LineStyle', '-', 'Marker', '^', 'Color', [0 0 0], 'LineWidth', 1.5);
xlim([350 650]);
ylim([-20 140]);
xlabel('RT (ms)')
ylabel('Congruency effects (ms)')
legend('+100 ms', '0 ms', '-50 ms', '-100 ms')
hold off

exportgraphics(t,'flankerDP_orientation.pdf')
