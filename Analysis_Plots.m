

% Author: Davide Ciccarese
% Date of creation: 22/06/2022
% Last modification: 19/05/2023
% License: CC BY 4.0 (https://creativecommons.org/licenses/by/4.0/)


%% --- LOAD ALL THE DATA ----

clc;
clear all;
close all;

brkn =2; %Intact == 1

if brkn ==1 %Intact
    rootDir = {'/INTACT_CA/tif_XY_05/'...
        '/INTACT_CA/tif_XY_06/'...
        '/INTACT_CA/tif_XY_07/'...
        '/INTACT_CA/tif_XY_08/'};
else
    % Broken diatom
    rootDir = {'/BROKEN_CA/tif_XY_05/'...
        '/BROKEN_CA//tif_XY_06/'...
        '/BROKEN_CA/tif_XY_07/'...
        '/BROKEN_CA/tif_XY_08/'};
end

phisV = 0.65;
Lngt = 9.73;
Hgt = 0.88;
A_area = 5
minutes = 60;
mltomm3 = 10;
tom = 10^-3;
SinkVel = 0.2

EdgeDID = {};
GrowthID = {};
DistEdgCh = {};
AreaCh_v = {};
OxIntCol = {};
OxSigDistParticleZones = {};
OxFract = {};

for k = 1:length(rootDir)
    
    cd(rootDir{1,k});
    load('MultipleROI_ID7.mat');
    
    EdgeDID{1,k} = StructMultipleROI_ID.EdgeDistances;
    GrowthID{1,k} = StructMultipleROI_ID.GrowthTime;
    DistEdgCh{1,k} = StructMultipleROI_ID.DistChannel;
    AreaCh_v{1,k} = StructMultipleROI_ID.AreaChannel;
    OxIntCol{1,k} = StructMultipleROI_ID.DynamicOx;
    OxSigDistParticleZones{1,k} = StructMultipleROI_ID.ParticleOX;
    OxFract{1,k}= StructMultipleROI_ID.DynamicFractOx;
    
end

Aqtime = 1;
limiT = 261*Aqtime*SinkVel;


if brkn ==1
    load '/INTACT_CA/OxygenCalibration2.mat'
else
    load '/BROKEN_CA/OxygenCalibration2.mat'
end

% Long vector
DistEdgChLong = {};
AreaChLong = {};
EdgeDIDLong = {};

for k = 1:length(rootDir)
    for jk = 1:length(DistEdgCh{1,k})
        DistEdgChLong{jk,k} = [DistEdgCh{1,k}{1,jk}{:}].*phisV;
        AreaChLong{jk,k} = [AreaCh_v{1,k}{1,jk}{:}].*phisV;
        EdgeDIDLong{jk,k} = [EdgeDID{1,k}{1,jk}{:}].*phisV;
        
    end
end

growthN = {};
BioCel = {};
BioCelox = {};
BioM = {};
BioOx = {};
BioOutOx = {};

for k = 1:length(rootDir)
    for jk = 1:length(EdgeDID{1,k})
        Ch = jk;
        if isempty(EdgeDID{1,k}{1,Ch})
        else
            [srtEdg,indx] = sort([EdgeDID{1,k}{1,Ch}{1,:}]*phisV);
        end
        for i = 1:length(indx)
            for jj = 1:size(GrowthID{1,k}{:,Ch},1)
                if ((jj>0) && (GrowthID{1,k}{Ch}{jj,indx(i)}<0))
                    growthN{jj,i}= NaN;
                    OxthN{jj,i} = NaN;
                    OxFractN{jj,i} = NaN;
                else
                    OxthN{jj,i} = OxIntCol{1,k}{:,Ch}{jj,indx(i)};
                    OxFractN40{jj,i} = OxFract{1,k}{:,Ch}{jj,indx(i)}{1,1};
                    OxFractN20{jj,i} = OxFract{1,k}{:,Ch}{jj,indx(i)}{1,2};
                    OxFractN10{jj,i} = OxFract{1,k}{:,Ch}{jj,indx(i)}{1,3};
                    OxFractN625{jj,i} = OxFract{1,k}{:,Ch}{jj,indx(i)}{1,4};
                    OxFractN5{jj,i} = OxFract{1,k}{:,Ch}{jj,indx(i)}{1,5};
                    OxFractN1{jj,i} = OxFract{1,k}{:,Ch}{jj,indx(i)}{1,6};
                    growthN{jj,i} = GrowthID{1,k}{:,Ch}{jj,indx(i)};
                end
            end
        end
        
        t = ([1:(size(growthN,1))])*Aqtime*SinkVel;
        for i = 1:length(indx)
            
            if ((sum(([growthN{:,indx(i)}])==0))> 1000)
            else
                yy2 = fillmissing([growthN{:,i}],'nearest');
                yy2_ox = fillmissing([OxthN{:,i}],'nearest');
                yy2_oxFr40 = fillmissing([OxFractN40{:,i}],'nearest');
                yy2_oxFr20 = fillmissing([OxFractN20{:,i}],'nearest');
                yy2_oxFr10 = fillmissing([OxFractN10{:,i}],'nearest');
                yy2_oxFr625 =fillmissing([OxFractN625{:,i}],'nearest');
                yy2_oxFr5 = fillmissing([OxFractN5{:,i}],'nearest');
                yy2_oxFr1 = fillmissing([OxFractN1{:,i}],'nearest');
            end
            BioOx{jk,i} = ([yy2_ox(:)]);
            BioM{jk,i} = ([yy2(:)]);
            BioOx40{jk,i} = ([yy2_oxFr40(:)]);
            BioOx20{jk,i} = ([yy2_oxFr20(:)]);
            BioOx10{jk,i} = ([yy2_oxFr10(:)]);
            BioOx625{jk,i} = ([yy2_oxFr625(:)]);
            BioOx5{jk,i} = ([yy2_oxFr5(:)]);
            BioOx1{jk,i} = ([yy2_oxFr1(:)]);
        end
    end
    BioCel{k} = ([BioM{:,:}])
    BioCelox{k} = ([BioOx{:,:}])
    BioCelox40{k} = ([BioOx40{:,:}])
    BioCelox20{k} = ([BioOx20{:,:}])
    BioCelox10{k} = ([BioOx10{:,:}])
    BioCelox625{k} = ([BioOx625{:,:}])
    BioCelox5{k} = ([BioOx5{:,:}])
    BioCelox1{k} = ([BioOx1{:,:}])
end

BioCeloxFract = {BioCelox40 BioCelox20 BioCelox10 BioCelox625 BioCelox5 BioCelox1};
%% Figure 2,a/e Broken and intact diatoms drives a spatially heterogeneous biomass accumulation

growthN = {};
BioM = {};
BioOx = {};
smothBiom= {};
celldaSub={};

for k = 1:length(rootDir)
    for jk = 1:length(EdgeDID{1,k})
        
        Ch = 1;
        if isempty(EdgeDID{1,k}{Ch})
        else
            [srtEdg,indx] = sort([EdgeDID{1,k}{1,Ch}{1,:}]*phisV);
        end
        for i = 1:length(indx)
            for jj = 1:size(GrowthID{1,k}{:,Ch},1)
                if ((jj>0) && (GrowthID{1,k}{Ch}{jj,indx(i)}<0))
                    growthN{jj,i}= NaN;
                    OxthN{jj,i} = NaN;
                else
                    OxthN{jj,i} = OxIntCol{1,k}{:,Ch}{jj,indx(i)};
                    growthN{jj,i} = GrowthID{1,k}{:,Ch}{jj,indx(i)};
                end
            end
        end
        t = ([1:(size(growthN,1))])*Aqtime*SinkVel;
        inDex = [1:length(indx)];
        indFx = {inDex(srtEdg<300) inDex(srtEdg>300 & srtEdg<500) inDex(srtEdg>500 & srtEdg<700) inDex(srtEdg>700 & srtEdg<900), inDex(srtEdg>900 & srtEdg<1100), inDex(srtEdg>1100 & srtEdg<1300) inDex(srtEdg>1300)};
        zon1 = pi*(1500^2-1300^2).*phisV;
        zon2 = pi*(1300^2-1100^2).*phisV;
        zon3 = pi*(1100^2-900^2).*phisV;
        zon4 = pi*(900^2-700^2).*phisV;
        zon5 = pi*(700^2-500^2).*phisV;
        zon6 = pi*(500^2-300^2).*phisV;
        zon7 = pi*(300^2).*phisV;
        indaRad = {zon1 zon2 zon3 zon4 zon5 zon6 zon7};
        for indaFux= 1:length(indFx)
            for i = 1:length(indx)
                yy2 = fillmissing([growthN{:,i}],'nearest');
                yy2_ox = fillmissing([OxthN{:,i}],'nearest');
                BioOx{jk,i} = ([yy2_ox(:)]);
                BioM{jk,i} = ([yy2(:)]);
            end
            yy2_b = smooth(sum([BioM{jk,[indFx{indaFux}]}]'),0.2,'rlowess');
            celldaSub{indaFux} = yy2_b;
        end
        smothBiom{jk,k} = celldaSub
    end
end

matsmothBiom = {};
maBdir = zeros([length(t)],1);
limitT = 261*Aqtime*SinkVel;
for inj= 1:length(indFx)
    for ik = 1:length(rootDir)
        for jki = 1:length(t)
            matsmothBiom{inj,ik} = smothBiom{1, ik}{inj};
        end
    end
end

for inj= 1:length(indFx)
    for jki = 1:length(t)
        for ik= 1:length(rootDir)
            if (matsmothBiom{inj,ik})==0
            else
                maBdirFour(ik) = (matsmothBiom{inj,ik}(jki))*phisV;
            end
        end
        maBdir(jki,1) = mean(maBdirFour(:)')
    end
    plot(t,(maBdir),'linewidth',2,'Color',[0 0 0]+0.15*inj/1.5);
    hold on
    set(gca,'View',[90 90])
    hold on
    set(gca,'LineWidth',2)
    set(gca,'FontSize',20)
    ylabel('Biomass (\mum^2)','FontSize',20)
    xlabel('Depth (m)','FontSize',20)
    xlim([0 limitT])
    if brkn ==2;
        ylim([0 250000])
    else
    end
    hold on
    legend('<300 \mum','300-500 \mum', '500-700 \mum', '700-900 \mum', '900-1100 \mum', '1100-1200 \mum', '>1200 \mum')
end
%% Figure 2 b/f | Broken and intact diatoms drives a spatially heterogeneous biomass accumulation

k = 1:length(rootDir)
if (size(DistEdgChLong,1))==1
    clr = [0 0 0];
    binSize = 100;
    TrimVal = 1500;
    distVec = [0:binSize:TrimVal];
    distTrim = 1500;
    myfac     = 4;
    Adist = [((DistEdgChLong{1,k(1)})) ((DistEdgChLong{1,k(2)})) ((DistEdgChLong{1,k(3)})) ((DistEdgChLong{1,k(4)}))];
    A = [((AreaChLong{1,k(1)})) ((AreaChLong{1,k(2)})) ((AreaChLong{1,k(3)})) ((AreaChLong{1,k(4)}))];
    DistA = {};Acells = {};
    for i = 1:length(0:binSize:TrimVal)
        DistA{i} = Adist(Adist >= distVec(i) & Adist <= distVec(i)+binSize);
        Acells{i} = A(Adist >= distVec(i) & Adist <= distVec(i)+binSize);
        hold on
        r = ((mean(Acells{i})))
        x = nanmedian(DistA{i});
        y = r;
        d = r*2;
        px = x-r;
        py = -r;
        h = rectangle('Position',[px py d d],'Curvature',[1,1]);
        set(h,'FaceColor',[clr,0.1],'EdgeColor', [0,0,0.1])
        daspect([1,1,1])
        ylabel('Mean colony biomass (\mum^2)','FontSize',10)
        xlabel('Distance to the edge (\mum)','FontSize',15)
        set(gca,'LineWidth',2)
        set(gca,'LineWidth',2, 'FontSize',15)
        yticklocs = get(gca, 'YTick');
        yticklabels(cellstr(num2str(yticklocs'*myfac)));
    end
end
%% Figure 2 c/g | Broken and intact diatoms drives a spatially heterogeneous biomass accumulation

k = size(rootDir,2)
sbpl =size(EdgeDID{1,k},2);
for k = 1:length(rootDir)
    for jk = 1:length(EdgeDID{1,k})
        subplot(1,sbpl,1)
        scatter((([DistEdgCh{1,k}{1,1}{:}]).*phisV),(([AreaCh_v{1,k}{1,1}{:}]).*phisV),(([AreaCh_v{1,k}{1,1}{:}]).*phisV),'MarkerFaceColor','k','MarkerEdgeColor','r',...
            'MarkerFaceAlpha',.2,'MarkerEdgeAlpha',.2);
        if brkn==1
            ylim([0 30000])
            xlim([0 1500])
        else
            ylim([0 3000])
            xlim([0 1500])
        end
        ylabel('Colony size (\mum^2)','FontSize',10, 'FontWeight','bold')
        xlabel('Distance to particle edge (\mum)','FontSize',10, 'FontWeight','bold')
        set(gca,'LineWidth',2, 'FontSize',15)
    end
end
%% Figure 2 d/h | Broken and intact diatoms drives a spatially heterogeneous biomass accumulation

k = 1:length(rootDir)
A = [((DistEdgChLong{1,k(1)})) ((DistEdgChLong{1,k(2)})) ((DistEdgChLong{1,k(3)})) ((DistEdgChLong{1,k(4)}))];
A = A(A>10);
A = A(A<1500);
[f,xi] = ksdensity((A))
yn =f./sum(f)
plot(xi,yn,'k','LineWidth',5)
hold on
a = area(xi,yn,'FaceColor',[0 0 0])
a.FaceAlpha = 0.2;
xlabel({'Distance to the edge (\mum)'},'FontSize',5)
set(gca,'LineWidth',2)
set(gca,'FontSize',10)
legend('Biomass (\mum^2)','Location','northeast')
set(gca,'LineWidth',2, 'FontSize',15)
xlim([0 1500])
ylim([0 0.03])
%% Figure 3 a/b | DOC depth profiles of particles seeded with broken diatoms and intact diatoms. a) DOC

cd '/TOC_BROKEN/'
load 'TOC.mat'
figure,
subplot(1,2,1)
xt = ([1:(size(growthN,1))])*Aqtime*SinkVel;
plot(1:4.8:max([xt]),TOC(1:3,:),'Color',[0, 1, 1, 0.2],'LineWidth',0.1)
hold on
plot(1:4.8:max([xt]),nanmedian(TOC(1:3,:)),'-c','LineWidth',3)
xlim([1 max(xt)-2])
ylim([1 max(nanmedian(TOC(1:3,:)))+2])
ylabel('TOC (mgC/L)')
xlabel('Depth (m)','FontSize',20)
set(gca,'LineWidth',2, 'FontSize',15)
set(gca,'View',[90 90])
title('Broken diatom')

subplot(1,2,2)
cd '/TOC_INTACT/'

load 'TOC.mat'
xt = ([1:12]*24)*0.2;

load 'TOC.mat'

subplot(1,2,2)
cd '/TOC_INTACT/'
load 'TOC.mat'
xt = ([1:12]*24)*0.2;
plot(xt,TOC(:,1:3)','Color',[0, 1, 1, 0.2],'LineWidth',0.1)
hold on
plot(xt,nanmedian(TOC(:,1:3)'),'-c','LineWidth',3)
xlim([min(xt) 55])
ylim([0 30])
ylabel('TOC (mgC/L)')
xlabel('Depth (m)','FontSize',20)
set(gca,'LineWidth',2, 'FontSize',15)
title('Intact diatom')
set(gca,'View',[90 90])
%% Figure 4 a,b/g,h | Oxygen, biomass accumulation and nitrate and nitrite depth profiles.

t = ([1:(size(growthN,1))])*Aqtime*SinkVel
subplot(1,2,1)
sumOx = {};
for k = 1:length(rootDir)
    MatBioOx = ([BioCelox{k}]');
    MatBioOx(MatBioOx==0) = NaN;
    yc = smooth(nanmean(MatBioOx),0.2,'rlowess');
    if brkn ==1
        h = plot(t,yc-min(yc),'-b')
    else
        h = plot(t,yc,'-b')
    end
    hold on
    h.Color(4) = 0.25;
    sumOx{k} = nanmean(MatBioOx);
end

matBiomOx = zeros([length(t),length(rootDir)]);
for ik = 1:length(rootDir)
    for jki = 1:length(t)
        matBiomOx(jki,ik) = ([sumOx{ik}(jki)]);
    end
end

hold on
yc = smooth((nanmean(matBiomOx')),0.2,'rlowess');

if brkn ==1
    h2 = plot(t,yc-min(yc),'-b','LineWidth',2);
else
    h2 = plot(t,yc,'-b','LineWidth',2);
    
end

h2 = gca;
set(h2, 'YDir', 'reverse');
set(gca,'View',[90 90])
hold on
set(gca,'LineWidth',2)
set(gca,'FontSize',20)
xt = get(gca, 'XTick');
ylabel('Oxygen (\muM)','FontSize',20)
xlabel('Depth (m)','FontSize',20)
xlim([0 max(t)])

if brkn ==1
    ylim([0 250])
else
    ylim([100 300])
end

subplot(1,2,2)
sumBiom = {};
for k = 1:length(rootDir)
    yy2 = smooth(t,sum([BioCel{k}]'),0.2,'rlowess');
    h = plot(t,yy2.*phisV,'-k')
    h.Color(4) = 0.25;
    hold on
    sumBiom{k} = sum([BioCel{k}]')*phisV;
end

matBiom = zeros([length(t),length(rootDir)]);
for ik = 1:length(rootDir)
    for jki = 1:length(t)
        matBiom(jki,ik) = ([sumBiom{ik}(jki)]);
    end
end

hold on
yy2 = smooth(t,nanmean(matBiom'),0.2,'rlowess');
hold on
h = plot(t,yy2,'-k','LineWidth',2);
set(gca,'View',[90 90])
hold on
set(gca,'LineWidth',2)
set(gca,'FontSize',20)
xt = get(gca, 'XTick');
ylabel('Biomass (\mum^2)','FontSize',20)
xlabel('Depth (m)','FontSize',20)
xlim([0 max(t)])
ylim([0 max(max(yy2'))+100000])
%% Figure 4 c/i | Oxygen, biomass accumulation and nitrate and nitrite depth profiles.

BioCeloxFract = {BioCelox40 BioCelox20 BioCelox10 BioCelox625 BioCelox5 BioCelox1};

sumBiomOX = {};
sumBiom = {};

CM = magma(6);
for jk =  1:length(BioCeloxFract)
    for k = 1:length(rootDir)
        yy2OX = smooth(t,sum([BioCeloxFract{jk}{k}]')*phisV,0.2,'rlowess');
        yy2 = smooth(t,sum([BioCel{k}]')*phisV,0.2,'rlowess');
        hold on
        sumBiomOX{k} = sum([BioCeloxFract{jk}{k}]')*phisV;
        sumBiom{k} = sum([BioCel{k}]')*phisV;
    end
    
    matBiomoX = zeros([length(t),length(rootDir)]);
    matBiom = zeros([length(t),length(rootDir)]);
    
    for ik = 1:length(rootDir)
        for jki = 1:length(t)
            matBiomoX(jki,ik) = ([sumBiomOX{ik}(jki)]);
            matBiom(jki,ik) = ([sumBiom{ik}(jki)]);
        end
    end
    
    hold on
    yy2Ox = smooth(t,nanmean(matBiomoX'),0.2,'rlowess');
    yy2 = smooth(t,nanmean(matBiom'),0.2,'rlowess');
    hold on
    h = plot(t,yy2Ox./yy2,'-','color',CM(jk,:),'LineWidth',5);
    set(gca,'View',[90 90])
    hold on
    set(gca,'LineWidth',2)
    set(gca,'FontSize',20)
    xt = get(gca, 'XTick');
    ylabel('Biomass (\mum^2)','FontSize',20)
    xlabel('Depth (m)','FontSize',20)
    xlim([0 max(t)])
    Number = {40 20 10 6.25 5 1}
    box on;
    
    if (jk > 1 && brkn==1)
        ylim([0 0.01])
    elseif brkn==2
        ylim([0 0.1])
    end
end
%% Figure 4 d/j | Oxygen, biomass accumulation and nitrate and nitrite depth profiles

% ---- Identify switching points ------
growthN= {};
growthNth = {}
SWdataN1 = {};

for k = 1:length(rootDir)
    for jk = 1:size(GrowthID{1,k},2)
        [srtEdg,indx] = sort([EdgeDID{1,k}{1,jk}{1,:}]);
        for i = 1:length(indx)
            for jj = 1:size(GrowthID{1,k}{:,jk},1)
                if ((jj>0) && (GrowthID{1,k}{jk}{jj,indx(i)}<0))
                    growthN{1,k}{1,jk}{jj,i}= NaN;
                else
                    growthN{1,k}{1,jk}{jj,i} = GrowthID{1,k}{:,jk}{jj,indx(i)};
                end
            end
        end
        
        t = [1:size(GrowthID{1,k}{:,jk},1)]*Aqtime*SinkVel;
        
        for i = 1:length(indx)
            if (((sum((GrowthID{1,k}{:,jk}{jj,indx(i)})))==0)> 1000)
            else
                tmp3 = fillmissing([growthN{1,k}{1,jk}{:,i}],'previous');
                yy3 = smooth(t,tmp3*phisV,0.2,'rlowess');
                if any(yy3 >1000)
                    [ipt,residual]=findchangepts(yy3,'MaxNumChanges',4,'Statistic','linear');
                    growthNth{1,k}{1,jk}{i}=(t(ipt))
                else
                    growthNth{1,k}{1,jk}{i}=NaN;
                end
            end
            SWdataN1{1,k}{i} = [growthNth{1,k}{1,jk}{1,i}];
        end
    end
end

% ---- LONG VECTOR----
swDN1 ={};
swDN2 ={};
swDN3 = {};

for swtch = 1:4
    for dr = 1:length(rootDir)
        for jk = 1:size(GrowthID{1,dr},2)
            [srtEdg,indx] = sort([EdgeDID{1,dr}{1,jk}{1,:}].*phisV);
            for kjk = 1:length(indx)
                if isempty(SWdataN1{1, dr}{1,kjk})
                    swDN1{dr}(swtch,kjk) = NaN;
                elseif (length(SWdataN1{1, dr}{1,kjk}))==1
                    swDN1{dr}(swtch,kjk) = NaN;
                elseif (length(SWdataN1{1, dr}{1,kjk}))==2
                    swDN1{dr}(swtch,kjk) = NaN;
                elseif (length(SWdataN1{1, dr}{1,kjk}))==3
                    swDN1{dr}(swtch,kjk) = NaN;
                else
                    swDN1{dr}(swtch,kjk) = [SWdataN1{1, dr}{1,kjk}(1,swtch)];
                end
            end
        end
    end
end

% Waiting time between 'Diauxic Lag Phase (h)' and 'End diauxic Lag Phase (h)'
figure,
longswDN1 = [swDN1{1,:}]
lag1 = longswDN1(2,:); lag1 = lag1(lag1>0);lag1(~isnan(lag1));
lag2 = longswDN1(3,:); lag2 = lag2(lag2>0);lag2(~isnan(lag2));
al_goodplot(median(lag1)+(lag2-lag1),1,0.5,'c');
hold on
plot(nanmedian(median(lag1)+(lag2-lag1)),'-ok','LineWidth',50)
grid off
set(gca,'LineWidth',2)
set(gca,'FontSize',20)
xt = get(gca, 'XTick');
ylabel('Depth (m)','FontSize',20)
ylim([0 max(t)])

x = (lag2-lag1);
mean(x)
std(x)

cd /Switch_intact/
load('Lag1Lag2.mat')

lag1Intact = lag1;
lag2Intact = lag2;


cd /Switch_broken/
load('Lag1Lag2.mat')

lag1broken = lag1;
lag2broken = lag2;
xintact = (lag2Intact-lag1Intact);
xbroken = (lag2broken-lag1broken);

[h,p] = ttest2(xintact,xbroken);
%% Figure 4 e/k | Oxygen, biomass accumulation and nitrate and nitrite depth profiles.

cd '/NitriteNitrateQuantification/'


load 'NND_Is5_n1.mat'
load 'MD_Is5_n1.mat'
load 'OD_Is5_n1.mat'
load 'OM_Is5_n1.mat'


intercept =[138.71 1244.7 2420.6 3024.4 717.87 1371.9 8876.2 121.68 421.69 706.02 2659.1 768.3];
slope = [5195.2 5185.1 5667.3 5481 9045.2 7337.3 5349.6 9326.1 9330.8 10160 8553.5 9653.5];


vecMD = zeros(3,12);
for i = 1:12
    
    TMD = MD_Is5_n1(:,1)==i;
    MD = MD_Is5_n1(:,3);
    vMD = (MD(TMD)-intercept(i))./slope(i);
    %     plot(i,vMD,'*')
    %     hold on
    
    for jk = 1:length(vMD)
        vecMD(jk,i) = vMD(jk);
    end
end

vecNND = zeros(3,12);
for i = 1:12
    TNND = NND_Is5_n1(:,1)==i;
    NND = NND_Is5_n1(:,3);
    vNND = (NND(TNND)-intercept(i))./slope(i);
    %     plot(i,vNND,'*')
    %     hold on
    
    for jk = 1:length(vNND)
        vecNND(jk,i) = vNND(jk);
    end
end

vecOD = zeros(3,12);
for i = 1:12
    TOD = OD_Is5_n1(:,1)==i;
    OD = OD_Is5_n1(:,3);
    vOD = (OD(TOD)-intercept(i))./slope(i);
    %     plot(i,vOD,'*')
    %     hold on
    
    for jk = 1:length(vOD)
        vecOD(jk,i) = vOD(jk);
    end
end

vecOM = zeros(3,12);
for i = 1:12
    
    TOM = OM_Is5_n1(:,1)==i;
    OM = OM_Is5_n1(:,3);
    vOM = (OM(TOM)-intercept(i))./slope(i);
    %     plot(i,vOM,'*')
    %     hold on
    for jk = 1:length(vOM)
        vecOM(jk,i) = vOM(jk);
    end
end



[SinkVel] = SinkingVelocity(1);

xt = [1:12]*24*SinkVel;
plot(xt,vecMD,'Color',[0, 1, 0, 0.2],'LineWidth',0.1) %[1, 0, 0, 0.2]
hold on
plot(xt,mean(vecMD),'-g','LineWidth',3)
hold on

ylim([0 27])
xlim([min(xt) (xt(12))])
ylabel('Nitrate [\muM]')
xlabel('Depth (m)')
set(gca,'LineWidth',2, 'FontSize',15)
set(gca,'FontSize',20)
set(gca,'View',[90 90])

cd '/NitriteNitrateQuantification/'

load 'NND_Is5_n2.mat'
load 'MD_Is5_n2.mat'
load 'OD_Is5_n2.mat'

intercept_T_1_3 = 4675;
slope_T_1_3 = 9406.1;

intercept_TNND3_4_5 = 6721.6;
slope_TNND3_4_5 = 10247;

intercept_T6_7 = 5389.4;
slope_T6_7 = 8107;

intercept_T8_9 = 5441.6;
slope_T8_9 = 9771.9;

intercept_T10 = 4591;
slope_T10 = 7229.6;

intercept_T11 = 2881.6;
slope_T11 = 7533.4;

intercept_T12 = 6661.4;
slope_T12 = 5813.7;

intercept_T13 = 3845;
slope_T13 = 9723.6;


vecMD =zeros(3,13);

for i = 1:12
    
    TMD = MD_Is5_n2(:,1)==i;
    MD = MD_Is5_n2(:,3);
    if (i == 1) || (i == 2) || (i == 3)
        vMD = MD(TMD)+intercept_T_1_3./slope_T_1_3;
    elseif (i == 4) || (i == 5)
        vMD = MD(TMD)+intercept_TNND3_4_5./slope_TNND3_4_5;
    elseif (i == 6) || (i == 7)
        vMD = MD(TMD)+intercept_T6_7./slope_T6_7;
    elseif (i == 8) || (i == 9)
        vMD = MD(TMD)+intercept_T8_9./slope_T8_9;
    elseif (i == 10)
        vMD = MD(TMD)+intercept_T10./slope_T10;
    elseif (i == 11)
        vMD = MD(TMD)+intercept_T11./slope_T11;
    elseif (i == 12)
        vMD = MD(TMD)+intercept_T12./slope_T12;
    elseif (i == 13)
        vMD = MD(TMD)+intercept_T13./slope_T13;
    end
    for jk = 1:length(vMD)
        vecMD(jk,i) = vMD(jk);
    end
end


[SinkVel] = SinkingVelocity(1);
xt = [1:13]*24*SinkVel;

figure,
plot(xt,vecMD./10000,'Color',[0, 1, 0, 0.2],'LineWidth',0.1) %[1, 0, 0, 0.2]
hold on
plot(xt,mean(vecMD./10000),'-g','LineWidth',3)
ylim([0 27])
xlim([min(xt) (xt(12))])
ylabel('Nitrate [\muM]')
xlabel('Depth (m)')
% legend('Marine isolate + diatom','PAO1 + DIATOM', 'ONLY + DIATOM')
set(gca,'LineWidth',2, 'FontSize',15)
set(gca,'FontSize',20)
set(gca,'View',[90 90])
%% Figure 4 f/l | Oxygen, biomass accumulation and nitrate and nitrite depth profiles.


figure,
cd '/Nitrite_Is5_NND/'
load 'Variables.mat'

slope = 0.0244;
intercept = 0.0437;
[SinkVel] = SinkingVelocity(1);
xt = [1:length(MD)]*24*SinkVel;
minT = min(min(MD(:,1:12).*slope+intercept));
MDT = (MD(:,1:12).*slope+intercept)-minT;

plot(xt,MDT*1000,'Color',[1, 0, 0, 0.2],'LineWidth',0.1)
hold on
plot(xt,mean(MDT*1000),'-r','LineWidth',3)
xlim([min(xt) max(xt)])
ax = gca;
ax.YAxis.Exponent = 0;
ylim([0 1.5])
ylabel('Nitrite (nM)')
xlabel('Depth (m)')
set(gca,'LineWidth',2, 'FontSize',15)
set(gca,'View',[90 90])

clear all
figure,
cd '/Nitrite_Is5_NND_2/'
load 'Variables.mat'

slope = 0.0244; %m
intercept = 0.0437; % q

[SinkVel] = SinkingVelocity(1);
xt = [1:length(MD(:,1:12))]*24*SinkVel;

minT = min(min(MD(:,1:12).*slope+intercept)); %baseline correction
MDT = (MD(:,1:12).*slope+intercept)-minT;

plot(xt,MDT*1000,'Color',[1, 0, 0, 0.2],'LineWidth',0.1)
hold on
plot(xt,mean(MDT*1000),'-r','LineWidth',3)
xlim([min(xt) max(xt)])
ax = gca;
ax.YAxis.Exponent = 0;
ylim([0 1.5])
ylabel('Nitrite (nM)')
xlabel('Depth (m)')
set(gca,'LineWidth',2, 'FontSize',15)
set(gca,'View',[90 90])
%% Extended Data Fig. 1 Nitrous oxide quantification

cd '/N2O_microsensor/'

load Ch1_Is5.mat
load time.mat

NarNir_N2O = Ch1_Is5;

slope = 1.143;
intercept = 2.5299;
slope2 = 1.0558;
intercept2 = 0.7148;

timeH = time/60/60;
plot(timeH,(NarNir_N2O.*slope+intercept)-(min(NarNir_N2O.*slope+intercept)),'-k','LineWidth',3)
xlim([0 34.5])
ylim([0 60])
ylabel('Nitrous oxide [\muM]')
xlabel('Time (h)')
set(gca,'LineWidth',2, 'FontSize',15)
%% Extended Data Fig. 2 Nitrite quantification with Griess method.


cd '/Nitrte_PT_CA_NN/'

load('Standards.mat')


kk = 5 %for Is5_NND

% Silent this for Is5_NND
NO2_02 = [Standards(:,1) Standards(:,1+kk) Standards(:,1+kk*2)...
    Standards(:,1+kk*3) Standards(:,1+kk*4)];

meanNO2_02 = smooth(mean(NO2_02'),0.1,'rloess');
stdNO2_02 = smooth(std(NO2_02'),0.1,'rloess');

NO2_01 = [Standards(:,2) Standards(:,2+kk) Standards(:,2+kk*2)...
    Standards(:,2+kk*3) Standards(:,2+kk*4)];

meanNO2_01 = smooth(mean(NO2_01'),0.1,'rloess');
stdNO2_01 = smooth(std(NO2_01'),0.1,'rloess');


NO2_005 = [Standards(:,3) Standards(:,3+kk) Standards(:,3+kk*2)...
    Standards(:,3+kk*3) Standards(:,3+kk*4)];

meanNO2_005 = smooth(mean(NO2_005'),0.1,'rloess');
stdNO2_005 = smooth(std(NO2_005'),0.1,'rloess');

NO2_002 = [Standards(:,4) Standards(:,4+kk) Standards(:,4+kk*2)...
    Standards(:,4+kk*3) Standards(:,4+kk*4)];

meanNO2_002 = smooth(mean(NO2_002'),0.1,'rloess');
stdNO2_002 = smooth(std(NO2_002'),0.1,'rloess');


NO2_001 = [Standards(:,5) Standards(:,5+kk) Standards(:,5+kk*2)...
    Standards(:,5+kk*3) Standards(:,5+kk*4)]

meanNO2_001 = smooth(mean(NO2_001'),0.1,'rloess');
stdNO2_001 = smooth(std(NO2_001'),0.1,'rloess');

nm543 = abs(400 - 543);
nm750 = abs(400 - 750);

x = [400:1:800];
figure,
plot(x,meanNO2_02,'b','LineWidth',2)
hold on
plot(x,meanNO2_01,'g','LineWidth',2)
hold on
plot(x,meanNO2_005,'k','LineWidth',2)
hold on
plot(x,meanNO2_002,'c','LineWidth',2)
hold on
plot(x,meanNO2_001,'r','LineWidth',2)
legend('0.2', '0.1', '0.05', '0.02', '0.01','0')
line([543 543],get(gca,'YLim'),'Color',[1 0 0],'LineWidth',2)
xlim([400 700])
set(gca,'LineWidth',2)
set(gca,'LineWidth',2, 'FontSize',15)


% ---Standard curve

%---Plot
figure,
plot(0.22,0.043)
hold on
pl= plot(0.2,meanNO2_02(nm543),'gs',...
    'LineWidth',2,...
    'MarkerSize',10,...
    'MarkerEdgeColor','b',...
    'MarkerFaceColor',[0.5,0.5,0.5])
hold on
errorbar(0.2,meanNO2_02(nm543),stdNO2_02(nm543),'b','LineWidth',2)
hold on
pl= plot(0.1,meanNO2_01(nm543),'gs',...
    'LineWidth',2,...
    'MarkerSize',10,...
    'MarkerEdgeColor','b',...
    'MarkerFaceColor',[0.5,0.5,0.5])
hold on
errorbar(0.1,meanNO2_01(nm543),stdNO2_01(nm543),'b','LineWidth',2)
hold on
pl = plot(0.05,meanNO2_005(nm543),'gs',...
    'LineWidth',2,...
    'MarkerSize',10,...
    'MarkerEdgeColor','b',...
    'MarkerFaceColor',[0.5,0.5,0.5])
hold on
errorbar(0.05,meanNO2_005(nm543),stdNO2_005(nm543),'b','LineWidth',2)
hold on
pl = plot(0.02,meanNO2_002(nm543),'gs',...
    'LineWidth',2,...
    'MarkerSize',10,...
    'MarkerEdgeColor','b',...
    'MarkerFaceColor',[0.5,0.5,0.5])
hold on
errorbar(0.02,meanNO2_002(nm543),stdNO2_002(nm543),'b','LineWidth',2)
hold on
pl = plot(0.01,meanNO2_001(nm543),'gs',...
    'LineWidth',2,...
    'MarkerSize',10,...
    'MarkerEdgeColor','b',...
    'MarkerFaceColor',[0.5,0.5,0.5])
hold on
errorbar(0.01,meanNO2_001(nm543),stdNO2_001(nm543),'b','LineWidth',2)
hold on
plot(0.005,0.043)
set(gca, 'XDir','reverse')
set(gca,'LineWidth',2)
set(gca,'LineWidth',2, 'FontSize',15)

hold on
y1 = flip([meanNO2_02(nm543) meanNO2_01(nm543) meanNO2_005(nm543) meanNO2_002(nm543) meanNO2_001(nm543)]);
x = flip([0.2 0.1 0.05 0.02 0.01]);

P = polyfit(x,y1,1);
yfit = polyval(P,x);
hold on;
plot(x,yfit,'r-.');
eqn = string('Linear: y =' + P(1)) + 'x + ' + string(P(2));
% text(min(x),max(y1),eqn,'HorizontalAlignment','left','VerticalAlignment','top')

[meanNO2_02(nm543)...
    meanNO2_01(nm543) ...
    meanNO2_005(nm543) ...
    meanNO2_002(nm543) ...
    meanNO2_001(nm543)]
slope = 0.0285; %m
intercept = 0.0442; % q
%% Extended Data Fig. 3 Aerobic and anaerobic calibration of in situ oxygen nanoprobes in the particles.

rootDirCal = {'/tiff_XY_06/tif_XY_06_Z2/'};
rootDirAnaerobic = {'/tiff_XY_06/tif_XY_06_Z2/'};

for i=1
    
    cd(rootDirCal{i});
    
    filCorrupted1 = dir('*._T*');
    nameCorrupted1 = {filCorrupted1(:).name};
    for kj = 1:length(nameCorrupted1)
        delete(char(nameCorrupted1{kj}));
    end
    
    %---Define path---
    file=dir([rootDirCal{i}, '/', '*.tif']);
    name = {file(:).name};
    NBrf = zeros(1,length(name));
    nameSort=natsort(name);
    [radii,centers,biggest] = EdgeDiatPart2(name);
    
    % ---- OXNANO AEROBIC----
    nameCharOX = char(nameSort{1});
    OxImAE = imread(nameCharOX,2);
    OxImAE_Orig = OxImAE;
    
    % ---- OXNANO ANAEROBIC ----
    cd(rootDirAnaerobic{i});
    fileAN=dir([rootDirAnaerobic{i}, '/', '*.tif']);
    nameAN = {fileAN(:).name};
    nameCharOXAN = char(nameAN{1});
    OxImAN = imread(nameCharOXAN,2);
    OxImAN_Orig = OxImAN;
    [radiiAN,centersAN,biggestAN] = EdgeDiatPart2(nameAN);
    cd(rootDirCal{i});
    
    %--------IDENTIFY DIATOM-----
    % identify diatoms in the oxnano image in the first image.
    Elong = 1.25;
    AreD = 100;
    AreDmax = 10000;
    fudgeFactor = 0.9; %0.8
    filt = 2;
    I_BF = nameSort{1};
    [output,nmbrD] = SelectCountDiatom(biggest,I_BF,AreD,Elong,AreDmax,fudgeFactor,filt);
    outputInit = output;
    for kk = 1
        % --- OXNANO IMAGE ----
        nameCharOX = char(nameSort{1});
        OxIm2 = imread(nameCharOX,2);
        % Substract the diatoms from the Oxnano image
        output2 = (~output);
        OxIm2(~output2) = NaN;
        % Substract diatom from the first image
        outputInit2 = (~outputInit);
        OxIm2(~outputInit2) = NaN;
        % Substract the diatoms from the ANAEROBIC Oxnano image
        OxImAN(~output2) = NaN;
        % Substract diatom from the first image
        OxImAN(~outputInit2) = NaN;
        % Substract the diatoms from the AEROBIC Oxnano image
        OxImAE(~output2) = NaN;
        % Substract diatom from the first image
        OxImAE(~outputInit2) = NaN;
        
        %---------OX. PARTICLE LEVEL---------
        % Measure oxnano median in 6 zones
        CM = magma(6);
        for ijx = 1:6
            radY = [1:500:radii];
            %--Mask of zone for calibration aerobic particle---
            imageSize2 = size(OxIm2);
            ci2 = [centers, radii-radY(ijx)];
            [xx2,yy2] = ndgrid((1:imageSize2(1))-ci2(1),(1:imageSize2(2))-ci2(2));
            mask2 = uint8((xx2.^2 + yy2.^2)<ci2(3)^2);
            %--Mask of zone for calibration anaerobic particle---
            radYAN = [1:500:radiiAN];
            imageSize3 = size(OxImAN);
            ci3 = [centersAN, radiiAN-radYAN(ijx)];
            [xx3,yy3] = ndgrid((1:imageSize3(1))-ci3(1),(1:imageSize3(2))-ci3(2));
            mask3 = uint8((xx3.^2 + yy3.^2)<ci3(3)^2);
            % Signal at particle level AEROBIC CALIBRATION
            croImage3 = OxImAE;
            croImage3(~mask2)= NaN;
            % Signal at particle level ANAEROBIC CALIBRATION
            % For the broken dataset, the mask is different.
            croImage4 = OxImAN;
            croImage4(~mask3)= NaN;
            %anaerobic
            ImVc3 = [croImage4(:)];
            ImVc3 = ImVc3(~isnan(ImVc3));
            vectIm3 =(ImVc3(ImVc3>0));
            minYlocal = nanmedian(vectIm3);
            %aerobic
            ImVc4 = [croImage3(:)];
            ImVc4 = ImVc4(~isnan(ImVc4));
            vectIm4 =(ImVc4(ImVc4>0));
            maxYlocal = nanmedian(vectIm4);
            
            %--- Store local calibration ----
            zoneScal{ijx} = {[minYlocal maxYlocal]};
            [cnt1,bns1z] = imhist([croImage3(croImage3>0)]); %Signal at colony level
            [cnt2,bns2z] = imhist([croImage4(croImage4>0)]); %Signal at particle level
            % Plot the histogram of the image
            subplot(1,2,1)
            plot(bns1z,cnt1,'-','color',CM(ijx,:),'LineWidth',2);
            xlim([0 26000])
            ax = gca; % axes handle
            ax.XAxis.Exponent = 0;
            ylabel('Frequency')
            xlabel('Fluorescence')
            set(gca,'LineWidth',2, 'FontSize',15)
            hold on
            subplot(1,2,2)
            plot(bns2z,cnt2,'-','color',CM(ijx,:),'LineWidth',2);
            xlim([0 26000])
            ax = gca; % axes handle
            ax.XAxis.Exponent = 0;
            ylabel('Frequency')
            xlabel('Fluorescence')
            set(gca,'LineWidth',2, 'FontSize',15)
            hold on
        end
    end
end
%% Extended Data Fig. 5-8 Probability density function of oxygen change at colony level of particles seeded with intact/broken diatoms

Ch = 1
if brkn ==1
    brk = 20
    brk2 = 40
else
    brk = 0
    brk2 = 0
end

for col=1:2 %1= colony level; 2= particle level
    figure,
    CM = magma(size(GrowthID{1,k}{:,Ch},1));
    for k = 1:length(rootDir)
        for jj = 1:50:size(GrowthID{1,k}{:,Ch},1)
            sbF = ([1:50:size(GrowthID{1,k}{:,Ch},1)]);
            
            for zon = 1:6
                idxF = find(sbF ==jj)
                if zon ==1
                    Tv = [1 7 13 19 25 31];
                    subplot(6,6,Tv(idxF))
                    A = {};
                    for jk =1:length(OxSigDistParticleZones{1,k}{1,jj}{1,1}{1,col}(:,1))
                        A{1,jk} = (ones(OxSigDistParticleZones{1,k}{1,jj}{1,1}{1,col}(jk,1),1)).*(OxSigDistParticleZones{1,k}{1,jj}{1,1}{1,col}(jk,2));
                        
                    end
                    A(~cellfun('isempty',A))
                    V = cell2mat(A');
                    [f,xi] = ksdensity((V));
                    yn =f./sum(f);
                    plot(xi,yn,'-','color',CM(jj,:),'LineWidth',2)
                    hold on
                    set(gca,'LineWidth',2)
                    set(gca,'FontSize',10)
                    xlim([-100 250])
                    ylim([0 1])
                    if Tv == 1
                        xt = get(gca, 'XTick');
                        ylabel('Frequency','FontSize',10)
                        xlabel('Oxygen (\mumM)','FontSize',10)
                        title('Zone 500 \mumm)')
                    else
                    end
                    
                    
                elseif zon ==2
                    Tv2 = [2 8 14 20 26 32];
                    subplot(6,6,Tv2(idxF))
                    A = {};
                    for jk =1:length(OxSigDistParticleZones{1,k}{1,jj}{1,1}{1,col}(:,1))
                        A{1,jk} = (ones(OxSigDistParticleZones{1,k}{1,jj}{1,2}{1,col}(jk,1),1)).*(OxSigDistParticleZones{1,k}{1,jj}{1,2}{1,col}(jk,2)+brk);
                    end
                    A(~cellfun('isempty',A))
                    V = cell2mat(A');
                    [f,xi] = ksdensity((V));
                    yn =f./sum(f);
                    plot(xi,yn,'-','color',CM(jj,:),'LineWidth',2)
                    hold on
                    set(gca,'LineWidth',2)
                    set(gca,'FontSize',10)
                    xlim([-100 250])
                    ylim([0 1])
                    if Tv2 == 2
                        xt = get(gca, 'XTick');
                        ylabel('Frequency','FontSize',10)
                        xlabel('Oxygen (\mumM)','FontSize',10)
                        title('Zone 500 \mumm)')
                    else
                    end
                    
                elseif zon ==3
                    Tv3 = [3 9 15 21 27 33];
                    subplot(6,6,Tv3(idxF))
                    A = {};
                    for jk =1:length(OxSigDistParticleZones{1,k}{1,jj}{1,1}{1,col}(:,1))
                        A{1,jk} = (ones(OxSigDistParticleZones{1,k}{1,jj}{1,3}{1,col}(jk,1),1)).*(OxSigDistParticleZones{1,k}{1,jj}{1,3}{1,col}(jk,2)+brk);
                    end
                    A(~cellfun('isempty',A))
                    V = cell2mat(A');
                    [f,xi] = ksdensity((V));
                    yn =f./sum(f);
                    plot(xi,yn,'-','color',CM(jj,:),'LineWidth',2)
                    hold on
                    set(gca,'LineWidth',2)
                    set(gca,'FontSize',10)
                    xlim([-100 250])
                    ylim([0 1])
                    if Tv3 == 3
                        xt = get(gca, 'XTick');
                        ylabel('Frequency','FontSize',10)
                        xlabel('Oxygen (\mumM)','FontSize',10)
                        title('Zone 500 \mumm)')
                    else
                    end
                    
                elseif zon ==4
                    Tv4 = [4 10 16 22 28 34];
                    subplot(6,6,Tv4(idxF))
                    A = {};
                    for jk =1:length(OxSigDistParticleZones{1,k}{1,jj}{1,1}{1,col}(:,1))
                        A{1,jk} = (ones(OxSigDistParticleZones{1,k}{1,jj}{1,4}{1,col}(jk,1),1)).*(OxSigDistParticleZones{1,k}{1,jj}{1,4}{1,col}(jk,2)+brk);
                    end
                    A(~cellfun('isempty',A))
                    V = cell2mat(A');
                    [f,xi] = ksdensity((V));
                    yn =f./sum(f);
                    plot(xi,yn,'-','color',CM(jj,:),'LineWidth',2)
                    hold on
                    set(gca,'LineWidth',2)
                    set(gca,'FontSize',10)
                    xlim([-100 250])
                    ylim([0 1])
                    if Tv4 == 4
                        xt = get(gca, 'XTick');
                        ylabel('Frequency','FontSize',10)
                        xlabel('Oxygen (\mumM)','FontSize',10)
                        title('Zone 500 \mumm)')
                    else
                    end
                    
                elseif zon ==5
                    Tv5 = [5 11 17 23 29 35];
                    subplot(6,6,Tv5(idxF))
                    A = {};
                    for jk =1:length(OxSigDistParticleZones{1,k}{1,jj}{1,1}{1,col}(:,1))
                        A{1,jk} = (ones(OxSigDistParticleZones{1,k}{1,jj}{1,5}{1,col}(jk,1),1)).*(OxSigDistParticleZones{1,k}{1,jj}{1,5}{1,col}(jk,2)+brk2);
                    end
                    A(~cellfun('isempty',A))
                    V = cell2mat(A');
                    [f,xi] = ksdensity((V));
                    yn =f./sum(f);
                    plot(xi,yn,'-','color',CM(jj,:),'LineWidth',2)
                    hold on
                    set(gca,'LineWidth',2)
                    set(gca,'FontSize',10)
                    xlim([-100 250])
                    ylim([0 1])
                    if Tv5 == 5
                        xt = get(gca, 'XTick');
                        ylabel('Frequency','FontSize',10)
                        xlabel('Oxygen (\mumM)','FontSize',10)
                        title('Zone 500 \mumm)')
                    else
                    end
                    
                elseif zon ==6
                    Tv6 = [6 12 18 24 30 36];
                    subplot(6,6,Tv6(idxF))
                    A = {};
                    for jk =1:length(OxSigDistParticleZones{1,k}{1,jj}{1,1}{1,col}(:,1))
                        A{1,jk} = (ones(OxSigDistParticleZones{1,k}{1,jj}{1,6}{1,col}(jk,1),1)).*(OxSigDistParticleZones{1,k}{1,jj}{1,6}{1,col}(jk,2)+brk2);
                    end
                    A(~cellfun('isempty',A))
                    V = cell2mat(A');
                    [f,xi] = ksdensity((V));
                    yn =f./sum(f);
                    plot(xi,yn,'-','color',CM(jj,:),'LineWidth',2)
                    hold on
                    set(gca,'LineWidth',2)
                    set(gca,'FontSize',10)
                    xlim([-100 250])
                    ylim([0 1])
                    if Tv6 == 6
                        xt = get(gca, 'XTick');
                        ylabel('Frequency','FontSize',10)
                        xlabel('Oxygen (\mumM)','FontSize',10)
                        title('Zone 500 \mumm)')
                    else
                    end
                end
            end
        end
    end
end

figure
imagesc(1:1:size(GrowthID{1,k}{:,Ch},1)*Aqtime*SinkVel);
colormap(magma)
set(gca,'FontSize',20)
%% Extended Data Fig. 9 Spatial expression of NarK and NirS genes for experiment with Phaeodactylum tricornutum

clc;
clear all;
close all;

whichDiatom = 2;
if (whichDiatom ==1)
    strgn_Mxy = '/PT_PAO1/';
    cd('/PT_PAO1/');
else
    strgn_Mxy = '/CA_diatom_PAO1/';
    cd('/CA_diatom_PAO1/');
end

% strgn_Mxy = pwd;
expname_Mxy = '\w*ShortChain\w*' ; %  begin with any number of alphanumeric or underscore carachter
nameF_Mxy = regexp(strgn_Mxy, expname_Mxy, 'match'); %match
if isempty(nameF_Mxy)
    load('WorkSpaceLong.mat');
else nameF_Mxy{:} == 'ShortChain'
    load('WorkSpaceShort.mat');
end

DistEdgChLong = {};
AreaChLong = {};
EdgeDIDLong = {};

for k = 1:length(rootDir)
    for jk = 1:length(DistEdgCh{1,k}) %P or WT or C
        
        DistEdgChLong{jk,k} = [DistEdgCh{1,k}{1,jk}{:}].*phisV;
        AreaChLong{jk,k} = [AreaCh_v{1,k}{1,jk}{:}].*phisV;
        EdgeDIDLong{jk,k} = [EdgeDID{1,k}{1,jk}{:}].*phisV;
        
    end
end


for k = 1:size(rootDir,2)
    
    A = (([DistEdgChLong{1,k}])); B = (([DistEdgChLong{2,k}]));
    
    
    A = A(A>0);
    A = A(A<1500)
    B = B(B>0);
    B = B(B<1500)
    
    Acell{k} = A;
    nanmedian([Acell{:}])
    std([Acell{:}])
    
    Bcell{k} = B;
    nanmean([Bcell{:}])
    std([Bcell{:}])
    
    
    [f,xi] = ksdensity((A))
    yn =f./sum(f)
    Acelly{k} = (yn);
    Acellx{k} = xi;
    
    p1 = plot(xi,yn,'g','LineWidth',1)
    Ac{k} = xi(find(yn==(max(yn))));
    
    p1.Color(4) = 0.2;
    hold on
    [f,xi] = ksdensity((B))
    yn =f./sum(f)
    h = plot(xi,yn,'r','LineWidth',1)
    h.Color(4) = 0.2;
    Bcelly{k} = yn;
    Bcellx{k} = xi;
    Bc{k} = xi(find(yn==(max(yn))));
    
    xlabel({'Distance to the edge (\mum)'},'FontSize',5)
    set(gca,'LineWidth',2)
    set(gca,'FontSize',10)
    legend('NarK-GFP','NirS-dsRED','Location','northeast')
    xlim([0 1500])
    hold on
end

hold on
for k = 1:size(rootDir,2)
    for kj= 1:100
        Amaty(k,kj) = Acelly{1,k}(kj);
        Amatx(k,kj) = Acellx{1,k}(kj);
    end
end

plot(mean(Amatx),mean(Amaty),'g','LineWidth',2)


hold on

for k = 1:size(rootDir,2)
    for kj= 1:100
        Bmaty(k,kj) = Bcelly{1,k}(kj);
        Bmatx(k,kj) = Bcellx{1,k}(kj);
    end
end

plot(mean(Bmatx),mean(Bmaty),'r','LineWidth',2)
xlabel({'Distance to the edge (\mum)'},'FontSize',5)

set(gca,'LineWidth',2)
set(gca,'FontSize',10)
legend('NarK-GFP','NirS-dsRED','Location','northeast')
xlim([0 1500])
hold on

% Figure settings
fig = gcf;
fig.Units               = 'centimeters';
fig.Position(3)         = 12;
fig.Position(4)         = 6;



mean([Ac{:}])
std([Ac{:}])

mean([Bc{:}])
std([Bc{:}])
%% Extended Data Fig. 10 Nitrite evolution in the water column of particles seeded with Pseudomonas aeruginosa and Phaeodactylum tricornutum or Chaetoceros affinis


cd '/Nitrte_PT_CA_NN/NN_PT/'
load 'Variables.mat'

figure,
% subplot(2,5,1)
subplot(1,2,1)

%y = mx +q
% Excell from this path
% %y = mx +q
slope = 0.0285; %m
intercept = 0.0442; % q

[SinkVel] = SinkingVelocity(1);

xt = [1:length(pTNO2)]*24*SinkVel;

SinkingVelocity(1)

xt = [1:length(pTNO2)]*24*SinkVel;

nitT = pTNO2.*slope+intercept;
nitT = nitT-0.0442;


plot(xt,nitT*1000,'Color',[1, 0, 0, 0.2],'LineWidth',0.1) %1000 transform to nM
hold on
plot(xt,mean(nitT)*1000,'-r','LineWidth',3)
ax = gca; % axes handle
ax.YAxis.Exponent = 0;

xlim([min(xt) max(xt)])
ylim([0 1])
% xticklabels({'24','48','72','96','120','144','168','192','216', '240', '264', '288', '312'})

% ylabel('Nitrite [\muM]')
ylabel('Nitrite (nM)')
xlabel('Depth (m)')
% legend('PAO1 + Phaeodactylum tricornutum')
set(gca,'LineWidth',2, 'FontSize',15)
set(gca,'View',[90 90])

hold on
% Nitrite Plotting CA


clear all

cd '/Nitrte_PT_CA_NN//NN_CA/'
load 'Variables.mat'

% subplot(2,5,2)
subplot(1,2,2)
%y = mx +q

% Excell from this path
% %y = mx +q
slope = 0.0285; %m
intercept = 0.0442; % q

[SinkVel] = SinkingVelocity(1);

xt = [1:length(pTNO2)]*24*SinkVel;

nitT = pTNO2.*slope+intercept;
nitT = nitT-0.0442;
format long
plot(xt,nitT*1000,'Color',[1, 0, 0, 0.2],'LineWidth',0.1)
hold on
plot(xt,mean(nitT)*1000,'-r','LineWidth',3)
ax = gca; % axes handle
ax.YAxis.Exponent = 0;

xlim([min(xt) max(xt)])
ylim([0 1])

% xticklabels({'24','48','72','96','120','144','168','192','216', '240', '264', '288', '312'})


% ylabel('Nitrite [\muM]')
% ylabel('Nitrite [\muM]')
ylabel('Nitrite (nM)')
xlabel('Depth (m)')
% legend('PAO1 + Chaetocerotaceae affinis')
set(gca,'LineWidth',2, 'FontSize',15)
set(gca,'View',[90 90])
%% Extended Data Fig. 11 Nitrate evolution in the water column.

cd '/NitriteNitrateQuantification/'

load 'NND_Is5_n1.mat'
load 'MD_Is5_n1.mat'
load 'OD_Is5_n1.mat'
load 'OM_Is5_n1.mat'

intercept =[138.71 1244.7 2420.6 3024.4 717.87 1371.9 8876.2 121.68 421.69 706.02 2659.1 768.3];
slope = [5195.2 5185.1 5667.3 5481 9045.2 7337.3 5349.6 9326.1 9330.8 10160 8553.5 9653.5];

vecMD = zeros(3,12);
for i = 1:12
    TMD = MD_Is5_n1(:,1)==i;
    MD = MD_Is5_n1(:,3);
    vMD = (MD(TMD)-intercept(i))./slope(i);
    for jk = 1:length(vMD)
        vecMD(jk,i) = vMD(jk);
    end
end

vecNND = zeros(3,12);
for i = 1:12
    TNND = NND_Is5_n1(:,1)==i;
    NND = NND_Is5_n1(:,3);
    vNND = (NND(TNND)-intercept(i))./slope(i);
    for jk = 1:length(vNND)
        vecNND(jk,i) = vNND(jk);
    end
end

vecOD = zeros(3,12);
for i = 1:12
    TOD = OD_Is5_n1(:,1)==i;
    OD = OD_Is5_n1(:,3);
    vOD = (OD(TOD)-intercept(i))./slope(i);
    for jk = 1:length(vOD)
        vecOD(jk,i) = vOD(jk);
    end
end

vecOM = zeros(3,12);
for i = 1:12
    TOM = OM_Is5_n1(:,1)==i;
    OM = OM_Is5_n1(:,3);
    vOM = (OM(TOM)-intercept(i))./slope(i);
    for jk = 1:length(vOM)
        vecOM(jk,i) = vOM(jk);
    end
end

[SinkVel] = SinkingVelocity(1);
xt = [1:12]*24*SinkVel;

subplot(1,5,1)
plot(xt,vecNND,'Color',[0, 1, 0, 0.2],'LineWidth',0.1)
hold on
plot(xt,mean(vecNND),'-g','LineWidth',3)
ylim([0 27])
xlim([min(xt) (xt(12))])
ylabel('Nitrate [\muM]')
xlabel('Depth (m)')
set(gca,'LineWidth',2, 'FontSize',15)
set(gca,'FontSize',20)
set(gca,'View',[90 90])
hold on
subplot(1,5,2)
plot(xt,vecOD,'Color',[0, 1, 0, 0.2],'LineWidth',0.1)
hold on
plot(xt,mean(vecOD),'-g','LineWidth',3)
ylim([0 27])
xlim([min(xt) (xt(12))])
ylabel('Nitrate [\muM]')
xlabel('Depth (m)')
set(gca,'LineWidth',2, 'FontSize',15)
set(gca,'FontSize',20)
set(gca,'View',[90 90])

% hold on
subplot(1,5,3)
plot(xt,vecOM,'Color',[0, 1, 0, 0.2],'LineWidth',0.1)
hold on
plot(xt,mean(vecOM),'-g','LineWidth',3)
ylim([0 27])
xlim([min(xt) (xt(12))])
ylabel('Nitrate [\muM]')
xlabel('Depth (m)')
set(gca,'LineWidth',2, 'FontSize',15)
set(gca,'FontSize',20)
set(gca,'View',[90 90])
hold on

% ---Nitrate Is5 #2 intact diatoms


cd '/NitriteNitrateQuantification/'

load 'NND_Is5_n2.mat'
load 'MD_Is5_n2.mat'
load 'OD_Is5_n2.mat'

intercept_T_1_3 = 4675;
slope_T_1_3 = 9406.1;

intercept_TNND3_4_5 = 6721.6;
slope_TNND3_4_5 = 10247;

intercept_T6_7 = 5389.4;
slope_T6_7 = 8107;

intercept_T8_9 = 5441.6;
slope_T8_9 = 9771.9;

intercept_T10 = 4591;
slope_T10 = 7229.6;

intercept_T11 = 2881.6;
slope_T11 = 7533.4;

intercept_T12 = 6661.4;
slope_T12 = 5813.7;

intercept_T13 = 3845;
slope_T13 = 9723.6;


vecMD =zeros(3,13);

for i = 1:12
    
    TMD = MD_Is5_n2(:,1)==i;
    MD = MD_Is5_n2(:,3);
    if (i == 1) || (i == 2) || (i == 3)
        vMD = MD(TMD)+intercept_T_1_3./slope_T_1_3;
    elseif (i == 4) || (i == 5)
        vMD = MD(TMD)+intercept_TNND3_4_5./slope_TNND3_4_5;
    elseif (i == 6) || (i == 7)
        vMD = MD(TMD)+intercept_T6_7./slope_T6_7;
    elseif (i == 8) || (i == 9)
        vMD = MD(TMD)+intercept_T8_9./slope_T8_9;
    elseif (i == 10)
        vMD = MD(TMD)+intercept_T10./slope_T10;
    elseif (i == 11)
        vMD = MD(TMD)+intercept_T11./slope_T11;
    elseif (i == 12)
        vMD = MD(TMD)+intercept_T12./slope_T12;
    elseif (i == 13)
        vMD = MD(TMD)+intercept_T13./slope_T13;
    end
    for jk = 1:length(vMD)
        vecMD(jk,i) = vMD(jk);
    end
end

vecNND =zeros(3,13);
for i = 1:12
    
    TNND = NND_Is5_n2(:,1)==i;
    NND = NND_Is5_n2(:,3);
    if (i == 1) || (i == 2)
        vNND = NND(TNND)+intercept_T_1_3./slope_T_1_3;
    elseif (i == 3) || (i == 4) || (i == 5)
        vNND = NND(TNND)+intercept_TNND3_4_5./slope_TNND3_4_5;
    elseif (i == 6) || (i == 7)
        vNND = NND(TNND)+intercept_T6_7./slope_T6_7;
    elseif (i == 8) || (i == 9)
        vNND = NND(TNND)+intercept_T8_9./slope_T8_9;
    elseif (i == 10)
        vNND = NND(TNND)+intercept_T10./slope_T10;
    elseif (i == 11)
        vNND = NND(TNND)+intercept_T11./slope_T11;
    elseif (i == 12)
        vNND = NND(TNND)+intercept_T12./slope_T12;
    elseif (i == 13)
        vNND = NND(TNND)+intercept_T13./slope_T13;
    end
    for jk = 1:length(vNND)
        vecNND(jk,i) = vNND(jk);
    end
end

vecOD =zeros(3,13);

for i = 1:12
    
    TOD = OD_Is5_n2(:,1)==i;
    OD = OD_Is5_n2(:,3);
    if (i == 1) || (i == 2) || (i == 3)
        vOD = OD(TOD)+intercept_T_1_3./slope_T_1_3;
    elseif  (i == 4) || (i == 5)
        vOD = OD(TOD)+intercept_TNND3_4_5./slope_TNND3_4_5;
    elseif (i == 6) || (i == 7)
        vOD = OD(TOD)+intercept_T6_7./slope_T6_7;
    elseif (i == 8) || (i == 9)
        vOD = OD(TOD)+intercept_T8_9./slope_T8_9;
    elseif (i == 10)
        vOD = OD(TOD)+intercept_T10./slope_T10;
    elseif (i == 11)
        vOD = OD(TOD)+intercept_T11./slope_T11;
    elseif (i == 12)
        vOD = OD(TOD)+intercept_T12./slope_T12;
    elseif (i == 13)
        vOD = OD(TOD)+intercept_T13./slope_T13;
    end
    for jk = 1:length(vOD)
        vecOD(jk,i) = vOD(jk);
    end
end

[SinkVel] = SinkingVelocity(1);

xt = [1:13]*24*SinkVel;

subplot(1,5,4)
plot(xt,vecNND./10000,'Color',[0, 1, 0, 0.2],'LineWidth',0.1)
hold on
plot(xt,nanmean(vecNND)./10000,'-g','LineWidth',3)
ylim([0 27])
xlim([min(xt) (xt(12))])
ylabel('Nitrate [\muM]')
xlabel('Depth (m)')
set(gca,'LineWidth',2, 'FontSize',15)
set(gca,'FontSize',20)
set(gca,'View',[90 90])

% hold on
subplot(1,5,5)
plot(xt,vecOD./10000,'Color',[0, 1, 0, 0.2],'LineWidth',0.1)
hold on
plot(xt,mean(vecOD)./10000,'-g','LineWidth',3)

ylim([0 27])
xlim([min(xt) (xt(12))])
ylabel('Nitrate [\muM]')
xlabel('Depth (m)')
set(gca,'LineWidth',2, 'FontSize',15)
set(gca,'FontSize',20)
set(gca,'View',[90 90])
%% Extended Data Fig. 12 Nitrite evolution in the water column.


cd '/Nitrite_Is5_NND/'
load 'Variables.mat'
figure,
subplot(1,5,1)
minT = min(min((D.*slope+intercept)));
nnT = (D.*slope+intercept)-minT;
plot(xt,((D.*slope+intercept)-min(D.*slope+intercept))*1000,'Color',[1, 0, 0, 0.2],'LineWidth',0.1)
hold on
plot(xt,(mean(D.*slope+intercept)-min(D.*slope+intercept))*1000,'-r','LineWidth',3)
xlim([min(xt) max(xt)])
ylim([0 1.5])
ax = gca;
ax.YAxis.Exponent = 0;
ylabel('Nitrite (nM)')
xlabel('Depth (m)')
set(gca,'LineWidth',2, 'FontSize',15)
set(gca,'View',[90 90])

subplot(1,5,3)
minT = min(min((M.*slope+intercept)));
nnT = (M.*slope+intercept)-minT;
plot(xt,((M.*slope+intercept)-min(M.*slope+intercept))*1000,'Color',[1, 0, 0, 0.2],'LineWidth',0.1)
hold on
plot(xt,(mean(M.*slope+intercept)-min(M.*slope+intercept))*1000,'-r','LineWidth',3)
xlim([min(xt) max(xt)])
ylim([0 1.5])
ax = gca;
ax.YAxis.Exponent = 0;
ylabel('Nitrite (nM)')
xlabel('Depth (m)')
set(gca,'LineWidth',2, 'FontSize',15)
set(gca,'View',[90 90])
hold on

% Nitrite Plotting intact diatom CA Is5_NND #2
cd '/Nitrite_Is5_NND_2/'
load 'Variables.mat'

slope = 0.0244; %m
intercept = 0.0437; % q
[SinkVel] = SinkingVelocity(1);
xt = [1:length(MD(:,1:12))]*24*SinkVel;

minT = min(min(NND(:,1:12).*slope+intercept)); %baseline correction
nnT = (NND(:,1:12).*slope+intercept)-minT;

subplot(1,5,4)
plot(xt,nnT*1000,'Color',[1, 0, 0, 0.2],'LineWidth',0.1)
hold on
plot(xt,mean(nnT*1000),'-r','LineWidth',3)
xlim([min(xt) max(xt)])
ax = gca; % axes handle
ax.YAxis.Exponent = 0;
ylim([0 1.5])
ylabel('Nitrite (nM)')
xlabel('Depth (m)')
set(gca,'LineWidth',2, 'FontSize',15)
set(gca,'View',[90 90])

minT = min(min(D(:,1:12).*slope+intercept));
nnT = (D(:,1:12).*slope+intercept)-minT;

subplot(1,5,5)
plot(xt,nnT*1000,'Color',[1, 0, 0, 0.2],'LineWidth',0.1)
hold on
plot(xt,mean(nnT)*1000,'-r','LineWidth',3)
xlim([min(xt) max(xt)])
ax = gca; % axes handle
ax.YAxis.Exponent = 0;
ylim([0 1.5])
% ylabel('Nitrite [\muM]')
ylabel('Nitrite (nM)')
xlabel('Depth (m)')
% legend('Marine isolate + diatom','PAO1 + DIATOM', 'ONLY + DIATOM')
set(gca,'LineWidth',2, 'FontSize',15)
set(gca,'View',[90 90])
%
subplot(1,5,2)
minT = min(min(M(:,1:12).*slope+intercept)); %baseline correction
nnT = (M(:,1:12).*slope+intercept)-minT;

plot(xt,nnT*1000,'Color',[1, 0, 0, 0.2],'LineWidth',0.1)
hold on
plot(xt,(mean(nnT)*1000),'-r','LineWidth',3)
xlim([min(xt) max(xt)])
ax = gca; % axes handle
ax.YAxis.Exponent = 0;
ylim([0 1.5])
% ylabel('Nitrite [\muM]')
ylabel('Nitrite (nM)')

xlabel('Depth (m)')
% legend('Marine isolate + diatom','PAO1 + DIATOM','Only Marine', 'ONLY + DIATOM')
set(gca,'LineWidth',2, 'FontSize',15)
set(gca,'View',[90 90])
%% Extended Data Fig. 13 DOC evolution in the water column

cd '/TOC_BROKEN/'
load 'TOC.mat'

TOCbroke = TOC;
[SinkVel] = SinkingVelocity(1);
xt = [1:12]*24*SinkVel;

figure,
subplot(1,5,1)
plot(xt,TOC(4:6,:),'Color',[0, 1, 1, 0.2],'LineWidth',0.1)
hold on
plot(xt,nanmedian(TOC(4:6,:)),'-c','LineWidth',3)
xlim([min(xt) max(xt)])
ylim([0 30])
ylabel('TOC (mgC/L)')
xlabel('Hours')
set(gca,'LineWidth',2, 'FontSize',15)
set(gca,'View',[90 90])

% hold on
subplot(1,5,2)
plot(xt,TOC(7:9,:),'Color',[0, 1, 1, 0.2],'LineWidth',0.1)
hold on
plot(xt,nanmedian(TOC(7:9,:)),'-c','LineWidth',3)
xlim([min(xt) max(xt)])
ylim([0 30])
ylabel('TOC (mgC/L)')
xlabel('Hours')
set(gca,'LineWidth',2, 'FontSize',15)
set(gca,'View',[90 90])

% hold on
subplot(1,5,3)
plot(xt,TOC(10:12,:),'Color',[0, 1, 1, 0.2],'LineWidth',0.1)
hold on
plot(xt,nanmedian(TOC(10:12,:)),'-c','LineWidth',3)
xlim([min(xt) max(xt)])
ylim([0 30])
ylabel('TOC (mgC/L)')
xlabel('Hours')
set(gca,'LineWidth',2, 'FontSize',15)
set(gca,'View',[90 90])
hold on

cd '/TOC_INTACT/'

load 'TOC.mat'
TOCintact = TOC;
[SinkVel] = SinkingVelocity(1);
xt = [1:12]*24*SinkVel;

[p,h] = ttest2(TOCintact(1,1:3)', TOCbroke(1:3,1))
mean([(TOCbroke(1:3,1))])
std([(TOCbroke(1:3,1))])
mean([(TOCintact(1,1:3)')])
std([(TOCintact(1,1:3)')])

subplot(1,5,4)
plot(xt,TOC(:,4:6)','Color',[0, 1, 1, 0.2],'LineWidth',0.1)
hold on
plot(xt,nanmedian(TOC(:,4:6)'),'-c','LineWidth',3)
xlim([min(xt) max(xt)])
ylim([0 30])
ylabel('TOC (mgC/L)')
xlabel('Hours')
set(gca,'LineWidth',2, 'FontSize',15)
set(gca,'View',[90 90])

subplot(1,5,5)
plot(xt,TOC(:,7:9)','Color',[0, 1, 1, 0.2],'LineWidth',0.1)
hold on
plot(xt,nanmedian(TOC(:,7:9)'),'-c','LineWidth',3)
xlim([min(xt) max(xt)])
ylim([0 30])
ylabel('TOC (mgC/L)')
xlabel('Hours')
set(gca,'LineWidth',2, 'FontSize',15)
set(gca,'View',[90 90])