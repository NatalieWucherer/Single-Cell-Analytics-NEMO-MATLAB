function nemo_umap_gene_browser(matfile)

if nargin==0
    [f,p] = uigetfile('*.mat');
    if isequal(f,0); return; end
    matfile = fullfile(p,f);
end

S = load(matfile);

%%  DATA LOAD 
if isfield(S,'X')
    X = double(S.X);
elseif isfield(S,'counts')
    X = double(S.counts);
elseif isfield(S,'expr')
    X = double(S.expr);
else
    error('No expression matrix found (X / counts / expr)');
end

umap = double(S.umap);
genes = upper(string(S.genes));

if isfield(S,'louvain')
    clusters = categorical(string(S.louvain));
elseif isfield(S,'clusters')
    clusters = categorical(string(S.clusters));
else
    error('No cluster field found');
end

nCells = size(umap,1);

%%  REGION 
[~, matname, ~] = fileparts(matfile);
prefix = upper(matname(1:2));
titlePrefix = upper(extractBefore(matname,7));   % first 6 letters

regionMap = containers.Map();
regionMap('TH') = "thalamus";
regionMap('CX') = "cortex";
regionMap('GE') = "ganglionic eminence";
regionMap('HP') = "hippocampus";

region = regionMap(prefix);

%%  ALIASES 
aliasMap = containers.Map();
aliasMap('COUPTF1') = 'NR2F1';
aliasMap('COUPTF2') = 'NR2F2';
aliasMap('NKX2-1') = 'NKX2.1';

priorityMap = containers.Map();

priorityMap('SST') = 3;
priorityMap('LHX6') = 2;
priorityMap('GAD1') = 1;
priorityMap('GAD2') = 1;
%%  MARKERS 
markerDB = containers.Map();

markerDB("thalamus") = {
"THALAMIC_NEURONS", ["SLC17A6","SLC17A7"];
"MIGRATING NEURONS", ["ROBO1","ROBO2","SLIT2"];
"GLUTAMATERGIC", ["SLC17A6"];
"GABAERGIC", ["GAD1","GAD2"];
"RADIAL GLIA AND ASTROCYTES", ["HES1","PAX6","GFAP"];
"OPC", ["OLIG2","PDGFRA"];
"ACTIVELY DIVIDING PROGENITORS",["MKI67"];
"MICROGLIA",["CCL3"];
"PERICYTES AND ENDOTHELIAL CELLS",["BGN"];
"RBC",["HBA1"];
};

markerDB("cortex") = {
"EXC L2/3", ["CUX2"];
"EXC L4",  ["RORB"];
"EXC L5",  ["CRYM"];
"EXC L6",  ["TLE4"];
"MIGRATING NEURONS", ["ROBO1","ROBO2","SLIT2"];
"GLUTAMATERGIC", ["SLC17A7"];
"GABAERGIC (LHX6)",["LHX6"];
"GABERGIC (NR2F2)",["NR2F2"];
"RADIAL GLIA AND ASTROCYTES", ["HES1","PAX6","GFAP"];
"OPC", ["OLIG2","PDGFRA"];
"ACTIVELY DIVIDING PROGENITORS",["MKI67"];
};

markerDB("ganglionic eminence") = {
"GABAERGIC GROUP 2", ["LHX6","SST","NPY"];
"GABAERGIC GROUP 1",["SP8"];
"PROLIFERATIVE CGE", ["NR2F1"];
"PALLIDAL PROJECTION NEURONS", ["NKX2-1"];
"VENTRAL CGE", ["COUPTF2"];
"MIGRATING GABAERGIC",["ARX"];
"MIGRATING NEURONS", ["ROBO1","ROBO2","SLIT2"];
"GLUTAMATERGIC", ["SLC17A6","SLC17A7"];
"GABAERGIC", ["GAD1","GAD2"];
"RADIAL GLIA AND ASTROCYTES", ["HES1","PAX6","GFAP"];
"OPC", ["OLIG2","PDGFRA"];
"ACTIVELY DIVIDING PROGENITORS",["MKI67"];
"MICROGLIA",["CCL3"];
"PERICYTES AND ENDOTHELIAL CELLS",["BGN"];
"RBC",["HBA1"];
};

markerSets = markerDB(region);

%%  GENE INDEX 
geneIndex = containers.Map(cellstr(genes),1:length(genes));

%%  CLUSTER LABELS 
uniqueClusters = unique(clusters);
clusterLabel = strings(length(uniqueClusters),1);

for c = 1:length(uniqueClusters)

    mask = clusters == uniqueClusters(c);
    scores = -inf(size(markerSets,1),1);

    for m = 1:size(markerSets,1)

        markers = upper(string(markerSets{m,2}));
        idxs = [];

        for k = 1:length(markers)
            g = markers(k);
            if isKey(geneIndex,g)
                idxs(end+1)=geneIndex(g);
            end
        end

        if ~isempty(idxs)
            scores(m)=mean(X(mask,idxs),'all');
        end
    end

    % TOP 2 MARKER SETS 
    [~,sortedIdx] = sort(scores,'descend');
    topN = min(2,numel(sortedIdx));

   topLabels = strings(1,topN);

for i = 1:topN
    topLabels(i) = string(markerSets{sortedIdx(i),1});
end

topLabels = sort(topLabels);

clusterLabel(c) = strjoin(topLabels," / ");
end

%%  COLOR MAP 
uniqueLabels = unique(clusterLabel);
colorMap = containers.Map();

for i=1:length(uniqueLabels)

    key = char(uniqueLabels(i));
    seed = sum(double(key));

    colorMap(key)=hsv2rgb([mod(seed,360)/360 0.85 0.95]);

end

%%  CELL COLORS 
cellColor = zeros(nCells,3);

for i=1:length(uniqueClusters)

    mask = clusters==uniqueClusters(i);
    col = colorMap(char(clusterLabel(i)));

    cellColor(mask,:) = repmat(col,sum(mask),1);

end

%%  PADDED LIMITS 
padX = 0.08 * range(umap(:,1));
padY = 0.08 * range(umap(:,2));

xlimUMAP = [min(umap(:,1))-padX, max(umap(:,1))+padX];
ylimUMAP = [min(umap(:,2))-padY, max(umap(:,2))+padY];

%%  UI 
fig = uifigure('Name',titlePrefix + " NeMO UMAP Browser",...
    'Position',[100 100 1250 760]);

ax = uiaxes(fig,'Position',[280 60 850 660]);
axis(ax,'equal')
ax.XTick=[]; ax.YTick=[];

axLegend = uiaxes(fig,'Position',[1140 120 300 520]);
axis(axLegend,'off')
title(axLegend,'Cluster Atlas','FontWeight','bold');

edit = uieditfield(fig,'text',...
    'Position',[20 650 240 30]);

uibutton(fig,'Text','UMAP',...
    'Position',[20 600 120 30],...
    'ButtonPushedFcn',@(s,e)showUMAP());

uibutton(fig,'Text','Clusters',...
    'Position',[150 600 120 30],...
    'ButtonPushedFcn',@(s,e)showClusters());

uibutton(fig,'Text','Plot Gene',...
    'Position',[20 560 120 30],...
    'ButtonPushedFcn',@(s,e)plotGene());

list = uilistbox(fig,...
    'Position',[20 250 240 280],...
    'Items',cellstr(genes),...
    'ValueChangedFcn',@(s,e)set(edit,'Value',s.Value));

showUMAP();
showClusterGrid();

%%  FUNCTIONS 
function showUMAP()
    cla(ax)
    scatter(ax,umap(:,1),umap(:,2),8,[0.75 0.75 0.75],'filled');
    title(ax,titlePrefix + " UMAP - " + region);
    xlim(ax,xlimUMAP)
    ylim(ax,ylimUMAP)
end

function showClusters()
    cla(ax)
    scatter(ax,umap(:,1),umap(:,2),10,cellColor,'filled');
    axis(ax,'equal')
    title(ax,titlePrefix + " Atlas Cluster Map")
    xlim(ax,xlimUMAP)
    ylim(ax,ylimUMAP)
    showClusterGrid();
end

function plotGene()

g = upper(strtrim(string(edit.Value)));
idx = find(strcmpi(genes,g),1);

cla(ax)
if isempty(idx); return; end

expr = X(:,idx);
[~,ord] = sort(expr,'ascend');

scatter(ax,umap(ord,1),umap(ord,2),8,expr(ord),'filled');

axis(ax,'equal')
title(ax,"Gene: " + g);

applyFeatureMap(ax,expr)
colorbar(ax)

xlim(ax,xlimUMAP)
ylim(ax,ylimUMAP)

end

function showClusterGrid()

cla(axLegend)
hold(axLegend,'on')

for i=1:length(uniqueLabels)

    lab = uniqueLabels(i);
    col = colorMap(char(lab));

    y = length(uniqueLabels)-i;

    plot(axLegend,1,y,'s',...
        'MarkerSize',14,...
        'MarkerFaceColor',col,...
        'MarkerEdgeColor',col);

    text(axLegend,1.2,y,lab,...
        'FontSize',10,...
        'FontWeight','bold',...
        'Color',col,...
        'Interpreter','none');

end

xlim(axLegend,[0.5 3])
ylim(axLegend,[0 length(uniqueLabels)+1])
axis(axLegend,'off')

end

end

%% COLOR MAP 
function applyFeatureMap(axh,expr)

n = 256;
x = linspace(0,1,n)';

c1 = [0.92 0.92 0.92];
c2 = [1.00 1.00 0.85];
c3 = [1.00 0.90 0.20];
c4 = [1.00 0.50 0.00];
c5 = [0.70 0.00 0.00];

cmap = zeros(n,3);

for k=1:3
    cmap(:,k)=interp1([0 0.20 0.45 0.75 1],...
                      [c1(k) c2(k) c3(k) c4(k) c5(k)],...
                      x,'pchip');
end

colormap(axh,cmap)

mx = prctile(expr,99);
if mx<=0; mx=max(expr); end
if mx>0
    caxis(axh,[0 mx]);
end

end
