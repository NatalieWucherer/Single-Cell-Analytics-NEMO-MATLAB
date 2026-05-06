function nemo_umap_multigene_panels(matfile, geneList)

if nargin < 2
    error('Usage: nemo_umap_multigene_panels(matfile, geneList)');
end

% LOAD DATA

S = load(matfile);

X        = double(S.X);
umap     = double(S.umap);
genes    = upper(string(S.genes));
clusters = categorical(string(S.louvain));

nCells = size(umap,1);

if size(X,1) ~= nCells
    error('Rows of X must match rows of umap');
end

% OUTPUT

[outdir, matname, ~] = fileparts(matfile);

if isempty(outdir)
    outdir = pwd;
end

prefix = upper(matname(1:2));

saveprefix= upper(matname(1:6));

% GENE LOOKUP

geneIndex = containers.Map(cellstr(genes),1:length(genes));

%%  REGION MAP
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

%% MARKERS 
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

% GENE INPUT PARSER

if ischar(geneList) || isstring(geneList)

    geneList = string(geneList);
    geneList = erase(geneList,"[");
    geneList = erase(geneList,"]");
    geneList = erase(geneList,"""");
    geneList = split(geneList,",");

elseif iscell(geneList)

    geneList = string(geneList);

end

geneList = upper(strtrim(geneList));
geneList = geneList(geneList ~= "");
geneList = geneList(:);


% CLUSTER LABELING (TOP 2 MARKERS)


uClust = categories(clusters);
nClust = numel(uClust);

clusterLabel = strings(nClust,1);

for c = 1:nClust

    mask = clusters == uClust{c};
    scores = -inf(size(markerSets,1),1);

    % score each marker set
    for m = 1:size(markerSets,1)

        markers = upper(string(markerSets{m,2}));
        idxs = [];

        for g = 1:length(markers)
            gg = markers(g);

            if isKey(aliasMap,gg)
                gg = upper(string(aliasMap(gg)));
            end

            if isKey(geneIndex,gg)
                idxs(end+1) = geneIndex(gg);
            end
        end

        if ~isempty(idxs)
            scores(m) = mean(X(mask,idxs),'all');
        end
    end

    % TOP 2 LABEL LOGIC

    [~,sortedIdx] = sort(scores,'descend');
    topN = min(2,numel(sortedIdx));

    topLabels = strings(1,topN);

    for i = 1:topN
        topLabels(i) = string(markerSets{sortedIdx(i),1});
    end

    % enforce stable ordering (prevents flip like A/B vs B/A)
    topLabels = sort(topLabels);

    clusterLabel(c) = strjoin(topLabels," / ");

end

% UNIQUE LABELS

atlasLabels = unique(clusterLabel);
nAtlas = numel(atlasLabels);

% COLOR MAP (UNCHANGED)

colorMap = containers.Map();

hues = linspace(0,1,nAtlas+1);
hues(end) = [];

for i = 1:nAtlas

    key = char(atlasLabels(i));

    seed = sum(double(key)) + 98765;
    rng(seed,'twister');

    h = hues(i) + (rand-0.5)*0.05;
    h = mod(h,1);

    colorMap(key) = hsv2rgb([h 0.85 0.95]);

end

% CELL COLORS (UNCHANGED)

cellColor = zeros(nCells,3);

for c = 1:nClust

    mask = clusters == uClust{c};

    lab = char(clusterLabel(c));
    col = colorMap(lab);

    cellColor(mask,:) = repmat(col,sum(mask),1);

end

% PADDED LIMITS (UNCHANGED)

padX = 0.08 * range(umap(:,1));
padY = 0.08 * range(umap(:,2));

xmin = min(umap(:,1)) - padX;
xmax = max(umap(:,1)) + padX;
ymin = min(umap(:,2)) - padY;
ymax = max(umap(:,2)) + padY;

% FIGURE 1 (UNCHANGED)

fig1 = figure('Color','w','Position',[100 100 1400 900]);

scatter(umap(:,1),umap(:,2),10,cellColor,'filled');
axis equal off
hold on

xlim([xmin xmax + range(umap(:,1))*0.45]);
ylim([ymin ymax]);

title(saveprefix + " Atlas Cluster Map",'FontSize',16,'FontWeight','bold');

sideX = xmax + range(umap(:,1))*0.10;
step  = range(umap(:,2)) / (nAtlas + 1);

for i = 1:nAtlas

    lab = atlasLabels(i);
    col = colorMap(char(lab));

    yPos = ymax - i*step;

    plot(sideX,yPos,'s','MarkerSize',12,...
        'MarkerFaceColor',col,'MarkerEdgeColor',col);

    text(sideX + range(umap(:,1))*0.02, yPos, lab,...
        'Color',col,'FontSize',11,'FontWeight','bold');

end

atlasFile = fullfile(outdir,saveprefix + "_AtlasClusterMap.tif");
exportgraphics(fig1,atlasFile,'Resolution',300);

disp("Saved: " + atlasFile);

% FIGURE 2 (UNCHANGED)

nGenes = numel(geneList);

nCols = 4;
nRows = ceil(nGenes / nCols);

fig2 = figure('Color','w',...
    'Position',[100 100 360*nCols 320*nRows]);

tiledlayout(nRows,nCols,'Padding','compact','TileSpacing','compact');

for i = 1:nGenes

    nexttile

    g = geneList(i);

    if ~isKey(geneIndex,g)

        scatter(umap(:,1),umap(:,2),8,[0.85 0.85 0.85],'filled');
        axis equal off
        title(g + " (NOT FOUND)");
        continue;

    end

    expr = X(:,geneIndex(g));

    [~,ord] = sort(expr,'ascend');

    scatter(umap(ord,1),umap(ord,2),4,expr(ord),'filled');

    axis equal off
    title(g,'Interpreter','none');

    n = 256;
    x = linspace(0,1,n)';

    c1 = [0.92 0.92 0.92];
    c2 = [1 0.98 0.85];
    c3 = [1 0.90 0.20];
    c4 = [1 0.45 0];
    c5 = [0.7 0 0];

    cmap = zeros(n,3);

    for k = 1:3
        cmap(:,k) = interp1([0 0.20 0.45 0.72 1],...
            [c1(k) c2(k) c3(k) c4(k) c5(k)],x,'pchip');
    end

    colormap(gca,cmap);
    colorbar;

end

geneFile = fullfile(outdir,...
    saveprefix + "_" + strjoin(geneList,"_") + "_MultiUMAP.tif");

exportgraphics(fig2,geneFile,'Resolution',300);

disp("Saved: " + geneFile);

%% FIGURE 3 — DOT PLOT (CLEAN SPACING + SHORT LABELS)

% --- Map each cell → atlas label ---
cellAtlasLabel = strings(nCells,1);

for c = 1:nClust
    mask = clusters == uClust{c};
    cellAtlasLabel(mask) = clusterLabel(c);
end

% --- Unique atlas labels ---
atlasLabels = unique(cellAtlasLabel);
nAtlas = numel(atlasLabels);
nGenes = numel(geneList);

% --- Compute metrics ---
pctExpr = zeros(nGenes,nAtlas);
avgExpr = zeros(nGenes,nAtlas);

for a = 1:nAtlas
    mask = cellAtlasLabel == atlasLabels(a);

    for g = 1:nGenes
        gene = geneList(g);

        if isKey(geneIndex,gene)
            vals = X(mask, geneIndex(gene));

            pctExpr(g,a) = mean(vals > 0) * 100;
            avgExpr(g,a) = mean(vals);
        end
    end
end

% --- Normalize expression ---
avgExprScaled = (avgExpr - min(avgExpr(:))) ./ ...
                (max(avgExpr(:)) - min(avgExpr(:)) + eps);

%% --- AUTO-SHORTEN LABELS ---
shortLabels = atlasLabels;

for i = 1:numel(shortLabels)

    lab = shortLabels(i);

    lab = replace(lab,"RADIAL GLIA AND ASTROCYTES","RG/AST");
    lab = replace(lab,"ACTIVELY DIVIDING PROGENITORS","DIV PROG");
    lab = replace(lab,"PERICYTES AND ENDOTHELIAL CELLS","PERI/ENDO");
    lab = replace(lab,"GLUTAMATERGIC","GLU");
    lab = replace(lab,"GABAERGIC","GABA");
    lab = replace(lab,"MIGRATING NEURONS","MIG");

    lab = replace(lab," / ","/");

    if strlength(lab) > 18
        lab = replace(lab,"/","/\newline");
    end

    shortLabels(i) = lab;
end

%% --- PLOT ---
fig3 = figure('Color','w','Position',[200 200 900 600]); % narrower width
hold on

for g = 1:nGenes
    for a = 1:nAtlas
        dotSize = 20 + pctExpr(g,a)*1.5;
        scatter(a, g, dotSize, avgExprScaled(g,a), 'filled');
    end
end

% --- SAME colormap as UMAP ---
colormap(parula);
colorbar;


%% --- AXIS FORMATTING ---
set(gca,'YDir','reverse',...
    'XTick',1:nAtlas,...
    'XTickLabel',shortLabels,...
    'YTick',1:nGenes,...
    'YTickLabel',geneList);

xtickangle(45);

xlim([0.5, nAtlas + 0.5]);

ylim([0.3, nGenes + 0.7]);

ax = gca;
ax.FontSize = 8;
ax.XAxis.FontSize = 7;
ax.YAxis.FontSize = 9;
ax.YAxis.TickLabelGapOffset = 6;  % <-- key line

xlabel('Cell Types');
ylabel('Genes');
title('Dot Plot');

% --- SAVE ---
dotFile = fullfile(outdir,...
    saveprefix + "_" + strjoin(geneList,"_") + "_DotPlot.tif");
exportgraphics(fig3,dotFile,'Resolution',300);

disp("Saved: " + dotFile);


end
