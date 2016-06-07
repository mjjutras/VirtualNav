function YC1_lassoReport(dataDir,ana_name,saveDir)

% find all subjects
subjFiles     = dir(fullfile(dataDir,['*.mat']));
getSubjString = @(x) x(1:end-10);
subjs         = cellfun(getSubjString,{subjFiles.name},'UniformOutput',false);

% process each subject
for s = 1:length(subjs);
    fprintf('Processing %s.\n',subjs{s});
    process_subj(subjs{s},dataDir,ana_name,saveDir)
end

function process_subj(subj,dataDir,ana_name,saveDir)

% load individual subject data
subjData       = load(fullfile(dataDir,[subj '_lasso.mat']));
powerAndParams = load(fullfile(dataDir,'power',[subj '_binnedPower.mat']));

% figure out unique object positions
[uniqueObjPos,~,objGroups] = unique(subjData.objLocs,'rows','stable');

% One row per timebin
% Y x Residuals, Y x Prediction, 2D error scatter, Reg. Info

figs = [];
nTimes = length(subjData.res);
for t = 1:nTimes
    
    % for each time bin, get timebin specific info
    time   = subjData.res(t).timeBin;
    yHat   = vertcat(subjData.res(t).yHat{:});
    mse    = vertcat(subjData.res(t).mse{:});
    lambda = subjData.res(t).lambda;
    df     = subjData.res(t).df;
    
    % mean yHat for each object position (bc two trials per object)
    yHat_meanByObj = grpstats(abs(abs(subjData.Y)-abs(yHat)),objGroups);
    
    % Make Figure
    close all
    figure('units','normalized','paperpositionmode','auto','position',[.2  .2  .5  .25]);
    titleStr = sprintf('%s: %s, Time = %d-%d, DF = %d, Lambda = %.4f',strrep(subj,'_',' '),ana_name,time(1),time(2),df,lambda);
    title(titleStr,'fontsize',14)
    axis off
        
    % Y x Prediction
    ax = axes('Position',[.08 .2 .25 .7]);
    h=scatter(ax,subjData.Y,yHat);
    ax.YLim = [0 1];
    ax.YTick = 0:.2:1;
    ax.XTick = 0:.2:1;
    grid on
    ax.GridLineStyle = ':';
    ax.GridColor = 'k';
    xlabel('Performance Factor','fontsize',22)
    ylabel('Prediction','fontsize',22) 
    ax.Clipping = 'off';
    ax.FontSize = 22;
    
    % Y x Residuals        
    ax = axes('Position',[.41 .2 .25 .7]);
    h = scatter(ax,subjData.Y,subjData.Y-yHat);
    ax.XTick = 0:.2:1;
    grid on
    ax.GridLineStyle = ':';
    ax.GridColor = 'k';
    xlabel('Performance Factor','fontsize',22)
    ylabel('Residual','fontsize',22) 
    ax.Clipping = 'off';
    ax.FontSize = 22;
    
    % 2D scatter prediction error
    ax = axes('Position',[.68 .56 .3 .4]);
    scatter(ax,uniqueObjPos(:,1),uniqueObjPos(:,2),100,yHat_meanByObj,'filled');
    ax.YLim = [-18 18];
    ax.XLim = [-32 32];   
    ax.XTick = [];
    ax.YTick = [];
    ax.Clipping = 'off';
    title('Pred. Error')
    h = colorbar;
    h.Location = 'northoutside';    
    box on
    
    % 2D scatter mse
    ax = axes('Position',[.68 .14 .3 .4]);
    scatter(ax,uniqueObjPos(:,1),uniqueObjPos(:,2),100,mse,'filled');
    ax.YLim = [-18 18];
    ax.XLim = [-32 32];    
    ax.XTick = [];
    ax.YTick = [];
    ax.Clipping = 'off';
    title('MSE')
    box on
    h = colorbar;
    h.Location = 'northoutside';
    
    figDir = fullfile(saveDir,strrep(ana_name,' ','_'));
    if ~exist(figDir,'dir')
        mkdir(figDir);
    end
    figName = sprintf('%s_timebin%d',subj,t);    
    figs(t).fname = fullfile(figDir,figName);
    print(figs(t).fname,'-depsc2','-loose')
end


texName = [subj '_' strrep(ana_name,' ','_') '_report.tex'];
write_texfile(figDir,texName,subj,figs)


curr_dir = pwd;
cd(figDir);
fprintf('Compiling pdf...\n');
unix(['pdflatex -shell-escape ' fullfile(figDir, texName)]);
unix(['rm ' texName(1:end-3) 'aux']);
unix(['rm ' texName(1:end-3) 'log']);
fprintf('Done!\n');
cd(curr_dir);


% Start making the tex file
function write_texfile(saveDir,texName, subj, figs)

% Write the document. If you do not have write permission, this will crash.
fid = fopen(fullfile(saveDir,texName),'w');

if fid==-1;
    error(sprintf('cannot open %s',texName))
end

% Write out the preamble to the tex doc. This is standard stuff and doesn't
% need to be changed
fprintf(fid,'\\documentclass[a4paper]{article} \n');
fprintf(fid,'\\usepackage[usenames,dvipsnames,svgnames,table]{xcolor}\n');
fprintf(fid,'\\usepackage{graphicx,multirow} \n');
fprintf(fid,'\\usepackage{epstopdf} \n');
fprintf(fid,'\\usepackage[small,bf,it]{caption}\n');
fprintf(fid,'\\usepackage{subfig,amsmath} \n');
fprintf(fid,'\\usepackage{wrapfig} \n');
fprintf(fid,'\\usepackage{longtable} \n');
fprintf(fid,'\\usepackage{pdfpages}\n');
fprintf(fid,'\\usepackage{mathtools}\n');
fprintf(fid,'\\usepackage{array}\n');
fprintf(fid,'\\usepackage{enumitem}\n');
fprintf(fid,'\\usepackage{sidecap} \\usepackage{soul}\n');

% fprintf(fid,'\\setlength\\belowcaptionskip{5pt}\n');
fprintf(fid,'\n');
fprintf(fid,'\\addtolength{\\oddsidemargin}{-.875in} \n');
fprintf(fid,'\\addtolength{\\evensidemargin}{-.875in} \n');
fprintf(fid,'\\addtolength{\\textwidth}{1.75in} \n');
fprintf(fid,'\\addtolength{\\topmargin}{-.75in} \n');
fprintf(fid,'\\addtolength{\\textheight}{1.75in} \n');
fprintf(fid,'\n');
fprintf(fid,'\\newcolumntype{C}[1]{>{\\centering\\let\\newline\\\\\\arraybackslash\\hspace{0pt}}m{#1}} \n');

fprintf(fid,'\\usepackage{fancyhdr}\n');
fprintf(fid,'\\pagestyle{fancy}\n');
fprintf(fid,'\\fancyhf{}\n');
fprintf(fid,'\\lhead{Report: %s }\n',strrep(subj,'_','\_'));
fprintf(fid,'\\rhead{Date created: %s}\n',date);

fprintf(fid,'\\usepackage{hyperref}\n');

% Start the document
fprintf(fid,'\\begin{document}\n\n\n');

% fprintf(fid,'\\hypertarget{%s}{}\n',region{1});

% This section writes the figures
for f = 1:length(figs)
    
    if f == 1
        fprintf(fid,'\\begin{figure}[!h]\n');
        fprintf(fid,'\\centering\n');
    end
    fprintf(fid,'\\includegraphics[width=0.99\\textwidth]{%s}\n',figs(f).fname);
    fprintf(fid,'\\vspace{.2in}\n')
    if f == length(figs)
%         fprintf(fid,'\\caption{\\textbf{test}.}\n');
        fprintf(fid,'\\end{figure}\n\n\n');
    
    end
end

fprintf(fid,'\\end{document}\n\n\n');










