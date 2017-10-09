function fig = plotIDAcurve(obj,results,plotMode)
% PLOTIDACURVE Plot results data in IDA format
%
%   `plotMode` can be 'single' or 'multiple', with 'single' plotting everything
%       in monochrome on a single plot and 'multiple' generating enough figures
%       to plot using the default color scheme. Thus, 'multiple' can actually
%       output a single plot, though not in monochrome.
%

nDistinctColors = size(unique(lines,'rows'), 1);    % Length of `lines` color palette

if nargin == 2
    if length(results.gm) > nDistinctColors
        plotMode = 'single';
    else
        plotMode = 'multiple';
    end
end

switch lower(plotMode)
case 'single'
    fig = figure;
    hold on
    for gmIndex = 1:obj.optionsIDA.nMotions
        plot([0 results.gm(gmIndex).maxDriftRatio*100],[0 results.gm(gmIndex).ST],'ko-');
        plot(results.gm(gmIndex).maxDriftRatio(end-1:end)*100,results.gm(gmIndex).ST(end-1:end),'k.','MarkerSize',20)
    end
    h(1) = plot(xlim,[results.SCT_hat results.SCT_hat],'--','LineWidth',2);
    h(2) = plot(xlim,[results.SMT     results.SMT],    '--','LineWidth',2);
    legend(h,{'$\hat{S}_{CT}$','$S_{MT}$'},'Location','southeast','Interpreter','latex')
    xlim([0 4*obj.optionsIDA.collapseDriftRatio*100])
    xlabel('Maximum Story Drift Ratio (\%)','Interpreter','latex')
    ylabel('Ground Motion Intensity $S_T$ (g)','Interpreter','latex')
    grid on

case 'multiple'
    nFigs = ceil(obj.optionsIDA.nMotions/nDistinctColors);
    fig = gobjects(nFigs,1);
    for iFig = 1:nFigs
        fig(iFig) = figure;
        hold on
        gmStart = (iFig-1)*nDistinctColors + 1;
        nLines = min(obj.optionsIDA.nMotions - (gmStart-1),nDistinctColors);
        gmEnd = gmStart + nLines - 1;
        GMs = gmStart:gmEnd;
        legendentries = cell(nLines+2,1);
        lineActuals = gobjects(nLines+2,1);
        for iLine = 1:nLines
            iGM = GMs(iLine);
            lineActuals(iLine) = plot([0 results.gm(iGM).maxDriftRatio*100],[0 results.gm(iGM).ST],'o-');
        end
        for iLine = 1:nLines
            iGM = GMs(iLine);
            plot(results.gm(iGM).maxDriftRatio(end-1:end)*100,results.gm(iGM).ST(end-1:end),'.','MarkerSize',20,'Color',lineActuals(iLine).Color)
            legendentries{iLine} = results.gm(iGM).ID;
        end
        lineActuals(end-1) = plot(xlim,[results.SCT_hat results.SCT_hat],'b--','LineWidth',2);
        lineActuals(end)   = plot(xlim,[results.SMT     results.SMT],    'r--','LineWidth',2);
        legendentries{end-1} = '$\hat{S}_{CT}$';
        legendentries{end}   = '$S_{MT}$';
        legend(lineActuals,legendentries,'Interpreter','latex')
        xlim([0 4*obj.optionsIDA.collapseDriftRatio*100])
        xlabel('Maximum Story Drift Ratio (\%)','Interpreter','latex')
        ylabel('Ground Motion Intensity $S_T$ (g)','Interpreter','latex')
        grid on;
    end

end

if nargout == 0
    clear fig
end

end
