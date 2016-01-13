function makeBarGut(xvals, yvals, titleString, xlabels, ylabelString, outputDir,labels,specialVals)
    % multiply by three to give some space between labels
    figure('Visible','off');
    %INSIDIOUS BUG, here and below, default xlim generates some
    %sort of overhang, appears constant wrt paper size, so could
    %really squeeze data columns, also appears bigger at right, add
    %one to length(xvals) here to compensate for tighg xlim below
    width=(length(xvals)+1)*10/72*3;
    height = (max(cellfun(@(x) length(x),xlabels)) + max(cellfun(@(x) length(x),values(specialVals))))*10/72;
    set(gcf,'PaperUnits','inches');
    papersize = get(gcf,'PaperSize');
    left = (papersize(1)-width)/2;
    bottom = (papersize(2)-height)/2;
    myfiguresize = [left,bottom,width,height];
    set(gcf,'PaperPosition',myfiguresize);
    
    origyvals = yvals;
    specialValsKeys = keys(specialVals);
    if isnumeric(specialValsKeys{1})
        specialValsKeys = cell2mat(specialValsKeys);
    end
    yvals(ismember(yvals,specialValsKeys))=0;

    barH=bar(xvals,yvals,'b','FaceColor','b');
    if size(yvals,2)==2
        set(barH(2),'FaceColor','r');
    end
    ylabel(ylabelString,'FontSize',20);
    % related to width setting above, set really tight bounds to
    % avoid large overhang distortion for few xvals
    xlim([0.5 length(xvals)+0.5]);
    if all(yvals==yvals(1))
        titleString = ['WARNING_UNIFOM_VALS_' titleString];
        ylim([yvals(1) yvals(1)+.01]);
    else
        ylim([-1*max(yvals(:)) 1.5*max(yvals(:))]);
    end
    temp = ylim(); yLimMin = temp(1); yLimMax = temp(2);
    title(titleString,'FontSize',20);
    for j=1:length(xlabels)
        text(j,-1*max(yvals(:)),xlabels{j},'Rotation',90,'FontSize',10);
    end

    if ~isempty(labels)
        prevLabel = labels{1};
        for j=1:length(labels)
            if ~strcmp(labels{j},prevLabel)
                prevLabel = labels{j};
                line([j+.5 j+.5], [yLimMin yLimMax], 'LineWidth', 2, 'Color', 'r'); 
            end
        end

        uniqLabels = {};
        for j=1:length(labels)
            if ~any(strcmp(labels{j},uniqLabels))
                uniqLabels{end+1}=labels{j};
            end
        end

        for j=1:length(uniqLabels)
            dispLabel = uniqLabels{j};
            temp = find(strcmp(labels,dispLabel));
            labelBegin = temp(1);
            labelEnd = temp(end);
            maxCharLength = (labelEnd-labelBegin+1)*3;
            labelOffset = 0;
            while length(dispLabel) > maxCharLength
                dispLabelPiece = dispLabel(1:maxCharLength);
                text(labelBegin,yLimMax*0.9-(10/72)/(height*yLimMax/(yLimMax-yLimMin))*yLimMax*labelOffset, dispLabelPiece,'FontSize',10);
                dispLabel = dispLabel(maxCharLength+1:end);
                labelOffset = labelOffset+1;
            end
            dispLabelPiece = dispLabel;
            text(labelBegin,yLimMax*0.9-(10/72)/(height*yLimMax/(yLimMax-yLimMin))*yLimMax*labelOffset, dispLabelPiece,'FontSize',10);
        end
    end

    for j=1:length(origyvals)
        if isKey(specialVals,origyvals(j))
            text(j,0,specialVals(origyvals(j)),'Rotation',90,'FontSize',10);
        end
    end
    
    set(gca,'XTickLabel',{});
    line([0 length(xlabels)], [0 0]);
    
    saveas(gcf,[outputDir filesep titleString '.png']);
    close(gcf);
end