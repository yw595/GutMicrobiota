function makeBarGut(xvals,yvals,titleString,outputDir,ylabelString,isBoxPlot,isScatter,xlabels,labels,specialVals,indLines,xlabelString,xlabelsFontSize,specialValsFontSize,labelsFontSize,ylabelFontSize,xlabelFontSize,titleFontSize)
    
    %p = inputParser;
    %addParameter(p,'isBoxPlot',0);
    %addParameter(p,'isScatter',0);
    %addParameter(p,'xlabels',{});
    %addParameter(p,'labels',{});
    %addParameter(p,'specialVals',{});
    %addParameter(p,'indLines',[]);
    %addParameter(p,'xlabelsFontSize',10);
    %addParameter(p,'specialValsFontSize',10);
    %addParameter(p,'labelsFontSize',20);
    
    % multiply by three to give some space between labels
    figure('Visible','off');
    %INSIDIOUS BUG, here and below, default xlim generates some
    %sort of overhang, appears constant wrt paper size, so could
    %really squeeze data columns, also appears bigger at right, add
    %one to length(xvals) here to compensate for tight xlim below
    width=(length(xvals)+1)*10/72*3;
    if strcmp(titleString,'Complex_Diff_Coverage')
        width=(length(xvals)+1)*10/72*30;
    end
    if width>42
        width=41;
    end
    height = 8.5;
    maxSpecialValLen = 0; maxXLabelLen = 0;
    if ~isempty(xlabels) || ~isempty(specialVals)
        if ~isempty(specialVals)
            maxSpecialValLen = max(cellfun(@(x) length(x),values(specialVals)));
            height = height+max(maxSpecialValLen*10/72*19/26-8.5,0);
        end
        if ~isempty(xlabels)
            set(gca,'XTickLabel',{});
            maxXLabelLen = max(cellfun(@(x) length(x),xlabels));
            height = height+maxXLabelLen*10/72*19/26;
        end
    end
    if ~isempty(labels)
        uniqLabels = {};
        for j=1:length(labels)
            if ~any(strcmp(labels{j},uniqLabels))
                uniqLabels{end+1}=labels{j};
            end
        end
        
        for j=1:length(uniqLabels)
            maxOffset = 0;
            dispLabel = uniqLabels{j};
            temp = find(strcmp(labels,dispLabel));
            labelBegin = temp(1);
            labelEnd = temp(end);
            maxCharLength = floor((labelEnd-labelBegin+1)*3/2);
            if ~isempty(regexp(titleString,'BLAST'))
                %disp(length(dispLabel));
                %disp(maxCharLength);
            end
            while length(dispLabel)/maxCharLength > 1;
                maxOffset =maxOffset+1;
                dispLabel = dispLabel(maxCharLength+1:end);
            end
            if ~isempty(regexp(titleString,'BLAST'))
                %disp(dispLabel)
                %disp(maxOffset);
            end
        end
        %disp(maxOffset)
        height = height+20/72*maxOffset;
    end
    set(gcf,'PaperUnits','inches');
    papersize = get(gcf,'PaperSize');
    left = (papersize(1)-width)/2;
    bottom = (papersize(2)-height)/2;
    myfiguresize = [left,bottom,width,height];
    set(gcf,'PaperPosition',myfiguresize);
    
    %addParameter(p,'ylabelFontSize',20);
    %addParameter(p,'titleFontSize',20);
    %addParameter(p,'xlabelFontSize',20);
    
    if ~isempty(specialVals)
        origyvals = yvals;
        specialValsKeys = keys(specialVals);
        if isnumeric(specialValsKeys{1})
            specialValsKeys = cell2mat(specialValsKeys);
        end
        yvals(ismember(yvals,specialValsKeys))=0;
    end

    if isScatter
        scatter(xvals,yvals,[],'b','MarkerFaceColor','b','SizeData',40);
    elseif isBoxPlot
        boxplot(yvals, xlabels);
    else
        barH=bar(xvals,yvals,'b','FaceColor','b');
        if size(yvals,2)==2
            set(barH(2),'FaceColor','r');
        end
    end
    ylabel(ylabelString,'FontSize',20,'FontName','FixedWidth');
    if exist('xlabelString','var')
        xlabel(xlabelString,'FontSize',20,'FontName','FixedWidth');
    end
    % related to width setting above, set really tight bounds to
    % avoid large overhang distortion for few xvals
    if isScatter
        xlim([0 max(xvals(:))]);
    elseif ~isBoxPlot
        xlim([0.5 length(xvals)+0.5]);
    end
    if all(yvals==yvals(1))
        titleString = ['WARNING_UNIFOM_VALS_' titleString];
        yLimData = [yvals(1) yvals(1)+.01];
    else
        yLimData = [min(yvals(:)) max(yvals(:))];
    end
    %disp(width)
    if isScatter
        ylim([0 yLimData(2)]);
    elseif ~isBoxPlot
        if yLimData(1)<0
            ylim([yLimData(1) yLimData(2)]);
        else
            ylim([0 yLimData(2)]);
        end
    end
    temp = ylim();
    yLimDataTop = temp(2)/(temp(2)-temp(1))*8.5;
    if temp(1)~=0
        yLimMax = max(temp(2),maxSpecialValLen*(10/72)*(19/26)/yLimDataTop*temp(2));
    else
        yLimMax = max(temp(2),maxSpecialValLen*10/72/yLimDataTop*temp(2));
    end
    if ~isempty(labels)
        yLimMax = yLimMax+maxOffset*20/72/yLimDataTop*temp(2);
    end
    yLimDataBottom = temp(1)/(temp(2)-temp(1))*8.5;
    if temp(1)~=0
        yLimMin = min(temp(1),-maxXLabelLen*(10/72)*(19/26)/yLimDataBottom*temp(1)+yLimDataBottom);
    else
        yLimMin = min(temp(1),-maxXLabelLen*10/72/yLimDataTop*temp(2)+yLimDataBottom);
    end
    if temp(1)~=0
        height = height-yLimDataBottom;
        set(gcf,'PaperUnits','inches');
        papersize = get(gcf,'PaperSize');
        left = (papersize(1)-width)/2;
        bottom = (papersize(2)-height)/2;
        myfiguresize = [left,bottom,width,height];
        set(gcf,'PaperPosition',myfiguresize);
    end
    if ~isempty(regexp(titleString,'_Diff_Coverage'))
    %if ~isempty(regexp(titleString,'Complex_I vs cytochrome_bd_I_oxidase'))
        disp(xvals)
        disp(yvals)
        %disp(maxOffset)
        %disp(maxSpecialValLen)
        disp(temp)
        disp(height)
        disp(yLimDataTop)
        disp(yLimMax)
        disp(maxXLabelLen)
        disp(yLimDataBottom)
        disp(yLimMin)
        disp(xlim)
        disp(width)
        %absurd = absurd+1;
    end
    ylim([yLimMin yLimMax]);
    title(titleString,'FontSize',20,'FontName','FixedWidth');
    if ~isempty(xlabels)
        for j=1:length(xlabels)
            tempfont = 10;
            if strcmp(titleString,'Complex_Diff_Coverage')
                tempfont = 25;
            end
            text(j,yLimMin,xlabels{j},'Rotation',90,'FontSize',tempfont,'FontName','FixedWidth');
        end
    end

    if ~isempty(labels)
        prevLabel = labels{1};
        for j=1:length(labels)
            if ~strcmp(labels{j},prevLabel)
                prevLabel = labels{j};
                line([j-.5 j-.5], [yLimMin yLimMax], 'LineWidth', 2, 'Color', 'r'); 
            end
        end

        for j=1:length(uniqLabels)
            dispLabel = uniqLabels{j};
            temp = find(strcmp(labels,dispLabel));
            labelBegin = temp(1);
            labelEnd = temp(end);
            maxCharLength = floor((labelEnd-labelBegin+1)*3/2);
            labelOffset = 0;
            while length(dispLabel) > maxCharLength
                dispLabelPiece = dispLabel(1:maxCharLength);
                text(labelBegin-.5,yLimMax-(20/72)/(height*yLimMax/(yLimMax-yLimMin))*yLimMax*labelOffset, dispLabelPiece,'FontSize',20,'FontName','FixedWidth');
                dispLabel = dispLabel(maxCharLength+1:end);
                labelOffset = labelOffset+1;
            end
            dispLabelPiece = dispLabel;
            text(labelBegin-.5,yLimMax-(20/72)/(height*yLimMax/(yLimMax-yLimMin))*yLimMax*labelOffset, dispLabelPiece,'FontSize',20,'FontName','FixedWidth');
        end
    end
    if ~isempty(specialVals)
        for j=1:length(origyvals)
            if isKey(specialVals,origyvals(j))
                text(j,0,specialVals(origyvals(j)),'Rotation',90,'FontSize',10,'FontName','Courier');
            end
        end
    end
    
    if ~isempty(indLines)
        for j=1:length(indLines)
            line([indLines(j) indLines(j)], [yLimMin yLimMax], 'LineWidth', 2, 'Color', 'g');
        end
    end
    
    line([0 length(xvals)], [0 0]);
    
    saveas(gcf,[outputDir filesep titleString '.png']);
    close(gcf);
end