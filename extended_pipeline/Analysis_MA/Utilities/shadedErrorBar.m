function varargout = shadedErrorBar(x, y, errBar, lineProps, transparent)
% Purpose 
% Makes a 2-d line plot with a pretty shaded error bar made
% using patch. Error bar plot is transparent by default.
%
% Inputs:
% x - vector of x values [optional, can be left empty]
% y - vector of y values or a matrix of n observations by m cases
%     where m has length(x);
% errBar - if a vector we draw symmetric errorbars. If it has a size
%          of [2,length(x)] then we draw asymmetric error bars with
%          row 1 being the upper bar and row 2 being the lower bar
% lineProps - [optional] if not empty determines line properties of
%             the main line. Default is 'b-'.
% transparent - [optional] if true make patch transparent. Default is true
%
% Outputs:
% H - structure with handles to plot objects
%
% Examples:
% y=randn(30,80); 
% x=1:size(y,2);
% shadedErrorBar(x,mean(y,1),std(y),'g-',0);
% shadedErrorBar([],mean(y,1),std(y),'g-',0);
% shadedErrorBar(x,mean(y,1),std(y),'g-');

    % Check input arguments
    if nargin<2
        error('Not enough input arguments')
    end
    
    if nargin<3
        errBar=zeros(size(y));
    end
    
    if nargin<4 || isempty(lineProps)
        lineProps='b-';
    end
    
    if nargin<5 || isempty(transparent)
        transparent=1;
    end
    
    % If no x vector, create one
    if isempty(x)
        x=1:length(y);
    end
    
    % If error bars are symmetrical, make them both positive and negative
    if size(errBar,1)==1
        errBar=[errBar;errBar];
    end
    
    % Plot main line
    H.mainLine=plot(x,y,lineProps);
    hold on
    
    % Plot error bar edges
    edgeColor=get(H.mainLine,'Color');
    patchSaturation=0.15; % How dark should the shaded area be?
    
    % Make patch
    yP=[y+errBar(1,:), fliplr(y-errBar(2,:))];
    xP=[x,fliplr(x)];
    
    % Remove NaN values
    xP = xP(~isnan(yP));
    yP = yP(~isnan(yP));
    
    H.patch=patch(xP,yP,1,'FaceColor',edgeColor,...
                 'EdgeColor','none',...
                 'FaceAlpha',patchSaturation);
    
    % Make pretty edges around the patch
    H.edge(1)=plot(x,y+errBar(1,:),'-','Color',edgeColor);
    H.edge(2)=plot(x,y-errBar(2,:),'-','Color',edgeColor);
    
    % Make transparent if requested
    if transparent
        set([H.edge,H.mainLine],'Color',edgeColor)
        alpha(H.patch,patchSaturation)
    end
    
    % Only output handles if requested
    if nargout>0
        varargout{1}=H;
    end
    
    hold off
end