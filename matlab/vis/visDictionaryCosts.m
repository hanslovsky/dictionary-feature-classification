function [ h, obsAxes, plotAxes, dictPatchAxes ] = visDictionaryCosts( xLR, D, szD, dsFactor, costs, indices, params  )
% Usage:
%   h = visDictionaryCosts( xLR, D, costs )
%
%   defaultParameters = visDictionaryCosts() 
%       returns default parameters
%  

%% parameter parsing

if( ~exist('xLR','var' ))
    xLR = [];
end
if( ~exist('D','var' ))
    D = [];
end
if( ~exist('szLR','var' ))
    szLR = [];
end
if( ~exist('dsFactor','var' ))
    dsFactor = 1;
end
if( ~exist('costs','var' ))
    cost = [];
end


if( ~exist('params','var' ))
    params = defaultVisDictionaryParameters();
end

if( isempty( xLR ))
    h = params;
    return;
end

if( isempty(params.fig_handle))
    h = figure( 'color', 'w' );
else
    h = params.fig_handle;
%     set( groot, 'CurrentFigure', params.fig_handle );
end

%% place observation axes

szIn = size( xLR );

sh = params.obsFrac - params.boundaryPad;
sw = 1 - params.boundaryPad;
px = params.boundaryPad./2;
py = 1 - params.boundaryPad./2 - sh;

obsAxes = axes('units','norm','pos',[px py sw sh]);


%% place plot axis
plotAxes = axes('units','norm','pos', ...
    [ params.boundaryPad./2 (py-params.imPad-params.costPlotFrac.*params.costFrac) ...
      (1-params.boundaryPad) (params.costPlotFrac.*params.costFrac) ]);

%% place dictionary patch axes

doPlotTwoDictRows = ( szIn(1) ~= szIn(2) && params.showHRdictWhenLR );
numDictPlotRows   = 1;

dictFrac = 1 - params.obsFrac - params.costFrac;
if( doPlotTwoDictRows )
    dictHeight = dictFrac / 2;
    numDictPlotRows = 2;
else
    dictHeight = dictFrac;
end

dictWidth = (( 1 - params.boundaryPad ) ./ params.numCostsToShow);
dictImPad = dictWidth .* params.imPad;
dictWidth = dictWidth - dictImPad;

dictPatchAxes = gobjects( params.numCostsToShow, numDictPlotRows );


px = params.boundaryPad./2;
for i = 1:params.numCostsToShow
    
    for j = 1:numDictPlotRows
        if( j == 1 )
            py = 0;
        else
            py = dictHeight;
        end
        
        dictPatchAxes(i,j) = axes('units','norm','pos', ...
            [ px py dictWidth (dictHeight-params.imPad./2) ] );
    end
    
    px = px + dictWidth + dictImPad;
end

%% plot everything

% plot cost
axes( plotAxes );
bar( costs( 1:params.numCostsToShow ), 'b' );
axbnds = axis;
axbnds(1) = 0.5;
axbnds(2) = params.numCostsToShow + 0.5;
textOffset = params.textOffsetMultiplier .* axbnds(4);
axis( axbnds );
for i = 1:params.numCostsToShow
    if( params.rotateDictIdx )
        text( i, (-2.25.*textOffset), sprintf('%04d', indices(i) ), ...
            'color', [0.3 0.3 1], 'fontsize', 12, 'rotation', 90); 
    else
        text( i-0.4, ( -textOffset-(textOffset*mod(i,2))), sprintf('%04d', indices(i) ), ...
       'color', [0.3 0.3 1], 'fontsize', 12); 
    end
   
end

% render observations
axes( obsAxes );

if( size(xLR,1) ~= size(xLR,2) )
    M = min( size( xLR ));
    upInds = reshape( repmat( 1:M, dsFactor, 1 ), [] , 1 );
    pVec = find( size(xLR) == M );
    pVec = [ pVec setdiff( 1:2, pVec ) ];
    patch = permute( xLR, pVec );
    patch = patch( upInds, : );
else
    patch = xLR;
end

imdisp( patch );


% % plot dictionary patches
for i = 1:params.numCostsToShow
    
    axes( dictPatchAxes(i,1) );
    dict_im = reshape( D(i,:), szD(1:2) );
    
    if( szIn(1) ~= szIn(2 ))
        
        if( doPlotTwoDictRows )
            % plot the high res patch if we ask for it
            imdisp( dict_im );
            axes( dictPatchAxes(i,2) );
%             pause;
        end
        
        % now plot the low res patch
        dict_im = PatchConstraints.downsamplePatch( dict_im, szD, dsFactor  );
        
        M = min( size( dict_im ));
        upInds = reshape( repmat( 1:M, dsFactor, 1 ), [] , 1 );
        pVec = find( size(xLR) == M );
        pVec = [ pVec setdiff( 1:2, pVec ) ];
        patch = permute( dict_im, pVec );
        dict_im = patch( upInds, : );
    end
    
    % if theres no downsampling, this will just plot the original
    % dictionary patch
    imdisp( dict_im );
%     pause;
end


end

function defaultParameters = defaultVisDictionaryParameters()

defaultParameters.boundaryPad = 0.05;
defaultParameters.imPad = 0.05;
defaultParameters.obsFrac = 0.5;
defaultParameters.costFrac = 0.25;
defaultParameters.costPlotFrac = 0.8;
defaultParameters.numCostsToShow = -1;
defaultParameters.textOffsetMultiplier = 0.15;
defaultParameters.showHRdictWhenLR = true;
defaultParameters.rotateDictIdx = true;

defaultParameters.fig_handle = [];

end
