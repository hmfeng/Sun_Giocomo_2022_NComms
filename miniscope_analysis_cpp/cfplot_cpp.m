function cfplot_cpp(x, y, condition, figname)
%cfplot_cpp plots curve fitting (linear regression) results of CPP scores
%vs. place field properties
%   Detailed explanation goes here
if ~exist('figname', 'var') || isempty(figname)
    figname = 'Figure'
end
mdl = fitlm(x, y);
switch condition
    case{'peak_shift', 'shift', 'Peak_shift'}
    %plot for correlation on test1
    [xData, yData] = prepareCurveData( x, y );
    % Set up fittype and options.
    ft = fittype( 'poly1' );
    % Fit model to data.
    [fitresult, gof] = fit( xData, yData, ft );
    % Plot fit with data.
    figure( 'Name', figname );
    h = plot( fitresult, xData, yData, 'o');
    legend( h, 'Raw data', 'Fit', 'Location', 'SouthEast' );
    text(-0.75,200, ['R2 = ', num2str(gof.rsquare)]);
    text(-0.75,220, ['P = ', num2str(mdl.Coefficients.pValue(end))])
    % Label axes
    xlabel('Delta Place field shift (cm)')
    ylabel('CPP score')
    axis square
    
    case{'correlation', 'corr', 'Corr', 'Correlation'}
    %plot for correlation on test1
    [xData, yData] = prepareCurveData( x, y );
    % Set up fittype and options.
    ft = fittype( 'poly1' );
    % Fit model to data.
    [fitresult, gof] = fit( xData, yData, ft );
    % Plot fit with data.
    figure( 'Name', figname );
    h = plot( fitresult, xData, yData, 'o');
    legend( h, 'Raw data', 'Fit', 'Location', 'SouthEast' );
    text(0,100, ['R2 = ', num2str(gof.rsquare)]);
    text(0,50, ['P = ', num2str(mdl.Coefficients.pValue(end))])
    % Label axes
    xlabel('Correlation')
    ylabel('CPP score')
    axis square

    case{'pv_corr', 'PV_corr', 'pv_correlation', 'PV_correlation'}
    %plot for correlation on test1
    [xData, yData] = prepareCurveData( x, y );
    % Set up fittype and options.
    ft = fittype( 'poly1' );
    % Fit model to data.
    [fitresult, gof] = fit( xData, yData, ft );
    % Plot fit with data.
    figure( 'Name', figname );
    h = plot( fitresult, xData, yData, 'o');
    legend( h, 'Raw data', 'Fit', 'Location', 'SouthEast' );
    text(0,100, ['R2 = ', num2str(gof.rsquare)]);
    text(0,50, ['P = ', num2str(mdl.Coefficients.pValue(end))])
    % Label axes
    xlabel('PV correlation')
    ylabel('CPP score')
    axis square

end

end

