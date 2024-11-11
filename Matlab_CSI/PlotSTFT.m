%ø… ”ªØ
function PlotSTFT(T,F,S,X)
    % Plots STFT
    plotOpts = struct();
    plotOpts.isFsnormalized = false;
    plotOpts.cblbl = getString(message('signal:dspdata:dspdata:MagnitudedB'));
    plotOpts.title = 'Short-time Fourier Transform';
    plotOpts.threshold = max(20*log10(abs(S(:))+eps))-X;
    signalwavelet.internal.convenienceplot.plotTFR(T,F,20*log10(abs(S)+eps),plotOpts);
end
%}
