function hfig = DisplayRgbStimulusProtocol(protocol, varargin)

hfig = myparse(varargin,'hfig',[]);

if isempty(hfig) || ~ishandle(hfig),
  hfig = figure('color', 'white') ;
end

% Compute the signals
%[Yrgbo, to] = ComputeRgbStimulusProtocol(protocol) ;
[Yrgb, t] = computeRgbStimulusProtocol(protocol) ;
Yr = Yrgb(:,1) ;
Yg = Yrgb(:,2) ;
Yb = Yrgb(:,3) ;

% Clear the figure
clf(hfig);

%start to plot
n = numel(t) ;
if n < numel(Yr),
  warning('Protocol pulses defined require longer than the specified duration, cutting off at duration');
  Yr = Yr(1:n);
  Yg = Yg(1:n);
  Yb = Yb(1:n);
end
haxr = axes('Parent', hfig) ;
subplot(3,1,1,haxr) ;
protocolRL = line('Xdata',t, 'Ydata', Yr,'color','r','LineStyle','-','Parent',haxr);  %#ok<NASGU> 

haxg = axes('Parent', hfig) ;
subplot(3,1,2,haxg) ;
protocolGL = line('Xdata',t, 'Ydata', Yg,'color','g','LineStyle','-','Parent',haxg);  %#ok<NASGU> 

haxb = axes('Parent', hfig) ;
subplot(3,1,3,haxb) ;
protocolBL = line('Xdata',t, 'Ydata', Yb,'color','b','LineStyle','-','Parent',haxb);  %#ok<NASGU> 

stepStartSec = 0;

%plot steps start and stop line
for stepIndex = 1:length(protocol.stepNum) 
    %plot steps sstart and stop line
    stepStartSec = stepStartSec + protocol.duration(stepIndex)/1000 ;
    line(haxr, 'XData', [stepStartSec stepStartSec], 'YData', [0 1], 'Color', 'c', 'LineStyle','--');
    line(haxg, 'XData', [stepStartSec stepStartSec], 'YData', [0 1], 'Color', 'c', 'LineStyle','--');
    line(haxb, 'XData', [stepStartSec stepStartSec], 'YData', [0 1], 'Color', 'c', 'LineStyle','--');
end

set(haxr,'XLim',[0,t(end)+1]);
set(haxg,'XLim',[0,t(end)+1]);
set(haxb,'XLim',[0,t(end)+1]);
xlabel(haxb,'Time (s)');
ylabel(haxr,'Red LED');
ylabel(haxg,'Green LED');
ylabel(haxb,'Blue LED');
