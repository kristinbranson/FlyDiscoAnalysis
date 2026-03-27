function	[dip, p_value, xlow,xup]=HartigansDipSignifTest(xpdf,nboot)

%  function	[dip,p_value, xlow,xup]=HartigansDipSignifTest(xpdf,nboot)
%
% Calculates Hartigan's DIP statistic and its significance for the empirical
% p.d.f  XPDF (vector of sample values).
%
% This routine calls the matlab routine 'HartigansDipTest' that actually
% calculates the DIP. NBOOT is the user-supplied sample size of boot-strap.
% Code by F. Mechler (27 August 2002)
%
% Source: VH-Lab/vhlab-toolbox-matlab (https://github.com/VH-Lab/vhlab-toolbox-matlab)
%   Original file: stats/hartigansdipsigniftest.m
%   License: GPL-3.0 (https://github.com/VH-Lab/vhlab-toolbox-matlab/blob/master/LICENSE)
%   Retrieved: 2026-03-05

% calculate the DIP statistic from the empirical pdf
[dip,xlow,xup, ifault, gcm, lcm, mn, mj]=HartigansDipTest(xpdf);
N=length(xpdf);

% calculate a bootstrap sample of size NBOOT of the dip statistic for a uniform pdf of sample size N (the same as empirical pdf)
boot_dip=zeros(nboot,1);
for i=1:nboot
   unifpdfboot=sort(unifrnd(0,1,1,N));
   [unif_dip]=HartigansDipTest(unifpdfboot);
   boot_dip(i)=unif_dip;
end;
boot_dip=sort(boot_dip);
p_value=sum(dip<boot_dip)/nboot;
