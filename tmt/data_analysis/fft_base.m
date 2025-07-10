function x = fft_base(N)

% Generates a frequency line to go with an N-point fft.  Frequencies are in
% cycles per sample, i.e. they go from about -1/2 to about 1/2.
% Multiply by fs to convert to frequency in Hz.

hi_x_sample_index=ceil(N/2);
x_pos=linspace(0,hi_x_sample_index-1,hi_x_sample_index)';
x_neg=linspace(-(N-hi_x_sample_index),...
               -1,...
               N-hi_x_sample_index)';
x=(1/N)*[x_pos ; x_neg];

end
