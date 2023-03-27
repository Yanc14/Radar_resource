function [out1] = myspecgramnew(data,window,nfft,shift)
%%对Data时频分析
%window为窗口长度
%shift为窗口移动步长
%nfft为FFT点数，一般等于窗口长度
%%
N = floor((length(data)-window-1)/shift); % 104
for i=1:N
    %每一列为一次窗口分析时频分析
    tmp =(fft(data((i-1)*shift+1:(i-1)*shift+window).'.*hann(window),nfft));
    out1(:,i) = tmp;
end
