% Author: Aayush Chhabra
% Time-Frequency Analysis

%% Analysis: Part I
clear; close all; clc;

load handel
v = y'/2; 

% Let's look at the song handel
fig = figure(1);
plot((1:length(v))/Fs,v); hold on;
xlabel('Time[sec]');
ylabel('Amplitude');
title('Signal of Interest v(n)');
print(fig, '-dpng', 'fig1')

% Now, let's set up some variables for 
% further analysis.

t2 = (1:length(v))/Fs; 
t = t2(1:end-1);
n = length(t2); L = t2(end);
k = (2*pi/L)*[0:n/2-1 -n/2:-1];
ks = fftshift(k);
S = v(1:end-1); 

% Now, let's look at the frequency content of this data.
fig = figure(2)
St = fft(S);
plot(ks, fftshift(abs(St)))
title("Frequency content of the song (Fourier Transform)")
xlabel('Fourier Modes');
print(fig, '-dpng', 'fig2')
% Gabor Transform - sliding window.
sampling_control = 0.1; % = 0.01(over), 1(under), 0.1(normal)
tslide = 0:sampling_control:9;
% Parameters for gabor window
a = 1;
b = 1; % higher b means thinner window and vice versa 

index = 1;
bvec = [1 10 20 40];
for b = bvec
    Sgt_spec = []; % Data collection for spectrogram.
    for center = tslide
        gabor = a*exp(-b*(t-center).^2);
        Sg = gabor.*S;
        Sgt = fft(Sg);
        Sgt_spec = [Sgt_spec; abs(fftshift(Sgt))];

        figure(3)
        subplot(3,1,1)
        plot(t, S,'b'); hold on;
        plot(t, gabor, 'k'); hold off;

        subplot(3,1,2)
        plot(t, Sg, 'r'); hold on;
        plot(t, gabor, 'k');hold off;

        subplot(3,1,3)
        plot(ks, abs(fftshift(Sgt)), 'b');
        pause(0.01);
    end
    % Let's now use the data we have collected to make a spectrogram.
    fig = figure(4);
    subplot(length(bvec)/2, 2, index);
    pcolor(tslide,ks,Sgt_spec.');
    shading interp;
    colormap(hot);
    title(strjoin(["b =",b]));
    xlabel('Time');
    ylabel('Frequency');
    index = index + 1;
end
print(fig, '-dpng', 'fig4')

%% Analysis: Part II
%% Piano
clear; close all; clc;

tr_piano=16; % record time in seconds 
y=audioread('music1.wav'); 
Fs=length(y)/tr_piano; 
plot((1:length(y))/Fs,y); 
xlabel('Time [sec]'); 
ylabel('Amplitude'); 
title('Mary had a little lamb (piano)');
%p8 = audioplayer(y,Fs); playblocking(p8);

% Let's set up some variables for analysis
n = length(y);
L = tr_piano;
t2 = linspace(0, L, n+1);
t=t2(1:n);
S = y';
St = fft(S);
k = (1/L)*[0:n/2-1 -n/2:-1];
ks = fftshift(k);

% Let's look at the frequency content
figure;
plot(ks, abs(fftshift(St)))

% Gabor Transform - sliding window.
sampling_control = .1; % = 0.01(over), 1(under), 0.1(normal)
tslide = 0:sampling_control:L

% Parameters for gabor window
a = 1;
b = 100; % higher b means thinner window and vice versa 
piano_notes=[]; % Data collection for piano notes.
Sgt_spec = []; % Data collection for spectrogram.
for center = tslide
    gabor = a*exp(-b*(t-center).^2);
    Sg = gabor.*S;
    Sgt = fft(Sg);
    [m, ind] = max(Sgt);
    piano_notes = [piano_notes; abs(k(ind))];
    Sgt_spec = [Sgt_spec; abs(fftshift(Sgt))];

    figure(3)
    subplot(3,1,1)
    plot(t, S,'b'); hold on;
    plot(t, gabor, 'k'); hold off;

    subplot(3,1,2)
    plot(t, Sg, 'r'); hold on;
    plot(t, gabor, 'k');hold off;

    subplot(3,1,3)
    plot(ks, abs(fftshift(Sgt)), 'b');
    pause(0.01);
end
% Let's now use the data we have collected to make a spectrogram.
fig = figure(4);
pcolor(tslide,ks,Sgt_spec.');
shading interp;
colormap(hot);
xlabel('Time');
ylabel('Frequency');
ylim([0 500]);
print(fig, '-dpng', 'piano_spect')
%%
% Let's generate the music score for the piano
fig = figure(5);
plot(tslide, piano_notes, '.');
yticks([220.00, 233.08, 246.94, 261.63, 277.18, 293.66, 311.13, 329.63, 349.23]);
yticklabels({'A','A#','B','C','C#','D', 'D#', 'E', 'F'});
ylim([200 350]);
title("Music Score (Piano)");
xlabel("Time"); 
ylabel("Notes");
print(fig, '-dpng', 'piano_score');
%% Recorder
clear; close all; clc;

tr_rec=14; 
% record time in seconds
y=audioread('music2.wav'); 
Fs=length(y)/tr_rec;
plot((1:length(y))/Fs,y); 
xlabel('Time [sec]'); 
ylabel('Amplitude'); 
title('Mary had a little lamb (recorder)'); 

% Let's set up some variables for analysis
n = length(y);
L = tr_rec;
t2 = linspace(0, L, n+1);
t=t2(1:n);
S = y';
St = fft(S);
k = (1/L)*[0:n/2-1 -n/2:-1];
ks = fftshift(k);

% Gabor Transform - sliding window.
sampling_control = .1; % = 0.01(over), 1(under), 0.1(normal)
tslide = 0:sampling_control:L;
% Parameters for gabor window
a = 1;
b = 100; % higher b means thinner window and vice versa 
recorder_notes = [];
Sgt_spec = []; % Data collection for spectrogram.
for center = tslide
    gabor = a*exp(-b*(t-center).^2);
    Sg = gabor.*S;
    Sgt = fft(Sg);
    Sgt_spec = [Sgt_spec; abs(fftshift(Sgt))];
    [m, ind] = max(Sgt);
    recorder_notes = [recorder_notes; abs(k(ind))];
    figure(3)
    subplot(3,1,1)
    plot(t, S,'b'); hold on;
    plot(t, gabor, 'k'); hold off;

    subplot(3,1,2)
    plot(t, Sg, 'r'); hold on;
    plot(t, gabor, 'k');hold off;

    subplot(3,1,3)
    plot(ks, abs(fftshift(Sgt)), 'b');
    pause(0.01);
end
% Let's now use the data we have collected to make a spectrogram.
fig = figure(6);
pcolor(tslide,ks,Sgt_spec.');
shading interp;
colormap(hot);
xlabel('Time');
ylabel('Frequency');
ylim([0 1500])
print(fig, '-dpng', 'recorder_spect')

%% Let's generate the music score for the recorder
fig = figure(7);
plot(tslide, recorder_notes, '.');
yticks([783.99, 830.61, 880.00, 932.33, 987.77, 1046.5, 1108.7, 1174.4]);
yticklabels({'G','G#','A','A#','B','C', 'C#', 'D'});
ylim([750, 1300]);
title("Music Score (Recorder)");
xlabel("Time"); 
ylabel("Notes");
print(fig, '-dpng', 'recorder_score');