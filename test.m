clear;
close all;
c = 340;                    % Sound velocity (m/s)
fs = 44100;                 % Sample frequency (samples/s)
r = [3 3 2];                % Receiver position [x y z] (m)
s = [-2 0.5 1];             % Source position [x y z] (m)
Ls = [4 4 6];                % Room dimensions [x y z] (m)
Lr = [4 4 6];                % Room dimensions [x y z] (m)
beta_s = 0.8;                 % Reverberation time (s)
beta_r = 0.8;                 % Reverberation time (s)
n = 4096;                   % Number of samples
mtype = 'omnidirectional';
order = -1;
dim =3;

r(1) = r(1)+Ls(1);
s(1) = s(1)+Ls(1);

mh = mexhost;
h_1 = StIM_rir_generator(c, fs, r, s, Lr, Ls, beta_r, beta_s, n, mtype, order, dim, 1);
plot(h_1);