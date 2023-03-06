% Optical simulation of dead pixel grids
%   Assumes scalar diffraction theory.
%   Uses 0.4mum pp for 640nm.

path = 'D:/1003_NRSH_CTC_Pipeline/nrsh/';
addpath(genpath(path))

dat = load('D:/04_HologramData/06_BCOM_REPOSingleChannel/s_bcom_dices4k_RNew.mat');
H = dat.bcom_dices4k_R;

pp = dat.pp;
wlen = dat.wlen;
z = dat.zrec;
N = size(H);
Npad = [2*N(1), N(2:end)];

p_period = 256;
p_inactive = 32;
titlestr = ['Used area: ' num2str((1-(p_inactive/p_period)).^2)];
disp(titlestr)

ftFun = @(x) fftshift(fft2(ifftshift(x)));
iftFun = @(x) ifftshift(ifft2(fftshift(x)));
cropFun = @(x) x(1:N(1), 1:N(2));
zoomFun = @(x) abs(x(1:3*p_period, 1:3*p_period));

bwFun = @(x) iftFun(padarray(ftFun(x), [size(x, 1), 0], 'post'));
binFun = @(x) real(x) > 0;
bbFun = @(x) binFun(bwFun(x));
ibbFun = @(x) iftFun(cropFun(ftFun(single(x))));


Hbin = bbFun(H);

nrsh_set = @(sstr) getSettings('ap_sizes', N, 'orthographic', false, 'h_pos', [0], 'v_pos', [0], 'dataset', 'mis_px_test', 'name_prefix', ['dice4k_mp_' sstr '_'], 'resize_fun', @(x) imresize(abs(x), 2048*[1,1], 'bilinear'), 'cfg_file', fullfile(path, 'config_files/bcom/', 'dices4kr_000.txt'));
% recFun = @(x, sstr) imresize(abs(ang_spectrum_nbl_fun3(x, pp, -z, wlen)), 2048*[1,1], "bilinear");
% recFunBin = @(x, sstr) imresize(abs(ang_spectrum_nbl_fun3(ibbFun(x), pp, -z, wlen)), 2048*[1,1], "bilinear");
recFun = @(x, sstr) nrsh(x, z, nrsh_set(sstr));
recFunBinRef = @(x, sstr) nrsh(ibbFun(x), z, nrsh_set(sstr));
recFunBin = @(x, sstr, cmin, cmax) nrsh(ibbFun(x), z, nrsh_set(sstr), 'exhaustive', {N}, 0, 0, cmin, cmax);
smFun = @(x) s_method(single(x(:, 1234)));

m = max(abs(H(:)))/2;
M0 = zeros(Npad);
M1 = ones(Npad);
MR = rand(Npad)>0.5; 
MC = M1;
for r=1:Npad(1)
    for c=mod(r+2,2)+1:2:N(2)
        MC(r,c)=0;
    end
end

[~, cmin, cmax] = recFunBinRef(Hbin, 'GT bin');

figure(1), clf
subplot(2,3,1), imagesc(recFun(H, 'GT')); title('GT'), RECformatFun()
subplot(2,3,2), imagesc(recFunBinRef(Hbin, 'GT bin')); title('GT bin'), RECformatFun()
subplot(2,3,3), imagesc(recFunBin(dh_pattern(Hbin, MR, p_period, p_inactive), 'Rand', cmin, cmax)); title('Rand'), RECformatFun()
subplot(2,3,4), imagesc(recFunBin(dh_pattern(Hbin, M0, p_period, p_inactive), 'Zero', cmin, cmax)); title('Zero'), RECformatFun()
subplot(2,3,5), imagesc(recFunBin(dh_pattern(Hbin, M1, p_period, p_inactive), 'Max', cmin, cmax)); title('Max'), RECformatFun()
subplot(2,3,6), imagesc(recFunBin(dh_pattern(Hbin, MC, p_period, p_inactive), 'Checker', cmin, cmax)); title('Checker'), RECformatFun()
colormap(gray(256))
sgtitle(titlestr)
print(gcf, 'Comparison_long.jpeg', '-djpeg', '-r300')

cl = [0, 1e3];
figure(2), clf
subplot(2,3,1), imagesc(smFun(iftFun(padarray(ftFun(H), [N(1), 0], 'post')))); title('GT'), clim([0, 1e7]), SMformatFun(N, Npad)
subplot(2,3,2), imagesc(smFun(Hbin)); title('GT bin'), clim(cl), SMformatFun(N, Npad)
subplot(2,3,3), imagesc(smFun(dh_pattern(Hbin, MR, p_period, p_inactive))); title('Rand'), clim(cl), SMformatFun(N, Npad)
subplot(2,3,4), imagesc(smFun(dh_pattern(Hbin, M0, p_period, p_inactive))); title('Zero'), clim(cl), SMformatFun(N, Npad)
subplot(2,3,5), imagesc(smFun(dh_pattern(Hbin, M1, p_period, p_inactive))); title('Max'), clim(cl), SMformatFun(N, Npad)
subplot(2,3,6), imagesc(smFun(dh_pattern(Hbin, MC, p_period, p_inactive))); title('Checker'), clim(cl), SMformatFun(N, Npad)
colormap parula
print(gcf, 'Comparison_long_sm.jpeg', '-djpeg', '-r300')


% figure(2), clf
% subplot(2,3,1), smFun(iftFun(padarray(ftFun(H), [N(1), 0], 'post'))); title('GT')
% subplot(2,3,2), smFun(Hbin); title('GT bin')
% subplot(2,3,3), smFun(dh_pattern(Hbin, MR, p_period, p_inactive)); title('Rand')
% subplot(2,3,4), smFun(dh_pattern(Hbin, M0, p_period, p_inactive)); title('Zero')
% subplot(2,3,5), smFun(dh_pattern(Hbin, M1, p_period, p_inactive)); title('Max')
% subplot(2,3,6), smFun(dh_pattern(Hbin, MC, p_period, p_inactive)); title('Checker')
% colormap parula
% print(gcf, 'Comparison_long_sm.jpeg', '-djpeg', '-r300')



figure(3), clf
subplot(2,3,1), imagesc(zoomFun(H)); title('GT'),ZOOMformatFun()
subplot(2,3,2), imagesc(zoomFun(Hbin)); title('GT bin'),ZOOMformatFun()
subplot(2,3,3), imagesc(zoomFun(dh_pattern(Hbin, MR, p_period, p_inactive))); title('Rand'),ZOOMformatFun()
subplot(2,3,4), imagesc(zoomFun(dh_pattern(Hbin, M0, p_period, p_inactive))); title('Zero'),ZOOMformatFun()
subplot(2,3,5), imagesc(zoomFun(dh_pattern(Hbin, M1, p_period, p_inactive))); title('Max'),ZOOMformatFun()
subplot(2,3,6), imagesc(zoomFun(dh_pattern(Hbin, MC, p_period, p_inactive))); title('Checker'),ZOOMformatFun()
colormap parula
print(gcf, 'Comparison_long_zoom.jpeg', '-djpeg', '-r300')

figure(4)
imagesc(abs(dh_pattern(Hbin, MR, p_period, p_inactive))), title('Random Mask'), RECformatFun()
colormap(gray(256))
print(gcf, 'Patterining_Random.jpeg', '-djpeg', '-r300')

% Replica of fig. 1 vor better visualization
figure(5)
subplot(2,2,3), imagesc(recFun(H, 'foo')); title('GT'), RECformatFun()
subplot(2,2,2), imagesc(recFunBinRef(Hbin, 'bar')); title('GT Bin'), RECformatFun()
subplot(2,2,1), imagesc(recFunBin(dh_pattern(Hbin, MC, p_period, p_inactive), 'baz', cmin, cmax)); title('Checker'), RECformatFun()
subplot(2,2,4), imagesc(recFunBin(dh_pattern(Hbin, MR, p_period, p_inactive), 'baz2', cmin, cmax)); title('Rand'), RECformatFun()
colormap(gray(256))
print(gcf, 'Comparison.jpeg', '-djpeg', '-r300')

function SMformatFun(N, Npad)
    xlabel('Samples')
    ylabel('Spat. Freq.')
    xl = xlim();
    xticks([1, mean(xl), xl(2)])
    xticklabels([1, N(1), Npad(1)])

    yl = ylim();
    yticks([1, mean(yl), yl(2)])
    yticklabels([-1, 0, 1])
    ax = gca();
    ax.FontSize = 16;
end

function ZOOMformatFun()
    xlabel('X')
    ylabel('Y')
    sgtitle('3x3 period zoom')

    xl = xlim();
    xticks([xl(1), mean(xl), xl(2)])

    yl = ylim();
    yticks([yl(1), mean(yl), yl(2)])

    ax = gca();
    ax.FontSize = 16;
end

function RECformatFun()
    axis off
    ax = gca();
    ax.FontSize = 16;
end