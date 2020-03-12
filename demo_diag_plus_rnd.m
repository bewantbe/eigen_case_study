% 
pic_common_include;

rand('state', 32466345);

n = 58;

A=[0.01*eye(n/2) zeros(n/2); 0.01*randn(n/2) eye(n/2)+0.01*randn(n/2)];
[eigvec, eigval] = eig(A);

pic_output_png = @(st) print('-dpng', ...
                   sprintf('%s%s.png', pic_prefix, st));

figure(1); imagesc(log10(abs(A))); colorbar();
pic_output_png('AQei0');

figure(2); imagesc(log10(abs(eigvec))); colorbar();
pic_output_png('AQei0-vec');

figure(3); plot(abs(diag(eigval)));
pic_output_png('AQei0-val');

