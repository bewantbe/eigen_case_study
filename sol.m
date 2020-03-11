% 
pic_common_include;

rand('state', 32466345);

id_case = 2;
switch id_case
case 0

  fc = load('connectivity.mat');
  A = fc.W_p;
%  A = A(7:end, 7:end);
  
  n = length(A);

case -2

  fc = load('connectivity.mat');
  A0 = fc.W_p;
  A = A0;
  n = length(A);
  % remove off-2diagonal
  dd = arrayfun(@(x)ones(2), 1:n/2, 'UniformOutput', false);
  A(blkdiag(dd{:})==0) = 0;

case -1
  % for matrix examination

  pic_output_png = @(st) print('-dpng', ...
                   sprintf('%s%s_orig.png', pic_prefix, st));

  fc = load('connectivity.mat');
  A0 = fc.W_p;
  n = length(A0);

  A3 = A0([1:2:n 2:2:n], [1:2:n 2:2:n]);
  A11 = A3(1:n/2, 1:n/2);
  A21 = A3(n/2+1:n, 1:n/2);

%  figure(45326); plot(A21(~eye(n/2)), A11(~eye(n/2)), '+')
  nondiag11 = A11(~eye(n/2));
  nondiag21 = A21(~eye(n/2));
  
  r11 = [min(nondiag11(nondiag11>0)), max(nondiag11)];
  r21 = [min(nondiag21(nondiag21>0)), max(nondiag21)];
  figure(45327); loglog(nondiag21, nondiag11, '+', r21, r11)
  xlabel('A21 diag');
  ylabel('A11 diag');
  legend('orig', 'fit', 'location', 'southeast');
  pic_output_png('A11diag_vs_A21diag_loglog');

  return

case 1

  n=28*2;

  % fill 2x2 blocks
  A = zeros(n, n);
  for k=1:n/2
    A(2*(k-1)+1:2*k, 2*(k-1)+1:2*k) = [0.09 -0.1; 1 -1] * (1+1e-2*k);
  end

  % add random interconnection
  mnz = floor(n*n*0.05);  % number of non-zero
  idnz = 2*randi(n*n/2,mnz,1);
  A(idnz) = A(idnz) + 0.01*randn(mnz,1);

  % random scale
  A(1:2:n-1, :) = 10.^(-3+2*rand(n/2,1)) .* A(1:2:n-1, :);

case 3
  % structured shuffle 'connectivity.mat'
  
  fc = load('connectivity.mat');
  A0 = fc.W_p;
  n = length(A0);

  A3 = A0([1:2:n 2:2:n], [1:2:n 2:2:n]);
  A11 = A3(1:n/2, 1:n/2);
  A21 = A3(n/2+1:n, 1:n/2);

  rid = randperm(n/2*(n/2-1));
  elem_nondiag = A11(eye(n/2)==0);
  A11(eye(n/2)==0) = elem_nondiag(rid);
  elem_nondiag = A21(eye(n/2)==0);
  A21(eye(n/2)==0) = elem_nondiag(rid);
  
  A3(1:n/2, 1:n/2)   = A11;
  A3(n/2+1:n, 1:n/2) = A21;
  
  A([1:2:n 2:2:n], [1:2:n 2:2:n]) = A3;

case 4
  % structured shuffle 'connectivity.mat'
  
  fc = load('connectivity.mat');
  A0 = fc.W_p;
  n = length(A0);

  A3 = A0([1:2:n 2:2:n], [1:2:n 2:2:n]);
  A21 = A3(n/2+1:n, 1:n/2);

  rid = randperm(n/2*(n/2-1));
  elem_nondiag = A21(eye(n/2)==0);
  A21(eye(n/2)==0) = elem_nondiag(rid);
  
  A3(n/2+1:n, 1:n/2) = A21;
  
  A([1:2:n 2:2:n], [1:2:n 2:2:n]) = A3;

case {2, 22}

  fc = load('connectivity.mat');
  A0 = fc.W_p;
  n = length(A0);

  figure(10086);
  imagesc(A0);
  title('A0');
  caxis([-1, 1])
  colorbar();
  
  A3 = A0([1:2:n 2:2:n], [1:2:n 2:2:n]);
  figure(10087);
  imagesc(A3);
  title('A3');
  caxis([-1, 1])
  colorbar();

  A1 = A0(7:end, 7:end);
  n0 = length(A1);
  A_local = [mean(A1(1:2*n0+2:n0*n0)) mean(A1(n0+1:2*n0+2:n0*n0));
             mean(A1(2:2*n0+2:n0*n0)) mean(A1(n0+2:2*n0+2:n0*n0))];

  % fill 2x2 blocks
  A = zeros(n, n);
  eigvalapp = zeros(2, n/2);
  for k=1:n/2
    twob = A_local .* [(0.83+0.6*k/n) 1; (0.89+0.4*k/n) 1];
    A(2*(k-1)+1:2*k, 2*(k-1)+1:2*k) = twob;
    eigvalapp(1, k) = twob(2,2) + twob(2,1)*twob(1,2)/twob(2,2);
    eigvalapp(2, k) = twob(1,1) - twob(2,1)*twob(1,2)/twob(2,2);
  end

  A_skel = A;
  [eigvec_skel, eigval_skel] = eig(A);

%  % add random interconnection
%  mnz = floor(n*n*0.02);  % number of non-zero
%  idnz = 2*randi(n*n/2,mnz,1);
%  A(idnz) = A(idnz) + 0.001*randn(mnz,1);

  % add random interconnection
  mnz = floor(n*n*0.02);  % number of non-zero
  id_col = 2*randi(n/2,mnz,1)-1;   % should be odd
  id_row = 2*randi(n/2,mnz,1);     % should be even
  % remove diagonal noise
  delid = find((id_col+1) == id_row);
  id_col(delid) = [];
  id_row(delid) = [];
  mnz = length(id_col);
  % convert to serial index
  idnz = id_row + (id_col-1) * n;
  idnz1 = id_row-1 + (id_col-1) * n;
  
  %scaled_rnds = rand(mnz,1) .* 10 .^ (-2+2*rand(mnz,1));
  scaled_rnds = rand(mnz,1) .* 10 .^ (-2.3+2*rand(mnz,1));
  A(idnz) = A(idnz) + scaled_rnds;
  A(idnz1) = A(idnz1) + scaled_rnds/8;

  if id_case == 22
    % diagonalize 2x2 blocks
    A = eigvec_skel \ A * eigvec_skel;
  end

%  % random scale
%  A(1:2:n-1, :) = 10.^(-3+2*rand(n/2,1)) .* A(1:2:n-1, :);
end

pic_output_png = @(st) print('-dpng', ...
                   sprintf('%s%s_c%d.png', pic_prefix, st, id_case));

if id_case
  % show diag
  figure(32490); plot(1:n/2, A0(1:2*n+2:end), '-+', 1:n/2, A(1:2*n+2:end), '-o')
  xlabel('index of 2x2 matrix');
  ylabel('value');
  legend('original', 'hand fit', 'location', 'southeast');
  pic_output_png('a11diag');
  figure(32491); plot(1:n/2, A0(2:2*n+2:end), '-+', 1:n/2, A(2:2*n+2:end), '-o')
  xlabel('index of 2x2 matrix');
  ylabel('value');
  legend('original', 'hand fit', 'location', 'southeast');
  pic_output_png('a21diag');
  figure(32492); plot(1:n/2, A0(n+1:2*n+2:end), '-+', 1:n/2, A(n+1:2*n+2:end), '-o')
  xlabel('index of 2x2 matrix');
  ylabel('value');
  legend('original', 'hand fit', 'location', 'southeast');
  pic_output_png('a12diag');
  figure(32493); plot(1:n/2, A0(n+2:2*n+2:end), '-+', 1:n/2, A(n+2:2*n+2:end), '-o')
  xlabel('index of 2x2 matrix');
  ylabel('value');
  legend('original', 'hand fit', 'location', 'southeast');
  pic_output_png('a22diag');
end

figure(10);
imagesc(A);
title('A');
caxis([-1, 1])
colorbar();
pic_output_png('A');

Aei = A([1:2:n 2:2:n], [1:2:n 2:2:n]);
figure(11);
imagesc(Aei);
title('Aei');
caxis([-1, 1])
colorbar();
pic_output_png('Aei');

figure(345);
imagesc(log10(abs(Aei)));
colorbar();
caxis([-6, 1])
pic_output_png('Aei_log10');

%% verification
%A615 = eigvec_skel \ A_skel * eigvec_skel;
%A615ei = A615([1:2:n 2:2:n], [1:2:n 2:2:n]);
%figure(346);
%imagesc(log10(abs(A615)));
%colorbar();
%caxis([-6, 1])
%%pic_output_png('Aei_log10');


[eigvec eigval] = eig(A);
eigval = diag(eigval);

%figure(1);
%plot(real(eigval));

%figure(2);
%imagesc(abs(eigvec));
%caxis([-1, 1])
%colorbar();

% sort by abs
[~, sidx] = sort(-eigval);
eigval = eigval(sidx);
eigvec = eigvec(:, sidx);

figure(3);
plot(real(eigval));
xlabel('idx');
ylabel('eig-val');
pic_output_png('eig-val');

figure(300);
plot(1:n, real(eigval), 1:n, sort(eigvalapp(:), 'descend'));
xlabel('idx');
ylabel('eig-val');
legend('eig', 'approx');
pic_output_png('eig-val-app');

figure(301);
semilogy(1:n, -real(eigval), 1:n, -sort(eigvalapp(:), 'descend'));
xlabel('idx');
ylabel('abs-eig-val');
legend('eig', 'approx');
pic_output_png('eig-val-app-log');

figure(4);
imagesc(abs(eigvec));
xlabel('idx');
%caxis([-1, 1])
title('eig-vec');
%colorbar();
pic_output_png('eig-vec');

figure(20);
semilogy(sort(sum(abs(A),2)))
ylabel('sum row');
xlabel('row id');
pic_output_png('sum-row');

