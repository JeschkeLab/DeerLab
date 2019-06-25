function void=make_tables(),
% Compute and save 2^n kernels for the DEER Transformation program
%
% G. Jeschke, 2002
%
make_bas;
[base,tnorm,ny,t,crosstalk]=get_pake_kernel(64);
save('kernel64.mat','base','tnorm','ny','t','crosstalk');
[base,tnorm,ny,t,crosstalk]=get_pake_kernel(128);
save('kernel128.mat','base','tnorm','ny','t','crosstalk');
[base,tnorm,ny,t,crosstalk]=get_pake_kernel(256);
save('kernel256.mat','base','tnorm','ny','t','crosstalk');
[base,tnorm,ny,t,crosstalk]=get_pake_kernel(512);
save('kernel512.mat','base','tnorm','ny','t','crosstalk');
[base,tnorm,ny,t,crosstalk]=get_pake_kernel(1024);
save('kernel1024.mat','base','tnorm','ny','t','crosstalk');
[base,tnorm,ny,t,crosstalk]=get_pake_kernel(2048);
save('kernel2048.mat','base','tnorm','ny','t','crosstalk');
