% Simulation reference estimate and model selection
clean;
pt = 4;           % patch type
Nc = 58;        % channels number
Nt = 5120;     % time length
Np = 1001;    % LMDs number
Ns = 89;         % subjects number
H0=Hsc(Nc);
re = zeros(Ns,10); rebm = zeros(Ns,22);
oli = zeros(Ns,5);   omi = zeros(Ns,3,6);
nsr = 0.01+0.30*rand(Ns,1);
kd = dir('.\subject\LF'); kd(1:2)=[];
Vi = zeros(Np,Nc,Nt,5);
msi = zeros(Np,6,Ns,5);
[~, ~, K3, K4] = trsK(kd, []);% trs avg and yao lf

for i=1:Ns
       
    [V0,V1, K1, K2] = rseegsim(kd, i, pt, nsr(i));
    Vr = H0*V1;
   % cmt line 44 and uncmt 43, 49
    [Vi(:,:,:,1),msi(:,:,i,1)] = refbys(Vr,[],'ar');      % AR
    [Vi(:,:,:,2),msi(:,:,i,2)] = refbys(Vr,K1,'rt');    % Individual LF
    [Vi(:,:,:,3),msi(:,:,i,3)] = refbys(Vr,K2,'rt');    % Sparse Ind LF
    [Vi(:,:,:,4),msi(:,:,i,4)] = refbys(Vr,K3,'rt');    % Trs avg LF
    [Vi(:,:,:,5),msi(:,:,i,5)] = refbys(Vr,K4,'rt');    % Yao LF
    
    M = squeeze(msi(:,:,i,:));
    
    [re(i,:), rebm(i,:), oli(i,:), omi(i,:,:)] = refpre(Vi,M,V0,Vr,K1, K2, K3, K4);
    
    display(strcat(num2str(i),'  subjects done!'));
end
save sim re rebm msi oli omi nsr;