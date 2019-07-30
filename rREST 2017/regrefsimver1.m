% Simulation reference estimate and model selection
clean;
pt = 4;           Nc = 58;       % patch type and channels number
Nt = 5120;    Np = 1001;  % time length and lmds number
Ns = 89;       Ne = 7;         % subjects and noise test number

re = zeros(Ns,10, Ne);   rebm = zeros(Ns,22,Ne);
oli = zeros(Ns,  5, Ne);     omi = zeros(Ns,3,6,Ne);

Vi = zeros(Np,Nc,Nt,5);

dbstop if error
for i=1:Ns
    for j=1:Ne
        
        [V0,V1, K1, K2, K3, K4] = rseegsim(i, pt, mksnr(j)); refvct;
        
        re(i,:,j) = refpre(Vi, msi, V0,Vr,K1, K2, K3, K4);
       
        display(strcat('subjects_', num2str(i), '  noise_', num2str(j) , '  done!'));
    end
    
end

save sim1 re;