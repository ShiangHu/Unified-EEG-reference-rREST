% Validation reference estimate and model selection
clean;
en=dir('.\subject\EEG');en(1:2)=[];
kd=dir('.\subject\LF');kd(1:2)=[];
Ns = size(en,1); % # of subjects
Np = 1001;        % # of lmds
msi = zeros(Np,6,Ns,4);

for i=1:Ns
    eeg = load(en(i).name);
    data=eeg.data; data(eeg.pc_malos,:)=[];
    data=data./norm(data,'fro');
    H0=Hsc(size(data,1)); data=H0*data;
    
    [K1, ~, K3, K4] = trsK( kd, eeg.pc_malos, i);

    % comment the line 43 and 49 in refbys to avoid exceeding the memory 
    [~, msi(:,:,i,1)]=refbys(data,[],'ar');
    [~, msi(:,:,i,2)]=refbys(data,K1,'rt'); % 'ind LF'
    [~, msi(:,:,i,3)]=refbys(data,K3,'rt'); % 'avg LF'
    [~, msi(:,:,i,4)]=refbys(data,K4,'rt'); % 'Yao LF'
    
    display(strcat(num2str(i),' subjects done!'));
end

save vld msi;