function refestplot(V0,nu_ar, nu_rest,nu_ar_invregc,nu_rest_invregc,L2)
% REFESTPLOT plot the potential estimate by AR and REST

nu_rest=pinv(L2)*nu_rest;
nu_rest_invregc=pinv(L2)*nu_rest_invregc;

for i=1:n
figure, 

subplot(1,2,1),
plot(V0,'r-');hold on
plot(nu_ar(:,i),'g-');hold on;plot(nu_rest(:,i),'b-');
xlabel('channel order');ylabel('Amplitude by evidencer.m');

subplot(1,2,2),
plot(V0,'r-');hold on
plot(nu_ar_invregc(:,i),'g-');hold on; plot(nu_rest_invregc(:,i),'b-');
xlabel('channel order');ylabel('Amplitude by evidencer_invregc.m');

legend('GR','AR','REST');
close all;

end