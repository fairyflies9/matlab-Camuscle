% compute maximum force vs frequency

nu_vec = logspace(log10(1),log10(1000),20);
force = zeros(size(nu_vec));

for nn = 1:length(nu_vec)
    nu = nu_vec(nn);
    ca_muscle_spikes_KT_v1;
    force(nn) = max(P);
end

figure(6);clf;
semilogx(nu_vec,force,'o-','linewidth',2)
set(gca,'fontsize',18)

load('force_square');
hold on;
semilogx(nu_vec,force,'v-','linewidth',2)