%Initializing
t = 1;
dt = [-1,-0.3,-0.15,0,0.15,0.5,1];
N = 5; %number of atoms
figure('Renderer', 'painters', 'Position', [200 200 1200 300]); %Position of output figure with plot
H = zeros(N,N); E = zeros(1,N); P = zeros(1,N);
%%
%Energy values on site
for dtID=1:length(dt)
	T = dt(dtID);
	for ii=1:N-1 %Defining hamiltonian
        H(ii,ii+1) = t + (-1)^(ii+1)*T;
        H(ii+1,ii) = H(ii,ii+1);
    end
	[V,D] = eig(H);
    for jj=1:N
        E(jj) = D(jj,jj); 
        %P(jj) = conj(V(jj,N/2))*V(jj,N/2);
    end
	%Plotting Energy on site
        ax = subplot(1,length(dt),dtID);
        plot(ax,1:N,E,'ko','MarkerSize',1);
        if(T>0) title(sprintf('Trivial\n \\deltat = %0.1f',T));
        elseif(T<0) title(sprintf('Non-trivial\n \\deltat = %0.1f',T));
        else title(sprintf('Topological phase transition\n\\deltat = %0.1f',T));
        end
        if(dtID==1) xlabel(ax,'Number of state\---->'); ylabel(ax,'E - energy'); end
end
suptitle(sprintf('\nEnergy for %d atoms\n',N));
%%
%Wavefunction density on site
N = 500;
dt = [-0.8,-0.4,-0.05,0];
%dt = [ -0.10, 0, 0.3];
P = zeros(length(dt),N);
edge_state = zeros(1,2);
tol = 1e-10; %tolerance for finding edge states
figure('Renderer', 'painters', 'Position', [200 200 1200 400]); %Position of output figure with plot
for dtID=1:length(dt)
    T = dt(dtID); 
    ID = 1;
    for ii=1:N-1 %Defining hamiltonian
        H(ii,ii+1) = t + (-1)^(ii+1)*T;
        H(ii+1,ii) = H(ii,ii+1);
    end
    [V,D] = eig(H);
    for jj=1:N 
        P(dtID,jj) = abs(V(jj,floor(N/2))); 
    end
end

%Plotting wavfunction density on site 5 different sites
plot(1:N,P(1,:),'k-',1:N,P(2,:),'r-',1:N,P(3,:),'b-',1:N,P(4,:),'g-','LineWidth',0.2);
title(sprintf('Wavefunction density for %d atoms - non trivial case',N));
legend(sprintf('\\deltat= %0.2f',dt(1)),sprintf('\\deltat= %0.2f',dt(2)),sprintf('\\deltat= %0.2f',dt(3)),sprintf('\\deltat= %0.2f',dt(4)));
xlabel('site'); ylabel('|\Psi|');



