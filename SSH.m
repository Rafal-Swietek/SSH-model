%%Initial parameters
t = 1; %hopping integral
%dt = [-0.9,-0.8,-0.7,-0.6,-0.5,-0.4,-0.3,-0.2,-0.1,0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9]; %different hop
dt = [-1,-0.6,-0.2,0,0.2,0.6,1];
fi_numeric1 = zeros(1,length(dt));% berry phase numeric calculation
fi_numeric2 = zeros(1,length(dt));
syms k;
n = 1000;
sig_x = [0 1;1 0]; sig_y = [0 -i;i 0]; sig_z = [1 0;0 -1]; % Pauli matrices
x = linspace(-pi,pi);
%d_z = 0; % 2D model
%%
%Numerical & analitycal eigenvalues
Enum_1 = zeros(1,2*n+1); %for numerical eigenvalues
Enum_2 = zeros(1,2*n+1);
d_x = zeros(1,2*n+1);
d_y = zeros(1,2*n+1);
y = -pi:pi/n:pi;
figure('Renderer', 'painters', 'Position', [300 200 1000 500]); %Position of output figure with plot
for dtID=1:length(dt) % dtID enumerates different hopping
    W1 = zeros(2,2*n+1); %temporary eigenvectors for each site
    W2 = zeros(2,2*n+1);
    for ii=1:2*n+1 %calculating eigenvalue and defining D-vector
        k = y(ii);
        T = dt(dtID);
        H = [0, t+T+(t-T)*exp(-i*k);t+T+(t-T)*exp(i*k), 0];
        [V,E] = eig(H);
        Enum_1(ii) = E(1,1);
        Enum_2(ii) = E(2,2);
        d_x(ii) = (t+T)+(t-T)*cos(k);
        d_y(ii) = (t-T)*sin(k); %d_z = 0;
        W1(:,ii) = V(:,1);
        W2(:,ii) = V(:,2);
    end
    %Berry Phase
    gamma1 = 1;
    for ii=1:2*n
        S(ii) = dot(W1(:,ii),W1(:,ii+1))/(abs(dot(W1(:,ii),W1(:,ii+1))));  %element of multiplication for Berry phase
        R(ii) = dot(W2(:,ii),W2(:,ii+1))/(abs(dot(W2(:,ii),W2(:,ii+1))));
    end
    gamma1 = prod(S); %capital PI symbol
    gamma2 = prod(R);
    fi_numeric1(dtID) = angle(gamma1); %argument of multiplication
    fi_numeric2(dtID) = angle(gamma2);
    
    %analitycal eigenvalue
    H_analitycal1 =  ( (t+T)^2 + (t-T)^2 + 2*(t+T)*(t-T)*cos(x) ).^0.5;
    H_analitycal2 = -( (t+T)^2 + (t-T)^2 + 2*(t+T)*(t-T)*cos(x) ).^0.5; 
%%
%Plotting energy spectrum and D vector
    % Energy spectrum - comparison between analitycal and numerical
    ax1 = subplot(2,length(dt),dtID);
    p = plot(ax1,y,Enum_1,'ro',y,Enum_2,'ro',x,H_analitycal1,'k-',x,H_analitycal2,'k-','MarkerIndices',1:n/10:length(Enum_1),'MarkerSize',4,'Linewidth',1.5);
    title(sprintf('\\deltat = %0.1f',T));
    axis(ax1,[-pi pi -2.5 2.5]);
    xL = xlim;  yL = ylim;  line([0 0], yL);  line(xL, [0 0]);  %x axis & y-axis
    xticks([-pi -pi/2 0 pi/2 pi]);   xticklabels({'-\pi','-\pi/2','0','\pi/2','\pi'});
    %D - vector
    ax2 = subplot(2,length(dt),dtID+length(dt));
    if(T==1) plot(ax2,d_x,d_y,'r*');
    else plot(ax2,d_x,d_y,'r-');
    end
    axis(ax2,[-2.5 2.5 -2.5 2.5]);
    xlabel(ax2,'D_x(k)'); 
    if(dtID==1) ylabel(ax2,'D_y(k)'); ylabel(ax1,'E(k)'); end
    xL = xlim;  yL = ylim;  line([0 0], yL);  line(xL, [0 0]);  %x axis & y-axis
    if(T>0) title('Trivial'); 
    elseif(T<0) title('Non-trivial');
    else
        title(sprintf('Topological phase\n    transition'));
    end 
end
suptitle(sprintf('Energy spectrum and corresponding D-vector for different hopping in infinite SSH system')); 
legend([p(1) p(3)], 'Numerical', 'Analitycal','Location','northwestoutside'); %no multiple legends
%plotting berry phase vs hopping dt
figure(2);
q = plot(dt,fi_numeric1,'kx', dt,fi_numeric2,'rs');
legend('\lambda_+','\lambda_-');
xlabel(sprintf('\\deltat - hopping parameter'));
ylabel(sprintf('\\phi - Berry phase'));
xL = xlim;  yL = ylim;  line(xL,[0 0]);  line(xL,[pi pi]); line(xL,[-pi -pi]); line(xL,[pi/2 pi/2]); line(xL,[-pi/2 -pi/2]);
yticks([-pi -pi/2 0 pi/2 pi]);   yticklabels({'-\pi','-\pi/2','0','\pi/2','\pi'});
title('Berry phase vs hopping');
