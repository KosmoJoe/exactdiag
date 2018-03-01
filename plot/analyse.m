%Auswertung

clear
format long

 %Where is the data?
 dir = 'data/runAI005/'; 
 
 
    N = 5; % Number of eigenfunctions to plot
    
numRow = 2;
numCol=ceil(N/numRow);
plotStyle = 0;
ratio = numRow/numCol;
redo = 0;

obdM = 0;
xmin = -2;
xmax = 2;
mMax = 10;
crange=[-2, 2];
% %=====================read in computed data================================
    

if exist([ dir 'rel.mat'], 'file') ~= 2 || redo
    
    filen = 'grid';
    file = load([dir filen '.txt']);
    x = file(:,1);
    dx  = file(:,2);
    nx = length(x); 

    filen = 'potN';
    Vpot = load([dir filen '.txt']);

    filen = 'epsi';
    epsilon = load([dir filen '.txt']);

    filen = 'eigen';
    energ = load([dir filen '.txt']);
    m=sqrt(length(energ));

    filen = 'eigenV1B';
    eV = load([dir filen '.txt']);

    %%Get coeff matrix
    for jj = 1:N
        filen = ['eigenV_' int2str(jj)];
        file = load([dir filen '.txt']);
        for j = 1:N
            Vec{jj} = file(:,:);
        end  
    end

    %%Get integrated coeff matrix
    for jj = 1:N
        vecT = Vec{jj};
        coeff = zeros(m,m);
        for kk1=1:m
            for kk2=1:m
                for qq=1:m
                    coeff(kk1,kk2) = coeff(kk1,kk2) + conj(vecT(kk1,qq))*vecT(kk2,qq);
                end
            end
        end
        coeff1B{jj} = coeff;
    end

    for jj = 1:N
        norm = 0;
        for kk=1:m
            norm = norm + coeff1B{jj}(kk,kk);
        end
        norm1B(jj)=norm;
    end
    % ***************************************************************
    % Derive Grid Rep.***********************************************
    % ***************************************************************

    for jj = 1:N
        vecT = Vec{jj};
        vecNew = zeros(nx,nx);
        for kk=1:m
            for qq=1:m
                vecNew = vecNew + vecT(kk,qq)*eV(kk,:)'*eV(qq,:);
            end
        end
        vecGrid{jj} = vecNew;
    end
    
    
    save([ dir 'rel.mat'],'x','vecGrid','energ','epsilon','eV','Vec','m','-v7.3')
else
    load([ dir 'rel.mat'])
end

    %% Test
    vecGround = Vec{1}(1,1)*eV(1,:)'*eV(1,:);% + Vec{1}(2,2)*eV(2,:)'*eV(2,:);
    
    %%One Body density matrix
    if obdM == 1
        for kk = 1:3
            obrdm = zeros(length(x),length(x));
            for ix=1:length(x)
                for jx=1:length(x)
                   obrdm(ix,jx) = sum(conj(vecGrid{kk}(ix,:)).*vecGrid{kk}(jx,:).*dx(1,:),2);
                end
            end
            densityMAtrix{kk} = obrdm;
        end
    end
% ***************************************************************
% Figure 1 ******************************************************
% ***************************************************************
%%
figure(1)
clf
nOfPoints=10000;
colors = colormap(jet(nOfPoints));

for jj=1:N+1
    
    hs=subplot(numRow,numCol,jj);

    if (jj <=N)
        imagesc(x,x,vecGrid{jj});   
        title([ '$\epsilon_{' num2str(jj) '}=' num2str(energ(jj,1)) '$'],'FontSize',12,'Interpreter','none')
    else
        imagesc(x,x,vecGround); 
        title([ 'Test $\epsilon = ' num2str(epsilon(1)+epsilon(1)) '$'],'FontSize',12,'Interpreter','none')
    end
    %colorbar

    
    if jj>numCol*(numRow-1)
        xlabel('$z$ (units of $R^*$)');
    end
    if mod(jj-1,numCol)==0
        ylabel('$z$ (units of $R^*$)');
    end
    set(gca,'YDir','normal')
    caxis(crange)
    xlim([xmin,xmax])
    ylim([xmin,xmax])
    axis equal
    box on

end    


hc = get(hs,'Position');
colorbar('Position', [hc(1)+hc(3)*1.1  hc(2)  0.025  hc(4)]) %[left bottom width height]
caxis(crange)


qfst(gcf, [ dir 'eigenfuncs' ],ratio,plotStyle)

% ***************************************************************
% Figure 2 ******************************************************
% ***************************************************************
%%
figure(2)
clf
nOfPoints=100;
colors = colormap(jet(nOfPoints));

for jj=1:N
    
    hs =subplot(numRow,numCol,jj);

    imagesc(1:m,1:m,Vec{jj});        
    %colorbar
    title([ '$\epsilon_{' num2str(jj) '}=' num2str(energ(jj,1)) '$'],'FontSize',12,'Interpreter','none')
    set(gca,'YDir','normal')
    %axis image
    caxis([-1,1])
    xlim([0.5,mMax+0.5])
    ylim([0.5,mMax+0.5])
    box on

end    

hc = get(hs,'Position');
colorbar('Position', [hc(1)+hc(3)*1.3  hc(2)  0.025  hc(4)]) %[left bottom width height]
caxis([-1,1])

qfst(gcf, [ dir 'coefficients' ],ratio,plotStyle)
%--------------------------------------
% end
%--------------------------------------