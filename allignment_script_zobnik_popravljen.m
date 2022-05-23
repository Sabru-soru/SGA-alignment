% SCRIPT FOR POINT CLOUD ALLIGNMENT
% NOTE: 3D points are given as rows of matrices R and S
% R = reference points
% S = scanned points that need to be alligned

% PODATKI so taki, da velja:
% R je centriran v izhodiscu, pri cemer je cilinder zobnika vzporeden z z-osjo 
% TODO: kaj ce je drugace? 

% IDEA:
% (1) Izracunaj tezisce S in premakni tocke v oblaku S v izhodisce, tako da
% premaknes za vektor tezisca. 
% (2) Izracunaj matriko dvojnih vztrajnostnih momentov.
% (3) Z uporabo rotacij unici (postavi na 0) izvendiagonalne elemente
%     matrike vztrajnostnih momentov za S
% (4) Poravnava notranjega cilindra
% (5) Poravnava visine (translacija v z - smeri)
% (6) Poravnava rotacije v ravnini x-y

tic; % start timer
draw_intermediate_plots = 1; % draw plots after each adjustment (1) - (6)

zobov = 20;     % stevilo zobov (kako to izracunati pri splosnem zobniku??)
kot = pi/zobov; % koliksen kot pripada enemu zobu (samo polovica!)

d = 20; %reference diameter
m = 1; %module
df = d - 2.5*m; %root diameter

% predpostavka, da je CAD model simetricen pri 0 stopinj, tj.
% polarni kot od prvega zoba je [-kot,kot]!!

% Referencne tocke
R = dlmread('Gear_3Dscan.asc',' ');
% potrebujemo samo point cloud (brez normal)
R = R(:,1:3);

% cilindricne koordinate za R
[theta_R,rho_R,z_R] = cart2pol(R(:,1),R(:,2),R(:,3));

% izracunaj podatke o zobniku: polmer cilindra, polmer zobov, polovicna
% visina zobnika
notr_polmer = 3.075; %min(rho_R);
zun_polmer = max(rho_R);
visina = max(abs(max(z_R)),abs(min(z_R))); % aplikate tock zobnika so [-visina,visina]

fprintf('notranji polmer cilindra: %f\n', notr_polmer);
fprintf('          polmer zobnika: %f\n', zun_polmer);
fprintf('          visina zobnika: %f\n', visina);



% preberi skenirane tocke, ki jih poravnavamo na R
ime_S = 'Gear_3Dscan_zamaknjena_luknja.asc';
S = dlmread(ime_S,' ');
% potrebujemo samo point cloud (brez normal)
S = S(:,1:3);

fprintf('\nBranje podatkov uspesno zakljuceno.\n');

if draw_intermediate_plots

    plot3(S(:,1),S(:,2),S(:,3),'.');
    hold on
    plot3(R(:,1),R(:,2),R(:,3),'.');
    axis equal
    xlabel('x'); ylabel('y'); zlabel('z'); 
    title('Slika 1 - zacetna pozicija obeh zobnikov','Color','k');

end

% stevilo tock v R in S
n = size(R,1);
m = size(S,1);

%% (1) Izracunaj tezisce S
S_center = sum(S)/m;

% od vsake tocke iz S odstej vektor tezisca, da centriras tocke v
% izhodisce.
S = S - repmat(S_center,m,1);


%% (2) Izracunaj matriko dvojnih vztrajnostnih momentov
fprintf('\nZacetni vztrajnostni matriki (dvojni momenti) za R in S:');
mi = 1/m; % mase

R_I_xx = norm(R(:,[2 3]),'fro')^2/n;
R_I_yy = norm(R(:,[1 3]),'fro')^2/n;
R_I_zz = norm(R(:,[1 2]),'fro')^2/n;

R_I = diag([R_I_xx, R_I_yy, R_I_zz])

S_I_xx = mi*norm(S(:,[2 3]),'fro')^2;
S_I_yy = mi*norm(S(:,[1 3]),'fro')^2;
S_I_zz = mi*norm(S(:,[1 2]),'fro')^2;
S_I_xy = -sum(prod(S(:,[1 2]),2))*mi;
S_I_xz = -sum(prod(S(:,[1 3]),2))*mi;
S_I_yz = -sum(prod(S(:,[2 3]),2))*mi;

S_I = [S_I_xx S_I_xy S_I_xz; S_I_xy S_I_yy S_I_yz; S_I_xz S_I_yz S_I_zz]



% iz diagonalnih elementov ugotovi kako je zobnik pozicioniran, tj.
% najvecji diagonalni element doloca os simetrije, tj. os vzporedno z
% notranjim valjem zobnika
[~,axis_symm] = max(diag(R_I));
names = ['x','y','z'];
fprintf('Os simetrije za referencne tocke R je %s.\n',names(axis_symm));

[lastni,lastne] = eig(S_I);
S_I = lastne

%izvedi transformacijo
S = (lastni'*S')';





% %% (3) z uporabo rotacij unici (postavi na 0) izvendiagonalne elemente matrike S_I
% iter = 0;
% prehodna = eye(3);
% % isci rotacije dokler absolutna vrednost najvecjega izvendiagonalnega elementa
% % v S_I ne pade pod toleranco 1e-10
% while max(max(abs(triu(S_I,1)))) > 1e-10
%     
%     iter = iter + 1;
%     
%     % poisci najvecji element po absolutni vrednosti
%     % in njegov indeks (doloca eno od treh moznih rotacij)
%     [~,ind] = max(abs(S_I(find(triu(abs(S_I),1)))));
%     
%     % uporabimo formule za transformacijo vztrajnostnih momentov
%     % https://calcresource.com/moment-of-inertia-rotation.html
%     
%     % resi nelinearno enacbo za iskani kot, da postavis izvendiagonalni
%     % element na 0
%     switch ind
%         case 1
%             f = @(x) cos(2*x)*S_I_xy + 0.5*sin(2*x)*(S_I_xx - S_I_yy);
%             fi = fzero(f,0);
%             Q = [cos(fi) -sin(fi) 0;sin(fi) cos(fi) 0; 0 0 1];
%         case 2
%             f = @(x) cos(2*x)*S_I_xz + 0.5*sin(2*x)*(S_I_xx - S_I_zz);
%             fi = fzero(f,0);
%             Q = [cos(fi) 0 -sin(fi);0 1 0;sin(fi) 0 cos(fi)];
%         case 3
%             f = @(x) cos(2*x)*S_I_yz + 0.5*sin(2*x)*(S_I_yy - S_I_zz);
%             fi = fzero(f,0);
%             Q = [1 0 0; 0 cos(fi) -sin(fi); 0 sin(fi) cos(fi)];
%     end
%     
%     % izvedi transformacijo
%     S = (Q*S')';
%     
%     prehodna = Q * prehodna;
%     
%     % ponovni izracun matrike momentov S_I
%     S_I_xx = mi*norm(S(:,[2 3]),'fro')^2;
%     S_I_yy = mi*norm(S(:,[1 3]),'fro')^2;
%     S_I_zz = mi*norm(S(:,[1 2]),'fro')^2;
%     S_I_xy = -sum(prod(S(:,[1 2]),2))*mi;
%     S_I_xz = -sum(prod(S(:,[1 3]),2))*mi;
%     S_I_yz = -sum(prod(S(:,[2 3]),2))*mi;
% 
%     S_I = [S_I_xx S_I_xy S_I_xz; S_I_xy S_I_yy S_I_yz; S_I_xz S_I_yz S_I_zz];
% end
% 
% fprintf('\nMatrika momentov po %d rotacijah:', iter);
% S_I
% prehodna

% poglej, ce moramo zobnik zarotirati za 90 stopinj, da poravnamo os
% simetrije z oblakom tock R: poisci najvecji diagonalni element
[~,ind_R] = max(diag(R_I)); % predpostavimo, da je z-os os simetrije za R
[~,ind_S] = max(diag(S_I));
if ind_S == ind_R
    fprintf('Zobnik S je pravilno rotiran glede na R.\n')
else
    % če imamo trenutno os x, hočemo pa imeti z, zavrtimo okrog y (vedno ta
    % tretjo!!
    fprintf('Zobnik S je potrebno zavrteti za 90 stopinj...');
    
    % zavrti okrog ta tretje osi
    if ind_S == 1 % rotiraj okrog y
        rotation = [0 0 -1;0 1 0;1 0 0];
    else % rotiraj okrog x
        rotation = [1 0 0; 0 0 -1; 0 1 0];
    end
   
    S = (rotation*S')';
    fprintf('narejeno\n');
end

if draw_intermediate_plots
    
    % plot
    figure
    plot3(S(:,1),S(:,2),S(:,3),'.');
    hold on
    plot3(R(:,1),R(:,2),R(:,3),'.');
    axis equal
    xlabel('x'); ylabel('y'); zlabel('z'); 
    title('Slika 2 - izvendiagonalni elementi postavljeni na 0','Color','k')

end

%% (4) Poravnava notranjega cilindra

% 1. korak: najdi tocke, ki so blizu cilindra, tj. s polarnim radijem blizu
% notranjemu polmeru R in vrednostmi aplikat na intervalu 
% abs(z) < 0.9 * visina
% POMEMBNO: 10% robnih tock odrezemo, ker so slabse pomerjene! 
z_S_approx = 0.9 * visina;
tol = 0.5;
[~,rho_S,z_S] = cart2pol(S(:,1),S(:,2),S(:,3));

% STARA KODA
%S_inner = S(abs(z_S) < z_S_approx & abs(rho_S - notr_polmer) < tol,:);

% NOVA KODA
S_inner = S(abs(z_S) < z_S_approx & abs(rho_S) < 1.2*notr_polmer,:);

% 2. korak: znotraj zanke za vsako od moznih rotacij minimiziraj funkcijo
% napake = sum_i (norm(Q*r - v) - notr_polmer)^2, tj. poisci premik v in
% kot fi, ki doloca rotacijo Q, tako da minimiziras funkcijo napake
options = optimoptions(@fminunc,'Display','off');

% x - rotacija
f_allign = @(x) std_napaka(S_inner,x(1),x(2:4),1,notr_polmer);
[x_allign, next_napaka] = fminunc(f_allign,zeros(4,1),options);
fi = x_allign(1);
v = x_allign(2:4);
Q = [cos(fi) -sin(fi) 0;sin(fi) cos(fi) 0; 0 0 1];
S = (Q*S' - repmat(v,1,size(S,1)))';

% izracun nove S_inner
[~,rho_S,z_S] = cart2pol(S(:,1),S(:,2),S(:,3));
%STARA KODA
%S_inner = S(abs(z_S) < z_S_approx & abs(rho_S - notr_polmer) < tol,:);
S_inner = S(abs(z_S) < z_S_approx & abs(rho_S) < 1.2*notr_polmer,:);

% dokler je najvecji element minimizatorja po absolutni vrednosti nad
% toleranco 1e-5, izvajaj zanko
while max(abs(x_allign)) > 1e-5
    
    % y-rotacija
    f_allign = @(x) std_napaka(S_inner,x(1),x(2:4),2,notr_polmer);
    x_allign = fminunc(f_allign,zeros(4,1),options);
    fi = x_allign(1);
    v = x_allign(2:4);
    Q = [cos(fi) 0 -sin(fi);0 1 0;sin(fi) 0 cos(fi)];
    S = (Q*S' - repmat(v,1,size(S,1)))';
    
    % izracun nove S_inner
    [~,rho_S,z_S] = cart2pol(S(:,1),S(:,2),S(:,3));
    %S_inner = S(abs(z_S) < z_S_approx & abs(rho_S - notr_polmer) < tol,:);
    S_inner = S(abs(z_S) < z_S_approx & abs(rho_S) < 1.2*notr_polmer,:);
    
    % z-rotacija
    f_allign = @(x) std_napaka(S_inner,x(1),x(2:4),3,notr_polmer);
    x_allign = fminunc(f_allign,zeros(4,1),options);
    fi = x_allign(1);
    v = x_allign(2:4);
    Q = [1 0 0; 0 cos(fi) -sin(fi); 0 sin(fi) cos(fi)];
    S = (Q*S' - repmat(v,1,size(S,1)))';
    
    % izracun nove S_inner
    [~,rho_S,z_S] = cart2pol(S(:,1),S(:,2),S(:,3));
    %S_inner = S(abs(z_S) < z_S_approx & abs(rho_S - notr_polmer) < tol,:);
    S_inner = S(abs(z_S) < z_S_approx & abs(rho_S) < 1.2*notr_polmer,:);
    
    % x-rotacija
    f_allign = @(x) std_napaka(S_inner,x(1),x(2:4),1,notr_polmer);
    [x_allign, napaka_allign] = fminunc(f_allign,zeros(4,1),options);
    fi = x_allign(1);
    v = x_allign(2:4);
    Q = [cos(fi) -sin(fi) 0;sin(fi) cos(fi) 0; 0 0 1];
    S = (Q*S' - repmat(v,1,size(S,1)))';
    
    % izracun nove S_inner
    [~,rho_S,z_S] = cart2pol(S(:,1),S(:,2),S(:,3));
    %S_inner = S(abs(z_S) < z_S_approx & abs(rho_S - notr_polmer) < tol,:);
    S_inner = S(abs(z_S) < z_S_approx & abs(rho_S) < 1.2*notr_polmer,:);
    
end

% TODO: spremeni formulo za napako
% napaka (vsota kvadratov) od radija: napake = sum_i (norma_(xy)(S_inner) - notr_polmer)^2
fprintf('Napaka po poravnavi cilindra je %f\n', sum((sqrt(sum(S_inner(:,1:2).^2,2)) - notr_polmer * ones(size(S_inner,1),1)).^2));

if draw_intermediate_plots
    
    figure
    plot3(S(:,1),S(:,2),S(:,3),'.');
    hold on
    plot3(R(:,1),R(:,2),R(:,3),'.');
    axis equal
    xlabel('x'); ylabel('y'); zlabel('z'); 
    title('Slika 3 - poravnava cilindra','Color','k')

end

%% (5) Poravnava visine (translacija v z - smeri)
% IDEJA: s pomocjo median odstrani osamelce in s pomocjo median 
% določi ravnino, ki aproksimira spodnjo in zgornjo osnovno ploskev zobnika
% source: https://eurekastatistics.com/using-the-median-absolute-deviation-to-find-outliers/

% cilindrične koordinate
[~,rho_S,z_S] = cart2pol(S(:,1),S(:,2),S(:,3));

% OPOMBA: na robu je par osamelcev (slaba natančnost meritev)
% -> za določitev ravnin, ki vsebuje osnovni ploskvi zobnika
% vzemi samo del točk, ki imajo kvečjemu radij (1-del) * zunanji polmer
% in vsaj radij (1+del) * notranji polmer
del = 0.3;

% zgornja osnovna ploskev
% STARA KODA
%S_inner = S(z_S > 0 & (rho_S < zun_polmer * (1-del)) & ((1+del) * notr_polmer < rho_S), :);

% NOVA KODA
S_inner = S(z_S > 0 & (rho_S < 0.45*df) & (rho_S >= 1.2*notr_polmer), :);

zg = median(S_inner(:,3));                      % mediana
med_povprecje = median(abs(S_inner(:,3) - zg)); % mediana absolutnih odmikov od mediane

% odstranimo osamelce, tj. vsi ki imajo relativne odmike od mediane več kot
% 1 se odstrani
S_final = S_inner( abs(S_inner(:,3) - zg)/med_povprecje < 1, :); 

% ravnina z = const, ki aproksimira zgornjo ploskev ima vrednost 
zg = median(S_final(:,3)); 


% ENAKO za spodnjo ploskev
% NOVA KODA
S_inner = S(z_S < 0 & (rho_S < 0.45*df) & (rho_S >= 1.2*notr_polmer), :);

sp = median(S_inner(:,3));
med_povprecje = median(abs(S_inner(:,3) - sp));

S_final = S_inner( abs(S_inner(:,3) - sp)/med_povprecje < 1, :);

% ravnina z = const, ki aproksimira spodnjo ploskev ima vrednost 
sp = median(S_final(:,3));

% doloci premik h, ki minimizira odmika zg in sp do visina in -visina
f = @(h) abs(sp + h + visina) + abs(zg + h - visina);

% zacetni priblizek
razlika = visina - zg;
premik = fminunc(f,razlika,options)

% popravi S
S(:,3) = S(:,3) + premik;

if draw_intermediate_plots

    figure
    plot3(S(:,1),S(:,2),S(:,3),'.');
    hold on
    plot3(R(:,1),R(:,2),R(:,3),'.');
    axis equal
    xlabel('x'); ylabel('y'); zlabel('z'); 
    title('Slika 4 - poravnava visine','Color','k')

end

% OPOMBA: min ||z - S_inner(:,3)||_1 vrne isto kot median(S_inner(:,3))


%% (6) poravnava rotacije v ravnini x-y
% za oba oblaka tock R in S izracunaj simetralo enega od zobov (zanima nas
% pri katerem kotu se to zgodi)

% R
tol = 1e-2;
ind = abs(zun_polmer - rho_R) < tol;
R_ind = R(ind,:);
theta_R_ind = theta_R(ind,:);
% doloci zob
ind = abs(theta_R_ind) < kot;
R_zob = R_ind(ind,:);
theta_zob = theta_R_ind(ind,:);

% izracunaj povprecje polarnih kotov
sredina_R = sum(theta_zob)/length(theta_zob);

% S
tol = 1e-1;
[theta,rho,z] = cart2pol(S(:,1),S(:,2),S(:,3));

% find angle where max(rho) is achieved
[~,ind_max] = max(rho);
theta_S = theta(ind_max);

ind = abs(max(rho) - rho) < tol;
S_ind = S(ind,:);
theta_ind = theta(ind);
% doloci zob
ind = abs(theta_ind - theta_S) < kot;
S_zob = S_ind(ind,:);
S_theta_zob = theta_ind(ind,:);

% izracunaj povprecje polarnih kotov
sredina_S = sum(S_theta_zob)/length(S_theta_zob);

% poisci kot rotacije
zavrti = sredina_R - sredina_S;
rotation = [cos(zavrti) -sin(zavrti) 0; sin(zavrti) cos(zavrti) 0; 0 0 1]; % okrog z
S = (rotation*S')';

if draw_intermediate_plots

    figure;
    plot3(S(:,1),S(:,2),S(:,3),'.');
    hold on
    plot3(R(:,1),R(:,2),R(:,3),'.');
    axis equal
    xlabel('x'); ylabel('y'); zlabel('z'); 
    lgd = legend('S - alligned','R - reference');
    lgd.FontSize = 12;
    title('Plot 5 - rotation adjustment (final allignment)','Color','k')

end


% save alignment
output_file = sprintf('poravnava - %s',ime_S);
dlmwrite(output_file, S, 'delimiter',' ','precision', 7);

% stop timer
toc;

%%
%%%%%% POMOZNE FUNKCIJE %%%%%%

% PORAVNAVA: oddaljenost od centralnega cilindra
% PAZI: funkcija je napisana tako, da predpostavlja, da je z os simetrije
% danega zobnika
% napaka = sum_i (norm(Q*r - v) - inner_rad)^2
function napaka = std_napaka(S,fi,v,os,inner_rad)

    switch os
        case 1
            Q = [cos(fi) -sin(fi) 0;sin(fi) cos(fi) 0; 0 0 1];
        case 2
            Q = [cos(fi) 0 -sin(fi);0 1 0;sin(fi) 0 cos(fi)];
        case 3
            Q = [1 0 0; 0 cos(fi) -sin(fi); 0 sin(fi) cos(fi)]; 
        otherwise
            error('Napacna os rotacije');
    end
    
    sos = (Q*S' - repmat(v,1,size(S,1)));
    % vzemi samo oddaljenost od centralnega radija. x in y coordinata
    sos = sos(1:2,:).^2;
    
    % 2. moznost: VSOTA KVADRATOV
    napaka = sum(sum(sos)) - 2*inner_rad*sum(sqrt(sum(sos))) + (inner_rad)^2*size(S,1);
    %napaka = napaka / size(S,1);

end 
            
           