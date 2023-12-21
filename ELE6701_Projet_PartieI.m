% H=[3 -2; 2 4]
% Htrans = transpose(H)
% I=eye(2);
% newmat = (I+Htrans*H)^-1*Htrans
% y=[6;4]
% newmat*y
clear all;
clc;
H=[1 0 0 0 0 0 0 0 0 0;
    -0.4 1 0 0 0 0 0 0 0 0;
    0 -0.4 1 0 0 0 0 0 0 0;
    0  0 -0.4 1 0 0 0 0 0 0;
    0 0  0 -0.4 1 0 0 0 0 0;
    0 0 0  0 -0.4 1 0 0 0 0;
    0 0 0 0  0 -0.4 1 0 0 0;
    0 0 0 0 0  0 -0.4 1 0 0;
    0 0 0 0 0 0  0 -0.4 1 0;
    0 0 0 0 0 0 0  0 -0.4 1];

%% G?n?ration des 1024 vecteurs hypoth?ses
vecteurs_cell = cell(1, 1024);

for i = 1:1024
    % Convertir l'indice en binaire sur 10 bits
    binaire = dec2bin(i-1, 10) - '0'; % Convertir en vecteur binaire
    
    % Remplacer les '0' par '-1' dans la repr?sentation binaire
    binaire(binaire==0) = -1;
    
    % Stocker le vecteur dans la cellule
    vecteurs_cell{i} = binaire';
end

% Convertir la cellule en une matrice de taille 10x1024
matrice_vecteurs_hypothese_s = cell2mat(vecteurs_cell);


%% G?n?ration du bruit n
% Variance d?sir?e
variance_desiree = (0.15);

% Nombre de vecteurs ? g?n?rer
nombre_paquets = 10000; %% VALEUR A MODIFIER
nombre_vecteurs = nombre_paquets;
taille_vecteur = 10;

% G?n?ration des 1024 vecteurs de bruit gaussien blanc avec moyenne nulle
vecteurs_bruit_gaussien = sqrt(variance_desiree) * randn(taille_vecteur, nombre_vecteurs);


%% G?n?ration de paquets s_i ? envoyer
% Nombre de vecteurs al?atoires ? g?n?rer
nombre_vecteurs = nombre_paquets;
taille_vecteur = 10;

% G?n?ration des vecteurs al?atoires
vecteurs_aleatoires_envoyes = randi([1, 2], taille_vecteur, nombre_vecteurs);
vecteurs_aleatoires_envoyes(vecteurs_aleatoires_envoyes == 2) = -1;

%% Vecteur y re?u
% Matrice des vecteurs y
y=vecteurs_aleatoires_envoyes + vecteurs_bruit_gaussien;


%% Detection ML
comparison_matrix=-1000000*ones(1,1024);
indice_guess = -1*ones(1,nombre_paquets);
vecteurs_s_guess = -3434314324*ones(10,nombre_paquets);
error_count_vector = 0;
error_count_symbol = 0;
for i=1:1:nombre_paquets
    
    for j=1:1:1024
        comparison_matrix(1,j) = real(transpose(matrice_vecteurs_hypothese_s(:,j))*y(:,i) - 0.5*transpose(matrice_vecteurs_hypothese_s(:,j))*matrice_vecteurs_hypothese_s(:,j));
    end
    [max_value, indice_guess(1,i)] = max(comparison_matrix);
    vecteurs_s_guess(:,i) = matrice_vecteurs_hypothese_s(:,indice_guess(1,i));
    comp = isequal(vecteurs_aleatoires_envoyes(:,i), vecteurs_s_guess(:,i));
    difference_per_symbol = nnz(vecteurs_s_guess(:,i) ~= vecteurs_aleatoires_envoyes(:,i));
    error_count_symbol = error_count_symbol + difference_per_symbol;
    if comp
        test = 0;
    else
        error_count_vector = error_count_vector+1;
        disp("error!");
    end
    
    
end
taux_erreur = error_count_vector/nombre_paquets
taux_erreur_par_symbole = error_count_symbol/(nombre_paquets*10)





