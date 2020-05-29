function [] = initializeModel()
    % Vytvori .m files s dosazen�mi parametry modelu pro urychlen� v�po?t?
    % pendCartC.m -> spojit�, neline�rn� model soustavy
    % AB.m -> symbolick� A a B matice -> pro linearizaci v libovoln�m
    % bod? stavov�ho prostoru
    
    % Definice parametr? soustavy
    p = getParameters();
    
    p_real = p;
    p_real ...
    % Dosazen� zn�m�ch parametr? do stavov�ch rovnic neline�rn�ho modelu
    % Nezn�m� "?�d�c�" prom?nn� jako stavov� vektor |X|, poruchy |d|, a ak?n�
    % veli?ina |u| se dosad� jako symbolick�.
    syms X [4 1] %   X [1 4] = stavov� vektor, sloupcov� 4x1 
    syms d [2 1]
    syms u
    % Vytvo?�me symbolick� rovnice s dosazen�mi parametry
    dXdt = pendCartC_symbolicPars(X,u,d,p);
    % Za poruchu dosazuji 0. Poruchy jsou dve, jedna pro s�lu a jedna pro
    % moment -> porucha d je vektor 2x1
    dXdt_incl_d = dXdt;
    dXdt = subs(dXdt, d, [0 0]');
    % Ze symbolick� rovnice vytvo?�me function handler, kter� ulo��me jako .m
    % file do slo�ky functions. Vytvo?en� funkce vrac� 4x1 vektor dY.
    matlabFunction(dXdt,...
        'file','functions/pendCartC',...
        'vars', {[X1; X2; X3; X4], u});
    
        matlabFunction(dXdt_incl_d,...
        'file','functions/pendCartC_d',...
        'vars', {[X1; X2; X3; X4], u, [d1; d2]});

    % Podobn? se postupuje p?i dosazov�n� parametr? do rovnice linearizovan�ho
    % popisu. Vytvo?en� funkce vrac� dv? matice, A a B.
    [A,B] = AB_symbolicPars(X,u,p);
    % AX + BU -> dXdt linearizovan�
    matlabFunction(A, B,...
        'file', 'functions/AB',...
        'vars', {[X1; X2; X3; X4], u});    
end