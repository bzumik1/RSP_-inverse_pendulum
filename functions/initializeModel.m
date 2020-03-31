function [] = initializeModel()
    % Vytvori .m files s dosazen�mi parametry modelu pro urychlen� v�po?t?
    % pendCartC.m -> spojit�, neline�rn� model soustavy
    % AB.m -> symbolick� A a B matice -> pro linearizaci v libovoln�m
    % bod? stavov�ho prostoru
    
    % Definice parametr? soustavy
    p = getParameters();

    % Dosazen� zn�m�ch parametr? do stavov�ch rovnic neline�rn�ho modelu
    % Nezn�m� "?�d�c�" prom?nn� jako stavov� vektor |Y|, poruchy |d|, a ak?n�
    % veli?ina |u| se dosad� jako symbolick�.
    syms Y [4 1] %   Y [1 4] = stavov� vektor, sloupcov� 4x1 
    syms d [2 1]
    syms u
    % Vytvo?�me symbolick� rovnice s dosazen�mi parametry
    dXdt = pendulumCart_symbolicPars(Y,u,d,p);
    % Za poruchu dosazuji 0. Poruchy jsou dve, jedna pro s�lu a jedna pro
    % moment -> porucha d je vektor 2x1
    dXdt_incl_d = dXdt;
    dXdt = subs(dXdt, d, [0 0]');
    % Ze symbolick� rovnice vytvo?�me function handler, kter� ulo��me jako .m
    % file do slo�ky functions. Vytvo?en� funkce vrac� 4x1 vektor dY.
    matlabFunction(dXdt,...
        'file','functions/pendCartC',...
        'vars', {[Y1; Y2; Y3; Y4], u});
    
        matlabFunction(dXdt_incl_d,...
        'file','functions/pendCartC_d',...
        'vars', {[Y1; Y2; Y3; Y4], u, [d1; d2]});

    % Podobn? se postupuje p?i dosazov�n� parametr? do rovnice linearizovan�ho
    % popisu. Vytvo?en� funkce vrac� dv? matice, A a B.
    [A,B] = AB_symbolicPars(Y,u,p);
    % AX + BU -> dXdt linearizovan�
    matlabFunction(A, B,...
        'file', 'functions/AB',...
        'vars', {[Y1; Y2; Y3; Y4], u});    
end