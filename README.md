# Inverzní kyvadlo

# mainScript.m
Tento skript je je hlavním skriptem celého projektu

# model_sym.mlx
Tento live script slouží k odvození symbolických nelineárních stavových rovnic soustavy pro funkci pendulumCart.m a symbolických Jakobiánů pro linearizovaný popis soustavy pro skript ABCD.m

# CoursewareResources
Tato složka obsahuje všechny soubory týkající se projektu

# Funkce
Všechny funkce použité v hlavním skriptu jsou ve složce functions. 

- ABCD.m
function [A,B,C,D] = ABCD(X, u, p)
Vrací linearizovaný popis soustavy ve stavové formulaci. Vstupem do funkce je stavový bod 'X', ve kterém chceme linearizovat; skalární akční veličina 'u', a struct 'p' s parametry soustavy.
- animRefresh.m
Vykreslí snímek animace. Na vstupu je stav 'Xs' a požadovaný stav 'Wx'.
- columnMatrixRearrange.m
Prohodí pořadí řádků ve sloupcovém vektoru.
- pendulumCart.m
Vrací časovou derivaci 'dX/dt' pro daný stavový bod 'X', vstupní veličinu 'u' a poruchovou veličinu 'd'. 
- plotRefresh.m
Vykreslí snímek grafu.
