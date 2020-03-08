# Inverzní kyvadlo - semestrální práce ŘSP
__Krystof Bystricky a Jakub Znamenáček__

## Odkazy na potřebné soubory
- [Domácí stránka projektu](https://www.quanser.com/products/linear-servo-base-unit-inverted-pendulum/)
- [SIMULINK COURSEWARE - odkaz přímo na soubory](https://quanserinc.box.com/shared/static/gu9ed72edso2r2bfbtlyi3k6m2kgq0ie.zip)

## mainScript.m
Tento skript je je hlavním skriptem celého projektu

## model_sym.mlx
Tento live script slouží k odvození symbolických nelineárních stavových rovnic soustavy pro funkci pendulumCart.m a symbolických Jakobiánů pro linearizovaný popis soustavy pro uživatelsky definovanou funkci ABCD.m

## CoursewareResources
Tato složka obsahuje všechny soubory týkající se projektu

## functions
- Obsahuje všechny uživatelsky vytvořené funkce použité v projektu. 
### columnMatrixRearrange.m
- Prohodí pořadí řádků ve sloupcovém vektoru.
### ABCD.m
- function [A,B,C,D] = ABCD(X, u, p)
- Vrací linearizovaný popis soustavy ve stavové formulaci. Vstupem do funkce je stavový bod 'X', ve kterém chceme linearizovat; skalární akční veličina 'u', a struct 'p' s parametry soustavy.
### animRefresh.m
- Vykreslí snímek animace. Na vstupu je stav 'Xs' a požadovaný stav 'Wx'.
### pendulumCart.m
- Vrací časovou derivaci 'dX/dt' pro daný stavový bod 'X', vstupní veličinu 'u' a poruchovou veličinu 'd'. 
### plotRefresh.m
- Vykreslí snímek grafu.
