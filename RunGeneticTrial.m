clear
clc
tic
pc_graph = [1000 975 980 990 995 985 925 875 810 785 735 690 680 655 650 625 625 650 650 400 375 250 10 0].*1.25; %chamber pressure history in psi
pc_graph2 = [600 650 675 700 710 715 720 725 730 735 740 745 750 750 750 750 750 750 750 750 600 100 10 0].*1.25; %chamber pressure history in psi
rng(5)
id = 1; % Combo id (previously iter)
numBest = 1; % Take best n from each generation
numGens = 10; % Generations to run for

import_combo = readmatrix("Design Matrix - Praeto.csv",'Range', 'A2:A127','OutputType','string');

for i = 1:length(import_combo)
    id = i;
    b = run(id,numBest,numGens,curve1=pc_graph,curve2=pc_graph2,stage=1);
end
toc