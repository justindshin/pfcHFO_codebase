clear all
close all
%%

animalprefix_list ={'ZT2','JS34','JS17','JS21','JS14','JS15','KL8','ER1'};

day_list = 1;
cellcountthresh = 5;
sleep = 1;
for theAnimal = 1:length(animalprefix_list)
    animalprefix = animalprefix_list{theAnimal};
    for day = day_list
        
        if sleep == 0
            eps = [2:2:16];
        elseif sleep == 1
            eps = [1:2:17];
        end
        
        for ep = eps
            jds_PCC_replay_decoding_CA1_allcells(animalprefix,day,ep,cellcountthresh)
        end
    end
end

