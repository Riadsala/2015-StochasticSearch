clear all
close all
people = [3 4 5 6 7 8 9];
ctr = 0 ;
for p = 1:length(people)
    ctr = ctr+1;
person = people(p);
load(['../results/personDat' int2str(person)]);

ages(ctr) = personDat.age
end

ages = [ages 22 23]; % add in Alex and Melissa
mean(ages)
min(ages)
max(ages)