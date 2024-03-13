%% Import the Old Faithful Dataset of Eruption Lenght vs. Waiting Time
data=readtable("faithful.csv");
X = normalize([data.waiting, data.eruptions]);

