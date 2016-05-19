format long;

nb_dataset = 5;
hand = 'r';
nb_trials_training = 27;
k = 7;

convertEEGLABdataIntoCsv(nb_dataset, ALLEEG, hand);
load('data_events.mat');

truc = kNN(final_mat_X, ex_events_Y, nb_trials_training, ALLEEG(nb_dataset).trials, k);