#!/bin/bash
#
# alias c2s='docker run -it --rm -v "$PWD:/data/workdir" jonasrauber/c2s'
# alias c2s >> ~/.bashrc
#
# execute c2s
for d in /gpfs01/nienborg/group/Katsuhisa/serotonin_project/LFP_project/Data/c2s/data_rest/*/;
  do cd "$d";
  # preprocess
  c2s preprocess data_base_cv10.mat preprocessed_base_cv10.mat &
  c2s preprocess data_drug_cv10.mat preprocessed_drug_cv10.mat &
  c2s preprocess data_lowFR_cv10.mat preprocessed_lowFR_cv10.mat &
  c2s preprocess data_highFR_cv10.mat preprocessed_highFR_cv10.mat &
  wait;
  # loo
  c2s leave-one-out preprocessed_base_cv10.mat predictions_base_cv10.mat &
  c2s leave-one-out preprocessed_drug_cv10.mat predictions_drug_cv10.mat &
  c2s leave-one-out preprocessed_lowFR_cv10.mat predictions_lowFR_cv10.mat &
  c2s leave-one-out preprocessed_highFR_cv10.mat predictions_highFR_cv10.mat &
  wait;
  # evaluate
  c2s evaluate -s 1 4 -o correlation_base_cv10.mat -m corr preprocessed_base_cv10.mat predictions_base_cv10.mat &
  c2s evaluate -s 1 4 -o MI_base_cv10.mat -m info preprocessed_base_cv10.mat predictions_base_cv10.mat &
  c2s evaluate -s 1 4 -o correlation_drug_cv10.mat -m corr preprocessed_drug_cv10.mat predictions_drug_cv10.mat &
  c2s evaluate -s 1 4 -o MI_drug_cv10.mat -m info preprocessed_drug_cv10.mat predictions_drug_cv10.mat &
  c2s evaluate -s 1 4 -o correlation_lowFR_cv10.mat -m corr preprocessed_lowFR_cv10.mat predictions_lowFR_cv10.mat &
  c2s evaluate -s 1 4 -o MI_lowFR_cv10.mat -m info preprocessed_lowFR_cv10.mat predictions_lowFR_cv10.mat &
  c2s evaluate -s 1 4 -o correlation_highFR_cv10.mat -m corr preprocessed_highFR_cv10.mat predictions_highFR_cv10.mat &
  c2s evaluate -s 1 4 -o MI_highFR_cv10.mat -m info preprocessed_highFR_cv10.mat predictions_highFR_cv10.mat &
  wait;
done
