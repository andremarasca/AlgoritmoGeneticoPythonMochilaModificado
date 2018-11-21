# -*- coding: utf-8 -*-
"""
Created on Tue Jan 30 14:29:38 2018

@author: andre
"""

import matplotlib.image as mpimg

import numpy as np
import random as rand
import matplotlib.pyplot as plt
import copy

import OperacoesAG as ag
rand.seed(None)

imagem = mpimg.imread('entrada.jpg')

forma = imagem.shape

objetivo = []
for i in range(forma[0]):
    for j in range(forma[1]):
        for k in range(forma[2]):
            objetivo.append(float(imagem[i][j][k])/255)


######################################## Inicializacao do Algoritmo Genetico

tam_dna = imagem.size
qtd_ind = 100
geracoes = 2000
prob_mutacao = 5
amp_mut = 0.5
individuos_salvos = 3
tipo = 'minimizacao' # maximizacao ou minimizacao

###################################### Definicao dos limites

lim_sup = []
lim_inf = []
for i in range(0, tam_dna + 1):
    lim_sup.append(1);
    lim_inf.append(0);

######################################## Vetor de probabilidade

# S_selection = ini_LinearRanking(0, qtd_ind); % inicializar vetor de probabilidades
S_selection = ag.ini_ExponentialRanking(0.97, qtd_ind)

######################################### Avaliacao da populacao inicial

org = ag.gera_populacao (qtd_ind, tam_dna, lim_sup, lim_inf)
retorno = ag.avaliacao (org, qtd_ind, tam_dna, 0, objetivo)
org = ag.ordena(org, tipo)
print('best fit: ', org['fitness'][-1])

def dna2imagem(dna, forma):
    imagem = np.zeros(forma)

    u = 0
    for i in range(forma[0]):
        for j in range(forma[1]):
            for k in range(forma[2]):
                imagem[i][j][k] = dna[u]
                u += 1

    return imagem


########################################## Loop do Genético

melhor = org['fitness'][-1]

for geracao in range(0, geracoes):
    org2 = copy.deepcopy(org)

    for filho in range(0, qtd_ind - individuos_salvos):
        pai = ag.RankingSelection(S_selection, qtd_ind)
        mae = ag.RankingSelection(S_selection, qtd_ind)
        ag.cruzamento(org2, org, pai, mae, filho, tam_dna)

        ag.mutacao(org, filho, tam_dna, lim_sup, lim_inf, geracao, geracoes, amp_mut, prob_mutacao)

    retorno = ag.avaliacao (org, qtd_ind, tam_dna, individuos_salvos, objetivo)
    org = ag.ordena(org, tipo)

    if geracao % 50 == 0 and melhor > org['fitness'][-1]:
        imagem2 = dna2imagem(org['dna'][-1], imagem.shape)
        mpimg.imsave('Saida/saida'+ str(geracao) + '.png', imagem2)
        melhor = org['fitness'][-1]
    print('G: ', geracao ,' best fit: ', org['fitness'][-1])

print("fim")