#!/bin/sh


python run.py data/stack.csv -a OpticsAnalysis -o ksi05m10e500 -p method="ksistep",ksi=0.05,M=10,eps=500
python run.py data/stack.csv -a OpticsAnalysis -o ksi05m20e500 -p method="ksistep",ksi=0.05,M=20,eps=500
python run.py data/stack.csv -a OpticsAnalysis -o ksi05m30e500 -p method="ksistep",ksi=0.05,M=30,eps=500

python run.py data/stack.csv -a OpticsAnalysis -o ksi025m10e500 -p method="ksistep",ksi=0.025,M=10,eps=500
python run.py data/stack.csv -a OpticsAnalysis -o ksi025m20e500 -p method="ksistep",ksi=0.025,M=20,eps=500
python run.py data/stack.csv -a OpticsAnalysis -o ksi025m30e500 -p method="ksistep",ksi=0.025,M=30,eps=500
             
    
