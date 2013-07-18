#!/bin/sh


python run.py data/stack.csv -a OpticsAnalysis -o ksi001m10 -p method="ksistep",ksi=0.01,M=10,eps=8000
python run.py data/stack.csv -a OpticsAnalysis -o ksi001m20 -p method="ksistep",ksi=0.01,M=20,eps=8000
python run.py data/stack.csv -a OpticsAnalysis -o ksi001m30 -p method="ksistep",ksi=0.01,M=30,eps=8000

python run.py data/stack.csv -a OpticsAnalysis -o ksi01m10 -p method="ksistep",ksi=0.1,M=10,eps=8000
python run.py data/stack.csv -a OpticsAnalysis -o ksi01m20 -p method="ksistep",ksi=0.1,M=20,eps=8000
python run.py data/stack.csv -a OpticsAnalysis -o ksi01m30 -p method="ksistep",ksi=0.1,M=30,eps=8000

python run.py data/stack.csv -a OpticsAnalysis -o ksi01m10e2000 -p method="ksistep",ksi=0.1,M=10,eps=2000
python run.py data/stack.csv -a OpticsAnalysis -o ksi01m20e2000 -p method="ksistep",ksi=0.1,M=20,eps=2000
python run.py data/stack.csv -a OpticsAnalysis -o ksi01m30e2000 -p method="ksistep",ksi=0.1,M=30,eps=2000

python run.py data/stack.csv -a OpticsAnalysis -o ksi001m10e2000 -p method="ksistep",ksi=0.01,M=10,eps=2000
python run.py data/stack.csv -a OpticsAnalysis -o ksi001m20e2000 -p method="ksistep",ksi=0.01,M=20,eps=2000
python run.py data/stack.csv -a OpticsAnalysis -o ksi001m30e2000 -p method="ksistep",ksi=0.01,M=30,eps=2000


python run.py data/stack.csv -a OpticsAnalysis -o ksi0001m10 -p method="ksistep",ksi=0.001,M=10,eps=8000
python run.py data/stack.csv -a OpticsAnalysis -o ksi0001m20 -p method="ksistep",ksi=0.001,M=20,eps=8000
python run.py data/stack.csv -a OpticsAnalysis -o ksi0001m30 -p method="ksistep",ksi=0.001,M=30,eps=8000

python run.py data/stack.csv -a OpticsAnalysis -o ksi0001m10e2000 -p method="ksistep",ksi=0.001,M=10,eps=2000
python run.py data/stack.csv -a OpticsAnalysis -o ksi0001m20e2000 -p method="ksistep",ksi=0.001,M=20,eps=2000
python run.py data/stack.csv -a OpticsAnalysis -o ksi0001m30e2000 -p method="ksistep",ksi=0.001,M=30,eps=2000

 
              
             
             
    
