# ========================================================
# script for microbial co-abundance network analysis
# uses SparCC network inference tool: https://github.com/JCSzamosi/SparCC3/archive/refs/heads/master.zip
# ========================================================


### Download the tools
# ===========================================
wget https://github.com/JCSzamosi/SparCC3/archive/refs/heads/master.zip
unzip master.zip
export PATH= home/spracc/SparCC3-master:$PATH


### Upload files
# ===========================================
mkdir /home/spracc/SparCC3-master/mydata
cd /home/spracc/SparCC3-master/mydata
mkdir basis_corr Resamplings Bootstraps pvalues

rz Species_data.txt


### Calculate Pearson r
# ===========================================
cd /home/spracc/SparCC3-master
python SparCC.py mydata/Species_data.txt -i 5 --cor_file=mydata/basis_corr/cor_sparcc.txt


### Bootstrap
# ===========================================
cd /home/spracc/SparCC3-master
python MakeBootstraps.py mydata/Species_data.txt -n 1000 -t permutation_#.txt -p mydata/Resamplings/


### Permutation
# ===========================================
cd /home/spracc/SparCC3-master
for i in 0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 40 41 42 43 44 45 46 47 48 49 50 51 52 53 54 55 56 57 58 59 60 61 62 63 64 65 66 67 68 69 70 71 72 73 74 75 76 77 78 79 80 81 82 83 84 85 86 87 88 89 90 91 92 93 94 95 96 97 98 99
do
python python SparCC.py mydata/Resamplings/permutation_$i.txt --cor_file=mydata/Bootstraps/sim_cor_$i.txt
done


### Two_sided p-value
# ===========================================
cd /home/spracc/SparCC3-master
python PseudoPvals.py mydata/basis_corr/cor_sparcc.txt mydata/Bootstraps/sim_cor_#.txt 100 -o mydata/pvalues/pvals_two_sided.txt -t two_sided