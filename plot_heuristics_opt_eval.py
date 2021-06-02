import pylab as py
import sys

py.rcParams['font.size'] = 15
filename = str(sys.argv[1])
filename2 = filename.replace("opt","optV4")
data = py.loadtxt("raw_results/"+filename,skiprows=1)
data2 = py.loadtxt("raw_results/"+filename2,skiprows=1)
data_max = py.loadtxt("raw_results/"+filename+"_MAX",skiprows=1)
data_min = py.loadtxt("raw_results/"+filename+"_MIN",skiprows=1)
data_max2 = py.loadtxt("raw_results/"+filename2+"_MAX",skiprows=1)
data_min2 = py.loadtxt("raw_results/"+filename2+"_MIN",skiprows=1)
data_headers = py.loadtxt("raw_results/"+filename,dtype=str,max_rows=1)
#data = py.loadtxt("ompss-arm-10M",skiprows=1)
#data_headers = py.loadtxt("ompss-arm-10M",dtype=str,max_rows=1)


py.figure(figsize=(20,5))

nb_algos = 10
nb_matrices = len(data_headers)

data_headers = py.append(data_headers, "Average")
posOpt = py.array([(nb_algos+1)*i-1 for i in range(nb_matrices)])
posOptC = py.array([(nb_algos+1)*i for i in range(nb_matrices)])
posNested = py.array([(nb_algos+1)*i+1 for i in range(nb_matrices)])
posD100 = py.array([(nb_algos+1)*i+2 for i in range(nb_matrices)])
posD300 = py.array([(nb_algos+1)*i+3 for i in range(nb_matrices)])
posD500 = py.array([(nb_algos+1)*i+4 for i in range(nb_matrices)])
posD700 = py.array([(nb_algos+1)*i+5 for i in range(nb_matrices)])
posNonested = py.array([(nb_algos+1)*i+6 for i in range(nb_matrices)])
posS50k = py.array([(nb_algos+1)*i+7 for i in range(nb_matrices)])
posS70k = py.array([(nb_algos+1)*i+8 for i in range(nb_matrices)])
posS90k = py.array([(nb_algos+1)*i+9 for i in range(nb_matrices)])

print(posNonested)

valNonested = data[8,:]
valNested = data[9,:]
valOpt = data[0,:]
valOptC = data2[0,:]
valD100 = data[1,:]
valD300 = data[2,:]
valD500 = data[3,:]
valD700 = data[4,:]
valS50k = data[5,:]
valS70k = data[6,:]
valS90k = data[7,:]


'''
minNonested = data_min[8,:]
minNested = data_min[9,:]
minOpt = data_min[0,:]
minOptC = data_min2[0,:]
minD100 = data_min[1,:]
minD300 = data_min[2,:]
minD500 = data_min[3,:]
minD700 = data_min[4,:]
minS50k = data_min[5,:]
minS70k = data_min[6,:]
minS90k = data_min[7,:]

maxNonested = data_max[8,:]
maxNested = data_max[9,:]
maxOpt = data_max[0,:]
maxOptC = data_max2[0,:]
maxD100 = data_max[1,:]
maxD300 = data_max[2,:]
maxD500 = data_max[3,:]
maxD700 = data_max[4,:]
maxS50k = data_max[5,:]
maxS70k = data_max[6,:]
maxS90k = data_max[7,:]

errorNonested = py.array([valNonested-minNonested,maxNonested-valNonested])
errorNested = py.array([valNested-minNested,maxNested-valNested])
errorOpt = py.array([valOpt-minOpt,maxOpt-valOpt])
errorOptC = py.array([valOptC-minOptC,maxOptC-valOptC])
errorD100 = py.array([valD100-minD100,maxD100-valD100])
errorD300 = py.array([valD300-minD300,maxD300-valD300])
errorD500 = py.array([valD500-minD500,maxD500-valD500])
errorD700 = py.array([valD700-minD700,maxD700-valD700])
errorS50k = py.array([valS50k-minS50k,maxS50k-valS50k])
errorS70k = py.array([valS70k-minS70k,maxS70k-valS70k])
errorS90k = py.array([valS90k-minS90k,maxS90k-valS90k])
'''

print(posOpt)
print(valNonested)
print(valOpt)

py.bar(posOpt+1,valNonested/valOpt,label='Opt-D',width=1.0,color='deepskyblue')
py.bar(posOptC+1,valNonested/valOptC,label='Opt-D-Cost',width=1.0,color='dodgerblue')
py.bar(posNested[valNested>0]+1,valNonested[valNested>0]/valNested[valNested>0],label=r'Nested',width=1.0,color='red')
py.bar(posD100+1,valNonested/valD100,label=r'Desc-100',width=1.0,color='orangered')
py.bar(posD300+1,valNonested/valD300,label=r'Desc-300',width=1.0,color='darkorange')
py.bar(posD500+1,valNonested/valD500,label=r'Desc-500',width=1.0,color='orange')
py.bar(posD700+1,valNonested/valD700,label=r'Desc-700',width=1.0,color='gold')
py.bar(posS50k,valNonested/valS50k,label=r'SN-50k',width=1.0,color='indigo')
py.bar(posS70k,valNonested/valS70k,label=r'SN-70k',width=1.0,color='darkviolet')
py.bar(posS90k,valNonested/valS90k,label=r'SN-90k',width=1.0,color='violet')

py.bar(posOpt[-1]+1+nb_algos+1, (valNonested/valOpt).mean(),width=1.0, color='deepskyblue')
py.bar(posOptC[-1]+1+nb_algos+1, (valNonested/valOptC).mean(),width=1.0, color='dodgerblue')
py.bar(posNested[-1]+1+nb_algos+1,(valNonested[valNested>0]/valNested[valNested>0]).mean(),width=1.0,color='red')
py.bar(posD100[-1]+1+nb_algos+1, (valNonested/valD100).mean(),width=1.0, color='orangered')
py.bar(posD300[-1]+1+nb_algos+1, (valNonested/valD300).mean(),width=1.0, color='darkorange')
py.bar(posD500[-1]+1+nb_algos+1, (valNonested/valD500).mean(),width=1.0, color='orange')
py.bar(posD700[-1]+1+nb_algos+1, (valNonested/valD700).mean(),width=1.0, color='gold')
py.bar(posS50k[-1]+nb_algos+1, (valNonested/valS50k).mean(),width=1.0, color='indigo')
py.bar(posS70k[-1]+nb_algos+1, (valNonested/valS70k).mean(),width=1.0, color='darkviolet')
py.bar(posS90k[-1]+nb_algos+1, (valNonested/valS90k).mean(),width=1.0, color='violet')

py.axhline(1.0, color='black',linestyle='--',linewidth=1)
py.ylim(bottom=0,top=2.5)
py.xlim([-2,(nb_algos+1)*(nb_matrices+1)])
py.xticks(py.append(posD500,posD500[-1]+1+nb_algos+1),data_headers,rotation='45')
py.xlabel("Matrices")
py.ylabel("Speedup")
py.legend(ncol=5)
py.tight_layout()
#py.title("Speedup of "+str(nb_matrices)+" matrices on MareNostrum4")
py.gcf().subplots_adjust(bottom=0.36,left=0.06)
py.savefig("figs/"+filename+"_S.pdf")
#py.title("Speedup of "+str(nb_matrices)+" matrices on CTE-ARM")
#py.savefig("ompss-arm_S.pdf")
py.clf()

'''
py.bar(posNested[valNested>0],valOpt[valNested>0]/valNested[valNested>0],label=r'Nested',width=1.0,color='red')
py.bar(posD100,valOpt/valD100,label=r'Desc-100',width=1.0,color='orangered')
py.bar(posD300,valOpt/valD300,label=r'Desc-300',width=1.0,color='darkorange')
py.bar(posD500,valOpt/valD500,label=r'Desc-500',width=1.0,color='orange')
py.bar(posD700,valOpt/valD700,label=r'Desc-700',width=1.0,color='gold')
py.bar(posNonested,valOpt/valNonested,label=r'Non-nested',width=1.0,color='limegreen')
py.bar(posS50k,valOpt/valS50k,label=r'SN-50k',width=1.0,color='indigo')
py.bar(posS70k,valOpt/valS70k,label=r'SN-70k',width=1.0,color='darkviolet')
py.bar(posS90k,valOpt/valS90k,label=r'SN-90k',width=1.0,color='violet')
py.axhline(1.0, color='black',linestyle='--',linewidth=1)
py.ylim(bottom=0.3,top=2)
py.yscale('log')
py.xlim([-2,(nb_algos+1)*(nb_matrices)])
py.xticks(posD500,data_headers,rotation='45')
py.xlabel("Matrices")
py.ylabel("Speedup")
py.legend(ncol=5)
py.tight_layout()
#py.title("Speedup of "+str(nb_matrices)+" matrices on MareNostrum4")
py.gcf().subplots_adjust(bottom=0.36,left=0.08)
py.savefig("figs/"+filename+"_Sopt.pdf")
#py.title("Speedup of "+str(nb_matrices)+" matrices on CTE-ARM")
#py.savefig("ompss-arm_S.pdf")
py.clf()

py.bar(posOpt,valOpt,label=r'Opt-D-Cost',width=1.0,yerr=errorOpt,color='dodgerblue')
py.bar(posNonested,valNonested,label=r'Non-nested',width=1.0,yerr=errorNonested,color='limegreen')
py.bar(posNested[valNested>0],valNested[valNested>0],label=r'Nested',width=1.0,yerr=errorNested[:,valNested>0],color='red')
py.bar(posD100,valD100,label=r'Desc-100',width=1.0,yerr=errorD100,color='orangered')
py.bar(posD300,valD300,label=r'Desc-300',width=1.0,yerr=errorD300,color='darkorange')
py.bar(posD500,valD500,label=r'Desc-500',width=1.0,yerr=errorD500,color='orange')
py.bar(posD700,valD700,label=r'Desc-700',width=1.0,yerr=errorD700,color='gold')
py.bar(posS50k,valS50k,label=r'SN-50k',width=1.0,yerr=errorS50k,color='indigo')
py.bar(posS70k,valS70k,label=r'SN-70k',width=1.0,yerr=errorS70k,color='darkviolet')
py.bar(posS90k,valS90k,label=r'SN-90k',width=1.0,yerr=errorS90k,color='violet')
py.xlim([-2,(nb_algos+1)*(nb_matrices)])
py.xticks(posD500,data_headers,rotation='45')
py.xlabel("Matrices")
py.ylabel("Time (s)")
py.legend(ncol=5)
#py.title("Time of "+str(nb_matrices)+" matrices on MareNostrum4")
py.gcf().subplots_adjust(bottom=0.36,left=0.06)
py.savefig("figs/"+filename+"_T.pdf")
#py.title("Time of "+str(nb_matrices)+" matrices on CTE-ARM")
#py.savefig("ompss-arm_T.pdf")
'''
