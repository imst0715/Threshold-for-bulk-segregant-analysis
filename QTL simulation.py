import numpy as np
import os

#Please check the path
path="F:\\myPhD_2015\\Materials and Methods\\QTL mapping\\simulation by Jenny"

os.chdir(path)

rounds = 10000 # number of rounds of simulation
pool1=16 #change to cells
pool2=20
cut1 = 0.9      # cut-off  p = np.percentile(a, cut*100)
cut2 = 0.95      # cut-off  p = np.percentile(a, cut*100)
cut3 = 0.99      # cut-off  p = np.percentile(a, cut*100)
cstart=20
cstop=1001

os.mkdir(str(pool1)+"_"+str(pool2)+" cells "+str(cstop-1)+"x")
os.mkdir(str(pool1)+"_"+str(pool2)+" cells "+str(cstop-1)+"x\\All")
os.mkdir(str(pool1)+"_"+str(pool2)+" cells "+str(cstop-1)+"x\\Select")
os.chdir(path+"\\"+str(pool1)+"_"+str(pool2)+" cells "+str(cstop-1)+"x")
#out1=open('all data.txt','w')
#out2=open('select data.txt','w')
out3=open('threshold_'+str(pool1)+"_"+str(pool2)+'.txt','w')#threshold for delta SNP
out4=open('threshold_snp_'+str(pool1)+"_"+str(pool2)+'.txt','w')# threshold for SNP  (1 pool)
out5=open('threshold_all_'+str(pool1)+"_"+str(pool2)+'.txt','w')# Not for analysis



allresultD={}
resultD={}
thrD={}
thrD_snp={}
thrD_all={}
for coverage in range (cstart,cstop):
    resultD[coverage]=[]
    allresultD[coverage]=[]
    ds=[]
    ds_all=[]
    allsnp1=[]
    allsnp2=[]
    for i in range (1,10001):
        #Simulate the allele frequency by random
        genotype1 = (np.random.choice([0,1], pool1))
        af1= list(genotype1).count(0)/float(pool1)#comparing to reference genotype
        genotype2= (np.random.choice([0,1], pool2))
        af2= list(genotype2).count(0)/float(pool2)#comparing to reference genotype
        #print (genotype1)
        #binomial distribution based on the coverage and allele frequency(af)
        s1 = np.random.binomial(1, af1, coverage)
        s2 = np.random.binomial(1, af2, coverage)
        #print (s1)
        snp1= list(s1).count(0)/float(coverage)#comparing to reference genotype
        snp2= list(s2).count(0)/float(coverage)#comparing to reference genotype
        dsnp=abs(snp1-snp2)
        ds_all.append(dsnp)
        allresultD[coverage].append([snp1,snp2,dsnp])
        allsnp1.append(snp1)
        allsnp2.append(snp2)
        if not ((snp1<=0.3) and (snp2<=0.3)):
            resultD[coverage].append([snp1,snp2,dsnp])
            ds.append(dsnp)
            
    data=np.array(ds)#data form ds, means the deltaSNP and the allele frequency from both bulk1 and bulk2 are larger than 0.3
    p1=np.percentile(data,cut1*100)
    p2=np.percentile(data,cut2*100)
    p3=np.percentile(data,cut3*100)
    thrD[coverage]=[p1,p2,p3]

    data2=np.array(allsnp1+allsnp2)#take both snp1 and snp2, and do not set the cutoff for SNP, since the simulation is 10000, I add 2 group, thus, the simulation time=20000
    ps1=np.percentile(data2,cut1*100)
    ps2=np.percentile(data2,cut2*100)
    ps3=np.percentile(data2,cut3*100)
    thrD_snp[coverage]=[ps1,ps2,ps3]

    data3=np.array(ds_all)#all delta SNP, do not use the allel frequency>0.3 as threshold 
    pa1=np.percentile(data3,cut1*100)
    pa2=np.percentile(data3,cut2*100)
    pa3=np.percentile(data3,cut3*100)
    thrD_all[coverage]=[pa1,pa2,pa3]        

for cov in allresultD.keys():
    out1=open("All\\"+str(cov)+"_all.txt",'w')
    for items in allresultD[cov]:
        out1.write(str(cov)+"\t")
        for item in items:
            out1.write(str(item)+"\t")
        out1.write("\n")
    out1.close()

for cov in resultD.keys():
    out2=open("Select\\"+str(cov)+"_select.txt",'w')
    for items in resultD[cov]:
        out2.write(str(cov)+"\t")
        for item in items:
            out2.write(str(item)+"\t")
        out2.write("\n")
    out2.close()

for cov in thrD:
    out3.write(str(cov)+"\t")
    for item in thrD[cov]:
        out3.write(str(item)+"\t")
    out3.write("\n")
        
out3.close()

for cov in thrD_snp:
    out4.write(str(cov)+"\t")
    for item in thrD_snp[cov]:
        out4.write(str(item)+"\t")
    out4.write("\n")    
out4.close()    

for cov in thrD_all:
    out5.write(str(cov)+"\t")
    for item in thrD_all[cov]:
        out5.write(str(item)+"\t")
    out5.write("\n")
out5.close()


