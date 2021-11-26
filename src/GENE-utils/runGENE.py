import numpy as np
import math
import pylab as py
import matplotlib.pyplot as plt
import matplotlib.ticker as mtick
from decimal import Decimal
from matplotlib.ticker import FormatStrFormatter
from matplotlib.ticker import ScalarFormatter
import sys
import glob
import os
from shutil import copyfile
import subprocess
import fnmatch
plt.switch_backend('agg')

#set up parameters for plotting:
plt.style.use('ggplot')
plt.rcParams['font.size'] = 20
plt.rcParams['axes.labelsize'] = 20
plt.rcParams['axes.labelweight'] = 'heavy'
plt.rcParams['xtick.labelsize'] = 15
plt.rcParams['ytick.labelsize'] = 15
plt.rcParams['legend.fontsize'] = 10
plt.rcParams['figure.titlesize'] = 12
plt.rcParams['text.usetex']=True

'''omega_star = (T*k_y)/(e*B*L_n) + (T*k_y)/(e*B*L_T) 
with B in Tesla, e in Coulomb, T in Joule, k_y in m^(-1), L_T, L_n in m.'''


src = "/global/homes/g/giannos/GENE/genecode/"
trgt = "/global/cscratch1/sd/giannos/gene/"
char = ["B","D","H"]
charn = ["I1","I2"]

fine = np.arange(0.05,1.00,0.05)
coarse = np.arange(1.00,2.1,0.1)
ky_list = list(fine)+list(coarse)

#ky_list = np.arange(0.05,2,0.1)

L_ref = 0.75434953
B_ref = 2.1650553
V_ref = 309499.0
#v/l = 410286
#In DIII-D shear is 8x10^5 rad/sec. Shear rate in kHz is 1.27x10^2. In this case, sh/(v/l)=0.31. In case we need the angular frequency of the shear, it turns out to be 1.95.
#shear =sh/(v/l) 

"""CMOD Values

L_ref = 0.269
B_ref = 5.291

"""
joule_conversion = 1.6022e-19
m_i = 2*1.67236232e-27
e = 1.6022e-19

omti_dict = {'B':10.42, 'D':20.84, 'H':5.21}#L_ref/L_Ti
omte_dict = {'B':28.21, 'D':56.42, 'H':14.1}#L_ref/L_Te
omn_dict = {'B':32.01, 'D':64.02, 'H':16.00}#L_ref/L_n

omti_ndict = {'B':10.42, 'I1':7.815, 'I2':15.63, 'D':20.84, 'H':5.21}#L_ref/L_Ti
omte_ndict = {'B':28.21, 'I1':21.16, 'I2':42.32, 'D':56.42, 'H':14.1}#L_ref/L_Te
omn_ndict = {'B':32.01, 'I1':24.00, 'I2':48.02, 'D':64.02, 'H':16.00}#L_ref/L_n


def example():
    char_list = ["A","B"]
    omti_dict = {'A':'10.0', 'B':'20.0'}
    omte_dict = {'A':'30.0', 'B':'60.0'}
    ky_list = [0.05,0.1]
    for x in range(1):
        for y in range(2):
            if not os.path.exists(trgt + "%s%s" % (char_list[x],char_list[y])):
                os.makedirs(trgt + "%s%s" % (char_list[x],char_list[y]))
            for ky in range(2):
                if not os.path.exists(trgt + "%s%s/%s" % (char_list[x],char_list[y],ky_list[ky])):
                    os.makedirs(trgt + "%s%s/%s" % (char_list[x],char_list[y],ky_list[ky]))
                    copyfile(src + "prob06/parameters", trgt + "%s%s/%s/parameters" % (char_list[x],char_list[y],ky_list[ky]))
                    copyfile(src + "prob06/gene_cori", trgt + "%s%s/%s/gene_cori" % (char_list[x],char_list[y],ky_list[ky]))
                    copyfile(src + "prob06/scanscript", trgt + "%s%s/%s/scanscript" % (char_list[x],char_list[y],ky_list[ky]))
                    copyfile(src + "prob06/launcher.cmd", trgt + "%s%s/%s/launcher.cmd" % (char_list[x],char_list[y],ky_list[ky]))
                    copyfile(src + "prob06/submit.cmd", trgt + "%s%s/%s/submit.cmd" % (char_list[x],char_list[y],ky_list[ky]))
                    copyfile(src + "prob06/submit_knl.cmd", trgt + "%s%s/%s/submit_knl.cmd" % (char_list[x],char_list[y],ky_list[ky]))
                    copyfile(src + "prob06/runit.sh", trgt + "%s%s/%s/runit.sh" % (char_list[x],char_list[y],ky_list[ky]))
                    f = open(trgt + "%s%s/%s/parameters" % (char_list[x],char_list[y],ky_list[ky]),'r')
                    lines=f.readlines()
                    lines[18] = 'kymin = '+ str(ky_list[ky])+'\n'
                    lines[27] = 'diagdir = ' + "\'" + trgt + "%s%s/%s/" % (char_list[x],char_list[y],ky_list[ky]) + "\'" + "\n"
                    lines[79] = 'omt = '+ str(omti_dict[char_list[x]]) + '\n'
                    lines[81] = 'omn = 100.00\n'
                    lines[90] = 'omt = '+ str(omte_dict[char_list[y]]) + '\n'
                    lines[92] = 'omn = 100.00\n'
                    f.close()
                    f = open(trgt + "%s%s/%s/parameters" % (char_list[x],char_list[y],ky_list[ky]),'w')
                    f.writelines(lines)
                    f.close()
                    #os.chdir(trgt + "%s%s/%s/" % (char_list[x],char_list[y],ky_list[ky]))
                    #os.system("sbatch runit.sh")
                    
                        
def base():
    ky_b = np.arange(0.1,3.7,0.1)
    if not os.path.exists(trgt + "ti262/base/"):
        os.makedirs(trgt + "ti262/base/")
    for ky in range(len(ky_b)):
        if not os.path.exists(trgt + "/ti262/base/%s" % (ky_b[ky])):
            os.makedirs(trgt + "ti262/base/%s" % (ky_b[ky]))
            copyfile(src + "prob06/parameters", trgt + "ti262/base/%s/parameters" % (ky_b[ky]))
            copyfile(src + "prob06/gene_cori", trgt + "ti262/base/%s/gene_cori" % (ky_b[ky]))
            copyfile(src + "prob06/scanscript", trgt + "ti262/base/%s/scanscript" % (ky_b[ky]))
            copyfile(src + "prob06/launcher.cmd", trgt + "ti262/base/%s/launcher.cmd" % (ky_b[ky]))
            copyfile(src + "prob06/submit.cmd", trgt + "ti262/base/%s/submit.cmd" % (ky_b[ky]))
            copyfile(src + "prob06/submit_knl.cmd", trgt + "ti262/base/%s/submit_knl.cmd" % (ky_b[ky]))
            copyfile(src + "prob06/runit.sh", trgt + "ti262/base/%s/runit.sh" % (ky_b[ky]))
            f = open(trgt + "ti262/base/%s/parameters" % (ky_b[ky]),'r')
            lines=f.readlines()
            lines[18] = 'kymin = '+ str(ky_b[ky])+'\n'
            lines[21] = 'x0 = 0.982\n'
            lines[25] = 'diagdir = ' + "\'" + trgt + "ti262/base/%s/" % (ky_b[ky]) + "\'" + "\n"
            lines[61] = 'geomdir = ' + "\'" + "/global/homes/g/giannos/GENE/genecode/prob02" + "\'" + "\n"
            lines[62] = 'geomfile = ' + "\'" + "g1100212023.01236_262" + "\'" + "\n"
            lines[73] = 'temp = 0.088\n'
            lines[74] = 'omt = 46.7\n'
            lines[76] = 'omn = 54.73\n'
            lines[84] = 'temp = 0.0335\n'
            lines[85] = 'omt = 97.15\n'
            lines[87] = 'omn = 54.73\n'
            lines[92] = 'nref = 1.9300000000000E+01\n'
            f.close()
            f = open(trgt + "ti262/base/%s/parameters" % (ky_b[ky]),'w')
            f.writelines(lines)
            f.close()
            os.chdir(trgt + "ti262/base/%s/" % (ky_b[ky]))
            os.system("sbatch runit.sh")
            
def base_shear():
    ky_b = ky_list#np.arange(0.1,3.7,0.1)
    if not os.path.exists(trgt + "ti255/base_shear_act/"):
        os.makedirs(trgt + "ti255/base_shear_act/")
    for ky in range(len(ky_b)):
        if not os.path.exists(trgt + "/ti255/base_shear_act/%s" % (ky_b[ky])):
            os.makedirs(trgt + "ti255/base_shear_act/%s" % (ky_b[ky]))
            copyfile(src + "prob06/parameters", trgt + "ti255/base_shear_act/%s/parameters" % (ky_b[ky]))
            copyfile(src + "prob06/gene_cori", trgt + "ti255/base_shear_act/%s/gene_cori" % (ky_b[ky]))
            copyfile(src + "prob06/scanscript", trgt + "ti255/base_shear_act/%s/scanscript" % (ky_b[ky]))
            copyfile(src + "prob06/launcher.cmd", trgt + "ti255/base_shear_act/%s/launcher.cmd" % (ky_b[ky]))
            copyfile(src + "prob06/submit.cmd", trgt + "ti255/base_shear_act/%s/submit.cmd" % (ky_b[ky]))
            copyfile(src + "prob06/submit_knl.cmd", trgt + "ti255/base_shear_act/%s/submit_knl.cmd" % (ky_b[ky]))
            copyfile(src + "prob06/runit.sh", trgt + "ti255/base_shear_act/%s/runit.sh" % (ky_b[ky]))
            f = open(trgt + "ti255/base_shear_act/%s/parameters" % (ky_b[ky]),'r')
            lines=f.readlines()
            lines[18] = 'kymin = '+ str(ky_b[ky])+'\n'
            lines[21] = 'x0 = 0.9673\n'
            lines[25] = 'diagdir = ' + "\'" + trgt + "ti255/base_shear_act/%s/" % (ky_b[ky]) + "\'" + "\n"
            lines[57] = 'ExBrate = 1.95\n'
            #lines[61] = 'geomdir = ' + "\'" + "/global/homes/g/giannos/GENE/genecode/prob02" + "\'" + "\n"
            #lines[62] = 'geomfile = ' + "\'" + "g1100212023.01236_262" + "\'" + "\n"
            lines[75] = 'temp = 0.475\n'
            lines[76] = 'omt = 10.42\n'
            lines[78] = 'omn = 32.01\n'
            lines[86] = 'temp = 0.4445\n'
            lines[87] = 'omt = 28.21\n'
            lines[89] = 'omn = 32.01\n'
            lines[94] = 'nref = 5.8000000000000E+00\n'
            f.close()
            f = open(trgt + "ti255/base_shear_act/%s/parameters" % (ky_b[ky]),'w')
            f.writelines(lines)
            f.close()
            os.chdir(trgt + "ti255/base_shear_act/%s/" % (ky_b[ky]))
            os.system("sbatch runit.sh")
            

def peak():
    ky_p1 = list(np.arange(0.01,0.1,0.01))
    ky_p2 = list(np.arange(0.1,5.1,0.1))
    ky_p = ky_p1+ky_p2
    if not os.path.exists(trgt + "peak/"):
        os.makedirs(trgt + "peak/")
    for ky in range(len(ky_p)):
        if not os.path.exists(trgt + "peak/%s" % (ky_p[ky])):
            os.makedirs(trgt + "peak/%s" % (ky_p[ky]))
            copyfile(src + "prob06/parameters", trgt + "peak/%s/parameters" % (ky_p[ky]))
            copyfile(src + "prob06/gene_cori", trgt + "peak/%s/gene_cori" % (ky_p[ky]))
            copyfile(src + "prob06/scanscript", trgt + "peak/%s/scanscript" % (ky_p[ky]))
            copyfile(src + "prob06/launcher.cmd", trgt + "peak/%s/launcher.cmd" % (ky_p[ky]))
            copyfile(src + "prob06/submit.cmd", trgt + "peak/%s/submit.cmd" % (ky_p[ky]))
            copyfile(src + "prob06/submit_knl.cmd", trgt + "peak/%s/submit_knl.cmd" % (ky_p[ky]))
            copyfile(src + "prob06/runit.sh", trgt + "peak/%s/runit.sh" % (ky_p[ky]))
            f = open(trgt + "peak/%s/parameters" % (ky_p[ky]),'r')
            lines=f.readlines()
            lines[18] = 'kymin = '+ str(ky_p[ky])+'\n'
            lines[21] = 'x0 = 0.85\n'
            lines[27] = 'diagdir = ' + "\'" + trgt + "peak/%s/" % (ky_p[ky]) + "\'" + "\n"
            lines[78] = 'temp = 0.746\n'
            lines[79] = 'omt = 4.573\n'
            lines[81] = 'omn = -0.553\n'
            lines[89] = 'temp = 0.6755\n'
            lines[90] = 'omt = 5.336\n'
            lines[92] = 'omn = -0.553\n'
            lines[97] = 'nref = 8.1060000000000000E+19\n'
            f.close()
            f = open(trgt + "peak/%s/parameters" % (ky_p[ky]),'w')
            f.writelines(lines)
            f.close()
            os.chdir(trgt + "peak/%s/" % (ky_p[ky]))
            os.system("sbatch runit.sh")

def peak_avg():
    ky_ph = ky_list
    
    if not os.path.exists(trgt + "peak_avg/"):
        os.makedirs(trgt + "peak_avg/")
    for ky in range(len(ky_ph)):
        if not os.path.exists(trgt + "peak_avg/%s" % (ky_ph[ky])):
            os.makedirs(trgt + "peak_avg/%s" % (ky_ph[ky]))
            copyfile(src + "prob06/parameters", trgt + "peak_avg/%s/parameters" % (ky_ph[ky]))
            copyfile(src + "prob06/gene_cori", trgt + "peak_avg/%s/gene_cori" % (ky_ph[ky]))
            copyfile(src + "prob06/scanscript", trgt + "peak_avg/%s/scanscript" % (ky_ph[ky]))
            copyfile(src + "prob06/launcher.cmd", trgt + "peak_avg/%s/launcher.cmd" % (ky_ph[ky]))
            copyfile(src + "prob06/submit.cmd", trgt + "peak_avg/%s/submit.cmd" % (ky_ph[ky]))
            copyfile(src + "prob06/submit_knl.cmd", trgt + "peak_avg/%s/submit_knl.cmd" % (ky_ph[ky]))
            copyfile(src + "prob06/runit.sh", trgt + "peak_avg/%s/runit.sh" % (ky_ph[ky]))
            f = open(trgt + "peak_avg/%s/parameters" % (ky_ph[ky]),'r')
            lines=f.readlines()
            lines[18] = 'kymin = '+ str(ky_ph[ky])+'\n'
            lines[21] = 'x0 = 0.93\n'
            lines[25] = 'diagdir = ' + "\'" + trgt + "peak_avg/%s/" % (ky_ph[ky]) + "\'" + "\n"
            lines[73] = 'temp = 0.475\n'
            lines[74] = 'omt = 10.42\n'
            lines[76] = 'omn = 32.01\n'
            lines[84] = 'temp = 0.4445\n'
            lines[85] = 'omt = 28.21\n'
            lines[87] = 'omn = 32.01\n'
            f.close()
            f = open(trgt + "peak_avg/%s/parameters" % (ky_ph[ky]),'w')
            f.writelines(lines)
            f.close()
            os.chdir(trgt + "peak_avg/%s/" % (ky_ph[ky]))
            os.system("sbatch runit.sh")

"""Scan of x_o values to locate proper value based on normalized psi_pol. In GENE x_o specifies rho_tor.
Only 5 jobs are allowed to be submitted in the debug partition every time."""          
def x_scan(): 
    x_list = np.linspace(0.91,0.935,5)
    #x_list = [0.967,0.9673,0.9676,0.968]
    
    if not os.path.exists(trgt + "ti262/x_scan/"):
        os.makedirs(trgt + "ti262/x_scan/")
    for x in x_list:
        if not os.path.exists(trgt + "ti262/x_scan/%s" %(x)):
            os.makedirs(trgt + "ti262/x_scan/%s" % (x))
            copyfile(src + "prob06/parameters", trgt + "ti262/x_scan/%s/parameters" % (x))
            copyfile(src + "prob06/gene_cori", trgt + "ti262/x_scan/%s/gene_cori" % (x))
            copyfile(src + "prob06/scanscript", trgt + "ti262/x_scan/%s/scanscript" % (x))
            copyfile(src + "prob06/launcher.cmd", trgt + "ti262/x_scan/%s/launcher.cmd" % (x))
            copyfile(src + "prob06/submit.cmd", trgt + "ti262/x_scan/%s/submit.cmd" % (x))
            copyfile(src + "prob06/submit_knl.cmd", trgt + "ti262/x_scan/%s/submit_knl.cmd" % (x))
            copyfile(src + "prob06/runit_debug.sh", trgt + "ti262/x_scan/%s/runit_debug.sh" % (x))
            f = open(trgt + "ti262/x_scan/%s/parameters" % (x),'r')
            lines=f.readlines()
            lines[1] = 'n_procs_s = 2\n'
            lines[2] = 'n_procs_v = 1\n'
            lines[3] = 'n_procs_w = 32\n'
            lines[4] = 'n_procs_x = 1\n'
            lines[5] = 'n_procs_y = 1\n'
            lines[6] = 'n_procs_z = 1\n'
            lines[7] = 'n_procs_sim = 64\n'         
            lines[18] = 'kymin = 0.25\n'
            lines[21] = 'x0 ='+ str(x)+ '\n'
            lines[25] = 'diagdir = ' + "\'" + trgt + "ti262/x_scan/%s/" % (x) + "\'" + "\n"
            lines[61] = 'geomdir = ' + "\'" + "/global/homes/g/giannos/GENE/genecode/prob02" + "\'" + "\n"
            lines[62] = 'geomfile = ' + "\'" + "g1100212023.01236_262" + "\'" + "\n"
            lines[73] = 'temp = 0.088\n'
            lines[74] = 'omt = 25.2\n'
            lines[76] = 'omn = 7.42\n'
            lines[84] = 'temp = 0.0335\n'
            lines[85] = 'omt = 32.8\n'
            lines[87] = 'omn = 7.42\n'
            lines[92] = 'nref = 1.9300000000000E+01\n'
            f.close()
            f = open(trgt + "ti262/x_scan/%s/parameters" % (x),'w')
            f.writelines(lines)
            f.close()
            os.chdir(trgt + "ti262/x_scan/%s/" % (x))
            os.system("sbatch runit_debug.sh")
            
def run():  
    for x in range(3):
        for y in range(3):
            for z in range(3):
                if not os.path.exists(trgt + "ti255/new/%s%s%s" % (char[x],char[y],char[z])):
                    os.makedirs(trgt + "ti255/new/%s%s%s" % (char[x],char[y],char[z]))
                for ky in range(len(ky_list)):
                    if not os.path.exists(trgt + "ti255/new/%s%s%s/%s" % (char[x],char[y],char[z],ky_list[ky])):
                        os.makedirs(trgt + "ti255/new/%s%s%s/%s" % (char[x],char[y],char[z],ky_list[ky]))
                    counter = 0
                    for file in os.listdir(trgt + "ti255/new/%s%s%s/%s" % (char[x],char[y],char[z],ky_list[ky])):
                        if fnmatch.fnmatch(file,'slurm*') == True:
                            counter = counter+1
                    if counter == 0:
                        os.remove(trgt + "ti255/new/%s%s%s/%s/runit.sh" % (char[x],char[y],char[z],ky_list[ky]))
                        copyfile(src + "prob06/runit.sh", trgt + "ti255/new/%s%s%s/%s/runit.sh" % (char[x],char[y],char[z],ky_list[ky]))
                        os.chdir(trgt + "ti255/new/%s%s%s/%s/" % (char[x],char[y],char[z],ky_list[ky]))
                        os.system("sbatch runit.sh")

                            #copyfile(src + "prob06/parameters", trgt + "ti255/new/%s%s%s/%s/parameters" % (char[x],char[y],char[z],ky_list[ky]))
                            #copyfile(src + "prob06/gene_cori", trgt + "ti255/new/%s%s%s/%s/gene_cori" % (char[x],char[y],char[z],ky_list[ky]))
                            #copyfile(src + "prob06/scanscript", trgt + "ti255/new/%s%s%s/%s/scanscript" % (char[x],char[y],char[z],ky_list[ky]))
                            #copyfile(src + "prob06/launcher.cmd", trgt + "ti255/new/%s%s%s/%s/launcher.cmd" % (char[x],char[y],char[z],ky_list[ky]))
                            #copyfile(src + "prob06/submit.cmd", trgt + "ti255/new/%s%s%s/%s/submit.cmd" % (char[x],char[y],char[z],ky_list[ky]))
                            #copyfile(src + "prob06/submit_knl.cmd", trgt + "ti255/new/%s%s%s/%s/submit_knl.cmd" % (char[x],char[y],char[z],ky_list[ky]))
                            #copyfile(src + "prob06/runit.sh", trgt + "ti255/new/%s%s%s/%s/runit.sh" % (char[x],char[y],char[z],ky_list[ky]))
                            #f = open(trgt + "ti255/new/%s%s%s/%s/parameters" % (char[x],char[y],char[z],ky_list[ky]),'r')
                            #lines=f.readlines()
                            #lines[18] = 'kymin = '+ str(ky_list[ky])+'\n'
                            #lines[21] = 'x0 = 0.9673\n'
                            #lines[25] = 'diagdir = ' + "\'" + trgt + "ti255/new/%s%s%s/%s/" % (char[x],char[y],char[z],ky_list[ky]) + "\'" + "\n"
                            #lines[73] = 'temp = 0.475\n'
                            #lines[74] = 'omt = '+ str(omti_dict[char[y]]) + '\n'
                            #lines[76] = 'omn = '+ str(omn_dict[char[x]]) + '\n'
                            #lines[84] = 'temp = 0.4445\n'
                            #lines[85] = 'omt = '+ str(omte_dict[char[z]]) + '\n'
                            #lines[87] = 'omn = '+ str(omn_dict[char[x]]) + '\n'
                            #f.close()
                            #f = open(trgt + "ti255/new/%s%s%s/%s/parameters" % (char[x],char[y],char[z],ky_list[ky]),'w')
                            #f.writelines(lines)
                            #f.close()
                        #os.chdir(trgt + "ti255/new/%s%s%s/%s/" % (char[x],char[y],char[z],ky_list[ky]))
                        #os.system("sbatch runit.sh")

                        
def run_inter_n():#gradn varies  
    for x in range(2):
        if not os.path.exists(trgt + "ti255/intermediate/%s%s%s" % (charn[x],"B","B")):
            os.makedirs(trgt + "ti255/intermediate/%s%s%s" % (charn[x],"B","B"))
            for ky in range(len(ky_list)):
                if not os.path.exists(trgt + "ti255/intermediate/%s%s%s/%s" % (charn[x],"B","B",ky_list[ky])):
                    os.makedirs(trgt + "ti255/intermediate/%s%s%s/%s" % (charn[x],"B","B",ky_list[ky]))
                    
                    copyfile(src + "prob06/parameters", trgt + "ti255/intermediate/%s%s%s/%s/parameters" % (charn[x],"B","B",ky_list[ky]))
                    copyfile(src + "prob06/gene_cori", trgt + "ti255/intermediate/%s%s%s/%s/gene_cori" % (charn[x],"B","B",ky_list[ky]))
                    copyfile(src + "prob06/scanscript", trgt + "ti255/intermediate/%s%s%s/%s/scanscript" % (charn[x],"B","B",ky_list[ky]))
                    copyfile(src + "prob06/launcher.cmd", trgt + "ti255/intermediate/%s%s%s/%s/launcher.cmd" % (charn[x],"B","B",ky_list[ky]))
                    copyfile(src + "prob06/submit.cmd", trgt + "ti255/intermediate/%s%s%s/%s/submit.cmd" % (charn[x],"B","B",ky_list[ky]))
                    copyfile(src + "prob06/submit_knl.cmd", trgt + "ti255/intermediate/%s%s%s/%s/submit_knl.cmd" % (charn[x],"B","B",ky_list[ky]))
                    copyfile(src + "prob06/runit.sh", trgt + "ti255/intermediate/%s%s%s/%s/runit.sh" % (charn[x],"B","B",ky_list[ky]))
                    f = open(trgt + "ti255/intermediate/%s%s%s/%s/parameters" % (charn[x],"B","B",ky_list[ky]),'r')
                    lines=f.readlines()
                    lines[18] = 'kymin = '+ str(ky_list[ky])+'\n'
                    lines[21] = 'x0 = 0.9673\n'
                    lines[25] = 'diagdir = ' + "\'" + trgt + "ti255/intermediate/%s%s%s/%s/" % (charn[x],"B","B",ky_list[ky]) + "\'" + "\n"
                    lines[73] = 'temp = 0.475\n'
                    lines[74] = 'omt = '+ str(omti_dict["B"]) + '\n'
                    lines[76] = 'omn = '+ str(omn_ndict[charn[x]]) + '\n'
                    lines[84] = 'temp = 0.4445\n'
                    lines[85] = 'omt = '+ str(omte_dict["B"]) + '\n'
                    lines[87] = 'omn = '+ str(omn_ndict[charn[x]]) + '\n'
                    f.close()
                    f = open(trgt + "ti255/intermediate/%s%s%s/%s/parameters" % (charn[x],"B","B",ky_list[ky]),'w')
                    f.writelines(lines)
                    f.close()
                    os.chdir(trgt + "ti255/intermediate/%s%s%s/%s/" % (charn[x],"B","B",ky_list[ky]))
                    os.system("sbatch runit.sh")

def run_inter_ti():#gradti varies  
    for x in range(2):
        if not os.path.exists(trgt + "ti255/intermediate/%s%s%s" % ("B",charn[x],"B")):
            os.makedirs(trgt + "ti255/intermediate/%s%s%s" % ("B",charn[x],"B"))
            for ky in range(len(ky_list)):
                if not os.path.exists(trgt + "ti255/intermediate/%s%s%s/%s" % ("B",charn[x],"B",ky_list[ky])):
                    os.makedirs(trgt + "ti255/intermediate/%s%s%s/%s" % ("B",charn[x],"B",ky_list[ky]))
                    
                    copyfile(src + "prob06/parameters", trgt + "ti255/intermediate/%s%s%s/%s/parameters" % ("B",charn[x],"B",ky_list[ky]))
                    copyfile(src + "prob06/gene_cori", trgt + "ti255/intermediate/%s%s%s/%s/gene_cori" % ("B",charn[x],"B",ky_list[ky]))
                    copyfile(src + "prob06/scanscript", trgt + "ti255/intermediate/%s%s%s/%s/scanscript" % ("B",charn[x],"B",ky_list[ky]))
                    copyfile(src + "prob06/launcher.cmd", trgt + "ti255/intermediate/%s%s%s/%s/launcher.cmd" % ("B",charn[x],"B",ky_list[ky]))
                    copyfile(src + "prob06/submit.cmd", trgt + "ti255/intermediate/%s%s%s/%s/submit.cmd" % ("B",charn[x],"B",ky_list[ky]))
                    copyfile(src + "prob06/submit_knl.cmd", trgt + "ti255/intermediate/%s%s%s/%s/submit_knl.cmd" % ("B",charn[x],"B",ky_list[ky]))
                    copyfile(src + "prob06/runit.sh", trgt + "ti255/intermediate/%s%s%s/%s/runit.sh" % ("B",charn[x],"B",ky_list[ky]))
                    f = open(trgt + "ti255/intermediate/%s%s%s/%s/parameters" % ("B",charn[x],"B",ky_list[ky]),'r')
                    lines=f.readlines()
                    lines[18] = 'kymin = '+ str(ky_list[ky])+'\n'
                    lines[21] = 'x0 = 0.9673\n'
                    lines[25] = 'diagdir = ' + "\'" + trgt + "ti255/intermediate/%s%s%s/%s/" % ("B",charn[x],"B",ky_list[ky]) + "\'" + "\n"
                    lines[73] = 'temp = 0.475\n'
                    lines[74] = 'omt = '+ str(omti_ndict[charn[x]]) + '\n'
                    lines[76] = 'omn = '+ str(omn_dict["B"]) + '\n'
                    lines[84] = 'temp = 0.4445\n'
                    lines[85] = 'omt = '+ str(omte_dict["B"]) + '\n'
                    lines[87] = 'omn = '+ str(omn_dict["B"]) + '\n'
                    f.close()
                    f = open(trgt + "ti255/intermediate/%s%s%s/%s/parameters" % ("B",charn[x],"B",ky_list[ky]),'w')
                    f.writelines(lines)
                    f.close()
                    os.chdir(trgt + "ti255/intermediate/%s%s%s/%s/" % ("B",charn[x],"B",ky_list[ky]))
                    os.system("sbatch runit.sh")

                    
def run_inter_te():#gradte varies  
    for x in range(2):
        if not os.path.exists(trgt + "ti255/intermediate/%s%s%s" % ("B","B",charn[x])):
            os.makedirs(trgt + "ti255/intermediate/%s%s%s" % ("B","B",charn[x]))
            for ky in range(len(ky_list)):
                if not os.path.exists(trgt + "ti255/intermediate/%s%s%s/%s" % ("B","B",charn[x],ky_list[ky])):
                    os.makedirs(trgt + "ti255/intermediate/%s%s%s/%s" % ("B","B",charn[x],ky_list[ky]))
                    
                    copyfile(src + "prob06/parameters", trgt + "ti255/intermediate/%s%s%s/%s/parameters" % ("B","B",charn[x],ky_list[ky]))
                    copyfile(src + "prob06/gene_cori", trgt + "ti255/intermediate/%s%s%s/%s/gene_cori" % ("B","B",charn[x],ky_list[ky]))
                    copyfile(src + "prob06/scanscript", trgt + "ti255/intermediate/%s%s%s/%s/scanscript" % ("B","B",charn[x],ky_list[ky]))
                    copyfile(src + "prob06/launcher.cmd", trgt + "ti255/intermediate/%s%s%s/%s/launcher.cmd" % ("B","B",charn[x],ky_list[ky]))
                    copyfile(src + "prob06/submit.cmd", trgt + "ti255/intermediate/%s%s%s/%s/submit.cmd" % ("B","B",charn[x],ky_list[ky]))
                    copyfile(src + "prob06/submit_knl.cmd", trgt + "ti255/intermediate/%s%s%s/%s/submit_knl.cmd" % ("B","B",charn[x],ky_list[ky]))
                    copyfile(src + "prob06/runit.sh", trgt + "ti255/intermediate/%s%s%s/%s/runit.sh" % ("B","B",charn[x],ky_list[ky]))
                    f = open(trgt + "ti255/intermediate/%s%s%s/%s/parameters" % ("B","B",charn[x],ky_list[ky]),'r')
                    lines=f.readlines()
                    lines[18] = 'kymin = '+ str(ky_list[ky])+'\n'
                    lines[21] = 'x0 = 0.9673\n'
                    lines[25] = 'diagdir = ' + "\'" + trgt + "ti255/intermediate/%s%s%s/%s/" % ("B","B",charn[x],ky_list[ky]) + "\'" + "\n"
                    lines[73] = 'temp = 0.475\n'
                    lines[74] = 'omt = '+ str(omti_dict["B"]) + '\n'
                    lines[76] = 'omn = '+ str(omn_dict["B"]) + '\n'
                    lines[84] = 'temp = 0.4445\n'
                    lines[85] = 'omt = '+ str(omte_ndict[charn[x]]) + '\n'
                    lines[87] = 'omn = '+ str(omn_dict["B"]) + '\n'
                    f.close()
                    f = open(trgt + "ti255/intermediate/%s%s%s/%s/parameters" % ("B","B",charn[x],ky_list[ky]),'w')
                    f.writelines(lines)
                    f.close()
                    os.chdir(trgt + "ti255/intermediate/%s%s%s/%s/" % ("B","B",charn[x],ky_list[ky]))
                    os.system("sbatch runit.sh")
                    
def base_an():
    rho_i = Larmor_radius(475)
    #om_norm = (V_ref/L_ref)*(1./(2.*math.pi))
    omg_norm = V_ref/L_ref
    krho_sim = 0.246
    #Doppler Shift to the plasma frame.
    omega_sim = -640000*(2.*math.pi)#(cycl/sec)
    k_perp = 122#(cycl/m)
    V_phase = omega_sim/k_perp
    V_pol = -1.25e4
    V_plasma_frame = V_phase-V_pol
    omega_plfr = V_plasma_frame*k_perp

    #case = "/global/cscratch1/sd/giannos/gene/base/"
    #case = "/global/cscratch1/sd/giannos/gene/peak/"
    #case = "/global/cscratch1/sd/giannos/gene/peak_high/"
    #case = "/global/cscratch1/sd/giannos/gene/ti255/base_shear_act"
    case = "/global/cscratch1/sd/giannos/gene/ti255/new/BBB/"
    sub_dir_list = os.listdir(case)
    sub_dir_list_num = []
    for sub_dir in sub_dir_list:
        sub_dir_list_num.append(float(sub_dir))

    ky_range = sorted(sub_dir_list_num)
    kys = []
    omegas = []
    gammas = []
    for ky in range(len(ky_range)):
        if os.path.exists(trgt + "ti255/new/BBB/%s/omega.dat" % (ky_range[ky])):
            if os.stat(trgt + "ti255/new/BBB/%s/omega.dat" % (ky_range[ky])).st_size > 0:
                f = open(trgt + "ti255/new/BBB/%s/omega.dat" % (ky_range[ky]),'r')
                for columns in (raw.strip().split() for raw in f):
                    kys.append(float(columns[0]))
                    gammas.append(float(columns[1])*omg_norm)
                    omegas.append(float(columns[2])*omg_norm)
                f.close()
            else:
                kys.append((ky_range[ky]))
                gammas.append(0)
                omegas.append(0)
        else:
            pass
    
    file = open("/global/cscratch1/sd/giannos/gene/ti255/new/BBB/base_gene(fig.6).txt",'a')
    file.write("ky"+"\t"+"gamma"+"\t"+"omega"+"\n")
    for i in range(len(kys)):
        file.write(str(kys[i])+"\t"+str(gammas[i])+"\t"+str(omegas[i])+"\n")
    file.close()

    
    fig,ax = plt.subplots()
    ax.set_title(r"Base case rerun, $\psi = 0.97$")    
    ax.plot(kys[:],gammas[:])
    ax.axvline(x=krho_sim)
    ax.set_xlim(0.0,2.0)
    ax.set_xlabel(r"$k_y\rho_i$")
    ax.set_ylabel(r"$\gamma$")
    ax.yaxis.set_major_formatter(FormatStrFormatter(('%.2e')))
    ax.yaxis.set_major_locator(plt.MaxNLocator(10))
    ax.xaxis.set_major_locator(plt.MaxNLocator(4))
    ax.yaxis.set_minor_locator(plt.MaxNLocator(50))
    ax.xaxis.set_minor_locator(plt.MaxNLocator(20))
    ax.grid(which = 'major',linestyle='-')
    #plt.show()
    plt.savefig(trgt+"ti255/new/BBB/test_gamma2.png", bbox_inches='tight')
    plt.close()
    fig,ax = plt.subplots()
    ax.set_title(r"Base case rerun, $\psi = 0.97$")
    ax.plot(kys[:],omegas[:])
    ax.plot(krho_sim, omega_plfr,'bo',markersize=14)
    ax.axvline(x=krho_sim)
    ax.set_xlim(0.0,2.0)
    ax.set_xlabel(r"$k_y\rho_i$")
    ax.set_ylabel(r"$\omega$")
    ax.set_ylim(-2.75e6,0.0)
    ax.yaxis.set_major_formatter(FormatStrFormatter(('%.2e')))
    ax.yaxis.set_major_locator(plt.MaxNLocator(10))
    ax.xaxis.set_major_locator(plt.MaxNLocator(4))
    ax.yaxis.set_minor_locator(plt.MaxNLocator(50))
    ax.xaxis.set_minor_locator(plt.MaxNLocator(20))
    ax.grid(which = 'major',linestyle='-')
    #plt.show()
    plt.savefig(trgt+"ti255/new/BBB/test_omega2.png", bbox_inches='tight')
    plt.close()
    #plt.xlim(0.0,4.0)
    #plt.xlabel(r"$k_y\rho_i$")
    #plt.ylabel(r"$\omega$")
    #plt.savefig(trgt+"peak_avg/test_omega.png", bbox_inches='tight')
    
def test():
    sub_dir_list = os.listdir("/global/cscratch1/sd/giannos/gene/ti255/new/BBB/0.05")
    sub_dir_list_num = []
    for sub_dir in sub_dir_list:
        sub_dir_list_num.append(float(sub_dir))
        
    return sorted(sub_dir_list_num)

def test2():
    counter = 0
    for file in os.listdir("/global/cscratch1/sd/giannos/gene/ti255/new/BBB/0.05/"):
        if fnmatch.fnmatch(file,'slurm*')==True:
            counter = counter + 1
        if counter == 1:
            print("File exists")
        
def analyze_all():
    rho_i = Larmor_radius(475)
    omg_norm = (V_ref/L_ref)
    krho_sim = 0.246
    f_sim = -640000 #in the direction of w_*e which is negative in GENE.
    for x in range(3):
        for y in range(3):
            for z in range(3):
                case = trgt +  "ti255/new/%s%s%s" % (char[x],char[y],char[z])
                sub_dir_list = os.listdir(case)
                sub_dir_list_num = []
                for sub_dir in sub_dir_list:
                    sub_dir_list_num.append(float(sub_dir))

                ky_range = sorted(sub_dir_list_num)
                kys = []
                omegas = []
                gammas = []
                unfinished = 0
                for ky in range(len(ky_list)):
                    if os.path.exists(trgt + "ti255/new/%s%s%s/%s/omega.dat" % (char[x],char[y],char[z],ky_range[ky])):
                        if os.stat(trgt + "ti255/new/%s%s%s/%s/omega.dat" % (char[x],char[y],char[z],ky_range[ky])).st_size > 0:
                            f = open(trgt + "ti255/new/%s%s%s/%s/omega.dat" % (char[x],char[y],char[z],ky_range[ky]),'r')
                            for columns in (raw.strip().split() for raw in f):
                                kys.append(float(columns[0]))
                                gammas.append(float(columns[1])*omg_norm)
                                omegas.append(float(columns[2])*omg_norm)
                            f.close()
                        else:
                            kys.append((ky_range[ky]))
                            gammas.append(0)
                            omegas.append(0)
                            unfinished = unfinished +1 
                    else:
                        pass
    
                '''new_file = open(trgt+"figures/%s%s%s" % (char[x],char[y],char[z]),'a')
                new_file.write("k_y"+"\t"+"gamma"+"\t"+"omega"+"\n")
                for i in range(0,len(kys)):
                    new_file.write(str(kys[i])+"\t"+str(gammas[i])+"\t"+str(omegas[i])+"\n")
                new_file.close()'''
                
                fig, ax = plt.subplots()
                ax.set_title(r"case %s%s%s, $\psi = 0.97$" %(char[x],char[y],char[z]))    
                ax.plot(kys[:],gammas[:])
                ax.axvline(x=krho_sim)
                ax.set_xlim(0.0,2.0)
                ax.set_xlabel(r"$k_y \rho_i$")
                ax.set_ylabel(r"$\gamma$")
                ax.yaxis.set_major_formatter(FormatStrFormatter(('%.02f')))
                ax.yaxis.set_major_locator(plt.MaxNLocator(10))
                ax.xaxis.set_major_locator(plt.MaxNLocator(4))
                ax.yaxis.set_minor_locator(plt.MaxNLocator(50))
                ax.xaxis.set_minor_locator(plt.MaxNLocator(20))
                ax.grid(which = 'major',linestyle='-')
                plt.savefig(trgt+"ti255/new/new_figures/norm/%s%s%s_gamma.png" %(char[x],char[y],char[z]))
                plt.close()
                fig, ax = plt.subplots()
                ax.set_title(r"case %s%s%s, $\psi = 0.97$" %(char[x],char[y],char[z])) 
                ax.plot(kys[:],omegas[:])
                ax.plot(krho_sim, f_sim,"*")
                ax.set_xlim(0.0,2.0)
                ax.set_xlabel(r"$k_y \rho_i$")
                ax.set_ylabel(r"$\omega$")
                ax.yaxis.set_major_formatter(FormatStrFormatter(('%.02f')))
                ax.yaxis.set_major_locator(plt.MaxNLocator(10))
                ax.xaxis.set_major_locator(plt.MaxNLocator(4))
                ax.yaxis.set_minor_locator(plt.MaxNLocator(50))
                ax.xaxis.set_minor_locator(plt.MaxNLocator(20))
                ax.grid(which = 'major',linestyle='-')
                plt.savefig(trgt+"ti255/new/new_figures/norm/%s%s%s_freq.png" %(char[x],char[y],char[z]))
                plt.close()      

def an_inter_n():
    rho_i = Larmor_radius(475)
    omg_norm = (V_ref/L_ref)
    krho_sim = 0.246
    f_sim = -640000 #in the direction of w_*e which is negative in GENE.
    for x in range(2):
        case = trgt +  "ti255/intermediate/%s%s%s" % (charn[x],"B","B")
        sub_dir_list = os.listdir(case)
        sub_dir_list_num = []
        for sub_dir in sub_dir_list:
            sub_dir_list_num.append(float(sub_dir))

        ky_range = sorted(sub_dir_list_num)
        kys = []
        omegas = []
        gammas = []
        unfinished = 0
        for ky in range(len(ky_list)):
            if os.path.exists(trgt + "ti255/intermediate/%s%s%s/%s/omega.dat" % (charn[x],"B","B",ky_range[ky])):
                if os.stat(trgt + "ti255/intermediate/%s%s%s/%s/omega.dat" % (charn[x],"B","B",ky_range[ky])).st_size > 0:
                    f = open(trgt + "ti255/intermediate/%s%s%s/%s/omega.dat" % (charn[x],"B","B",ky_range[ky]),'r')
                    for columns in (raw.strip().split() for raw in f):
                        kys.append(float(columns[0]))
                        #gammas.append(float(columns[1])*omg_norm)
                        #omegas.append(float(columns[2])*omg_norm)
                        gammas.append(float(columns[1]))
                        omegas.append(float(columns[2])*omg_norm)
                    f.close()
                else:
                    kys.append((ky_range[ky]))
                    gammas.append(0)
                    omegas.append(0)
                    unfinished = unfinished +1 
            else:
                pass
    
                '''new_file = open(trgt+"figures/%s%s%s" % (char[x],char[y],char[z]),'a')
                new_file.write("k_y"+"\t"+"gamma"+"\t"+"omega"+"\n")
                for i in range(0,len(kys)):
                    new_file.write(str(kys[i])+"\t"+str(gammas[i])+"\t"+str(omegas[i])+"\n")
                new_file.close()'''
                
        fig, ax = plt.subplots()
        ax.set_title(r"case %s%s%s, $\psi = 0.97$" %(charn[x],"B","B"))    
        ax.plot(kys[:],gammas[:])
        ax.axvline(x=krho_sim)
        ax.set_xlim(0.0,2.0)
        ax.set_xlabel(r"$k_y \rho_i$")
        ax.set_ylabel(r"$\gamma$")
        ax.yaxis.set_major_formatter(FormatStrFormatter(('%.02f')))
        ax.yaxis.set_major_locator(plt.MaxNLocator(10))
        ax.xaxis.set_major_locator(plt.MaxNLocator(4))
        ax.yaxis.set_minor_locator(plt.MaxNLocator(50))
        ax.xaxis.set_minor_locator(plt.MaxNLocator(20))
        ax.grid(which = 'major',linestyle='-')
        plt.savefig(trgt+"ti255/intermediate/new_figures/norm/%s%s%s_gamma.png" %(charn[x],"B","B"), bbox_inches='tight')
        plt.close()
        fig, ax = plt.subplots()
        ax.set_title(r"case %s%s%s, $\psi = 0.97$" %(charn[x],"B","B")) 
        ax.plot(kys[:],omegas[:])
        ax.plot(krho_sim, f_sim,"*")
        ax.set_xlim(0.0,2.0)
        ax.set_xlabel(r"$k_y \rho_i$")
        ax.set_ylabel(r"$\omega$")
        ax.yaxis.set_major_formatter(FormatStrFormatter(('%.02f')))
        ax.yaxis.set_major_locator(plt.MaxNLocator(10))
        ax.xaxis.set_major_locator(plt.MaxNLocator(4))
        ax.yaxis.set_minor_locator(plt.MaxNLocator(50))
        ax.xaxis.set_minor_locator(plt.MaxNLocator(20))
        ax.grid(which = 'major',linestyle='-')
        plt.savefig(trgt+"ti255/intermediate/new_figures/norm/%s%s%s_freq.png" %(charn[x],"B","B"), bbox_inches='tight')
        plt.close()      

        
def an_inter_Ti():
    rho_i = Larmor_radius(475)
    omg_norm = (V_ref/L_ref)
    krho_sim = 0.246
    f_sim = -640000 #in the direction of w_*e which is negative in GENE.
    for x in range(2):
        case = trgt +  "ti255/intermediate/%s%s%s" % ("B",charn[x],"B")
        sub_dir_list = os.listdir(case)
        sub_dir_list_num = []
        for sub_dir in sub_dir_list:
            sub_dir_list_num.append(float(sub_dir))

        ky_range = sorted(sub_dir_list_num)
        kys = []
        omegas = []
        gammas = []
        unfinished = 0
        for ky in range(len(ky_list)):
            if os.path.exists(trgt + "ti255/intermediate/%s%s%s/%s/omega.dat" % ("B",charn[x],"B",ky_range[ky])):
                if os.stat(trgt + "ti255/intermediate/%s%s%s/%s/omega.dat" % ("B",charn[x],"B",ky_range[ky])).st_size > 0:
                    f = open(trgt + "ti255/intermediate/%s%s%s/%s/omega.dat" % ("B",charn[x],"B",ky_range[ky]),'r')
                    for columns in (raw.strip().split() for raw in f):
                        kys.append(float(columns[0]))
                        #gammas.append(float(columns[1])*omg_norm)
                        #omegas.append(float(columns[2])*omg_norm)
                        gammas.append(float(columns[1]))
                        omegas.append(float(columns[2])*omg_norm)
                    f.close()
                else:
                    kys.append((ky_range[ky]))
                    gammas.append(0)
                    omegas.append(0)
                    unfinished = unfinished +1 
            else:
                pass
    
                '''new_file = open(trgt+"figures/%s%s%s" % (char[x],char[y],char[z]),'a')
                new_file.write("k_y"+"\t"+"gamma"+"\t"+"omega"+"\n")
                for i in range(0,len(kys)):
                    new_file.write(str(kys[i])+"\t"+str(gammas[i])+"\t"+str(omegas[i])+"\n")
                new_file.close()'''
                
        fig, ax = plt.subplots()
        ax.set_title(r"case %s%s%s, $\psi = 0.97$" %("B",charn[x],"B"))    
        ax.plot(kys[:],gammas[:])
        ax.axvline(x=krho_sim)
        ax.set_xlim(0.0,2.0)
        ax.set_xlabel(r"$k_y \rho_i$")
        ax.set_ylabel(r"$\gamma$")
        ax.yaxis.set_major_formatter(FormatStrFormatter(('%.02f')))
        ax.yaxis.set_major_locator(plt.MaxNLocator(10))
        ax.xaxis.set_major_locator(plt.MaxNLocator(4))
        ax.yaxis.set_minor_locator(plt.MaxNLocator(50))
        ax.xaxis.set_minor_locator(plt.MaxNLocator(20))
        ax.grid(which = 'major',linestyle='-')
        plt.savefig(trgt+"ti255/intermediate/new_figures/norm/%s%s%s_gamma.png" %("B",charn[x],"B"), bbox_inches='tight')
        plt.close()
        fig, ax = plt.subplots()
        ax.set_title(r"case %s%s%s, $\psi = 0.97$" %("B",charn[x],"B")) 
        ax.plot(kys[:],omegas[:])
        ax.plot(krho_sim, f_sim,"*")
        ax.set_xlim(0.0,2.0)
        ax.set_xlabel(r"$k_y \rho_i$")
        ax.set_ylabel(r"$\omega$")
        ax.yaxis.set_major_formatter(FormatStrFormatter(('%.02f')))
        ax.yaxis.set_major_locator(plt.MaxNLocator(10))
        ax.xaxis.set_major_locator(plt.MaxNLocator(4))
        ax.yaxis.set_minor_locator(plt.MaxNLocator(50))
        ax.xaxis.set_minor_locator(plt.MaxNLocator(20))
        ax.grid(which = 'major',linestyle='-')
        plt.savefig(trgt+"ti255/intermediate/new_figures/norm/%s%s%s_freq.png" %("B",charn[x],"B"), bbox_inches='tight')
        plt.close()      
        
        
def an_inter_Te():
    rho_i = Larmor_radius(475)
    omg_norm = (V_ref/L_ref)
    krho_sim = 0.246
    f_sim = -640000 #in the direction of w_*e which is negative in GENE.
    for x in range(2):
        case = trgt +  "ti255/intermediate/%s%s%s" % ("B","B",charn[x])
        sub_dir_list = os.listdir(case)
        sub_dir_list_num = []
        for sub_dir in sub_dir_list:
            sub_dir_list_num.append(float(sub_dir))

        ky_range = sorted(sub_dir_list_num)
        kys = []
        omegas = []
        gammas = []
        unfinished = 0
        for ky in range(len(ky_list)):
            if os.path.exists(trgt + "ti255/intermediate/%s%s%s/%s/omega.dat" % ("B","B",charn[x],ky_range[ky])):
                if os.stat(trgt + "ti255/intermediate/%s%s%s/%s/omega.dat" % ("B","B",charn[x],ky_range[ky])).st_size > 0:
                    f = open(trgt + "ti255/intermediate/%s%s%s/%s/omega.dat" % ("B","B",charn[x],ky_range[ky]),'r')
                    for columns in (raw.strip().split() for raw in f):
                        kys.append(float(columns[0]))
                        #gammas.append(float(columns[1])*omg_norm)
                        #omegas.append(float(columns[2])*omg_norm)
                        gammas.append(float(columns[1]))
                        omegas.append(float(columns[2])*omg_norm)
                    f.close()
                else:
                    kys.append((ky_range[ky]))
                    gammas.append(0)
                    omegas.append(0)
                    unfinished = unfinished +1 
            else:
                pass
    
                '''new_file = open(trgt+"figures/%s%s%s" % (char[x],char[y],char[z]),'a')
                new_file.write("k_y"+"\t"+"gamma"+"\t"+"omega"+"\n")
                for i in range(0,len(kys)):
                    new_file.write(str(kys[i])+"\t"+str(gammas[i])+"\t"+str(omegas[i])+"\n")
                new_file.close()'''
                
        fig, ax = plt.subplots()
        ax.set_title(r"case %s%s%s, $\psi = 0.97$" %("B","B",charn[x]))    
        ax.plot(kys[:],gammas[:])
        ax.axvline(x=krho_sim)
        ax.set_xlim(0.0,2.0)
        ax.set_xlabel(r"$k_y \rho_i$")
        ax.set_ylabel(r"$\gamma$")
        ax.yaxis.set_major_formatter(FormatStrFormatter(('%.02f')))
        ax.yaxis.set_major_locator(plt.MaxNLocator(10))
        ax.xaxis.set_major_locator(plt.MaxNLocator(4))
        ax.yaxis.set_minor_locator(plt.MaxNLocator(50))
        ax.xaxis.set_minor_locator(plt.MaxNLocator(20))
        ax.grid(which = 'major',linestyle='-')
        plt.savefig(trgt+"ti255/intermediate/new_figures/norm/%s%s%s_gamma.png" %("B","B",charn[x]), bbox_inches='tight')
        plt.close()
        fig, ax = plt.subplots()
        ax.set_title(r"case %s%s%s, $\psi = 0.97$" %("B","B",charn[x])) 
        ax.plot(kys[:],omegas[:])
        ax.plot(krho_sim, f_sim,"*")
        ax.set_xlim(0.0,2.0)
        ax.set_xlabel(r"$k_y \rho_i$")
        ax.set_ylabel(r"$\omega$")
        ax.yaxis.set_major_formatter(FormatStrFormatter(('%.02f')))
        ax.yaxis.set_major_locator(plt.MaxNLocator(10))
        ax.xaxis.set_major_locator(plt.MaxNLocator(4))
        ax.yaxis.set_minor_locator(plt.MaxNLocator(50))
        ax.xaxis.set_minor_locator(plt.MaxNLocator(20))
        ax.grid(which = 'major',linestyle='-')
        plt.savefig(trgt+"ti255/intermediate/new_figures/norm/%s%s%s_freq.png" %("B","B",charn[x]), bbox_inches='tight')
        plt.close()      
        
        
def Larmor_radius(T_i):
    om_cycl = (e*B_ref)/m_i
    v_ti = np.sqrt((2.*joule_conversion*T_i)/m_i)
    rho_i = v_ti/om_cycl
    return rho_i

def omega_star(T,k_y,L_n,L_T):
    omega = (T*joule_conversion*k_y)/(e*B_ref*L_n) + (T*joule_conversion*k_y)/(e*B_ref*L_T)
    freq = (1./(2.*math.pi))*omega
    return freq

