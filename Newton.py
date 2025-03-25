#from numba import jit
import equations
import math
import numpy as np

#@jit(nopython=True) # Set "nopython" mode for best performance, equivalent to @njit
def Jac(func,x,k):

    epsfcn = 1.0e-3
    epsmch = 1.0e-3
    eps = math.sqrt(max(epsfcn,epsmch))

    Jfun=[[0 for i in range(len(x))] for i in range (len(x))]


    wa1=func(x,k)
    for i in range(len(x)):
        temp=x[i]
        h=eps*abs(temp)
        if (h==0):
            h=eps
        x[i]=temp+h
        fvec=func(x,k)
        x[i]=temp
        for j in range(len(x)):
            Jfun[j][i]=(-wa1[j]+fvec[j])/h

    return Jfun

#@jit(nopython=True) # Set "nopython" mode for best performance, equivalent to @njit
def newton(func,x,k,type,spec = 'rat',gender = 'male',sup_jux = 'sup',diabete='Non',inhib=None):

    fun=equations.conservation_eqs
    f = np.matrix(fun(x,k))
    # print(np.linalg.norm(f))
    TOLpcn = 1
    i = 1
    iter=0
    while(np.linalg.norm(f) > 0.0001) and (iter<150): #(iter<300): #male: (iter<300)  female: (iter<100)
#        print("Iteration Times: " + str(i) + " with TOL " + str(TOLpcn) + "%")
        i += 1
        J = np.matrix(Jac(fun,x,k))
        IJ = J.I
        F = np.matrix(fun(x,k))

        if spec == 'rat':
            if type=='DCT':
                amp = 1
            elif type == 'CNT':
                if gender =='female':
                    if sup_jux =='sup':
                        if np.linalg.norm(f)>100: # sup: 100
                            amp = 0.17# sup: 0.5
                        else:
                            amp=1.0
                    elif sup_jux == 'jux1':
                        if np.linalg.norm(f)>5000:
                            if k == 0:
                                amp = 0.5
                            else:
                                amp = 0.17
                        else:
                            amp = 1.0
                    elif sup_jux == 'jux2':
                        if np.linalg.norm(f)>5000:
                            amp = 0.13
                        else:
                            amp = 1.0
                    elif sup_jux == 'jux3':
                        if np.linalg.norm(f)>5000:
                            amp = 0.13
                        else:
                            amp = 1.0
                    elif sup_jux == 'jux4':
                        if np.linalg.norm(f)>5000:
                            amp = 0.13
                        else:
                            amp = 1.0
                    elif sup_jux == 'jux5':
                        if np.linalg.norm(f)>5000:
                            amp = 0.13
                        else:
                            amp = 1.0
                    else:
                        if np.linalg.norm(f)>5000:
                            amp = 0.15
                        else:
                            amp = 1.0
                elif gender =='male':
                    if sup_jux =='sup':
                        if np.linalg.norm(f)>5000: # sup: 100
                            amp = 0.5# sup: 0.5; saline: 0.5;
                        else:
                            amp=0.8
                    elif sup_jux == 'jux1':
                        if np.linalg.norm(f)>5000:
                            amp = 0.13 # nhe50: 0.3;
                        else:
                            amp = 1.0
                    elif sup_jux == 'jux2':
                        if np.linalg.norm(f)>5000:
                            amp = 0.13 # nhe50: 0.15 ncc: 0.05
                        else:
                            amp = 1.0
                    elif sup_jux == 'jux3':
                        if np.linalg.norm(f)>5000:
                            amp = 0.13 # nhe50: 0.15
                        else:
                            amp = 1.0
                    elif sup_jux == 'jux4':
                        if np.linalg.norm(f)>5000:
                            amp = 0.13 # nhe50: 0.15
                        else:
                            amp = 1.0
                    elif sup_jux == 'jux5':
                        if np.linalg.norm(f)>5000:
                            amp = 0.13 # nhe50: 0.15 nhe80:0.1
                        else:
                            amp = 1.0
                    else:
                        if np.linalg.norm(f)>5000:
                            amp = 0.15
                        else:
                            amp = 1.0
            elif type == 'SDL':
                amp = 0.2
            elif type == 'IMCD':
                if gender == 'female':
                    if np.linalg.norm(f)>1000: # 100
                        if k==0:
                            amp = 0.2
                        else:
                            amp = 0.2
                    else:
                        if k==0:
                            amp = 0.8#1.0 male:0.7 female:0.8
                        else:
                            amp = 0.8
                elif gender == 'male':
                    if np.linalg.norm(f)>5000:
                        if k==0:
                            amp = 0.1 # saline: 0.17
                        else:
                            amp = 0.2
                    else:
                        if k==0:
                            amp = 1.0
                        else:
                            amp = 1.0
            elif type == 'CCD':
                if np.linalg.norm(f)>1000:
                    if k == 0:
                        amp = 0.5#0.005 male:0.5 female:0.005
                    else:
                        amp = 0.1
                else:
                    amp = 0.8
            elif type == 'OMCD':
                if np.linalg.norm(f)>100:
                    amp = 0.8#0.8 #normal male/female: 1.0, diabetic male: 0.8
                else:
                    amp = 1.0# normal male: 1.0, diabetic male: 0.8
            elif type == 'cTAL' or type == 'MD':
                if np.linalg.norm(f)>100:
                    amp = 0.2
                else:
                    amp = 0.8 # normal male: 1.0
            elif type == 'mTAL':
                if np.linalg.norm(f)>100:
                    amp = 0.2 #0.2
                else:
                    amp = 1.0
            elif type == 'LDL':
                if np.linalg.norm(f)>5000:
                    amp = 0.5
                else:
                    amp = 1.0
            elif type == 'LAL':
                if np.linalg.norm(f)>5000:
                    amp = 0.5
                else:
                    amp = 1.0
            else:
                amp = 1
        elif spec == 'hum' and inhib == 'ACE':
            if type == 'S3':
                amp = 1.0
            elif type == 'SDL':
                amp = 1.0
            elif type == 'mTAL':
                if np.linalg.norm(f)>100:
                    amp = 0.2
                else:
                    amp = 1.0
            elif type == 'cTAL' or type == 'MD':
                if np.linalg.norm(f)>100:
                    amp = 0.2
                else:
                    amp = 0.8
            elif type=='DCT':
                if gender == 'male':
                    amp = 0.5
                elif gender == 'female':
                    amp = 1.0
            elif type == 'CNT':
                if gender == 'male':
                    if np.linalg.norm(f)>100:
                        amp = 0.5
                    else:
                        amp=0.8
                elif gender == 'female':
                    if sup_jux == 'sup':
                        if np.linalg.norm(f)>100:
                            amp = 0.5
                        else:
                            amp = 0.8
                    elif sup_jux == 'jux1':
                        if np.linalg.norm(f)>1000:
                            if k==0:
                                amp = 0.09
                            else:
                                amp = 0.13
                        else:
                            amp = 0.8
                    elif sup_jux == 'jux2':
                        if np.linalg.norm(f)>1000:
                            if k==0:
                                amp = 0.5
                            else:
                                amp = 0.5
                        else:
                            amp = 0.8
                    elif sup_jux == 'jux3':
                        if np.linalg.norm(f)>1000:
                            if k==0:
                                amp = 0.5
                            else:
                                amp = 0.5
                        else:
                            amp = 0.8
                    elif sup_jux == 'jux4':
                        if np.linalg.norm(f)>1000:
                            if k==0:
                                amp = 0.5
                            else:
                                amp = 0.5
                        else:
                            amp = 0.8
                    elif sup_jux == 'jux5':
                        if np.linalg.norm(f)>1000:
                            if k==0:
                                amp = 0.1
                            else:
                                amp = 0.17
                        else:
                            amp = 1.0
            elif type == 'CCD':
                if gender == 'male':
                    if np.linalg.norm(f)>100:
                        if k == 0:
                            amp = 0.8
                        else:
                            amp = 0.2
                    else:
                        amp = 0.8
                elif gender == 'female':
                    if np.linalg.norm(f)>100:
                        if k == 0:
                            amp = 0.4
                        else:
                            amp = 0.2
                    else:
                        amp = 0.8
            elif type == 'OMCD':
                if np.linalg.norm(f)>100:
                    amp = 0.5 #male: 0.5 female:0.8 (0.5 works for male and female)
                else:
                    amp = 0.8#male:0.8 female:
            elif type == 'IMCD':
                if np.linalg.norm(f)>100:
                    if k==0:
                        amp = 0.2 # male:0.1 female:0.2
                    else:
                        amp = 0.1 # male:0.2 female:0.2
                else:
                    if k==0:
                        amp = 1.5# male:0.5 female:0.5
                    else:
                        amp = 0.9 # male:0.5 female:0.5
            else:
                amp = 1.0
        elif spec == 'hum' and diabete == 'Non':
            if type == 'S3':
                amp = 1.0
            elif type == 'SDL':
                amp = 1.0
            elif type == 'mTAL':
                if np.linalg.norm(f)>100:
                    amp = 0.2
                else:
                    amp = 1.0
            elif type == 'cTAL' or type == 'MD':
                if np.linalg.norm(f)>100:
                    amp = 0.2
                else:
                    amp = 0.8
            elif type=='DCT':
                if gender == 'female' and sup_jux == 'jux1':
                    amp = 0.9
                elif gender == 'female' and sup_jux == 'jux3':
                    amp = 0.7
                elif gender == 'female' and sup_jux == 'jux2':
                    amp = 0.7
                elif gender == 'female' and sup_jux == 'sup':
                    if np.linalg.norm(f)>2000:
                        amp = 0.8#0.5
                    else:
                        amp = 0.5
                else:
                    amp = 0.5
            elif type == 'CNT':
                if gender == 'male':
                    if np.linalg.norm(f)>100:
                        amp = 0.5
                    else:
                        amp=0.8
                elif gender == 'female':
                    if sup_jux == 'sup':
                        if np.linalg.norm(f)>100:
                            amp = 0.5
                        else:
                            amp = 0.8
                    elif sup_jux == 'jux1':
                        if np.linalg.norm(f)>1000:
                            if k==0:
                                amp = 0.5
                            else:
                                amp = 0.13
                        else:
                            amp = 0.8
                    elif sup_jux == 'jux2':
                        if np.linalg.norm(f)>1000:
                            if k==0:
                                amp = 0.5
                            else:
                                amp = 0.5
                        else:
                            amp = 0.8
                    elif sup_jux == 'jux3':
                        if np.linalg.norm(f)>1000:
                            if k==0:
                                amp = 0.5
                            else:
                                amp = 0.5
                        else:
                            amp = 0.8
                    elif sup_jux == 'jux4':
                        if np.linalg.norm(f)>1000:
                            if k==0:
                                amp = 0.3
                            else:
                                amp = 0.5
                        else:
                            amp = 0.8
                    elif sup_jux == 'jux5':
                        if np.linalg.norm(f)>1000:
                            if k==0:
                                amp = 0.1
                            else:
                                amp = 0.17
                        else:
                            amp = 1.0
            elif type == 'CCD':
                if gender == 'male':
                    if np.linalg.norm(f)>100:
                        if k == 0:
                            amp = 0.5
                        else:
                            amp = 0.2
                    else:
                        amp = 0.8
                elif gender == 'female':
                    if np.linalg.norm(f)>100:
                        if k == 0:
                            amp = 0.2
                        else:
                            amp = 0.2
                    else:
                        amp = 0.8
            elif type == 'OMCD':
                if gender == 'male':
                    if np.linalg.norm(f)>100:
                        amp = 0.5 #male: 0.5 female:0.8 (0.5 works for male and female)
                    else:
                        amp = 0.8#male:0.8 female:
                elif gender == 'female':
                    if np.linalg.norm(f)>100:
                        amp = 0.5 #male: 0.5 female:0.8 (0.5 works for male and female)
                    else:
                        amp = 0.8#male:0.8 female:
            elif type == 'IMCD':
                if gender == 'female':
                    if np.linalg.norm(f)>100:
                        if k==0:
                            amp = 0.2 # male:0.1 female:0.2
                        else:
                            amp = 0.1 # male:0.2 female:0.2
                    else:
                        if k==0:
                            amp = 0.5# male:0.5 female:0.5
                        else:
                            amp = 0.5 # male:0.5 female:0.5
                elif gender == 'male':
                    if np.linalg.norm(f)>100:
                        if k==0:
                            amp = 0.27#0.19 # male:0.1 female:0.2
                        else:
                            amp = 0.2 # male:0.2 female:0.2
                    else:
                        if k==0:
                            amp = 0.8# male:0.5 female:0.5
                        else:
                            amp = 0.8 # male:0.5 female:0.5
            else:
                amp = 1.0
        elif spec == 'hum' and diabete == 'Moderate' and inhib != 'SGLT2':
            if type == 'S3':
                amp = 1.0
            elif type == 'SDL':
                amp = 1.0
            elif type == 'mTAL':
                if np.linalg.norm(f)>100:
                    amp = 0.2
                else:
                    amp = 1.0
            elif type == 'cTAL' or type == 'MD':
                if np.linalg.norm(f)>100:
                    amp = 0.2
                else:
                    amp = 0.8
            elif type=='DCT':
                if gender == 'female' and sup_jux == 'jux1':
                    amp = 0.9
                elif gender == 'female' and sup_jux == 'jux3':
                    amp = 0.7
                elif gender == 'female' and sup_jux == 'jux2':
                    amp = 0.7
                elif gender == 'female' and sup_jux == 'sup':
                    if np.linalg.norm(f)>2000:
                        amp = 0.8#0.5
                    else:
                        amp = 0.5
                else:
                    amp = 0.5
            elif type == 'CNT':
                if gender == 'male':
                    if np.linalg.norm(f)>100:
                        amp = 0.5
                    else:
                        amp=0.8
                elif gender == 'female':
                    if sup_jux == 'sup':
                        if np.linalg.norm(f)>100:
                            amp = 0.5
                        else:
                            amp = 0.8
                    elif sup_jux == 'jux1':
                        if np.linalg.norm(f)>1000:
                            if k==0:
                                amp = 0.5
                            else:
                                amp = 0.13
                        else:
                            amp = 0.8
                    elif sup_jux == 'jux2':
                        if np.linalg.norm(f)>1000:
                            if k==0:
                                amp = 0.5
                            else:
                                amp = 0.5
                        else:
                            amp = 0.8
                    elif sup_jux == 'jux3':
                        if np.linalg.norm(f)>1000:
                            if k==0:
                                amp = 0.5
                            else:
                                amp = 0.5
                        else:
                            amp = 0.8
                    elif sup_jux == 'jux4':
                        if np.linalg.norm(f)>1000:
                            if k==0:
                                amp = 0.3
                            else:
                                amp = 0.5
                        else:
                            amp = 0.8
                    elif sup_jux == 'jux5':
                        if np.linalg.norm(f)>1000:
                            if k==0:
                                amp = 0.1
                            else:
                                amp = 0.17
                        else:
                            amp = 1.0
            elif type == 'CCD':
                if gender == 'male':
                    if np.linalg.norm(f)>100:
                        if k == 0:
                            amp = 0.5
                        else:
                            amp = 0.2
                    else:
                        amp = 0.8
                elif gender == 'female':
                    if np.linalg.norm(f)>100:
                        if k == 0:
                            amp = 0.2
                        else:
                            amp = 0.2
                    else:
                        amp = 0.8
            elif type == 'OMCD':
                if gender == 'male':
                    if np.linalg.norm(f)>100:
                        amp = 0.5 #male: 0.5 female:0.8 (0.5 works for male and female)
                    else:
                        amp = 0.8#male:0.8 female:
                elif gender == 'female':
                    if np.linalg.norm(f)>100:
                        amp = 0.5 #male: 0.5 female:0.8 (0.5 works for male and female)
                    else:
                        amp = 0.8#male:0.8 female:
            elif type == 'IMCD':
                if gender == 'female':
                    if np.linalg.norm(f)>100:
                        if k==0:
                            amp = 0.2 # male:0.1 female:0.2
                        else:
                            amp = 0.1 # male:0.2 female:0.2
                    else:
                        if k==0:
                            amp = 0.5# male:0.5 female:0.5
                        else:
                            amp = 0.5 # male:0.5 female:0.5
                elif gender == 'male':
                    if np.linalg.norm(f)>100:
                        if k==0:
                            amp = 0.19#0.19 # male:0.1 female:0.2
                        else:
                            amp = 0.17 # male:0.2 female:0.2
                    else:
                        if k==0:
                            amp = 0.8# male:0.5 female:0.5
                        else:
                            amp = 0.8 # male:0.5 female:0.5
            else:
                amp = 1.0
        elif spec == 'hum' and diabete == 'Moderate' and inhib == 'SGLT2':
            if type == 'S3':
                amp = 1.0
            elif type == 'SDL':
                amp = 1.0
            elif type == 'mTAL':
                if np.linalg.norm(f)>100:
                    amp = 0.2
                else:
                    amp = 1.0
            elif type == 'cTAL' or type == 'MD':
                if np.linalg.norm(f)>100:
                    amp = 0.2
                else:
                    amp = 0.8
            elif type=='DCT':
                if gender == 'female' and sup_jux == 'jux1':
                    amp = 0.9
                elif gender == 'female' and sup_jux == 'jux3':
                    amp = 0.7
                elif gender == 'female' and sup_jux == 'jux2':
                    amp = 0.7
                elif gender == 'female' and sup_jux == 'sup':
                    if np.linalg.norm(f)>2000:
                        amp = 0.8#0.5
                    else:
                        amp = 0.5
                else:
                    amp = 0.5
            elif type == 'CNT':
                if gender == 'male':
                    if np.linalg.norm(f)>100:
                        amp = 0.5
                    else:
                        amp=0.8
                elif gender == 'female':
                    if sup_jux == 'sup':
                        if np.linalg.norm(f)>100:
                            amp = 0.5
                        else:
                            amp = 0.8
                    elif sup_jux == 'jux1':
                        if np.linalg.norm(f)>1000:
                            if k==0:
                                amp = 0.5
                            else:
                                amp = 0.13
                        else:
                            amp = 0.8
                    elif sup_jux == 'jux2':
                        if np.linalg.norm(f)>1000:
                            if k==0:
                                amp = 0.5
                            else:
                                amp = 0.5
                        else:
                            amp = 0.8
                    elif sup_jux == 'jux3':
                        if np.linalg.norm(f)>1000:
                            if k==0:
                                amp = 0.5
                            else:
                                amp = 0.5
                        else:
                            amp = 0.8
                    elif sup_jux == 'jux4':
                        if np.linalg.norm(f)>1000:
                            if k==0:
                                amp = 0.3
                            else:
                                amp = 0.5
                        else:
                            amp = 0.8
                    elif sup_jux == 'jux5':
                        if np.linalg.norm(f)>1000:
                            if k==0:
                                amp = 0.1
                            else:
                                amp = 0.17
                        else:
                            amp = 1.0
            elif type == 'CCD':
                if gender == 'male':
                    if np.linalg.norm(f)>100:
                        if k == 0:
                            amp = 0.5
                        else:
                            amp = 0.2
                    else:
                        amp = 0.8
                elif gender == 'female':
                    if np.linalg.norm(f)>100:
                        if k == 0:
                            amp = 0.2
                        else:
                            amp = 0.2
                    else:
                        amp = 0.8
            elif type == 'OMCD':
                if gender == 'male':
                    if np.linalg.norm(f)>100:
                        amp = 0.5 #male: 0.5 female:0.8 (0.5 works for male and female)
                    else:
                        amp = 0.8#male:0.8 female:
                elif gender == 'female':
                    if np.linalg.norm(f)>100:
                        amp = 0.5 #male: 0.5 female:0.8 (0.5 works for male and female)
                    else:
                        amp = 0.8#male:0.8 female:
            elif type == 'IMCD':
                if gender == 'female':
                    if np.linalg.norm(f)>100:
                        if k==0:
                            amp = 0.2 # male:0.1 female:0.2
                        else:
                            amp = 0.1 # male:0.2 female:0.2
                    else:
                        if k==0:
                            amp = 0.5# male:0.5 female:0.5
                        else:
                            amp = 0.5 # male:0.5 female:0.5
                elif gender == 'male':
                    if np.linalg.norm(f)>100:
                        if k==0:
                            amp = 0.29#0.19 # male:0.1 female:0.2
                        else:
                            amp = 0.17 # male:0.2 female:0.2
                    else:
                        if k==0:
                            amp = 0.8# male:0.5 female:0.5
                        else:
                            amp = 0.8 # male:0.5 female:0.5
            else:
                amp = 1.0
        elif spec == 'hum' and diabete == 'Severe' and inhib !='SGLT2':
            if type == 'S3':
                amp = 1.0
            elif type == 'SDL':
                amp = 1.0
            elif type == 'mTAL':
                if np.linalg.norm(f)>100:
                    amp = 0.2
                else:
                    amp = 1.0
            elif type == 'cTAL' or type == 'MD':
                if np.linalg.norm(f)>100:
                    amp = 0.2
                else:
                    amp = 0.8
            elif type=='DCT':
                if gender == 'female' and sup_jux == 'jux1':
                    amp = 0.9
                elif gender == 'female' and sup_jux == 'jux3':
                    amp = 0.7
                elif gender == 'female' and sup_jux == 'jux2':
                    amp = 0.7
                elif gender == 'female' and sup_jux == 'sup':
                    if np.linalg.norm(f)>2000:
                        amp = 0.8#0.5
                    else:
                        amp = 0.5
                else:
                    amp = 0.5
            elif type == 'CNT':
                if gender == 'male':
                    if np.linalg.norm(f)>100:
                        amp = 0.5
                    else:
                        amp=0.8
                elif gender == 'female':
                    if sup_jux == 'sup':
                        if np.linalg.norm(f)>100:
                            amp = 0.5
                        else:
                            amp = 0.8
                    elif sup_jux == 'jux1':
                        if np.linalg.norm(f)>1000:
                            if k==0:
                                amp = 0.5
                            else:
                                amp = 0.13
                        else:
                            amp = 0.8
                    elif sup_jux == 'jux2':
                        if np.linalg.norm(f)>1000:
                            if k==0:
                                amp = 0.5
                            else:
                                amp = 0.5
                        else:
                            amp = 0.8
                    elif sup_jux == 'jux3':
                        if np.linalg.norm(f)>1000:
                            if k==0:
                                amp = 0.5
                            else:
                                amp = 0.5
                        else:
                            amp = 0.8
                    elif sup_jux == 'jux4':
                        if np.linalg.norm(f)>1000:
                            if k==0:
                                amp = 0.3
                            else:
                                amp = 0.5
                        else:
                            amp = 0.8
                    elif sup_jux == 'jux5':
                        if np.linalg.norm(f)>1000:
                            if k==0:
                                amp = 0.1
                            else:
                                amp = 0.17
                        else:
                            amp = 1.0
            elif type == 'CCD':
                if gender == 'male':
                    if np.linalg.norm(f)>100:
                        if k == 0:
                            amp = 0.5
                        else:
                            amp = 0.2
                    else:
                        amp = 0.8
                elif gender == 'female':
                    if np.linalg.norm(f)>100:
                        if k == 0:
                            amp = 0.2
                        else:
                            amp = 0.2
                    else:
                        amp = 0.8
            elif type == 'OMCD':
                if gender == 'male':
                    if np.linalg.norm(f)>100:
                        amp = 0.5 #male: 0.5 female:0.8 (0.5 works for male and female)
                    else:
                        amp = 0.8#male:0.8 female:
                elif gender == 'female':
                    if np.linalg.norm(f)>100:
                        amp = 0.5 #male: 0.5 female:0.8 (0.5 works for male and female)
                    else:
                        amp = 0.8#male:0.8 female:
            elif type == 'IMCD':
                if gender == 'female':
                    if np.linalg.norm(f)>100:
                        if k==0:
                            amp = 0.2 # male:0.1 female:0.2
                        else:
                            amp = 0.1 # male:0.2 female:0.2
                    else:
                        if k==0:
                            amp = 0.5# male:0.5 female:0.5
                        else:
                            amp = 0.5 # male:0.5 female:0.5
                elif gender == 'male':
                    if np.linalg.norm(f)>100:
                        if k==0:
                            amp = 0.13 # male:0.1 female:0.2
                        else:
                            amp = 0.17 # male:0.2 female:0.2
                    else:
                        if k==0:
                            amp = 0.8# male:0.5 female:0.5
                        else:
                            amp = 0.8 # male:0.5 female:0.5
            else:
                amp = 1.0
        elif spec == 'hum' and diabete == 'Severe' and inhib == 'SGLT2':
            if type == 'S3':
                amp = 1.0
            elif type == 'SDL':
                amp = 1.0
            elif type == 'mTAL':
                if np.linalg.norm(f)>100:
                    amp = 0.2
                else:
                    amp = 1.0
            elif type == 'cTAL' or type == 'MD':
                if np.linalg.norm(f)>100:
                    amp = 0.2
                else:
                    amp = 0.8
            elif type=='DCT':
                if gender == 'female' and sup_jux == 'jux1':
                    amp = 0.9
                elif gender == 'female' and sup_jux == 'jux3':
                    amp = 0.7
                elif gender == 'female' and sup_jux == 'jux2':
                    amp = 0.7
                elif gender == 'female' and sup_jux == 'sup':
                    if np.linalg.norm(f)>2000:
                        amp = 0.8#0.5
                    else:
                        amp = 0.5
                else:
                    amp = 0.5
            elif type == 'CNT':
                if gender == 'male':
                    if np.linalg.norm(f)>100:
                        amp = 0.5
                    else:
                        amp=0.8
                elif gender == 'female':
                    if sup_jux == 'sup':
                        if np.linalg.norm(f)>100:
                            amp = 0.5
                        else:
                            amp = 0.8
                    elif sup_jux == 'jux1':
                        if np.linalg.norm(f)>1000:
                            if k==0:
                                amp = 0.5
                            else:
                                amp = 0.13
                        else:
                            amp = 0.8
                    elif sup_jux == 'jux2':
                        if np.linalg.norm(f)>1000:
                            if k==0:
                                amp = 0.5
                            else:
                                amp = 0.5
                        else:
                            amp = 0.8
                    elif sup_jux == 'jux3':
                        if np.linalg.norm(f)>1000:
                            if k==0:
                                amp = 0.5
                            else:
                                amp = 0.5
                        else:
                            amp = 0.8
                    elif sup_jux == 'jux4':
                        if np.linalg.norm(f)>1000:
                            if k==0:
                                amp = 0.3
                            else:
                                amp = 0.5
                        else:
                            amp = 0.8
                    elif sup_jux == 'jux5':
                        if np.linalg.norm(f)>1000:
                            if k==0:
                                amp = 0.1
                            else:
                                amp = 0.17
                        else:
                            amp = 1.0
            elif type == 'CCD':
                if gender == 'male':
                    if np.linalg.norm(f)>100:
                        if k == 0:
                            amp = 0.5
                        else:
                            amp = 0.2
                    else:
                        amp = 0.8
                elif gender == 'female':
                    if np.linalg.norm(f)>100:
                        if k == 0:
                            amp = 0.2
                        else:
                            amp = 0.2
                    else:
                        amp = 0.8
            elif type == 'OMCD':
                if gender == 'male':
                    if np.linalg.norm(f)>100:
                        amp = 0.5 #male: 0.5 female:0.8 (0.5 works for male and female)
                    else:
                        amp = 0.8#male:0.8 female:
                elif gender == 'female':
                    if np.linalg.norm(f)>100:
                        amp = 0.5 #male: 0.5 female:0.8 (0.5 works for male and female)
                    else:
                        amp = 0.8#male:0.8 female:
            elif type == 'IMCD':
                if gender == 'female':
                    if np.linalg.norm(f)>100:
                        if k==0:
                            amp = 0.2 # male:0.1 female:0.2
                        else:
                            amp = 0.1 # male:0.2 female:0.2
                    else:
                        if k==0:
                            amp = 0.5# male:0.5 female:0.5
                        else:
                            amp = 0.5 # male:0.5 female:0.5
                elif gender == 'male':
                    if np.linalg.norm(f)>100:
                        if k==0:
                            amp = 0.23 # male:0.1 female:0.2
                        else:
                            amp = 0.17 # male:0.2 female:0.2
                    else:
                        if k==0:
                            amp = 0.8# male:0.5 female:0.5
                        else:
                            amp = 0.8 # male:0.5 female:0.5
            else:
                amp = 1.0
        delta = amp*np.array(F * IJ.T)[0]
        x -= delta
        f = np.matrix(fun(x,k))
        iter+=1
        print(iter,np.linalg.norm(f))
        TOLpcn = np.max(delta / x)
        #print(i)

         # Pause: Added by Dania
#        input("Pausing! Press Enter to continue...")
#    print("Iteration Times: " + str(i) + " with TOL " + str(TOLpcn) + "%")

    return x

#@jit(nopython=True) # Set "nopython" mode for best performance, equivalent to @njit
def broyden(func,x,k,type):
    fun=equations.conservation_eqs
    f=np.matrix(fun(x,k))
    J=np.matrix(Jac(fun,x,k))
#    IJ=np.linalg.inv(J)
    IJ=J.I
    dx=np.ones(x.shape)
    i=0
    iter=0
    #while(np.max(dx)>0.0001):
    while(np.linalg.norm(f)>0.0001) and (iter<500):

        f_previous=f
        x_previous=x

        x=x-np.array(f*IJ.T)[0]

        f=np.matrix(fun(x,k))
#                print('x',x)
#                print('f',f)

        df=f-f_previous
        dx=x-x_previous
        # #
        # #-------------------------------------------------------
        # #using good broyden
        # dx=np.array([dx])
        # df=np.array([df])
        # dx=dx.T
        # df=df.T

        # IJ = IJ+(dx-IJ*df)*(dx.T*IJ)/np.inner(dx.T*IJ,df)
        # #-------------------------------------------------------

        IJ=IJ-np.outer(IJ*f.T,dx)*IJ/np.inner(dx,dx+(IJ*f.T).T)
        #print(i)
        iter+=1
        print(iter,np.linalg.norm(f))#SS-deletechange "Newton2" also Newton above somewhere

        #Pause: Added by Dania
#        input("Pausing! Press Enter to continue...")
        #J=J+np.outer((df-dx*J.T),dx)/np.linalg.norm(x)**2
        #IJ=J.I

    return x
