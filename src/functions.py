import csv
import os
import constants as const
import numpy as np
import math as m
from scipy.optimize import leastsq
import matplotlib.pyplot as plt
from tqdm import tqdm

#takes path of csv file, gaugelength of tensile test (in m), and cross sectional area of films (in mm^2)
#returns 2 arrays _ stress, strain
def load_csv(path, guageLength, a_crossSectional):
    strain = []
    stress = []
    with open(path, newline= '') as csvfile:
        spamreader = csv.reader(csvfile, delimiter=',', quotechar='|')
        initialPos = 0
        i = 0
        for row in spamreader:
            if i == 1:
                initialPos = float(row[2])
            if i != 0:
                strain.append((float(row[2]) - float(initialPos)) / guageLength)
                stress.append(float(row[3]) * const.GRAVITY / a_crossSectional)
            i = i + 1
        csvfile.close()
    return np.array(stress), np.array(strain)

#get finds the index of the maximum value of the array
#retunrs the maximum value, and the index of that value
def get_max_index(array):
    maxVal = array[0]
    index = 0
    for i in range(len(array)):
        if maxVal < array[i]:
            maxVal = array[i]
            index = i
    return maxVal, index

#partitions data for EAS tensile Averager
#assumes strain and stress are same length
#@strain input strain array
#@stress input stress array
#@return strainLPartition, strainRPartition, stressLPartition, stressRPartition, splitIndex, buffRange
def partition(strain, stress):
    _ , splitIndex = get_max_index(stress)
    strainL = []
    strainR = []
    stressL = []
    stressR = []
    buffRange = int(splitIndex / 5) #collect extra data from each side for better fit
    for i in range(len(strain)):
        if i > splitIndex - buffRange:
            strainR.append(strain[i])
            stressR.append(stress[i])
        if i <= splitIndex + buffRange:
            strainL.append(strain[i])
            stressL.append(stress[i])
    return np.array(strainL), np.array(strainR), np.array(stressL), np.array(stressR), splitIndex, buffRange

#least squares regression model for line fits (4th order polynomial)
#@coeffs coefficients for model ax^4 + bx^3 + cx^2 + dx + e
#@eps (epsilon) strain input
def model(coeffs, eps):
    a = coeffs[0]
    b = coeffs[1]
    c = coeffs[2]
    d = coeffs[3]
    e = coeffs[4]

    return (a * (eps)**2) + (b * ((eps)**3)) + (c * ((eps)**2)) + (d * eps) + e

#residuals used for least squares
#@coeffs coefficients passed to function
#@function function model for fit
#@eps strain values for function
#@exp experimental values for fit
def residuals(coeffs, function, eps, exp):
    return exp - abs(function(coeffs, eps))

#takes a stress and strain array input and outputs the coresponding
#stress and strain points calculated based on equally partitioned arcs
#of pre and post maximum partitions of stress strain curve
#@stress values loaded from csv assumed to be MPa
#@strain values loaded from csv assumed to be % and equal to @stress
#@return stress_arc, strain_arc
def get_arc_points(stress, strain):
    #partitions the data into pre and post max stress
    strain_pre, strain_post, stress_pre, stress_post, splitIndex, buffRange = partition(strain, stress)

    # fit pre and post partitions with 4th order polynomial using least squares regression
    # and generate fitted data brocken into N segments where Npre and Npost are >> than 
    # the number of experimental data points in those repective partitions
    coeff_guess = [0, 0, 0, 0, 0]
    coeff_pre, _ = leastsq(residuals, coeff_guess, args=(model, strain_pre, stress_pre))
    fitStrain_pre = np.arange(strain_pre[0], strain_pre[splitIndex], (strain_pre[splitIndex] - strain_pre[0])/ const.N_PRE)
    fit_pre = model(coeff_pre, fitStrain_pre)

    coeff_post, _ = leastsq(residuals, coeff_guess, args=(model, strain_post, stress_post))
    fitStrain_post = np.arange(strain_post[buffRange], strain_post[len(strain_post) - 1], (strain_post[len(strain_post) - 1] - strain_post[buffRange]) / const.N_POST)
    fit_post = model(coeff_post, fitStrain_post)

    # calculate arc length of each point int fit
    d_pre = []
    for i in range(1, len(fit_pre)):
        d_pre.append(m.sqrt(((fitStrain_pre[i] - fitStrain_pre[i - 1]) ** 2) + m.sqrt(((fit_pre[i] - fit_pre[i - 1]) ** 2))))

    d_post = []
    for i in range(1, len(fit_post)):
        d_post.append(m.sqrt(((fitStrain_post[i] - fitStrain_post[i - 1]) ** 2) + m.sqrt(((fit_post[i] - fit_post[i - 1]) ** 2))))

    #calculate total arc length from segmenteded arc lengths       
    Spre = sum(d_pre)
    Spost = sum(d_post)

    #segment partitions into equal arc lengths and find corresponding x and y coordinates
    #and add to xarc and y arc for both partitions
    spre_seg = Spre / (const.N_PRE / const.n_pre_factor)
    spost_seg = Spost / (const.N_POST / const.n_post_factor)
    xarc = []
    yarc = []
    sgs_pre = []
    sgs_post = []

    for i in range(1, int(const.N_PRE/const.n_pre_factor) + 1):
        seg = 0
        j = 0
        while seg < spre_seg * i and j < len(d_pre):
            seg = seg + d_pre[j]
            j = j + 1
        if i == 1:
            sgs_pre.append(seg)
        else:
            sgs_pre.append(seg - sum(sgs_pre))
        xarc.append(fitStrain_pre[j])
        yarc.append(fit_pre[j])
    
    for i in range(1, int(const.N_POST/const.n_post_factor) + 1):
        seg = 0
        j = 0
        while seg < spost_seg * i and j < len(d_post):
            seg = seg + d_post[j]
            j = j + 1
        if i == 1:
            sgs_post.append(seg)
        else:
            sgs_post.append(seg - sum(sgs_post))
        xarc.append(fitStrain_post[j])
        yarc.append(fit_post[j])
    
    #check to make sure that summed arc lengths equal total arc length
    if not np.isclose(Spre, sum(sgs_pre), 0.01):
        diff = Spre - sum(sgs_pre)
        print("Spre and sum of s_pre_i off by ", diff, " adjust N constants")
    elif not np.isclose(Spost, sum(sgs_post), 0.01):
        diff = Spost - sum(sgs_post)
        print("Spre and sum of s_pre_i off by ", diff, " adjust N constants")
    
    return yarc, xarc

#takes an array of plot data as an input and and outputs the averaged data
#@stress_plots array of stress values to be averaged
#@strain_plots array of streain values to be aveeraged (length assumed to be equal to stress_plots)
#@return averaged values stress_ave, strain_ave
def average_plots(stress_plots, strain_plots):
    strains = []
    stresses = []

    for i in tqdm(range(len(stress_plots)), desc="Processing..."):
        stress, strain = get_arc_points(stress_plots[i], strain_plots[i])
        stresses.append(stress)
        strains.append(strain)
        pass

    strain_ave = []
    stress_ave = []
    for i in range(len(strains[0])):
        strainVal = 0
        stressVal = 0
        for ars in  strains:
            strainVal = strainVal + ars[i]
        for ars in stresses:
            stressVal = stressVal + ars[i]
        strain_ave.append(strainVal / len(strains))
        stress_ave.append(stressVal / len(stresses))
    
    return stress_ave, strain_ave


# takes input files and outputs averaged data into csv file in output folder
# and will also output plots if plot is true
# @files csv files with stress strain data
# @thicknesses txt file with width (m) and cross sectional area (mm^2)
# @plot boolean value, diplays plot if true
def outputAverageFile(files, thicknesses):
    stresses = []
    strains = []
    outputFileName = files[0].split("/")
    print("sample " + outputFileName[1], " starting")
    outputFileName = "output/" + outputFileName[1]

    # for i in tqdm (range (25), desc="Loading..."):
    # pass
    for i in range(len(files)):
        stress, strain = load_csv(files[i], const.GAUGE_LENGTH, thicknesses[i] * const.WIDTH)
        stresses.append(stress)
        strains.append(strain)

    stress_ave, strain_ave = average_plots(stresses, strains)
    trimIndex = trimOscilationZero(stress_ave)
    if not trimIndex == -1:
        stress_ave = stress_ave[0:trimIndex]
        strain_ave = strain_ave[0:trimIndex]

    #write data to save file
    with open((outputFileName + ".csv"), "w") as f:
        for i in range(len(stress_ave)):
            string = "" + str(strain_ave[i]) + "," + str(stress_ave[i]) +"\n"
            f.write(string)
        f.close

    #save plots
    if const.PLOTS:
        for i in range(len(stresses)):
            plt.scatter(strains[i], stresses[i], s=const.MARKER_SIZE, color="blue", label="exp")

        plt.plot(strain_ave, stress_ave, color="red", label="ave")
        plt.legend()
        plt.xlabel("Strain (%")
        plt.ylabel("Stress (MPa)")
        plt.savefig(outputFileName + ".png")
        plt.close()
    print("sample complete and written to: ", outputFileName)
    print("..............")
    

def trimOscilationZero(stress):
    for i in range(1, len(stress)):
        if stress[i - 1] > stress[i] and stress[i] <= 0:
            return i
    return -1