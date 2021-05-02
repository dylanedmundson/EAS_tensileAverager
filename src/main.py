import functions as fun
import matplotlib.pyplot as plt
import constants as c
from scipy.optimize import leastsq
import numpy as np
import os

print("\nEAS averaging initialization")
print("............................")

directory = "input"
inputs = np.loadtxt("input/thickness.txt")
file1 = open("input/sampleLabels.txt")
labels = file1.read()
labels = labels.split('\n')
for subDir in os.listdir(directory):
    if not subDir.endswith(".txt"):
        files = []
        thicknesses = []
        for fileName in os.listdir(os.path.join(directory, subDir)):
            if fileName.endswith(".csv"):
                files.append(directory + "/" + subDir + "/" + fileName)
                for i in range(len(labels)):
                    if fileName.startswith(labels[i]):
                        thicknesses.append(inputs[i])
        fun.outputAverageFile(files, thicknesses)
print("..............\nall files complete\n")