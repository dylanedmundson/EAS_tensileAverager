A Program to average multiple stress strain curves and produce the averaged data as well
as plots of the experimental data with the average. Averaging is done using the equal
arc segmentation method described by Zhong and Wille:

    R. Zhong and K. Wille, “Equal Arc Segment Method for Averaging Data Plots Exemplified for 
    Averaging Stress versus Strain Curves of Pervious Concrete,” Journal of Materials in Civil 
    Engineering, vol. 28, no. 1, p. 04015071, Jan. 2016, doi: 10.1061/(ASCE)MT.1943-5533.0001345.




# INSTRUCTIONS:

input the sample labels into the sampleLabels.txt file in input
input the thickness (mm) associated with each sample lable into thicknesses.txt in the same exact order
      (row 1 sample lable corresponds to row 1 thickness)
      (easiest to organize in excel then copy an paste columns)
      (make sure there are no left over lines)
add folders with csv files of stress and strain data to input
     --> sample root folder 1
         -replicate1.csv
         -replicate2.csv
         -etc....
     --> sample root folder 2
     --> etc...........
NOTE: that the sample root folder name will become this output name of the averaged data
NOTE: input file should have position(col 3 or C) in m and load (col 4 or D) in kg
go into:
  -->EAS_tensileAverager
     -src
          -constants.py
              change constant GAUGE_LENGTH to your samples gauge length
              change constant WIDTH to your samples width
              change constant PLOTS to True if you want plots and False if not

data is saved into the output folder.
any subsequent runs of this program will overwirte previous data