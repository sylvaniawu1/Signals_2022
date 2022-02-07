import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.signal import fftconvolve
from scipy.signal import find_peaks

seq_csv = pd.read_csv("data/Staph_termseq_linz2_plus_strand_allCoverage.gr",sep="\t",names=["genome","time","signal"])

seq_depths = seq_csv['signal'].to_numpy()
win = 5000
win_avg = np.ones( (1,win) )/win
win_avg = win_avg[0, :]

print("Check average filter:\n")
print(win_avg[0:20])

seq_pd = pd.DataFrame(seq_depths,columns=['term-seq depth'])
#seq_pd['Moving average convolved'] = convl.tolist()
print(seq_pd)



# Fast Fourier Transform convolves moving average with term-seq depths
my_avg_convl = fftconvolve(in1=seq_depths,in2=win_avg,mode="same")


# Need for "computational formula for variance/stdev": Moving avg convolution of square of the data
my_avg_sqrd = np.square(my_avg_convl)
my_seq_sqrd = np.square(seq_depths)
my_seq_sqr_convl = fftconvolve(in1=my_seq_sqrd,in2=win_avg,mode="same")

print(f'Size of array: my_avg_convl {my_avg_convl.size}')
print(f'Size of array: my_avg_sqrd {my_avg_sqrd.size}')
print(f'Size of array: my_seq_sqrd {my_seq_sqrd.size}')
print(f'Size of array: my_seq_sqr_convl {my_seq_sqr_convl.size}')

my_stdev = np.sqrt( abs( my_seq_sqr_convl - my_avg_sqrd) )
z_score = abs( (my_avg_convl - seq_depths)/(0.0000000000000001+my_stdev) )

from scipy.signal import argrelextrema

print("Z-score table size:\n")
print(z_score.size)
# for local maxima
peaks_idx = argrelextrema(z_score, np.greater, order=500)
peaks_vals = z_score[peaks_idx[0]]

pd_peaks = pd.DataFrame(peaks_idx[0],columns=["peaks_idx"])
pd_peaks["peaks_vals"]=peaks_vals
pd_peaks[pd_peaks["peaks_vals"]>1.5].to_csv('/Users/sylvaniawu/PycharmProjects/termseqFFT/results/FFT_peaks_minZscore1.5_linezolid2_plus_strand_2022_02_07.tsv',sep='\t')

print("Number of peaks z-score > 1.5:\n")
print(peaks_vals[peaks_vals>1.5].size)
print("Number of peaks z-score < 1.5:\n")
print(peaks_vals[peaks_vals<1.5].size)
# plt.plot(my_stdev)
pd_stdev = pd.DataFrame(my_stdev,columns=["SD from scipy.fftconvolve"])
pd_stdev["time"]=seq_csv["time"]
pd_stdev["Z-score"]=pd.DataFrame(z_score,columns=["Z-score"])
#pd_stdev["Maxima_of_Z-score"]=pd.DataFrame(z_score[peaks],columns=["Maxima_of_Z-score"])


sd_fig = plt.figure(figsize=(20,15),dpi=100)
plt.ylim(0,12)

ax = plt.gca()
ax.scatter(x=pd_stdev["time"],y=pd_stdev["Z-score"],color="b")
ax.scatter(x=peaks_idx,y=peaks_vals,color="r")


plt.show()


#x = [1,2,3,4,5,6,7,8,9,10]
#y#=[11,12,13,14,15,16,17,18,19,20]
#plt.plot(x,y)

# See PyCharm help at https://www.jetbrains.com/help/pycharm/
#plt.show()
