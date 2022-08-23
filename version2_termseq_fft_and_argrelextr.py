import numpy as np
import pandas as pd
from scipy.signal import fftconvolve
from optparse import OptionParser
import matplotlib.pyplot as plt

parser = OptionParser()
parser.add_option("-q", "--earlyQuit", action="store", type="int", dest="quit_point")
# 1 2 3 4 5
parser.add_option("-g", "--genomeName", action="store", type="string", dest="genome_name")
parser.add_option("-t", "--termSeq", action="store", type="string", dest="termseq_covgr")
parser.add_option("-z", "--zScore", action="store", type="float", dest="z_threshold")
parser.add_option("-n", "--label", action="store", type="string", dest="output_label")
parser.add_option("-s", "--strand", action="store", type="string", dest="termseq_strand")
parser.add_option("-p", "--showGraph", action="store", type="string", dest="show_graph")
parser.add_option("-c", "--cutZeroCov", action="store", type="string", dest="filter_zero_cov")
parser.add_option("-f", "--movingAvgFilter", action="store", type="string", dest="avg_filter_type")

(options,args) = parser.parse_args()


import find_plateaus

find_plateaus.check_findp()


if options.quit_point == 0:
    quit()

if options.termseq_strand == "plus":
    strand = "+"
elif options.termseq_strand == "minus":
    strand = "-"
else:
    strand = "."
# RUNFFT="/Users/sylvaniawu/PycharmProjects/termseqFFT/termseq_fft_and_argrelextr.py"
# python ${RUNFFT} --genomeName="CP002120.1" --zScore=0.0 --strand="+" --termSeq=${COVGR} --label=${LABEL} --showGraph=="peaks"

#seq_csv = pd.read_csv("/Users/sylvaniawu/Desktop/Chapter5_May_June_July_2022/data/termseq_covgr/depth_gr/termseq_3end_single_base_vancomycin_1_plus_strand_flipped_sense_novo_aligned_25_06_2022.gr",sep="\t",names=["genome","time","signal"])
seq_csv = pd.read_csv(options.termseq_covgr,sep="\t",names=["genome","time","signal"])


seq_csv_nonzero = seq_csv[seq_csv['signal'] != 0].reset_index(drop=True)
gpos_nonzero = seq_csv[seq_csv['signal'] != 0]['time'].reset_index(drop=True)
gpos_nonzero.columns=["time"]

print("nonzero coverage:\n")
print(seq_csv_nonzero)


# LOOK HERE
# All convl's, i.e. standard deviation and z-score calcs, only done on reads > 1 (seq_csv_nonzero)
seq_depths = seq_csv_nonzero['signal'].to_numpy()
print("seq_depths:\n")
print(seq_depths)
print("gpos_nonzero:\n")
print(gpos_nonzero)

if options.quit_point == 1:
    quit()

win = 5000

if options.avg_filter_type == "simple":
    win_avg = np.ones( (1,win) )/win
    win_avg = win_avg[0, :]

if options.avg_filter_type == "gapped":
    from math import ceil
    length_data = len( gpos_nonzero )
    print(f" length_data = {length_data}\n")
    half_length_data = ceil(length_data/2)
    gap = 3
    half_width = ceil(win/2)
    pitted_avg_filter = np.zeros(length_data)
    pitted_avg_filter[half_length_data - half_width - gap : half_length_data - gap ] = np.ones( (1,half_width) )/win
    pitted_avg_filter[half_length_data + gap : half_length_data + half_width + gap  ] = np.ones( (1,half_width) )/win

#print("Check average filter:\n")
#print(win_avg[0:20])

seq_pd = pd.DataFrame(seq_depths,columns=['term-seq depth'])
#seq_pd['Moving average convolved'] = convl.tolist()

print(seq_pd)



# Fast Fourier Transform convolves moving average with term-seq depths
if options.avg_filter_type == "simple":
    my_avg_convl = fftconvolve(in1=seq_depths,in2=win_avg,mode="same")

if options.avg_filter_type == "gapped":
    my_avg_convl = fftconvolve(in1=seq_depths,in2=pitted_avg_filter,mode="same")

# Need for "computational formula for variance/stdev": Moving avg convolution of square of the data
my_avg_sqrd = np.square(my_avg_convl)
my_seq_sqrd = np.square(seq_depths)

if options.avg_filter_type == "simple":
    my_seq_sqr_convl = fftconvolve(in1=my_seq_sqrd,in2=win_avg,mode="same")

if options.avg_filter_type == "gapped":
    my_seq_sqr_convl = fftconvolve(in1=my_seq_sqrd,in2=pitted_avg_filter,mode="same")

print(f'Size of seq_csv_nonzero["signal"]:{len(seq_csv_nonzero["signal"])}')
print(f'Size of array: my_avg_convl {my_avg_convl.size}')
print(f'Size of array: my_avg_sqrd {my_avg_sqrd.size}')
print(f'Size of array: my_seq_sqrd {my_seq_sqrd.size}')
print(f'Size of array: my_seq_sqr_convl {my_seq_sqr_convl.size}')

my_stdev = np.sqrt( abs( my_seq_sqr_convl - my_avg_sqrd) )
z_score = abs( (my_avg_convl - seq_depths)/(0.0000000000000001+my_stdev) )

from scipy.signal import argrelextrema

print("Z-score table size:\n")
print(z_score.size)
print("Z-scores > 0.5 table size (before argrelextrema):\n")
print(z_score[z_score>0.5].size)
print("Z-scores > 1.0 table size (before argrelextrema):\n")
print(z_score[z_score>1.0].size)
print("Z-scores > 60 table size (before argrelextrema):\n")
print(z_score[z_score>60].size)

print("Z-scores > 60  (before argrelextrema):\n")
print(np.argwhere(z_score>60))

print("Test result of argwhere:\n")
print(f" >60 at : {z_score[42840]}, {z_score[42841]}, {z_score[51339]},{z_score[71271]},{z_score[81392]},{z_score[81393]}\n")

gpos_z_gt_60 = gpos_nonzero[ np.argwhere(z_score>60).flatten() ]
print("Genome positions at Z-scores > 60  (before argrelextrema):\n")

print(gpos_z_gt_60)
print("From seq_csv_nonzero ==> Raw depth at Z-scores > 60  (before argrelextrema):\n")
print( seq_csv_nonzero.merge(right=gpos_z_gt_60,on="time",how="inner") )
print("From seq_csv  ==> Raw depth at Z-scores > 60  (before argrelextrema):\n")
print( seq_csv.merge(right=gpos_z_gt_60,on="time",how="inner") )
#print("From seq_csv ==> Raw depth at Z-scores > 60  (before argrelextrema):\n")
#print(seq_csv[seq_csv["time"] == gpos_z_gt_60["time"][0] ])

if options.quit_point == 2:
    quit()

filterby="zscore"

# for local maxima
if filterby == "localextrema":
    peaks_idx = argrelextrema(z_score, np.greater, order=500)
    peaks_vals = z_score[peaks_idx[0]]
    pd_peaks = pd.DataFrame(peaks_idx[0],columns=["peaks_idx"])

if filterby == "zscore":
    peaks_idx = np.argwhere(z_score > options.z_threshold)
    gpos_z_gt = gpos_nonzero[ peaks_idx.flatten() ]
    print(f"filter by zscore.  this is argwhere(zscore > {str(options.z_threshold)}):\n")
    print(peaks_idx)
    print(f"This is genome positions where zscore > {str(options.z_threshold)}:\n")
    print(gpos_z_gt)
    peaks_vals = z_score[peaks_idx]
    print("These are the passing Z=scores:")
    print(peaks_vals[0:20])
    print("From seq_csv_nonzero ==>\n")
    df_nonzero_Zgt = seq_csv_nonzero.merge(right=gpos_z_gt,on="time",how="inner")
    print(df_nonzero_Zgt)
    print("From seq_csv ==>\n")
    df_Zgt = seq_csv.merge(right=gpos_z_gt,on="time",how="inner")
    df_Zgt["peaks_vals"] = pd.Series(z_score[peaks_idx].flatten())
    print(df_Zgt)
    pd_peaks = pd.DataFrame(peaks_idx,columns=["peaks_idx"])

if options.quit_point == 3:
    quit()


#pd_peaks[pd_peaks["peaks_vals"]>1.5].to_csv('/Users/sylvaniawu/PycharmProjects/termseqFFT/results/FFT_peaks_minZscore1.5_vancomycin_2_plus_strand_flipped_sense_2022_06_28.tsv',sep='\t')
#pd_peaks.to_csv('/Users/sylvaniawu/PycharmProjects/termseqFFT/results/FFT_peaks_noMinScore_vancomycin_2_plus_strand_flipped_sense_2022_06_28.tsv',sep='\t')


pd_bed = pd.DataFrame(columns=["chrom","chromStart","chromEnd","name","score","strand","signalValue","p-value","q-value","summit"])

pd_bed["chromStart"]=df_Zgt["time"]-5
pd_bed["chromEnd"]=df_Zgt["time"]+5
pd_bed["name"]="."
pd_bed["score"]=0
pd_bed["strand"]=str(strand)
pd_bed["signalValue"]=df_Zgt["peaks_vals"]
pd_bed["p-value"]=-1
pd_bed["q-value"]=4.5777559727315
pd_bed["summit"]=5
pd_bed["chrom"]=str(options.genome_name)
print(pd_bed)

top_ntile_z = pd_bed[ pd_bed["signalValue"] >= pd_bed.groupby('strand')["signalValue"].transform('quantile',0.999) ].sort_values(by="signalValue",ascending=False)

print("Here's the top 0.1 percent of peaks by highest Z-score:\n")
print(top_ntile_z)
if options.quit_point == 4:
    quit()
pd_bed.to_csv(f'/Users/sylvaniawu/Desktop/fft_gapped_central_filter/fft_peaks_2022_08_23/FFT_peaks_window{str(win)}_minZscore{str(options.z_threshold)}_{str(options.avg_filter_type)}MovingAvg_gap{str(gap)}bp_{str(options.output_label)}.bed',index=False,sep='\t') # label was vancomycin_1_plus_strand_flipped_sense_2022_06_28
top_ntile_z.to_csv(f'/Users/sylvaniawu/Desktop/fft_gapped_central_filter/fft_peaks_2022_08_23/Top.999_percentile_FFT_peaks_window{str(win)}_minZscore{str(options.z_threshold)}_{str(options.avg_filter_type)}MovingAvg_gap{str(gap)}bp_{str(options.output_label)}.bed',index=False,sep='\t') # label was vancomycin_1_plus_strand_flipped_sense_2022_06_28



# BED .oracle.NarrowPeak format for "IDR" (github nboley) or "termseq-peaks" (github nichd)
# (1) chrom string (2) chromStart int (3) chromEnd int (4) (use '.') (5) use 0 (IDR val)
# (6) strand [+-.] (7) signalValue float (8) -1 (9) q-value float
# (10) summit int (distance from chromStart to centre of the peak)

print(f"Number of peaks z-score > {str(options.z_threshold)}:\n")
print(peaks_vals[peaks_vals > options.z_threshold].size)
print("Number of peaks z-score < {str(options.z_threshold)}:\n")
print(peaks_vals[peaks_vals < options.z_threshold].size)
# plt.plot(my_stdev)
pd_stdev = pd.DataFrame(my_stdev,columns=["SD from scipy.fftconvolve"])
pd_stdev["time"]=seq_csv["time"]
pd_stdev["Z-score"]=pd.DataFrame(z_score,columns=["Z-score"])
#pd_stdev["Maxima_of_Z-score"]=pd.DataFrame(z_score[peaks],columns=["Maxima_of_Z-score"])



import math

if options.show_graph == "peaks":
    sd_fig = plt.figure(figsize=(20,15),dpi=100)
    plt.ylim(0,12)

    ax = plt.gca()
    #ax.scatter(x=pd_stdev["time"],y=pd_stdev["Z-score"],color="b")
    ax.scatter(x=seq_csv["time"],y=seq_csv["signal"]/500,color="b")
    ax.scatter(x=df_Zgt["time"],y=df_Zgt["peaks_vals"],color="r")
    start, end = ax.get_xlim()
    ax.xaxis.set_ticks(np.arange(start, end, 100000))
    plt.grid(True, which="both")
    plt.show()

if options.show_graph == "cdf":
    count, bins_count = np.histogram(z_score[z_score < 1] , bins=10000) #[z_score < 0.01]
    # finding the PDF of the histogram using count values
    pdf = count / sum(count)

    # using numpy np.cumsum to calculate the CDF
    # We can also find using the PDF values by looping and adding
    cdf = np.cumsum(pdf)

    plt.rcParams.update({'font.size': 52})

    # plotting PDF and CDF
    plt.plot(bins_count[1:], pdf, color="red", label="PDF")
    plt.plot(bins_count[1:], 1-cdf, label="1-CDF")
    plt.legend()
    plt.show()

#x = [1,2,3,4,5,6,7,8,9,10]
#y#=[11,12,13,14,15,16,17,18,19,20]
#plt.plot(x,y)

# See PyCharm help at https://www.jetbrains.com/help/pycharm/
#plt.show()
