import numpy as np
import pandas as pd
from optparse import OptionParser
import matplotlib.pyplot as plt

parser = OptionParser()
parser.add_option("-T", "--TEX", action="store", type="string", dest="term_exonucl") # plus or minus
parser.add_option("-C", "--Condition", action="store", type="string", dest="condition") # vancomycin or control
parser.add_option("-s", "--strand", action="store", type="string", dest="strand") # vancomycin or control
parser.add_option("-b", "--bed", action="store", type="string", dest="bed_peaks") #
parser.add_option("-g", "--gff", action="store", type="string", dest="gff_peaks") #
parser.add_option("-t", "--tsv", action="store", type="string", dest="tsv_peaks") #
parser.add_option("-w", "--window", action="store", type="int", dest="window_size") # nucleotides around peak
(options,args) = parser.parse_args()

dir_drnaseq_vanco = "/Users/sylvaniawu/PycharmProjects/termseqFFT/data/5end_only/Vanco/"
str_vanco_tex_minus_fwd = "Vanco_TEXminus_SYY4838A7-9_S6_R1_001.trimmed.depl.sorted_fwd.segemehl.5end_only.gr"
str_vanco_tex_minus_rev = "Vanco_TEXminus_SYY4838A7-9_S6_R1_001.trimmed.depl.sorted_rev.segemehl.5end_only.gr"
str_vanco_tex_plus_fwd = "Vanco_TEXplus_SYY4838A7-9-T_S3_R1_001.trimmed.depl.sorted_fwd.segemehl.5end_only.gr"
str_vanco_tex_plus_rev = "Vanco_TEXplus_SYY4838A7-9-T_S3_R1_001.trimmed.depl.sorted_rev.segemehl.5end_only.gr"

gr_header=["genome","position","depth"]
pd.set_option('display.max_columns', None)

if options.strand == "plus":
    if options.term_exonucl == "plus":
        gr_path = dir_drnaseq_vanco + str_vanco_tex_plus_fwd
    if options.term_exonucl == "minus":
        gr_path = dir_drnaseq_vanco + str_vanco_tex_minus_fwd
if options.strand == "minus":
    if options.term_exonucl == "plus":
        gr_path = dir_drnaseq_vanco + str_vanco_tex_plus_rev
    if options.term_exonucl == "minus":
        gr_path = dir_drnaseq_vanco + str_vanco_tex_minus_rev

cov_gr=pd.read_csv(gr_path, sep='\t', header=None,names=gr_header)

print(cov_gr.head())

def calc_total_upstream_cov(x):
    w=options.window_size
    total_upstream = cov_gr[x-w:x]["depth"].sum()
    return total_upstream

def calc_total_downstream_cov(x):
    W=options.window_size
    total_downstream = cov_gr[x:x+W]["depth"].sum()
    return total_downstream

peaks_df = pd.DataFrame()

if options.bed_peaks != None :
    print("HELLO\n")
    bed_header=["chrom","start","end","dot1","dot2","strand","dot3","dot4","dot5","distStartToPeakCentre"]
    bed_df = pd.read_csv(options.bed_peaks, sep='\t', header=None,names=bed_header)
    peaks_df["peakCentre"] = bed_df["start"] + bed_df["distStartToPeakCentre"]
    #print(bed_df.head())
    #print(peak_centres.head())

col_title_upstream = "total_upstream_cov_"+str(options.window_size)+"bp"
col_title_downstream = "total_downstream_cov_"+str(options.window_size)+"bp"
peaks_df[col_title_upstream] = peaks_df["peakCentre"].map(calc_total_upstream_cov)
peaks_df[col_title_downstream] = peaks_df["peakCentre"].map(calc_total_downstream_cov)


print(peaks_df)

fig = plt.figure(figsize =(10, 7))

data_1 = np.array(peaks_df[col_title_downstream])
data_2 = np.array(peaks_df[col_title_upstream])
data_3 = np.array(peaks_df[col_title_downstream])
data_4 = np.array(peaks_df[col_title_upstream])
data = [ data_1, data_2, data_3, data_4 ]

ax = fig.add_axes([0, 0, 1, 1])

# Creating plot
bp = ax.boxplot(data)

# show plot
plt.show()
