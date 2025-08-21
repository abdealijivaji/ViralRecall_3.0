
import matplotlib.pyplot as plt
import pandas
import scipy.signal as ssig
from scipy.signal import find_peaks

#annot_tbl = pandas.read_csv("../Chlamy_punui/Chlamy_punui.tsv", sep= "\t", header= 0)

plt.style.use('_mpl-gallery')




# y_peak, _ = find_peaks(y.array)

# #roll_welch, Pow_den  = scipy.signal.welch(y.array, 5e+4 )
# fs = 100

# sos = ssig.iirfilter(4, Wn=[0.5, 2.5], fs=fs, btype="bandpass",
#                              ftype="butter", output="sos")
# y_filt = ssig.sosfilt(sos, y.array)
# filt_peaks = find_peaks(y_filt, distance = fs, height = 10)



def plot_vreg(annot_tbl : pandas.DataFrame) -> None :
    x = annot_tbl['pstart']
    y = annot_tbl['rollscore']
    y_peak = ssig.find_peaks_cwt(y.array, widths= 20)

    fig, ax = plt.subplots()

    ax.plot(x, y)
    #ax.plot(x[y_peak], y[y_peak], "x")
    ax.plot(x[y_peak], y[y_peak], "x")

    plt.show()
    return None