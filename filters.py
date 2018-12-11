from scipy.signal import butter, filtfilt, iirnotch, freqz, lfilter, cheby2


def _butter_lowpass(cutoff, fs, order=5):
    nyq = 0.5 * fs
    normal_cutoff = cutoff / nyq
    b, a = butter(order, normal_cutoff, btype='low', analog=False)
    return b, a


def butter_lowpass_filter(data, cutoff, fs, order=5):
    b, a = _butter_lowpass(cutoff, fs, order=order)
    y = filtfilt(b, a, data)
    return y


def _butter_highpass(cutoff, fs, order=5):
    nyq = 0.5 * fs
    normal_cutoff = cutoff / nyq
    b, a = butter(order, normal_cutoff, btype='high', analog=False)
    return b, a


def butter_highpass_filter(data, cutoff, fs, order=5):
    b, a = _butter_highpass(cutoff, fs, order=order)
    y = filtfilt(b, a, data)
    return y


def notch_filter(data, f0, fs, Q=30):
    nyq = 0.5 * fs
    w0 = f0/nyq
    b, a = iirnotch(w0, Q)
    y = filtfilt(b, a, data)
    return y
  
  
def _cheb2_notch(cutoff, fs, order=5, rs=3, width=.1):
    nq = fs/2
    Wn_min, Wn_max = (cutoff - width) / nq, (cutoff + width) / nq
    Wn = [Wn_min, Wn_max]
    b, a = cheby2(N=order, rs=rs, Wn=Wn, btype='bandstop', analog=False, output='ba')
    return b, a
