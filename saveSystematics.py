import numpy as np

systemDict = {
    #"Normal": "Normal",
    #"rvtx_low": "rvtx_low",
    #"rvtx_high": "rvtx_high",
    #"zvtx_low": "zvtx_low",
    #"zvtx_high": "zvtx_high",
    "dca_low": "dca_low",
    "dca_high": "dca_high",
    "nhits_low": "nhits_low",
    "nhits_high": "nhits_high",
    #"nhitsdedx_low": "nhitsdedx_low",
    "nhitsdedx_high": "nhitsdedx_high",
    "nhitsratio_low": "nhitsratio_low",
    "nhitsratio_high": "nhitsratio_high",
    #"nSigPi_low": "nSigPi_low",
    #"nSigPi_high": "nSigPi_high",
    #"nSigKa_low": "nSigKa_low",
    #"nSigKa_high": "nSigKa_high",
    #"nSigPr_low": "nSigPr_low",
    #"nSigPr_high": "nSigPr_high",
    #"m2Pi_low": "m2Pi_low",
    #"m2Pi_high": "m2Pi_high",
    #"m2Ka_low": "m2Ka_low",
    #"m2Ka_high": "m2Ka_high",
    #"epd_low": "epd_low",
    #"epd_high": "epd_high"
}


np.save('dict_systematics.npy', systemDict)

print("Dictionary of current job IDs for systematics saved to dict_systematics.npy.")
