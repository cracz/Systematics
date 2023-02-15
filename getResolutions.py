import subprocess
import numpy as np

systemDict = np.load('dict_systematics.npy', allow_pickle='TRUE').item()

order_n = "3"

for key in systemDict.keys():
    cmd1 = "root -l -b -q resolutions.cxx\(\\\""+key+"\\\",\\\""+order_n+"\\\"\)"
    cmd2 = "mv resolutionInfo_INPUT.root resolutionInfo_INPUT_"+key+".root"

    sp = subprocess.Popen(cmd1, shell=True)
    sp.communicate()
    sp = subprocess.Popen(cmd2, shell=True)
    sp.communicate()

print("Done!")
