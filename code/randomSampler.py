"""
Script for random sampling 2,000 participants from a total of 3,821 participants with high 
myopia in at least one eye in the UK Biobank, as detailed in Figure 1 of the manuscript.

Author : Fabian Yii
Email  : fabian.yii@ed.ac.uk or fslyii@hotmail.com

2024
"""

import os
import numpy as np
import pandas as pd
import cv2 as cv
import matplotlib.pyplot as plt
from tqdm import tqdm

join = os.path.join
df   = pd.read_csv("cleaned_data_long_PM_cohort.csv")

# Only include eyes with high myopia (SER <= -5)
HM    = df[df.SER <= -5]
HM    = HM.sample(frac=1) # shuffle rows
HMids = np.unique(df[df.SER <= -5].id)
print(len(HM), "eyes of", len(HMids), "unique participants had high myopia") # 3821 participants had high myopia in at least one eye

# Randomly sample 2000 participants from the 3821 unique participants
np.random.seed(10)
randomHMids = np.random.choice(HMids, 2000, replace=False)

# 3024 eyes of 2000 participants were included
sampledHM = HM[HM.id.isin(randomHMids)]
print(len(sampledHM), "eyes randomly sampled")
bothEyesN = len(randomHMids) - (len(randomHMids)*2 - len(sampledHM))
print(bothEyesN, "participants had both eyes")
print(len(randomHMids)-bothEyesN, "participants only had one eligible eye")


