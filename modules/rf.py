# Differences in RandomForestRegressor compared to ranger:
# min_samples_split (=min.node.size), 2 (vs 5)

import os
import random
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from sklearn.utils import shuffle
from sklearn.ensemble import RandomForestRegressor
from scipy.stats import pearsonr, scoreatpercentile
from sklearn.inspection import permutation_importance
from sklearn.metrics import r2_score, explained_variance_score

from skranger.ensemble import RangerForestRegressor

from app.task_worker import progress

################################################################################
### Parameter optimization
################################################################################

def optim_rf(
  task_id,
  gt,
  phe,
  mtryvals,
  ntreevals):

    #repnum = 1
    #std = np.zeros((len(ntreevals), len(mtryvals)))
    pves = np.zeros((len(ntreevals), len(mtryvals)))
    for i, ntree in enumerate(ntreevals):
        for j, mtry in enumerate(mtryvals):
            progress(
              task_id,
              f"\rTesting {i*len(mtryvals)+j+1} of {len(ntreevals) * len(mtryvals)}",
              newline=False)
            #run = []
            #for _ in range(repnum):
            rf = RandomForestRegressor(
              n_estimators=ntree,
              max_features=mtry,
              oob_score=True,
              #max_depth=5,
              #n_jobs=os.cpu_count()-1
              n_jobs=-1
              )
            rf.fit(gt, phe)
            #run.append(rf.oob_score_)
            #pves[i, j] = np.mean(run)
            #std[i, j] = np.std(run)
            pves[i, j] = rf.oob_score_
    progress(task_id, "")
    return pves.tolist()

################################################################################
### Main RF function
################################################################################

def run_rf(
  task_id,
  gt,
  phe,
  mtry,
  ntree,
  niter,
  nperm):

################################################################################
### Run RandomForest with iterations based on correlation
################################################################################
    """
    ntree = 30
    mtry = round(0.1 * gt.shape[1])

    niter = 20
    iterdone = 0
    thresh = 0.000001
    rescor = 0

    incmse = np.zeros((gt.shape[1], niter))

    while rescor < thresh:
        for idx in range(iterdone, niter):
            print(f"Iteration {idx + 1} of {niter}", end="\r")
    #        rf = RandomForestRegressor(
    #            n_estimators=ntree,
    #            max_features=mtry,
    #            oob_score=True,
    #            n_jobs=-1)
            rf = RangerForestRegressor(
              n_estimators=ntree,
              mtry=mtry,
              oob_error=True,
              importance="permutation",
              scale_permutation_importance=True,
              seed=random.randint(0, 2**32 - 1))
            rf.fit(gt, phe)
            incmse[:, idx] = rf.feature_importances_

        rescor = np.corrcoef(
          np.mean(incmse[:, :niter // 2], axis=1),
          np.mean(incmse[:, niter // 2:], axis=1))[0,1]

        if rescor < thresh:
            incmse = np.hstack((incmse, np.zeros((gt.shape[1], 100))))
            iterdone = niter + 1
            niter += 100

        print(f"Correlation: {round(rescor, 2)}")
    """
################################################################################
### Run RandomForest with given n of iterations
################################################################################

    mtry = min(mtry, gt.shape[1])
    real_incmse = np.zeros((gt.shape[1], niter))
    for idx in range(niter):
        progress(
          task_id,
          f"\rIteration {idx + 1} of {niter}",
          newline=False)

        rf = RangerForestRegressor(
          n_estimators=ntree,
          mtry=mtry,
          oob_error=True,
          #min_node_size=5,
          importance="permutation",
          scale_permutation_importance=True,
          seed=random.randint(0, 2**32 - 1),
          #n_jobs=os.cpu_count()-1
          )
        rf.fit(gt, phe)
        real_incmse[:, idx] = rf.feature_importances_
    progress(task_id, "")

################################################################################
### Run permutation by correlation
################################################################################
    """
    niter = 20
    iterdone = 1
    thresh = 0.0000001
    permcor = 0

    perm_incmse = np.zeros((gt.shape[1], niter))

    while permcor < thresh:
        for idx in range(iterdone - 1, niter):
            print(f"Iteration {idx + 1} of {niter}", end="\r")
    #        perm_rf = RandomForestRegressor(
    #            n_estimators=ntree,
    #            max_features=mtry,
    #            oob_score=True)
            perm_rf = RangerForestRegressor(
              n_estimators=ntree,
              mtry=mtry,
              oob_error=True,
              importance="permutation",
              scale_permutation_importance=True,
              seed=random.randint(0, 2**32 - 1))
            shuffled_phe = shuffle(phe)
            perm_rf.fit(gt, shuffled_phe)
            perm_incmse[:, idx] = perm_rf.feature_importances_

        permcor = np.corrcoef(
            np.percentile(perm_incmse[:, :niter // 2], 99, axis=1),
            np.percentile(perm_incmse[:, niter // 2:], 99, axis=1))[0,1]

        if permcor < thresh:
            perm_incmse = np.hstack((perm_incmse, np.zeros((gt.shape[1], 100))))
            iterdone = niter + 1
            niter += 100

        print(f"Correlation: {round(permcor, 2)}")
    """
################################################################################
### Run permutation for given n of iterations
################################################################################

    nperm = nperm or 0
    perm_incmse = np.zeros((gt.shape[1], nperm))
    for idx in range(nperm):
        progress(
          task_id,
          f"\rPermutation {idx + 1} of {nperm}",
          newline=False)
        perm_rf = RangerForestRegressor(
          n_estimators=ntree,
          mtry=mtry,
          oob_error=True,
          importance="permutation",
          scale_permutation_importance=True,
          seed=random.randint(0, 2**32 - 1),
          #n_jobs=os.cpu_count()-1
          )
        perm_rf.fit(gt, shuffle(phe))
        perm_incmse[:, idx] = perm_rf.feature_importances_
    progress(task_id, "")

################################################################################
### Calculate importance thresholds
################################################################################

    meanincmse = np.mean(real_incmse, axis=1)
    if nperm == 0:
        incmse_final = meanincmse
    else:
        permquant = np.percentile(perm_incmse, 99, axis=1)
        #incmse_thresh = meanincmse > permquant
        #incmse_ratio = meanincmse / permquant
        incmse_final = meanincmse - permquant

#        print(meanincmse.tolist())
#        print(incmse_final.tolist())
#        print(np.corrcoef(meanincmse, incmse_final))

    # Save results
    output_df = pd.DataFrame({
#      "meanincmse": meanincmse,
#      "permquant": permquant,
#      "incmse_thresh": incmse_thresh,
#      "incmse_ratio": incmse_ratio,
      "Imp": incmse_final
    })
    return output_df