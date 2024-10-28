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
### Load data
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
### Parameter optimization
################################################################################
    """
    mtryvals = [
      10,
      round(np.sqrt(gt.shape[1])),
      2 * round(np.sqrt(gt.shape[1])),
      round(0.1 * gt.shape[1]),
      round(0.2 * gt.shape[1]),
      round(gt.shape[1] / 3)
    ]

    ntreevals = np.arange(10, 101, 10)
    repnum = 1

    opt = np.zeros((len(ntreevals), len(mtryvals)))
    std = np.zeros((len(ntreevals), len(mtryvals)))

    for i, ntree in enumerate(ntreevals):
        print(f"\r{ntree}", end="")
        for j, mtry in enumerate(mtryvals):
            run = []
            for _ in range(repnum):
                rf = RandomForestRegressor(
                  n_estimators=ntree,
                  max_features=mtry,
                  oob_score=True,
                  #max_depth=5,
                  n_jobs=-1)
                rf.fit(
                  gt,
                  phe)
                run.append(rf.oob_score_)
            opt[i, j] = np.mean(run)
            std[i, j] = np.std(run)

    # Plot the results
    colors = plt.cm.rainbow(np.linspace(0, 1, len(mtryvals)))
    for i in range(len(mtryvals)):
        plt.plot(
          ntreevals,
          opt[:, i],
          label=mtryvals[i],
          color=colors[i])

    if repnum > 1:
        for i in range(len(mtryvals)):
            plt.errorbar(
              ntreevals,
              opt[:, i],
              yerr=std[:, i],
              color=colors[i],
              fmt='-')

    plt.xlabel("ntree")
    plt.ylabel("PVE")
    plt.ylim([0, 1])
    plt.legend(title="mtry", loc="upper right")
    plt.show()
    """
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
        progress(task_id, f"Iteration {idx + 1} of {niter}")
        rf = RangerForestRegressor(
          n_estimators=ntree,
          mtry=mtry,
          oob_error=True,
          importance="permutation",
          scale_permutation_importance=True,
          seed=random.randint(0, 2**32 - 1))
        rf.fit(gt, phe)
        real_incmse[:, idx] = rf.feature_importances_

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

    perm_incmse = np.zeros((gt.shape[1], nperm))
    for idx in range(nperm):
        progress(task_id, f"Permutation {idx + 1} of {nperm}")
        perm_rf = RangerForestRegressor(
          n_estimators=ntree,
          mtry=mtry,
          oob_error=True,
          importance="permutation",
          scale_permutation_importance=True,
          seed=random.randint(0, 2**32 - 1))
        perm_rf.fit(gt, shuffle(phe))
        perm_incmse[:, idx] = perm_rf.feature_importances_

################################################################################
### Calculate importance thresholds
################################################################################

    meanincmse = np.mean(real_incmse, axis=1)
    permquant = np.percentile(perm_incmse, 99, axis=1)
    #incmse_thresh = meanincmse > permquant
    #incmse_ratio = meanincmse / permquant
    incmse_corrected = meanincmse - permquant

    # Save results
    output_df = pd.DataFrame({
#      "meanincmse": meanincmse,
#      "permquant": permquant,
#      "incmse_thresh": incmse_thresh,
#      "incmse_ratio": incmse_ratio,
      "Imp": incmse_corrected
    })
    return output_df