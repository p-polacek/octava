import io
import os
import sys
import json
import random
import shutil
import bed_reader
import subprocess
import portalocker
import numpy as np
import pandas as pd
from skranger.ensemble import RangerForestRegressor

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..")))

from modules.load_input import (
  load_gdata,
  load_config)

from app.dirs import (
  tmp_dir,
  data_dir,
  tool_dir,
  plink2_path)

from modules.assoc.kq import (
  kq_correct,
  kq_process)

from modules.mem import get_mem
from modules.assoc_plot import mqplot
from modules.forms import models as models_dict
from modules.assoc.assoc_process import assoc_process

models_dict = dict(models_dict)

def run_assoc(genome, params, task_id):
    config  = load_config(genome)
    gdata   = load_gdata(genome)
    markers = params["markers_sgl"]
    model   = params["model_sgl"]

    print("<b>Preparing input files</b>")
    assoc_dir = os.path.join(
      data_dir,
      genome,
      "assoc",
      markers)
    out_dir = os.path.join(
      assoc_dir,
      "res")
    in_root = os.path.join(
      assoc_dir,
      "gt",
      "gt")

    # Load samples and phenotype data
    fam_file = f"{in_root}.fam"
    fam_df = pd.read_csv(
      fam_file,
      sep="\t",
      header=None,
      index_col=None)
    samples = fam_df[1]

    phe_path = os.path.join(
      data_dir,
      genome,
      config["phe"])
    with portalocker.Lock(phe_path, "r", flags=portalocker.LOCK_SH) as locked_f:
        phe_df = (pd
          .read_csv(
            locked_f,
            sep="\t",
            header=0,
            index_col=0)
         .loc[samples])

    if model == "blink":
        # BLINK analyses either one or all traits
        # phenotype file rewritten each time to keep only selected traits
        (phe_df
          .fillna("NaN")
          .loc[:, params["trait"]]
          .rename_axis("Acc")
          .to_csv(
            f"{in_root}.txt",
            sep="\t"))

        params["k"] = True
        q_matr, k_matr = kq_process(
          genome,
          params["kq_from"] or markers,
          samples,
          False,
          int(params["q_num"]) if params["q"] else 0,
          params["cov"])

        if not q_matr.empty:
          q_matr.to_csv(
            f"{in_root}.cov",
            sep="\t")
        elif os.path.exists(f"{in_root}.cov"):
            os.remove(f"{in_root}.cov")

        blink_path = os.path.join(
          tool_dir,
          "BLINK",
          "blink_linux")
        print("<b>Running analysis</b>")
        command = [
          blink_path,
          "--gwas",
          "--file", in_root,
          #"--bin_size", "5", "10", "2", "0.3",
          "--out", os.path.join(out_dir, "BLINK"),
          #"--parallel", "1",
          #"--threads", "6",
          "--trait", "0"]

#        if os.path.exists(f"{in_root}.dat"):
#            command += ["--numeric"]
#        else:
        command += ["--plink"]
        subprocess.run(command, check=True)

        print("<b>Processing raw output</b>")
        num_traits = len(params["trait"])
        for trait_idx, trait in enumerate(params["trait"], start=1):
            print(f"{trait} ({trait_idx}/{num_traits})")
            if not os.path.exists(os.path.join(out_dir, f"BLINK_{trait}_GWAS_result.txt")):
                continue
            # Save parameters
            run_params = [
              ["Markers", "Samples", "PCs", "Kinship", "K/Q source", "Covariate files"],
              [sum(1 for l in open(os.path.join(out_dir, f"BLINK_{trait}_GWAS_result.txt"))) - 1,
               phe_df.loc[samples, trait].dropna().size,
               params["q_num"] if params["q"] else 0,
               "Yes",
               params["kq_from"] or markers,
               "<br>".join(params["cov"])]]
            with open(
              os.path.join(
                out_dir,
                f"{model}.{trait}.params"),
              "w") as params_out:
                json.dump(run_params, params_out)

            mqplot_title = ", ".join((markers, models_dict[model], trait))
            assoc_process(
              genome,
              assoc_dir,
              f"BLINK_{trait}_GWAS_result.txt",
              mqplot_title,
              f"{model}.{trait}",
              out_dir,
              ch_col="chr",
              pos_col="pos",
              sig_col="p_value")

        if os.path.exists(os.path.join(out_dir, "BLINK.log")):
            os.remove(os.path.join(out_dir, "BLINK.log"))

    elif model == "rf":
        # Correct for kinship and population structure
        phe_df = phe_df.loc[samples, params["trait"]]

        gt = bed_reader.open_bed(f"{in_root}.bed").read()

        print("<b>Running analysis</b>")
        num_traits = len(params["trait"])
        for trait_idx, trait in enumerate(params["trait"], start=1):
            print(f"{trait} ({trait_idx}/{num_traits})")
            raw_phe = (phe_df
              .loc[samples, trait]
              .dropna())
            phe_samples = raw_phe.index
            phe_gt = gt[np.isin(samples, phe_samples), :]

            phe_corr = kq_correct(
              genome,
              raw_phe,
              phe_samples,
              params["k"],
              int(params["q_num"]) if params["q"] else 0,
              params["kq_from"] or markers,
              params["cov"])

            if 1:
                mtry = min(params["mtry"], phe_gt.shape[1])

                rf = RangerForestRegressor(
                  n_estimators=params["ntree"],
                  mtry=mtry,
                  oob_error=True,
                  min_node_size=5,
                  importance="permutation",
                  scale_permutation_importance=True,
                  seed=random.randint(0, 2**32 - 1),
                  n_jobs=os.cpu_count() - 1,
                  )
                rf.fit(phe_gt, phe_corr)
                res = pd.DataFrame({"Imp": rf.feature_importances_})

            else:
                try:
                    res = pd.read_csv(f"../rf_raw/uncorr/{trait}.txt", sep="\t", index_col=0)
                except:
                    continue

            # Output files
            processed_out = os.path.join(out_dir, f"{model}.{trait}")
            mqplot_fn     = f"{processed_out}.mqplot.png"
            out_feather   = f"{processed_out}.feather"
            out_gz        = f"{processed_out}.tsv.gz"

            # Save raw results, rename numeric chromosomes back to original
            marker_map = (pd
              .read_csv(
                f"{in_root}.bim",
                sep="\t",
                header=None,
                index_col=None,
                usecols=[0,3])
              .set_axis(
                ["Chr", "Pos"],
                axis=1))

            chr_rename = dict(enumerate(gdata["chr_names"], start=1))
            res["Chr"]  = marker_map["Chr"].map(chr_rename)
            res["Pos"] = marker_map["Pos"]
            res = res.loc[:, ["Chr", "Pos", "Imp"]]
            res.to_csv(
              out_gz,
              sep="\t",
              index=False)

            with open(os.path.join(assoc_dir, "stats", "params")) as params_in:
                is_hap = json.load(params_in)[1][10] == "Yes"
            if is_hap:
                blockmap = (pd
                  .read_csv(
                    os.path.join(assoc_dir, "gt", "gt.blockmap"),
                    sep="\t",
                    header=0,
                    index_col=None))
                blockmap["chr"] = blockmap["chr"].astype(str)
                res = (pd
                  .merge(
                    res,
                    blockmap,
                    left_on=["Chr", "Pos"],
                    right_on=["chr", "start"],
                    how="left")
                  .drop(
                    ["start"],
                    axis=1))
                res["idx"] = res.groupby(["Chr", "Pos"]).cumcount() + 1
                c = res.columns.tolist()
                (res[c[2:3] + c[:2] + c[-2:]]
                  .sort_values(
                    by="Imp",
                    ascending=False)
                  .head(5000)
                  .reset_index(drop=True)
                  .to_feather(out_feather))
            else:
                c = res.columns.tolist()
                (res[c[2:3] + c[:2]]
                  .sort_values(
                    by="Imp",
                    ascending=False)
                  .head(5000)
                  .reset_index(drop=True)
                  .to_feather(out_feather))

            mqplot_title = ", ".join((markers, models_dict[model], trait))
            mqplot_fig = mqplot(
              genome,
              res,
              mqplot_fn,
              mqplot_title,
              json.loads(config["concat"]),
              ch_col="Chr",
              pos_col="Pos",
              sig_col="Imp",
              thresh=False,
              qqplot=False,
              end_col="end" if is_hap else None)

            # Save parameters
            run_params = [
              ["Markers", "Samples", "PCs", "Kinship", "K/Q source",
               "Covariate files", "Test markers", "Trees"],
              [res.shape[0],
               phe_samples.shape[0],
               params["q_num"] if params["q"] else 0,
               "Yes" if params["k"] else "No",
               params["kq_from"] or markers,
               "<br>".join(params["cov"]),
               params["mtry"],
               params["ntree"]]]
            with open(
              os.path.join(
                out_dir,
                f"{model}.{trait}.params"),
              "w") as params_out:
                json.dump(run_params, params_out)

    elif model in ("B", "C"):
        # GCTB requires space-separated .fam and cofactor/kinship files
        for ext in ("bed", "bim"):
            if not os.path.exists(f"{in_root}.gctb.{ext}"):
                os.link(f"{in_root}.{ext}", f"{in_root}.gctb.{ext}")

        (fam_df
          .assign(IID = samples.tolist())
          .loc[:, ["IID"] + fam_df.columns.tolist()[1:]]
          .to_csv(
            f"{in_root}.gctb.fam",
            sep=" ",
            header=False,
            index=False))
        (phe_df
          .assign(IID = phe_df.index.tolist())
          .loc[:, ["IID"] + params["trait"]]
          .to_csv(
            f"{in_root}.txt",
            header=False,
            sep=" "))

        q_matr, k_matr = kq_process(
          genome,
          params["kq_from"] or markers,
          samples,
          params["k"],
          int(params["q_num"]) if params["q"] else 0,
          params["cov"])
        if not q_matr.empty:
            (q_matr
              .assign(IID = samples.tolist())
              .loc[
                :,
                ["IID"] + [c for c in q_matr.columns.tolist()]]
              .to_csv(
                f"{in_root}.cov",
                header=False,
                sep=" "))
        elif os.path.exists(f"{in_root}.cov"):
            os.remove(f"{in_root}.cov")

        if params["k"]:
            (k_matr
              .assign(IID = samples.tolist())
              .loc[:, ["IID"] + k_matr.columns.tolist()]
              .to_csv(
                f"{in_root}.kin",
                sep=" ",
                header=False,
                index=True))

        gctb_path = os.path.join(
          tool_dir,
          "gctb_2.5.2_Linux",
          "gctb")
        print("<b>Running analysis</b>")

        num_traits = len(params["trait"])
        for trait_idx, trait in enumerate(params["trait"], start=1):
            print(f"{trait} ({trait_idx}/{num_traits})")
            cmd = [
              gctb_path,
              "--bfile", f"{in_root}.gctb",
              "--pheno", f"{in_root}.txt",
              "--mpheno", str(trait_idx),
              "--bayes", model,
              "--maf", "0",
              "--p-value", "1",
              "--out", os.path.join(out_dir, f"{model}.{trait}")]
            cmd = cmd + ["--burn-in", str(params["b_burn"])] if params["b_burn"] else cmd
            cmd = cmd + ["--chain-length", str(params["b_chain"])] if params["b_chain"] else cmd
            cmd = cmd + ["--thin", str(params["b_thin"])] if params["b_thin"] else cmd
            cmd = cmd + ["--covar", f"{in_root}.cov"] if not q_matr.empty else cmd
            cmd = cmd + ["--random-covar", f"{in_root}.kin"] if params["k"] else cmd
            pr = subprocess.run(cmd, check=True)

            # Save parameters
            run_params = [
              ["Markers", "Samples", "PCs", "Kinship", "K/Q source", "Covariate files",
               "Burn-in iterations", "Total iterations", "Thinning"],
              [sum(1 for l in open(os.path.join(out_dir, f"{model}.{trait}.snpRes"))) - 1,
               phe_df.loc[samples, trait].dropna().shape[0],
               params["q_num"] if params["q"] else 0,
               "No",
               params["kq_from"] or markers,
               "<br>".join(params["cov"]),
               params["b_burn"] if params["b_burn"] else 1000,
               params["b_chain"] if params["b_chain"] else 3000,
               params["b_thin"] if params["b_thin"] else 10]]
            with open(
              os.path.join(
                out_dir,
                f"{model}.{trait}.params"),
              "w") as params_out:
                json.dump(run_params, params_out)

            mqplot_title = ", ".join((markers, models_dict[model], trait))
            assoc_process(
              genome,
              assoc_dir,
              f"{model}.{trait}.snpRes",
              mqplot_title,
              f"{model}.{trait}",
              out_dir,
              ch_col="Chrom",
              pos_col="Position",
              sig_col="PIP",
              sig_log=False,
              eff_col=["A1Effect"],
              sep="\s+")

        for ext in ("covRes", "parRes", "randCovRes"):
            fn = os.path.join(out_dir, f"{model}.{trait}.{ext}")
            if os.path.exists(fn):
                os.remove(fn)

    elif model in ("bt", "qt"):
        (phe_df
          .fillna("NA")
          .assign(FID=phe_df.index)
          .assign(IID=phe_df.index)
          [["FID", "IID"] + params["trait"]]
          .to_csv(
            f"{in_root}.txt",
            sep="\t",
            index=False))

        q_matr, k_matr = kq_process(
          genome,
          params["kq_from"] or markers,
          samples,
          False,
          int(params["q_num"]) if params["q"] else 0,
          params["cov"])
        if not q_matr.empty:
            (q_matr
              .assign(FID=q_matr.index)
              .assign(IID=q_matr.index)
              .loc[:, ["FID", "IID"] + q_matr.columns.tolist()]
              .to_csv(
                f"{in_root}.cov",
                sep="\t",
                index=False))

        concat = json.loads(config["concat"])
        # copy/link bed
        if not os.path.exists(f"{in_root}.regenie.bed"):
            if concat:
                shutil.copy(f"{in_root}.bed", f"{in_root}.regenie.bed")
            else:
                os.link(f"{in_root}.bed", f"{in_root}.regenie.bed")
        # add FID to fam
        if not os.path.exists(f"{in_root}.regenie.fam"):
            fam = pd.read_csv(
              f"{in_root}.fam",
              sep="\t",
              header=None,
              usecols=[1,2,3,4,5],
              index_col=0)
            (fam
              .assign(FID=fam.index)
              .assign(IID=fam.index)
              .loc[:, ["FID", "IID"] + fam.columns.tolist()]
              .to_csv(
                f"{in_root}.regenie.fam",
                sep="\t",
                header=False,
                index=False))
        # concat contigs in bim
        if not os.path.exists(f"{in_root}.regenie.bim"):
            if concat:
                bim = pd.read_csv(
                  f"{in_root}.bim",
                  sep="\t",
                  header=None)
                bim_orig = bim.copy()

                concat_idxs = [
                  i+1 for i, c in enumerate(gdata["chr_names"])
                  if c in concat]
                zeroed_offsets = [
                  o - gdata["chr_offsets"][concat_idxs[0]-1]
                  for o in gdata["chr_offsets"]]

                mask = bim[0].isin(concat_idxs)
                bim.loc[mask, 3] += np.array(zeroed_offsets)[bim.loc[mask, 0]-1]
                bim.loc[mask, 0] = concat_idxs[0]
                bim.loc[mask, 1] = bim.loc[mask, 0].astype(str) + ":" + bim.loc[mask, 3].astype(str)
                bim.to_csv(
                  f"{in_root}.regenie.bim",
                  sep="\t",
                  index=False,
                  header=False)

                # Create a conversion dict for mapping to original chr and pos
                bim_map = pd.DataFrame({
                    "chr_new": bim[0],
                    "pos_new": bim[3],
                    "chr_old": bim_orig[0],
                    "pos_old": bim_orig[3]
                })
                bim_map.to_feather(f"{in_root}.bim_map.feather")

                cmd = [
                  plink2_path,
                  "--memory", f"{get_mem()}",
                  "--bfile", f"{in_root}.regenie",
                  "--allow-extra-chr",
                  "--sort-vars",
                  "--make-pgen",
                  "--out", f"{in_root}.regenie"]
                subprocess.run(cmd, check=True)

                cmd = [
                  plink2_path,
                  "--memory", f"{get_mem()}",
                  "--pfile", f"{in_root}.regenie",
                  "--allow-extra-chr",
                  "--make-bed",
                  "--out", f"{in_root}.regenie"]
                subprocess.run(cmd, check=True)

            else:
                os.link(f"{in_root}.bim", f"{in_root}.regenie.bim")

        ## thin down markers to 500,000
        if params["k"]:
            with open(os.path.join(assoc_dir, "stats", "cnts")) as cnts_in:
                num_markers = json.load(cnts_in)[1]
            if int(num_markers) > 500000:
                if any(
                  not os.path.exists(f"{in_root}.regenie.step1.{ext}")
                  for ext in ("bed", "bim", "fam")):
                    open(f"{in_root}.regenie.step1.bed", "a").close()
                    cmd = [
                      plink2_path,
                      "--memory", f"{get_mem()}",
                      "--bfile", f"{in_root}.regenie",
                      "--thin-count", "500000",
                      "--allow-extra-chr",
                      "--make-bed",
                      "--out", f"{in_root}.regenie.step1"]
                    subprocess.run(cmd, check=True)
                s1_path = f"{in_root}.regenie.step1"
            else:
                s1_path = f"{in_root}.regenie"

        regenie_path = os.path.join(
          tool_dir,
          "regenie")
        print("<b>Running analysis</b>")

        if params["k"]:
            command = [
              regenie_path,
              "--step", "1",
              "--bed", s1_path,
              "-p", f"{in_root}.txt",
              "--phenoCol", ",".join(params["trait"]),
              f"--{model}",
              "--bsize", "1000",
              "--out", os.path.join(tmp_dir, "regenie")]
            if not q_matr.empty:
                command += [
                  "--covarFile", f"{in_root}.cov",
                  "--covarColList", ",".join(q_matr.columns.tolist())]
            pr = subprocess.run(
              command,
              check=True)

        command = [
          regenie_path,
          "--step", "2",
          "--bed", f"{in_root}.regenie",
          "-p", f"{in_root}.txt",
          "--phenoCol", ",".join(params["trait"]),
          f"--{model}",
          "--bsize", "1000",
          "--out", os.path.join(out_dir, "regenie")]
        if not q_matr.empty:
            command += [
              "--covarFile", f"{in_root}.cov",
              "--covarColList", ",".join(q_matr.columns.tolist())]
        if params["k"]:
            command += ["--pred", os.path.join(tmp_dir, "regenie_pred.list")]
        else:
            command += ["--ignore-pred"]
        pr = subprocess.run(
          command,
          check=True)

        # Save parameters
        num_traits = len(params["trait"])
        for trait_idx, trait in enumerate(params["trait"], start=1):
            print(f"{trait} ({trait_idx}/{num_traits})")
            run_params = [
              ["Markers", "Samples", "PCs", "Kinship", "K/Q source", "Covariate files"],
              #[sum(1 for l in open(os.path.join(out_dir, f"regenie_{trait}.regenie"))) - 1,
              [sum(1 for l in open(f"modules/tools/regenie_test/bin_all/out/regenie_{trait}.regenie")) - 1, # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
               phe_df.loc[samples, trait].dropna().shape[0],
               params["q_num"] if params["q"] else 0,
               "Yes" if params["k"] else "No",
               params["kq_from"] or markers,
               "<br>".join(params["cov"])]]
            with open(
              os.path.join(
                out_dir,
                f"{model}.{trait}.params"),
              "w") as params_out:
                json.dump(run_params, params_out)

            # Output files
            processed_out = os.path.join(out_dir, f"{model}.{trait}")
            mqplot_fn     = f"{processed_out}.mqplot.png"
            out_feather   = f"{processed_out}.feather"
            out_gz        = f"{processed_out}.tsv.gz"

            bim_map = pd.read_feather(f"{in_root}.bim_map.feather")
            bim_map = bim_map.reset_index()

            res = (pd
              .read_csv(
                f"modules/tools/regenie_test/bin_all/out/regenie_{trait}.regenie",
                sep=" ",
                header=0,
                index_col=None)
              .loc[:, ["CHROM","GENPOS","BETA","SE","LOG10P"]])

            res = res.merge(
                bim_map,
                left_on=["CHROM", "GENPOS"],
                right_on=["chr_new", "pos_new"],
                how="left")
            res = (res
              .sort_values(by="index")
              .drop(["index"], axis=1))

            # Save raw results, rename numeric chromosomes back to original
            chr_rename = dict(enumerate(gdata["chr_names"], start=1))

            res["CHROM"] = res["chr_old"].map(chr_rename)
            res["GENPOS"] = res["pos_old"]
            res = res.drop(["chr_old", "pos_old", "chr_new", "pos_new"], axis=1)

            res.to_csv(
              out_gz,
              sep="\t",
              index=False)

            # Save results for MP visualisation
            with open(os.path.join(assoc_dir, "stats", "params")) as params_in:
                is_hap = json.load(params_in)[1][10] == "Yes"
            if is_hap:
                blockmap = (pd
                  .read_csv(
                    os.path.join(assoc_dir, "gt", "gt.blockmap"),
                    sep="\t",
                    header=0,
                    index_col=None))
                blockmap["chr"] = blockmap["chr"].astype(str)
                res = (pd
                  .merge(
                    res,
                    blockmap,
                    left_on=["CHROM", "GENPOS"],
                    right_on=["chr", "start"],
                    how="left")
                  .drop(
                    ["start"],
                    axis=1))
                res["idx"] = res.groupby(["CHROM", "GENPOS"]).cumcount() + 1
                c = res.columns.tolist()
                (res[["LOG10P", "CHROM", "GENPOS", "end", "BETA", "SE"]]
                  .sort_values(
                    by="LOG10P",
                    ascending=False)
                  .head(5000)
                  .reset_index(drop=True)
                  .to_feather(out_feather))
            else:
                c = res.columns.tolist()
                (res[["LOG10P", "CHROM", "GENPOS", "BETA", "SE"]]
                  .sort_values(
                    by="LOG10P",
                    ascending=False)
                  .head(5000)
                  .reset_index(drop=True)
                  .to_feather(out_feather))
            mqplot_title = ", ".join((markers, models_dict[model], trait))
            mqplot(
              genome,
              res,
              mqplot_fn,
              mqplot_title,
              json.loads(config["concat"]),
              ch_col  = "CHROM",
              pos_col = "GENPOS",
              sig_col = "LOG10P",
              qqplot  = True,
              thresh  = True,
              end_col = "end" if is_hap else None)

        log_file = os.path.join(out_dir, "regenie.log")
        if os.path.exists(log_file):
            os.remove(log_file)

    else:
        # GEMMA cannot accept row or column names in phenotype file
        # phenotype file rewritten each time to keep only selected traits
        (phe_df
          .fillna("NA")
          # GEMMA treats -9 as a missing value                                    #################################
          .replace(-9, -9.000000000001)
          .replace("-9", "-9.000000000001")
          .to_csv(
            f"{in_root}.txt",
            sep="\t",
            header=False,
            index=False))

        q_matr, k_matr = kq_process(
          genome,
          params["kq_from"] or markers,
          samples,
          model != "glm",
          int(params["q_num"]) if params["q"] else 0,
          params["cov"])
        if not q_matr.empty:
            (q_matr
              .assign(I=1)
              .loc[:, ["I"] + q_matr.columns.tolist()]
              .to_csv(
                f"{in_root}.cov",
                sep="\t",
                header=False,
                index=False))
        elif os.path.exists(f"{in_root}.cov"):
            os.remove(f"{in_root}.cov")

        if not k_matr.empty:
            (k_matr
              .to_csv(
                f"{in_root}.kin",
                sep="\t",
                header=False,
                index=False))

#        k_matr = k_matr.loc[k_lab.isin(samples), k_lab.isin(samples)]
#        eigenvalues, eigenvectors = np.linalg.eigh(k_matr)
#        y = phe_df.loc[samples, params["trait"][0]]
#        sigma_g2 = np.var(eigenvalues)
#        sigma_e2 = np.var(y.dropna()) - sigma_g2

        print("<b>Running analysis</b>")
        if model == "mtmlm":
            traits = params["trait_group"]
        else:
            traits = [([t], t) for t in params["trait"]]

        num_traits = len(params["trait"])
        for trait_idx, (trait, trait_name) in enumerate(traits, start=1):
            print(f"{trait_name} ({trait_idx}/{num_traits})")

            n = [str(phe_df.columns.get_loc(t) + 1) for t in trait]
            command = [
              "gemma",
              "-p", f"{in_root}.txt",
              "-n", *n,
              "-o", f"{model}.{trait_name}",
              "-outdir", out_dir,
              "-miss", "1",
              "-maf", "0", #0.04",
              "-hwe", "0"]
            if model == "glm":
                command += ["-lm", "1"]
            else:
                if model == "bslmm":
                    command += ["-bslmm", "1"]
                elif model == "rr":
                    command += ["-bslmm", "2"]
                elif model == "bin":
                    command += ["-bslmm", "3"]
                else:
                    command += ["-lmm", "1"]
                command += ["-k", f"{in_root}.kin"]
            #if not q_matr.empty and model not in ("bslmm", "rr", "bin"):
            #    command += ["-c", f"{in_root}.cov"]

            if os.path.exists(f"{in_root}.bimbam"):
                command += [
                  "-g", f"{in_root}.bimbam",
                  "-a", f"{in_root}.ann"]
            else:
                command += [
                  "-bfile", in_root]

            subprocess.run(command, check=True)

            # Save parameters
            run_params = [
              ["Markers", "Samples", "PCs", "Kinship", "K/Q source", "Covariate files"],
              [sum(1 for l in open(os.path.join(out_dir, f"{model}.{trait_name}.assoc.txt"))) - 1,
               phe_df.loc[samples, trait].dropna().shape[0],
               params["q_num"] if params["q"] else 0,
               "No" if model == "glm" else "Yes",
               params["kq_from"] or markers,
               "<br>".join(params["cov"])]]
            if model == "mtmlm":
                run_params[0].append("Traits")
                run_params[1].append("<br>".join(trait))
            with open(
              os.path.join(
                out_dir,
                f"{model}.{trait_name}.params"),
              "w") as params_out:
                json.dump(run_params, params_out)

            mqplot_title = ", ".join((markers, models_dict[model], trait_name))
            assoc_process(
              genome,
              assoc_dir,
              f"{model}.{trait_name}.assoc.txt",
              mqplot_title,
              f"{model}.{trait_name}",
              out_dir,
              ch_col="chr",
              pos_col="ps",
              sig_col="p_wald",
              eff_col=(
                ["beta"] if model in ("mlm", "glm")
                else ["beta_" + str(i) for i in range(1, len(n)+1)] if model == "mtmlm"
                else []),
              se_col="se" if model == "mlm" else None)

            log_file = os.path.join(out_dir, f"{model}.{trait_name}.log.txt")
            if os.path.exists(log_file):
                os.remove(log_file)

if __name__ == "__main__":
    genome, params, task_id = sys.argv[1:4]
    params = json.loads(params)
    run_assoc(genome, params, task_id)