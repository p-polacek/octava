import os
import sys

from wtforms import (
    Form,
    FloatField,
    RadioField,
    SelectField,
    StringField,
    BooleanField,
    IntegerField,
)

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..")))

from modules.forms import models as models_dict

models_dict = dict(models_dict)

from modules.assoc.kq import kq_correct, kq_process

from modules.mem import get_mem
from modules.assoc_plot import mqplot

###############################################################################
### Initialize your model
###############################################################################

additional_models = [
    # Add a model id and display name,
    # in the format ("id": "Display name"),

    # Reserved ids:

    # B
    # bin
    # blink
    # bslmm
    # bt
    # C
    # glm
    # mlm
    # mtmlm
    # qt
    # rf
    # rfc
    # rr

]

class ParametersForm(Form):
    # Add any additional parameters required by your function
    # as flask WTForms fields, in the format:

    # field_id: YourFlaskField(field_label, **{field_option: field_value})

    # For IntegerField and FloatField, it's recommended to add render_kw:
    #    IntegerField(render_kw={"class": "small-num-input"})
    #    FloatField(render_kw={"class": "small-num-input"})

    # For SelectField and StringField, it's recommended to add render_kw:
    #    SelectField(render_kw={"class": "default-input"})
    #    StringField(render_kw={"class": "default-input"})

    # Reserved field_ids:
	# b_burn
	# b_chain
	# b_thin
	# by_phe
	# chr
	# chrom
	# coord
	# cov
	# end
	# fmt
	# gene
	# hgl_coord_group
	# hgl_gene_group
	# k
	# kq_from
	# markers
	# markers_sgl
	# model_sgl
	# models
	# mtry
	# ntree
	# overwrite
	# pc_nums
	# plot
	# plot_thr
	# plothap
	# q
	# q_num
	# qtl_runs
	# r2_ld
	# run_ext
	# sample_group
	# start
	# trait
	# trait_group
	# window_ld


    pass

parameter_groups = {
    # For each model, associate model_id with its parameter fields
    # in the format: "model_id": [field_ids],

}

###############################################################################
### Setup your model
###############################################################################

# Use this template to define your own pipeline

def template_for_additional(genome, params, task_id, mem=0):
    #################################
    ### Better keep this part
    #################################

    config = load_config(genome, load_ext=params["run_ext"])
    gdata = load_gdata(genome, load_ext=params["run_ext"])

    markers = params["markers_sgl"]
    model = params["model_sgl"]
    # Specify memory allowance for your model
    mem = mem or get_mem()

    if params["run_ext"]:
        ext_data_dir = data_dir
        fn_data_dir = tmp_dir
    else:
        fn_data_dir = data_dir

    #################################
    ### 
    #################################

    assoc_dir = os.path.join(fn_data_dir, genome, "assoc", markers)
    out_dir = os.path.join(assoc_dir, "res")
    in_root = os.path.join(assoc_dir, "gt", "gt")

    # Load samples and phenotype data (example)
    fam_file = f"{in_root}.fam"
    fam_df = pd.read_csv(fam_file, sep="\t", header=None, index_col=None)
    samples = fam_df[1]

    phe_path = os.path.join(fn_data_dir, genome, config["phe"])
    with portalocker.Lock(phe_path, "r", flags=portalocker.LOCK_SH) as locked_f:
        phe_df = pd.read_csv(locked_f, sep="\t", header=0, index_col=0).loc[samples]

    # Process and save genotype and phenotype in required format


    # Process and save correction data (you can use the kq_process() function)
    q_matr, k_matr = kq_process(
        genome,
        params["kq_from"] or markers,
        samples,
        params["k"],
        int(params["q_num"]) if params["q"] else 0,
        params["cov"],
    )

    # Run your association analysis model
    for trait in traits:

    # Load results into pandas DataFrame


    # Save raw output as <model_id>.<trait>.tsv.gz


    # Save top 5,000 markers under <model_id>.<trait>.feather
    # Columns = [-log10(p-value), chrom, pos, *[additional columns]]
    # Optional: save top 5,000 marker effect sizes under <model_id>.<trait>.eff.feather
    # Columns = [abs(effect_size), chrom, pos, *[additional columns]]
    
    # For haplotype markers:
    # Columns = [-log10(p-value) OR abs(effect_size), chrom, pos, end, idx, *[additional columns]]

        if is_hap:
            blockmap = pd.read_csv(
                os.path.join(assoc_dir, "gt", "gt.blockmap"),
                sep="\t",
                header=0,
                index_col=None,
            )
            blockmap["chr"] = blockmap["chr"].astype(str)
            res = pd.merge(
                res,
                blockmap,
                left_on=[ch_col, pos_col],
                right_on=["chr", "start"],
                how="left",
            ).drop(["start"], axis=1)
            res["idx"] = res.groupby([ch_col, pos_col]).cumcount() + 1
        else:
            # your processing
            pass

        (
            res[sig_col, chr_col, pos_col], # [end_col, idx_col]
            .sort_values(by=sig_col, ascending=False)
            .head(5000)
            .reset_index(drop=True)
            .to_feather(out_feather)
        )


        # To generate static Manhattan and QQ-plots,
        # p-value column of res must be -log10-scaled already
        # mqplot_fn = 
        # qqplot and thresh only if applicable to your model

        qqplot = True or False
        thresh = True or False or float()

###############################################################################
### Standard process of generating static plots and saving parameters
###############################################################################

        mqplot(
            genome,
            res,
            f"{model}.{trait}.mqplot.png",
            ", ".join((markers, models_dict[model], trait_name))
            json.loads(config["concat"]),
            ch_col=ch_col,
            pos_col=pos_col,
            sig_col=sig_col,
            qqplot=qqplot,
            thresh=thresh,
            end_col="end" if is_hap else None,
        )

        form = ParametersForm()

        # Save parameters
        run_params = [
            [
                "Markers",
                "Samples",
                "PCs",
                "Kinship",
                "K/Q source",
                "Covariate files",
                *[getattr(form, field).label.text.strip(":")
                  for field in parameter_groups[model]]
            ],
            [
                res.shape[0],
                phe_samples.shape[0],
                params["q_num"] if params["q"] else 0,
                "Yes" if params["k"] else "No",
                params["kq_from"] or markers,
                "<br>".join(params["cov"]),
                *[params[field] for field in parameter_groups[model]]
            ],
        ]
        with open(
            os.path.join(out_dir, f"{model}.{trait}.params"), "w"
        ) as params_out:
            json.dump(run_params, params_out)


function_dict = {
    # Associate your model_id with function name
    # in the format "model_id": function_name,

}


def run_additional(genome, params, task_id, mem=0):
    model = function_dict[params["model_sgl"]]
    model(genome, params, task_id, mem=mem)

