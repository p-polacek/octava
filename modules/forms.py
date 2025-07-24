#!/usr/bin/env python

from markupsafe import Markup
from flask_wtf import FlaskForm

from wtforms import (
  widgets,
  FileField,
  FloatField,
  RadioField,
  HiddenField,
  SelectField,
  StringField,
  BooleanField,
  IntegerField,
  SelectMultipleField)

from modules.csq import var_csq

models = [
  ("glm", "GLM"),
  ("mlm", "MLM"),
  #("qt", "Fast MLM"),
  ("bt", "Binary trait MLM"),
  ("mtmlm", "Multi-trait MLM"),
  #("bslmm", "Bayesian sparse LMM"),
  #("rr", "Ridge regression / GBLUP"),
  #("bin", "Binary trait MLM"),
  ("blink", "BLINK"),
  ("rf", "Random forest"),
  ("B", "BayesB"),
  ("C", "BayesC")]

types = [
  ("snp",
   "Single nucleotide variant"),
  ("mnp",
   "Multi-nucleotide variant"),
  ("ins",
   "Insertion"),
  ("del",
   "Deletion"),
  ("dup",
   "Duplication"),
  ("inv",
   "Inversion"),
  ("cnv",
   "Copy number variant"),
  ("bnd",
   "Breakend variant")]

csqs = var_csq()

class MultiCheckboxField(SelectMultipleField):
    widget = widgets.ListWidget(prefix_label=False)
    option_widget = widgets.CheckboxInput()

class AddForm(FlaskForm):
    name = SelectField("Name:", choices=[])

    setup_from = SelectField(
      "Dataset:",
      choices=[])
    vcf = SelectField(
      "Variants:",
      choices=[])
    ann = SelectField(
      "Annotation:",
      choices=[])
    seq = SelectField(
      "Sequence:",
      choices=[])
    phe = SelectField(
      "Phenotype:",
      choices=[])
    exp = SelectField(
      "Gene expression:",
      choices=[])
    run = SelectField(
      "Run info:",
      choices=[])
    #proj = SelectField(
    #  "Study info:",
    #  choices=[])
    samp = SelectField(
      "Sample info:",
      choices=[])
    info = SelectField(
      "Gene info:",
      choices=[])
    #enr = MultiCheckboxField(
    #  "Enrichment for:",
    #  choices=[])
    bam = MultiCheckboxField(
      "BAM file(s):",
      choices=[])
    snpeff = BooleanField(
      "Annotate variants with SnpEff:")
    cds = SelectField(
      "CDS:",
      choices=[])
    pep = SelectField(
      "Proteins:",
      choices=[])
    sev = BooleanField(
      "Enable severity prediction:")
    #db_fasta = SelectField(
    #  "DB FASTA file:",
    #  choices=[])
    db_prefix = SelectField(
      "DB prefix:",
      choices=[])
    concat = MultiCheckboxField(
      "Concatenate contigs:",
      choices=[])

class EditForm(FlaskForm):
    phe = FileField(
      "File:")
    cof = FileField(
      "File:")

class RegionForm(FlaskForm):
    fmt = SelectField(
      "Format:",
      choices = [
        ("region",
         "Region"),
        ("gene",
         "Gene ID"),
        ("coord",
         "Coordinates")])
    chrom = SelectField(
      "Chrom:",
      choices = [])
    start = IntegerField(
      "Start:")
    end = IntegerField(
      "End:")
    coord = StringField(
      "")
    gene = StringField(
      "Gene:")

class FuncAnalysisForm(RegionForm):
    dataset = SelectField(
      "Markers:",
      choices = [("", "raw VCF")])
    target = SelectField(
      "Target:",
      choices = [
        ("full",
         "Entire region"),
         ("cds",
         "Coding regions"),
        ("prom",
         "Promoter regions")])
    gene_group = MultiCheckboxField(
      "Gene group:",
      choices = [])
    coord_group = MultiCheckboxField(
      "Region group:",
      choices = [])
    prom_size = IntegerField(
      "Size (bp):",
      default = 2000)
    impact = MultiCheckboxField(
      "Variant impact:",
      choices=[
        ("HIGH",
         "High"),
        ("MODERATE",
         "Moderate"),
        ("LOW",
         "Low"),
        ("MODIFIER",
         "Modifier"),
        ("geneonly",
         "Don't show variants"),
       ])
    ann_sev = HiddenField("")

    sample = MultiCheckboxField(
      "Sample:",
      choices = [])
    sample_group = MultiCheckboxField(
      "Sample group:",
      choices = [])
    types = MultiCheckboxField(
      "Type:",
      choices = types)
    csqs = MultiCheckboxField(
      "Consequence:",
      choices = csqs)

    mis = FloatField(
      "Missingness:")
    het = FloatField(
      "Heterozygosity:")
    maf = FloatField(
      "MAF:")

    skip_samp = MultiCheckboxField(
      "",
      choices=[
        ("·",
         "Reference-only samples"),
        ("—",
         "Absent-only samples")],
     default = [])
    skip_var = MultiCheckboxField(
      "",
      choices=[
        ("·",
         "Reference-only variants"),
        ("—",
         "Absent-only variants")],
      default = [])
    skip_no_var = BooleanField(
      "Variants without variation")

    exprun_group = MultiCheckboxField(
      "Run group:",
      choices = [])
    exp_col = MultiCheckboxField(
      "Column:",
      choices=[])

    #enr_col = RadioField(
    #  "Attribute:",
    #  choices=[])

class TraitViewForm(FlaskForm):
    trait_sgl = SelectField(
      "Trait:",
      choices=[])
    trait = MultiCheckboxField(
      "Trait:",
      choices=[])
    trait_group = MultiCheckboxField(
      "Trait group:",
      choices = [])
    sample = MultiCheckboxField(
      "Sample:",
      choices = [])
    sample_group = MultiCheckboxField(
      "Sample group:",
      choices = [])
    transform = RadioField(
      "Method:",
      choices=[
        ("", "(none)"),
        ("log", Markup("Log<sub><i>e</i></sub> transform")),
        ("bc", "Box-Cox transform"),
        ("pca", "PCA transform"),
        ("zs", "Z-score normalization"),
        ("corr", "Relatedness correction"),
        ("thr", "Threshold cutoff"),
        ("outl", "Outlier removal"),
        ("bin", "Make binary")],
      default="")
    k = BooleanField(
      "Kinship (K) correction",
      default = True)
    q = BooleanField(
      "PC (Q) correction",
      default = True)
    q_num = IntegerField(
      "Number of PCs:",
      default = 3)
    kq_from = RadioField(
      "K/Q from:",
      choices=[])
    cov = MultiCheckboxField(
      "Additional covariates:",
      choices = [])
    num_pc = IntegerField(
      "Number of PCs:",
      default=5)
    num_sd = IntegerField(
      "±SD from the mean >",
      default=3)
    binthr = FloatField(
      "Threshold:")
    thr_min = FloatField(
      "Minimum:")
    thr_max = FloatField(
      "Maximum:")
    sample_sel = StringField("")
    t_name = StringField("")

class PlotForm(FlaskForm):
    plot = SelectField(
      "Plot:",
      choices = [
        ("dens",
         "Marker density"),
        ("ld",
        "Linkage disequilibrium"),
        ("smiss",
         "Sample missingness"),
        ("vmiss",
         "Marker missingness"),
        ("shet",
         "Sample heterozygosity"),
        ("vhet",
         "Marker heterozygosity"),
        ("maf",
         "Minor allele frequency"),
        ("hwe",
         "-log10(p) HW equilibrium")
        ])
    plothap = SelectField(
      "Plot:",
      choices = [
        ("bdist",
         "Block distribution"),
        ("ld",
        "Linkage disequilibrium"),
        ("smiss",
         "Sample missingness"),
        ("vmiss",
         "Haplotype missingness"),
        ("shet",
         "Sample heterozygosity"),
        ("vhet",
         "Haplotype heterozygosity"),
        ("bsp",
         "Block size"),
        ("bmc",
         "Block marker count"),
        ("bhc",
         "Block total haplotype count"),
        ("bfhc",
         "Block MAF-filtered haplotype count"),
        ("baf",
         "Block allele frequency"),
        ("kin",
         "Kinship matrix"),
        ("pca",
         "Principal components")
        ])
    chr = SelectField(
      "Chrom:",
      choices=[])

    sample_group = MultiCheckboxField(
      "By sample group:",
      choices = [])
    by_phe = RadioField(
      "By phenotype:",
      choices = [])
    pc_nums = MultiCheckboxField(
      "PCs to plot:",
      choices = [(i, f"PC {i+1}") for i in range(15)],
      default = [0,1,2])

class AssocQCform(PlotForm):
    filt_from = RadioField(
      "From marker set:",
      choices=[])
    sample_group = MultiCheckboxField(
      "Limit to samples:",
      choices = [])
    coord = MultiCheckboxField(
      "Limit to regions:",
      choices = [])
    smiss = FloatField(
      "sample ≤")
    vmiss = FloatField(
      "marker ≤")
    shet = FloatField(
      "sample ≤")
    vhet = FloatField(
      "marker ≤")
    maf = FloatField(
      "MAF ≥")
    hwe = FloatField(
      Markup("-log<sub>10</sub>(p) HWE ≤"))
    prune = FloatField(
      "Threshold <")
    window = IntegerField(
      "Window (kb):")
    imp = BooleanField(
      "Impute")
    hap = BooleanField(
      "Make Haplotypes")
    btype = SelectField(
      "Method:",
      choices=[
        ("ld", "LD-based"),
        ("snpwin", "SNP number window"),
        ("bpwin", "Base-pair window")])
    bsize = IntegerField(
      "Max size (kb):",
      default=1000)
    bspan = IntegerField(
      "Block span:")
    bstep = IntegerField(
      "Block step:")
    name = StringField(
      "Set name:")

class AssocResForm(RegionForm, PlotForm):
    markers_sgl = RadioField(
      "Marker set:",
      choices = [])
    markers = MultiCheckboxField(
      "Marker set:",
      choices = [])
    model_sgl = RadioField(
      "Model:",
      choices = models)
    models = MultiCheckboxField(
      "Model:",
      choices = models)
    trait = MultiCheckboxField(
      "Trait:",
      choices = [])
    trait_group = MultiCheckboxField(
      "Trait group:",
      choices = [])
    mtry = IntegerField(
      "Test markers:")
    ntree = IntegerField(
      "Trees per forest:")

    b_burn = IntegerField(
      "Burn-in iterations:")
    b_chain = IntegerField(
      "Total iterations:")
    b_thin = IntegerField(
      "Thinning:")

    hgl_coord_group = MultiCheckboxField(
      "Saved regions:",
      choices=[])
    hgl_gene_group = MultiCheckboxField(
      "Gene groups:",
      choices = [])

    k = BooleanField(
      "Kinship (K) correction",
      default = True)
    q = BooleanField(
      "PC (Q) correction",
      default = True)
    q_num = IntegerField(
      "Number of PCs:",
      default = 3)
    kq_from = RadioField(
      "K/Q from:",
      choices=[])
    cov = MultiCheckboxField(
      "Additional covariates:",
      choices = [])
    overwrite = BooleanField(
      "Overwrite existing")

    plot_thr = FloatField(
      "Plot threshold:",
      default=8)
    window_ld = IntegerField(
      "Window (kb):",
      default=5000)
    r2_ld = FloatField(
      "r<sup>2</sup> to peak ≥",
      default=0.2)
    qtl_runs = StringField()