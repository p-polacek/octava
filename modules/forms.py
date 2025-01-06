#!/usr/bin/env python

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
  DecimalField,
  IntegerField,
  SelectMultipleField)

from modules.csq import var_csq

models = [
  ("glm", "GLM"),
  ("mlm", "MLM"),
  ("mtmlm", "Multi-trait MLM"),
  ("blink", "BLINK"),
  ("rf", "Random forest")]

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
    name = SelectField(
      "Name:",
      choices=[])
    icon = FileField(
      "Icon:")
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
    proj = SelectField(
      "Study info:",
      choices=[])
    samp = SelectField(
      "Sample info:",
      choices=[])
    info = SelectField(
      "Gene info:",
      choices=[])
    enr = MultiCheckboxField(
      "Enrichment for:",
      choices=[])
    bam = MultiCheckboxField(
      "BAM file(s):",
      choices=[])
    snpeff = BooleanField(
      "Variant effect:")
    cds = SelectField(
      "CDS:",
      choices=[])
    pep = SelectField(
      "Proteins:",
      choices=[])
    sev = BooleanField(
      "Variant severity:")
    blastp = SelectField(
      "Blast results:",
      choices=[])
    concat = MultiCheckboxField(
      "Concatenate contigs:",
      choices=[])
    edit = BooleanField(
      "Enable editing dataset later:")
    igv_vis = IntegerField(
      "Visibility window:",
      default=1000000)
    igv_h = IntegerField(
      "Track height:",
      default=150)

class EditForm(FlaskForm):
    icon = FileField(
      "Icon:")
    phe = FileField(
      "Phenotype:")
    run = FileField(
      "Run info:")
    proj = FileField(
      "Project info:")
    samp = FileField(
      "Sample names:")

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

class FuncAnalysisForm(RegionForm):
    dataset = SelectField(
      "Markers:",
      choices = [("", "raw VCF")])
    target = SelectField(
      "Target:",
      choices = [
        ("cds",
         "Coding regions"),
        ("prom",
         "Promoter regions"),
        ("full",
         "Entire region")])
    gene = StringField(
      "Gene:")
    gene_groups = MultiCheckboxField(
      "Gene group:",
      choices = [])
    coord_groups = MultiCheckboxField(
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
#        ("SEV",
#         "Only deleterious variants")
       ])
    ann_sev = HiddenField("")

    samples = MultiCheckboxField(
      "Sample:",
      choices = [])
    samples_group = MultiCheckboxField(
      "Sample group:",
      choices = [])
    types = MultiCheckboxField(
      "Type:",
      choices = types)
    csqs = MultiCheckboxField(
      "Consequence:",
      choices = csqs)

    mis = DecimalField(
      "Missingness:")
    het = DecimalField(
      "Heterozygosity:")
    maf = DecimalField(
      "MAF:")

    skip_samples = MultiCheckboxField(
      "",
      choices=[
        ("·",
         "Reference-only samples"),
        ("—",
         "Absent-only samples")],
     default = [])
    skip_variants = MultiCheckboxField(
      "",
      choices=[
        ("·",
         "Reference-only variants"),
        ("—",
         "Absent-only variants")],
      default = [])
    skip_no_var = BooleanField(
      "Variants without variation")

    exprun_groups = MultiCheckboxField(
      "Run group:",
      choices = [])
    exp_col = RadioField(
      "Column:",
      choices=[])

    enr_col = RadioField(
      "Attribute:",
      choices=[])

class TraitViewForm(FlaskForm):
    trait_sgl = SelectField(
      "Trait:",
      choices=[])
    traits = MultiCheckboxField(
      "Trait:",
      choices=[])
    traits_group = MultiCheckboxField(
      "Trait group:",
      choices = [])
    samples = MultiCheckboxField(
      "Sample:",
      choices = [])
    samples_group = MultiCheckboxField(
      "Sample group:",
      choices = [])
    transform = RadioField(
      "Method:",
      choices=[
        ("log", "Log(e) transform"),
        #("bc", "Box-cox transform"),
        ("pca", "PCA transform"),
        ("bin", "Make binary")])
    num_pc = IntegerField(
      "Number of PCs:",
      default=5)
    binthr = IntegerField(
      "Threshold:")
    samples_sel = StringField("")
    t_name = StringField("")
    t_type = StringField("")

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
         "-log10(P) HW equilibrium")
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

class AssocQCform(PlotForm):
    filt_ext = BooleanField(
      "Apply to existing:")
    samples_group = MultiCheckboxField(
      "By sample group:",
      choices = [])
    filt_from = RadioField(
      "Markers:",
      choices=[])
    smiss = DecimalField(
      "sample ≤",
      places  = None,
      #default = 0.1
      )
    vmiss = DecimalField(
      "marker ≤",
      places  = None,
      #default = 0.1
      )
    shet = DecimalField(
      "sample ≤",
      places  = None,
      #default = 0.1
      )
    vhet = DecimalField(
      "marker ≤",
      places  = None,
      #default = 0.1
      )
    maf = DecimalField(
      "Minor AF ≥",
      places  = None,
      #default = 0.05
      )
    hwe = DecimalField(
      "-log10(P) HWE ≤",
      places  = None)
    prune = DecimalField(
      "Threshold <",
      places  = None,
      #default = 0.95
      )
    window = IntegerField(
      "Window (kb):",
      #default = 200
      )
    imp = BooleanField(
      "Impute",
      default="checked")
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
    samples_group = MultiCheckboxField(
      "By sample group:",
      choices = [])
    traits = MultiCheckboxField(
      "Trait:",
      choices = [])
    traits_group = MultiCheckboxField(
      "Trait group:",
      choices = [])
    new_phe = FileField()
    mtryvals = StringField(
      "Test markers:")
    ntreevals = StringField(
      "Trees per forest:")
    optimize = BooleanField()
    mtry = IntegerField(
      "Test markers:")
    ntree = IntegerField(
      "Trees per forest:",
      default = 100)
    niter = IntegerField(
      "Iterations:",
      default = 20)
    nperm = IntegerField(
      "Permutations:",
      default = 100)

    pm_w = IntegerField(
      "Burn-in steps:")
    pm_s = IntegerField(
      "Sampling steps:")
    pm_num = IntegerField(
      "Thinning:")
    pm_hmin = FloatField(
      "Min heritability:")
    pm_hmax = FloatField(
      "Max heritability:")
    pm_pmin = IntegerField(
      "Min causal SNPs:")
    pm_pmax = IntegerField(
      "Max causal SNPs:")
    pm_nstart = IntegerField(
      "Init causal SNPs:")
    pm_smin = IntegerField(
      "Model min SNPs:")
    pm_smax = IntegerField(
      "Model max SNPs:")

    fm_ids = MultiCheckboxField(
      "",
      choices=[])
    fm_runs = MultiCheckboxField(
      "Saved analyses:",
      choices=[])

    assoc_loci = MultiCheckboxField(
      "Saved loci:",
      choices=[])

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
    overwrite = BooleanField(
      "Overwrite existing")
    fm_method = RadioField(
      "Method:",
      choices=[
        ("", "(none)"),
        ("pm", "PiMASS"),
        ("rf", "Random forest")])
    mrk  = StringField(
      "")
    mod  = StringField(
      "")
    trt  = StringField(
      "")
    pval = HiddenField(
      "")
