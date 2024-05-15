#!/usr/bin/env python

from flask_wtf import FlaskForm

from wtforms import (widgets,
                     RadioField,
                     HiddenField,
                     SelectField,
                     StringField,
                     BooleanField,
                     DecimalField,
                     IntegerField,
                     MultipleFileField,
                     SelectMultipleField,
                    )

from modules.csq import var_csq


models = [
  "BLINK",
  "CMLM",
  "FarmCPU",
  "GLM",
  "MMLM",
  "MLM",
  "SUPER"]

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

class AlleleMinerForm(RegionForm):
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
    prom_size = IntegerField(
      "Size:",
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

    samples = MultiCheckboxField(
      "Sample:",
      choices = [])
    sample_groups = MultiCheckboxField(
      "Sample group:",
      choices = [])
    types = MultiCheckboxField(
      "Type:",
      choices = types)
    csqs = MultiCheckboxField(
      "Consequence:",
      choices = csqs)

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
         "-log10(P) HW equilibrium")])
    ld_chr = SelectField(
      "Chrom:",
      choices=[])

class AssocQCform(PlotForm):
    smiss = DecimalField(
      "sample ≤",
      places  = None,
      default = 0.1)
    vmiss = DecimalField(
      "marker ≤",
      places  = None,
      default = 0.1)
    shet = DecimalField(
      "sample ≤",
      places  = None,
      default = 0.1)
    vhet = DecimalField(
      "marker ≤",
      places  = None,
      default = 0.1)
    maf = DecimalField(
      "Minor AF ≥",
      places  = None,
      default = 0.05)
    hwe = DecimalField(
      "-log10(P) HWE ≤",
      places  = None)
    prune = DecimalField(
      "LD ≤",
      places  = None,
      default = 0.95)
    window = IntegerField(
      "Window (kb):",
      default = 200)
    imp = BooleanField(
      "Impute",
      default="checked")
    hap = BooleanField(
      "Make Haplotypes")
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
    traits = MultiCheckboxField(
      "Trait:",
      choices = [])
    trait_groups = MultiCheckboxField(
      "Trait group:",
      choices = [])
    new_phe = MultipleFileField()
    k = BooleanField(
      "Kinship (K) correction",
      default = True)
    q = BooleanField(
      "PCA (Q) correction",
      default = True)
    q_num = IntegerField(
      "Number of PCs:",
      default = 3)
    kq_from = RadioField(
      "K/Q from:",
      choices=[])
    mrk  = StringField(
      "")
    mod  = StringField(
      "")
    trt  = StringField(
      "")
    pval = HiddenField(
      "")
