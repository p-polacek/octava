#!/usr/bin/env python

from flask_wtf import FlaskForm
from wtforms import (widgets,
                     FileField,
                     FloatField,
                     RadioField,
                     HiddenField,
                     SelectField,
                     StringField,
                     SubmitField,
                     BooleanField,
                     DecimalField,
                     IntegerField,
                     SelectMultipleField,
                     #, Form, FieldList, BooleanField, HiddenField, FormField
                    )
from wtforms.validators import (InputRequired,
                                DataRequired,
                                ValidationError,
                                Email)

from modules.getQueryOptions import getVariantConsequences

models = ["BLINK", "CMLM", "FarmCPU",
          "GLM", "MMLM", "MLM", "SUPER"]

types_ch = [("snp", "Single nucleotide variant"),
            ("mnp", "Multi-nucleotide variant"),
            ("ins", "Insertion"),
            ("del", "Deletion"),
            ("dup", "Duplication"),
            ("inv", "Inversion"),
            ("cnv", "Copy number variant"),
            ("bnd", "Breakend variant")]
            
csqs_ch = getVariantConsequences()


class MultiCheckboxField(SelectMultipleField):
    widget = widgets.ListWidget(prefix_label=False)
    option_widget = widgets.CheckboxInput()


class EmailForm(FlaskForm):
    email = StringField("Email", validators=[DataRequired(), Email()])
    submit = SubmitField("Login")


class HaploSmasherForm(FlaskForm):
    fmt             = SelectField("Format:",
                                  choices=[("region", "Region"),
                                           ("gene",   "Gene ID")])
    target          = SelectField("Target:",
                                  choices=[("cds",    "Coding regions"),
                                           ("prom",   "Promoter regions"),
                                           ("full",   "Entire region")])
    prom_size       = IntegerField("Size:", default=2000)

    gene            = StringField("Gene:")
    gene_groups     = MultiCheckboxField("Gene group:",  choices=[])
    chrom           = SelectField("Chrom:", choices=[])
    start           = IntegerField("Start:")
    end             = IntegerField("End:")

    impact          = MultiCheckboxField("Variant impact:",
                                         choices=[("HIGH", "High"),
                                                  ("MODERATE", "Moderate"),
                                                  ("LOW", "Low"),
                                                  ("MODIFIER", "Modifier"),
                                                  ("geneonly", "Don't show variants")])

    samples         = MultiCheckboxField("Sample:",        choices=[])
    sample_groups   = MultiCheckboxField("Sample group:",  choices=[])
    types           = MultiCheckboxField("Type:",          choices=types_ch)
    csqs            = MultiCheckboxField("Consequence:",   choices=csqs_ch)

    skip_samples    = MultiCheckboxField("",
                                         choices=[("·", "Reference-only samples"),
                                                  ("—", "Absent-only samples")],
                                         default=["·", "—"])
    skip_variants   = MultiCheckboxField("",
                                         choices=[("·", "Reference-only variants"),
                                                  ("—", "Absent-only variants")],
                                         default=["·", "—"])
    skip_no_var     = BooleanField("Variants without variation")

    exprun_groups   = MultiCheckboxField("Run group:",  choices=[])


class GWASQCform(FlaskForm):
    plot   = SelectField("Plot:",
                         choices=[("smiss", "Sample missingness"),
                                  ("vmiss", "Marker missingness"),
                                  ("shet",  "Sample heterozygosity"),
                                  ("vhet",  "Marker heterozygosity"),
                                  ("maf",   "Minor allele frequency"),
                                  ("hwe",   "-log10(P) HW equilibrium")])

    smiss  = DecimalField("sample ≤", places=None, default=0.1)
    vmiss  = DecimalField("marker ≤", places=None, default=0.1)
    shet   = DecimalField("sample ≤", places=None, default=0.1)
    vhet   = DecimalField("marker ≤", places=None, default=0.1)
    maf    = DecimalField("Minor AF ≥", places=None, default=0.05)
    hwe    = DecimalField("-log10(P) HWE ≤", places=None)
    prune  = DecimalField("LD ≤", places=None, default=0.95)
    window = IntegerField("Window (kb):", default=200)

    name   = StringField("Set name:")


class GWASForm(FlaskForm):
    markers = MultiCheckboxField("Marker set:", choices=[])
    models  = MultiCheckboxField("Model:", 
                                 choices=models)
    traits  = MultiCheckboxField("Trait:",  choices=[])

    trait_groups  = MultiCheckboxField("Trait group:",  choices=[])

    chrom   = SelectField("Chrom:", choices=[])
    start   = IntegerField("Start:")
    end     = IntegerField("End:")

    pval    = HiddenField("")
    fm_tool = RadioField("Fine-map tool:", choices=[("", "(none)"),
                                                    ("susie", "SuSiE")])

    submit = SubmitField("Submit")


if __name__ == "__main__":
    pass