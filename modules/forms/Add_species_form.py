#!/usr/bin/env python

from flask_wtf import FlaskForm
from wtforms import SelectField, SubmitField, StringField, FileField, IntegerField#, Form, FieldList, BooleanField, HiddenField, FormField
from wtforms.validators import InputRequired#, DataRequired


class SpeciesForm(FlaskForm):
    taxid = IntegerField("Taxon ID", validators=[InputRequired()])
    disp_name = StringField("Display name", render_kw={"disabled": True})
    search = SubmitField("Search")
    save = SubmitField("Save")
    image = FileField("Upload image")

    def enable_disp_name(self, disp_name):
        self.disp_name.render_kw = {"disabled": False}
        self.disp_name.data = disp_name
        self.disp_name.validators = [InputRequired()]

    def disable_disp_name(self):
        self.disp_name.render_kw = {"disabled": True}
        self.disp_name.data = ""
        self.disp_name.validators = []


class AssemblyForm(FlaskForm):
    assemblies = SelectField("Select assembly:", choices=[], validators=[InputRequired()])
    select = SubmitField("Select")


if __name__ == "__main__":
    pass