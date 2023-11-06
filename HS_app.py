#!/usr/bin/env python

import os
import sys
import json
import redis
import shutil
import secrets
import subprocess

from Bio import Entrez
from datetime import date
from zipfile import ZipFile

from flask import (Flask,
                   render_template,
                   request, redirect,
                   url_for,
                   send_from_directory)

from flask_wtf import FlaskForm
from wtforms import SelectField, SubmitField, StringField, SelectMultipleField#, Form, FieldList, BooleanField, HiddenField, FormField
from wtforms.validators import InputRequired#, DataRequired

from modules.forms.Add_species_form import SpeciesForm, AssemblyForm

from ncbi.datasets.package import dataset
from ncbi.datasets import GenomeApi as DatasetsGenomeApi
from ncbi.datasets.openapi import ApiClient as DatasetsApiClient
from ncbi.datasets.openapi import ApiException as DatasetsApiException
from ncbi.datasets.metadata.genome import (get_assembly_metadata_by_taxon,
                                           get_assembly_metadata_by_asm_accessions)


##### IMPLEMENT LOADING SCREEN: https://stackoverflow.com/questions/14525029/display-a-loading-message-while-a-time-consuming-function-is-executed-in-flask

    #allowed_extensions = ['.jpg', '.jpeg', '.png', '.gif']

app = Flask(__name__)

app.secret_key = secrets.token_hex(16)
app.config["UPLOAD_FOLDER"] = "static/icons"

Entrez.email = "peter.polacek@wur.nl"


#redis_client = redis.StrictRedis(host="localhost", port=6379, db=0)
##### REDIS HSET HGET (hash, key, value) and JSON.DUMPS/LOADS


@app.route("/image/<path:filename>")
def serve_image(filename):
    image_directory = os.getcwd()
    return send_from_directory(image_directory, filename)

##### IMAGE SANITIZATION        https://blog.miguelgrinberg.com/post/handling-file-uploads-with-flask
##### IMAGE CONVERSION          https://www.google.com/search?client=firefox-b-d&q=python+flask+convert+uploaded+image
##### MONGODB ANITIZATION       https://stackoverflow.com/questions/43925397/what-is-the-best-way-to-sanitize-inputs-with-flask-and-when-using-mongodb
#####                           https://www.google.com/search?client=firefox-b-d&q=flask+input+sanitization
##### FLASK OFFLOAD             https://www.google.com/search?client=firefox-b-d&q=flask+offload
##### POPEN                     https://stackoverflow.com/questions/7681715/whats-the-difference-between-subprocess-popen-and-call-how-can-i-use-them

##### ENTREZ EFETCH FOR BIOSAMPLE, LINK TO BIOPROJECT VIA ELINK

##### IS THERE ANYTHING IN FLASK TO DO ABOUT HTTPS?

@app.route("/")
def index():
    user_id = secrets.token_urlsafe(16)
    genomes = list()

    for taxid in sorted(os.listdir("data/")):
        icon = (f"static/icons/{taxid}.png"
                if os.path.exists(f"static/icons/{taxid}.png")
                else "static/icons/Icon_default.png")
        with open(f"data/{taxid}/metadata.json", "r") as json_file:
            metadata = json.load(json_file) 

        genomes.append({"taxid":      taxid,
                        "icon":       icon,
                        "metadata":   metadata})

    return render_template("index.html",
                           user_id=user_id,
                           genomes=genomes)


@app.route("/<user_id>/add_species", methods=["GET", "POST"])
def add_species(user_id):
    form = SpeciesForm()
    taxid = ""
    disp_name = ""
    image_path = "static/icons/Icon_default.png"
    image_alt = "Default icon"

    if form.validate_on_submit():
        taxid = request.form.get("taxid")

        if "search" in request.form:
            try:
                handle = Entrez.efetch(db="taxonomy", id=taxid, retmode="xml")
                taxon_record = Entrez.read(handle)[0]
                form.enable_disp_name(taxon_record["ScientificName"])
                handle.close()

            except Exception as e:
                form.disable_disp_name()

        elif "save" in request.form:
            os.mkdir(f"data/{taxid}/")
            metadata = {"disp_name":  request.form.get("disp_name"),
                        "created_on": str(date.today()),
                        "created_by": Entrez.email
                       }
            with open(f"data/{taxid}/metadata.json", "w") as out:
                out.write(json.dumps(metadata,
                                     indent=4)
                          )
            if "image" in request.files:
                image = request.files["image"]

                if image.filename != "":
                    image_path = f"{app.config['UPLOAD_FOLDER']}/{taxid}.png"
                    image.save(image_path)

            return redirect(url_for("index"))

    return render_template("add_species.html",
                           user_id=user_id,
                           taxid=taxid,
                           form=form,
                           image_path=image_path,
                           image_alt=image_alt)


##### FOR SELECT AND ADD ASSEMBLY
##### ADD SUMMARY INFO, SWITCHABLE WITH FULL INFO
##### SUMMARY INFO SHOWS IF ANNOTATED, IF REFSEQ, ORGANISM NAME, OTHER INFO IF USEFUL

@app.route("/<user_id>/<taxid>", methods=["GET", "POST"])
def select_assembly(user_id, taxid):
    form = AssemblyForm()

    if request.method == "GET":
        assembly_metadata = dict()
        assembly_accessions = sorted([path for path in os.listdir(f"data/{taxid}")
                                      if os.path.isdir(f"data/{taxid}/{path}")])

        for assembly in get_assembly_metadata_by_asm_accessions(assembly_accessions):
            assembly = assembly["assembly"]
            assembly_id = assembly["assembly_accession"]
            assembly_name = f"{assembly['display_name']} ({assembly_id})"

            assembly_metadata[assembly_id] = assembly.to_dict()
            form.assemblies.choices.append((assembly_id, assembly_name))

        return render_template("add_select_assembly.html",
                               user_id=user_id,
                               taxid=taxid,
                               form=form,
                               assembly_metadata=assembly_metadata,
                               add_new=False)

    elif request.method == "POST":
        assembly_id = request.form.get("assemblies")

        return redirect(url_for("select_datasets",
                                user_id=user_id,
                                taxid=taxid,
                                assembly_id=assembly_id))


@app.route("/<user_id>/<taxid>/add_assembly", methods=["GET", "POST"])
def add_assembly(user_id, taxid):

    if request.method == "GET":

        form = AssemblyForm()

        assembly_metadata = dict()
        for assembly in get_assembly_metadata_by_taxon(taxid):
            assembly = assembly["assembly"]
            assembly_id = assembly["assembly_accession"]
            assembly_name = f"{assembly['display_name']} ({assembly_id})"

            assembly_metadata[assembly_id] = assembly.to_dict()
            form.assemblies.choices.append((assembly_id, assembly_name))

        return render_template("add_select_assembly.html",
                               user_id=user_id,
                               taxid=taxid,
                               form=form,
                               assembly_metadata=assembly_metadata,
                               add_new=True)

    elif request.method == "POST":

        ### USE DATASETS COMMAND LINE TOOL IF TOTAL DATA SIZE > 15GB

        assembly_filetypes =    ["CDS_FASTA", "PROT_FASTA", "GENOME_GFF", "GENOME_FASTA"]
        downloaded_filetypes =  ["CDS_NUCLEOTIDE_FASTA", "PROTEIN_FASTA", "GFF3", "GENOMIC_NUCLEOTIDE_FASTA"]
        assembly_filenames =    ["cds.fa", "pep.fa", "ann.gff", "seq.fa"]

        assembly_id = request.form.get("assemblies")
        assembly_dir = f"data/{taxid}/{assembly_id}"
        os.mkdir(assembly_dir)

        zip_fn = f"{assembly_dir}/dataset.zip"

        with DatasetsApiClient() as api_client:
            genome_api = DatasetsGenomeApi(api_client)
            try:
                genome_download = genome_api.download_assembly_package(
                    [assembly_id],
                    include_annotation_type=assembly_filetypes,
                    _preload_content=False
                    #hydrated="DATA_REPORT_ONLY"
                    )

                with open(zip_fn, "wb") as f:
                    f.write(genome_download.data)

            except DatasetsApiException as e:
                sys.exit(f"Exception when calling download_assembly_package: {e}\n")

        package = dataset.AssemblyDataset(zip_fn)
        root_dir = package.get_file_root_dir() # ncbi_dataset/data

        with ZipFile(zip_fn, "r") as archive:
            for (downloaded_filetype, assembly_filename) in zip(downloaded_filetypes, assembly_filenames):
                try:
                    downloaded_filename = package.get_file_names_by_type(downloaded_filetype)[0]
                except IndexError:
                    continue

                downloaded_filename = f"{root_dir}/{downloaded_filename}"
                archive.getinfo(downloaded_filename).filename = f"{assembly_dir}/{assembly_filename}"
                archive.extract(downloaded_filename)

        os.remove(zip_fn)

        assembly = next(get_assembly_metadata_by_asm_accessions([assembly_id]))
        assembly = assembly["assembly"]
        assembly_name = assembly["display_name"]
        with open(f"{assembly_dir}/metadata.json", "w") as out:
            out.write(json.dumps({"display_name":  assembly_name,
                                   "metadata":     assembly.to_dict()},
                                 indent=4))

        return redirect(url_for("select_assembly",
                                user_id=user_id,
                                taxid=taxid))


@app.route("/<user_id>/<taxid>/<assembly_id>/select_datasets", methods=["GET", "POST"])
def select_datasets(user_id, taxid, assembly_id):

    form = AssemblyForm()
    if request.method == "GET":

        bp_ids = Entrez.read(Entrez.esearch(db="bioproject",
                                            term=f"txid{taxid}[Organism:exp]",
                                            retmode="xml",
                                            retmax=20))["IdList"]

        bp_data = Entrez.read(Entrez.esummary(db="bioproject",
                                              id=bp_ids,
                                              retmode="xml"))
        bp_data = bp_data["DocumentSummarySet"]["DocumentSummary"]

        dataset_metadata = dict()

        for bp_id, data in zip(bp_ids, bp_data):
            bp_title = f"{data['Project_Title']} ({bp_id})"
            data = {key: value for key, value in data.items() if value}
            dataset_metadata[bp_id] = data
            form.assemblies.choices.append((bp_id, bp_title))

        return render_template("select_datasets.html",
                               user_id=user_id,
                               taxid=taxid,
                               form=form,
                               assembly_metadata=dataset_metadata,
                               add_new=False)

        """
        for bp_id in bp_ids:

            bs_ids = [bs_link["Id"] 
                      for bs_link in Entrez.read(Entrez.elink(db="biosample",
                                                              dbfrom="bioproject",
                                                              id=bp_id,
                                                              retmode="xml"))[0]["LinkSetDb"][0]["Link"]]

            bs_data = Entrez.read(Entrez.esummary(db="biosample",
                                                  id=bs_ids,
                                                  retmode="xml"))
            bs_data = bp_data["DocumentSummarySet"]["DocumentSummary"]
        """
        form = AssemblyForm()
        return render_template("select_datasets.html",
                               user_id=user_id,
                               taxid=taxid,
                               assembly_id=assembly_id,
                               form=form,
                           add_new=False)

    elif request.method == "POST":
        bp_id = request.form.get("assemblies")
        bp_id = "123456"
        bs_id = "35882832"

        #bs_data = Entrez.read(Entrez.esummary(db="biosample",
        #                                      id="35882832",
        #                                      retmode="xml"))
        #bs_data = bs_data["DocumentSummarySet"]["DocumentSummary"][0]
        sra_id = "SRS18057999"

        FETCHNGS_PATH = "nf_core/nf-core-fetchngs_1.11.0/1_11_0/"
        OUT_DIR = f"data/{taxid}/{assembly_id}/{bp_id}"
        PARAMS_FILE = f"{OUT_DIR}/{sra_id}.csv"

        # Start the Sarek pipeline as a subprocess
        process = subprocess.Popen(['nextflow', 'run', FETCHNGS_PATH, "--input", PARAMS_FILE, "--outdir", OUT_DIR], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        stdout, stderr = process.communicate()
        return_code = process.wait()
        if return_code == 0:
            print("Workflow completed successfully.")
        print(stdout)
        print(stderr)

        return #redirect(url_for('index'))


if __name__ == "__main__":
    app.run(debug=True, port=8888)
