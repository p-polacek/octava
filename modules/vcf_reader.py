import os
import sys
import gzip
import re
import warnings
from collections import OrderedDict
from stat import S_ISREG

try:
    import pysam
except ImportError:
    pysam = None

"""
VCF parser created by David A. Parry
Obtained from: https://github.com/david-a-parry/parse_vcf.py

COPYRIGHT AND LICENSE

MIT License

Copyright (c) 2018 David A. Parry

Permission is hereby granted, free of charge, to any person obtaining
a copy of this software and associated documentation files (the "Software"),
to deal in the Software without restriction, including without limitation
the rights to use, copy, modify, merge, publish, distribute, sublicense,
and/or sell copies of the Software, and to permit persons to whom
the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included
in all copies or substantial portions of the Software.


Modified by Peter Poláček
for the use with HaploSmasher v2.0

Modifications:

Edited handling of "ANN" field and the "CSQ" property (lines 393-395)
Added "VAR_TYPE" property to the "Record" object,
edited from "var_type" property of "AltAllele" object (lines 691-717)
"""

# common INFO fields and their types in case absent from header
COMMON_INFO = {
    "1000G": {"Type": "Flag", "Class": None, "Number": 1, "Split": False},
    "AA": {"Type": "String", "Class": str, "Number": 1, "Split": False},
    "AC": {"Type": "Integer", "Class": int, "Number": "A", "Split": True},
    "AF": {"Type": "Float", "Class": float, "Number": "A", "Split": True},
    "AN": {"Type": "Integer", "Class": int, "Number": 1, "Split": False},
    "BQ": {"Type": "Float", "Class": float, "Number": 1, "Split": False},
    "CIGAR": {"Type": "String", "Class": str, "Number": 1, "Split": False},
    "DB": {"Type": "Flag", "Class": None, "Number": 1, "Split": False},
    "DP": {"Type": "Integer", "Class": int, "Number": 1, "Split": False},
    "END": {"Type": "Integer", "Class": int, "Number": 1, "Split": False},
    "H2": {"Type": "Flag", "Class": None, "Number": 1, "Split": False},
    "H3": {"Type": "Flag", "Class": None, "Number": 1, "Split": False},
    "MQ": {"Type": "Float", "Class": float, "Number": 1, "Split": False},
    "MQ0": {"Type": "Integer", "Class": int, "Number": 1, "Split": False},
    "NS": {"Type": "Integer", "Class": int, "Number": 1, "Split": False},
    "SB": {"Type": "String", "Class": str, "Number": 1, "Split": False},
    "SOMATIC": {"Type": "Flag", "Class": None, "Number": 1, "Split": False},
    "VALIDATED": {"Type": "Flag", "Class": None, "Number": 1, "Split": False},
    # for SVs
    "BKPTID": {"Type": "String", "Class": str, "Number": ".", "Split": False},
    "CICN": {"Type": "Integer", "Class": int, "Number": 2, "Split": True},
    "CICNADJ": {"Type": "Integer", "Class": int, "Number": 1, "Split": False},
    "CIEND": {"Type": "Integer", "Class": int, "Number": 2, "Split": True},
    "CILEN": {"Type": "Integer", "Class": int, "Number": 2, "Split": True},
    "CIPOS": {"Type": "Integer", "Class": int, "Number": 2, "Split": True},
    "CN": {"Type": "Integer", "Class": int, "Number": 1, "Split": False},
    "CNADJ": {"Type": "Integer", "Class": int, "Number": ".", "Split": False},
    "DBRIPID": {"Type": "String", "Class": str, "Number": 1, "Split": False},
    "DBVARID": {"Type": "String", "Class": str, "Number": 1, "Split": False},
    "DGVID": {"Type": "String", "Class": str, "Number": 1, "Split": False},
    "DPADJ": {"Type": "Integer", "Class": int, "Number": 1, "Split": False},
    "EVENT": {"Type": "String", "Class": str, "Number": 1, "Split": False},
    "HOMLEN": {"Type": "Integer", "Class": int, "Number": ".", "Split": False},
    "HOMSEQ": {"Type": "String", "Class": str, "Number": ".", "Split": False},
    "IMPRECISE": {"Type": "Flag", "Class": None, "Number": 1, "Split": False},
    "MATEID": {"Type": "String", "Class": str, "Number": 1, "Split": False},
    "MEINFO": {"Type": "String", "Class": str, "Number": 4, "Split": True},
    "METRANS": {"Type": "String", "Class": str, "Number": 4, "Split": True},
    "NOVEL": {"Type": "Flag", "Class": None, "Number": 1, "Split": False},
    "PARID": {"Type": "String", "Class": str, "Number": 1, "Split": False},
    "SVLEN": {"Type": "Integer", "Class": int, "Number": 1, "Split": False},
    "SVTYPE": {"Type": "String", "Class": str, "Number": 1, "Split": False},
}

# common genotype FORMAT Fields and their types in case absent from header
COMMON_FORMAT = {
    "DP": {"Type": "Integer", "Class": int, "Number": 1, "Split": False},
    "EC": {"Type": "Integer", "Class": int, "Number": "A", "Split": True},
    "FT": {"Type": "String", "Class": str, "Number": 1, "Split": False},
    "GL": {"Type": "Float", "Class": float, "Number": "G", "Split": True},
    "GLE": {"Type": "String", "Class": str, "Number": ".", "Split": False},
    "GP": {"Type": "Float", "Class": float, "Number": "G", "Split": True},
    "GQ": {"Type": "Integer", "Class": int, "Number": 1, "Split": False},
    "GT": {"Type": "String", "Class": str, "Number": 1, "Split": False},
    "HQ": {"Type": "Integer", "Class": int, "Number": 2, "Split": True},
    "MQ": {"Type": "Integer", "Class": int, "Number": 1, "Split": False},
    "PL": {"Type": "Integer", "Class": int, "Number": 1, "Split": False},
    "PQ": {"Type": "Integer", "Class": int, "Number": 1, "Split": False},
    "PS": {"Type": "Integer", "Class": int, "Number": 1, "Split": False},
    # for SVs
    "CN": {"Type": "Integer", "Class": int, "Number": 1, "Split": False},
    "CNQ": {"Type": "Float", "Class": float, "Number": 1, "Split": False},
    "CNL": {"Type": "Float", "Class": float, "Number": ".", "Split": False},
    "NQ": {"Type": "Integer", "Class": int, "Number": 1, "Split": False},
    "HAP": {"Type": "Integer", "Class": int, "Number": 1, "Split": False},
    "AHAP": {"Type": "Integer", "Class": int, "Number": 1, "Split": False},
}
_svalt_re = re.compile(r"""<(\w+)(:\w+)*>""")  # group 1 gives SV type
_bnd_re = re.compile(r"""^(([ACTGN]*)[\[\]]\w+):\d+[\]\[]([ACGTN]*)$""")
# group 1 gives VEP CSQ allele


class VcfReader(object):
    """
    A class for parsing Variant Call Format (VCF) files. Stores
    header information as a VcfHeader object in self.header. Is an
    iterable variant parser, returning VcfRecord objects.
    If your input is compressed and indexed, you may iterate over
    genomic regions instead of processing through the file line by
    line.
    Examples:
        #process line by line
        >>> v = VcfReader('path/to/my/input.vcf')
        >>> for record in v
        ... #do something with each record
        #process variants overlapping a given region
        >>> v = VcfReader('path/to/my/indexed.vcf.gz')
        >>> v.set_region(chrom='chr1', start=899999 end=1000000)
        >>> for record in v:
        ... #do something with each record
    """

    def __init__(self, filename, compressed=None, bcf=None, encoding="utf-8"):
        """
        Create a new VcfReader object
        Opens the given VCF file and stores the metaheader and
        sample information
        Args:
            filename:   VCF file to open. Required.
            compressed: Boolean indicating whether input should be treated
                        as bgzip compressed data. If not provided, this
                        will be inferred from the file extension.
            bcf:        Boolean indicating whether the input is from a BCF
                        file. If not provided, this will be inferred from
                        the file extension.
            encoding:   Encoding of input data. Default = 'utf-8'.
        """

        self.filename = filename
        self.compressed = compressed
        self.bcf = bcf
        self.encoding = encoding
        self._tabix = None
        if self.bcf is None:
            self.bcf = filename.endswith(".bcf")
        if self.compressed is None:
            self.compressed = filename.endswith((".gz", ".bgz"))
        if self.bcf:
            if pysam is None:
                raise ImportError(
                    "pysam not available. Please install (e.g. "
                    + "via 'pip install pysam' to parse bcf files"
                )
            self.file = pysam.VariantFile(filename)
            self.reader = (
                rec.__str__().rstrip()
                for rec in self.file.fetch()
                if rec.__str__().rstrip()
            )
            head = self.file.header.__str__().split("\n")
            head.pop()
            cols = head.pop().split("\t")
            self.header = VcfHeader(head, cols)
        else:
            if self.compressed:
                if filename == "-":
                    gz = sys.stdin
                else:
                    gz = filename
                # self.file = gzip.open(gz, encoding=encoding, mode='rt')
                self.file = gzip.open(
                    gz, encoding=encoding, errors="replace", mode="rt"
                )
            else:
                if filename == "-":
                    self.file = sys.stdin
                else:
                    self.file = open(
                        filename, encoding=encoding, errors="replace", mode="r"
                    )
            self.reader = (line.rstrip() for line in self.file if line.rstrip())
            self.header = self._read_header()
        if filename == "-":
            self._is_reg_file = False
        else:
            self._is_reg_file = S_ISREG(os.stat(self.filename).st_mode)
        # read header information
        # set some header values for convenience
        self.metadata = self.header.metadata
        self.col_header = self.header.col_header
        self.parser = (
            VcfRecord(
                line,
                self,
            )
            for line in self.reader
        )

    def __iter__(self):
        return self

    def __next__(self):
        return next(self.parser)

    def _read_header(self):
        """
        Called after opening VCF. This reads the meta header lines
        into a list, gets columns names and sample names
        """

        meta_header = []
        coln_header = []
        for line in self.reader:
            if line.startswith("##"):
                meta_header += [line]
            elif line.startswith("#CHROM"):
                coln_header = line.split("\t")
                break
            else:
                raise HeaderError(
                    "No column header found for VCF {}".format(self.filename)
                )
        return VcfHeader(meta_header, coln_header)

    def set_region(self, chrom, start=None, end=None):
        """
        Retrieve records by genomic location rather than reading
        records line-by-line.
        Sets self.reader and self.parser to iterators over the
        records retrieved.
        Args:
            chrom: chromosome or contig name. Required. You may
                   instead provide a region in format
                   'chr:start-end' if preferred instead of using
                   start/end arguments below.
            start: start position on chromosome/contig. 0-based
                   Default = None
            end:   end position on chromosome/contig.
                   Default = None
        >>> v = VcfReader(my_vcf)
        >>> v.set_region('chr1') #get all variants on chr1
        >>> for record in v:
        ... #do something with each record
        Because start coordinates are 0-based, to retrieve a variant
        at (or overlapping) chr1:1000000 the two examples below are
        equivalent:
        >>> v.set_region(chrom='chr1', start=999999 end=1000000)
        >>> v.set_region('chr1:1000000-1000000')
        """
        if not self._is_reg_file:
            raise ParseError("Cannot run set_region() on a non-regular file")
        if self.compressed or self.bcf:
            if not self._tabix:
                if not pysam:
                    raise ParseError(
                        "pysam not available. Please install "
                        + "(e.g. via 'pip install pysam' to "
                        + "search by location on bgzip compressed"
                        + "VCFs."
                    )
                if self.bcf:
                    idx = self.filename + ".csi"
                    preset = "bcf"
                else:
                    idx = self.filename + ".tbi"
                    preset = "vcf"
                if not os.path.isfile(idx):  # create index if it doesn't exist
                    pysam.tabix_index(self.filename, preset=preset)
                if self.bcf:
                    self._tabix = pysam.VariantFile(self.filename)
                else:
                    self._tabix = pysam.Tabixfile(
                        self.filename, index=idx, encoding=self.encoding
                    )
            try:
                if self.bcf:
                    self.reader = (
                        rec.__str__().rstrip()
                        for rec in self.file.fetch(str(chrom), start, end)
                        if rec.__str__().rstrip()
                    )
                else:
                    self.reader = self._tabix.fetch(str(chrom), start, end)
                self.parser = (
                    VcfRecord(
                        line,
                        self,
                    )
                    for line in self.reader
                )
            except ValueError:
                self.reader = iter([])
                self.parser = iter([])  # ignore missing contigs
        else:
            # easy solution - compress and index with pysam if not compressed,
            # but will be slow...
            # less easy solution, implement a custom index to seek to and from
            raise NotImplementedError(
                "Searching by location is not yet "
                + "implemented for non-bgzip compressed"
                + " VCFs."
            )


class VcfHeader(object):
    """Header class storing metadata and sample information for a vcf"""

    _meta_re = re.compile(r"""\#\#(\S+?)=(.*)""")  # should capture any metadata

    _dict_re = re.compile(
        r"""                 #for capturing dict-key metadata
                          \#\#(\S+)            #captures field name (e.g. INFO)
                          =<ID=(\S+?)          #captures metadata ID
                          (,(.*))*             #capture remaining keys/values
                          >""",  # dict line should end with a >
        re.X,
    )
    _subd_re = re.compile(
        r"""                 #for extracting keys/values from
                                               #group(3) of _dict_re.match()
                          ,(\S+?)=             #get key
                          (".+?"|[^\s,]+)      #value can either be quoted or
                                               #else should be all non-comma
                                               #and non-whitespace chars
                         """,
        re.X,
    )

    _required_keys = {
        "INFO": ["Number", "Type", "Description"],
        "FORMAT": ["Number", "Type", "Description"],
        "FILTER": ["Description"],
        "ALT": ["Description"],
    }
    __slots__ = [
        "meta_header",
        "col_header",
        "samples",
        "sample_cols",
        "metadata",
        "fileformat",
        "__csq_label",
        "__csq_fields",
        "_info_field_translater",
        "_format_field_translater",
    ]

    def __init__(self, meta_header, col_header):
        """
        Requires a list of metaheader lines and a list of column
        names
        """

        self.meta_header = meta_header
        self.col_header = col_header
        self.samples = col_header[9:] or None
        self.sample_cols = (
            dict([(x, i) for (i, x) in enumerate(col_header[9:], start=9)]) or None
        )
        for i, c in enumerate("#CHROM POS ID REF ALT QUAL FILTER INFO".split()):
            # check essential column names and order
            if self.col_header[i] != c:
                raise HeaderError(
                    "Invalid column name. Expected {}, got {}".format(
                        c, self.col_header[i]
                    )
                )
        if len(self.col_header) > 8:
            # 9th column must be FORMAT if present, but is optional
            if self.col_header[8] != "FORMAT":
                raise HeaderError(
                    "Invalid column name. Expected {}, got {}".format(
                        "FORMAT", self.col_header[8]
                    )
                )
        self.metadata = {"INFO": dict(), "FORMAT": dict()}
        self.csq_label = None
        self.csq_fields = None
        self._info_field_translater = {}
        self._format_field_translater = {}
        self._parse_metadata()

    def __str__(self):
        return (
            str.join("\n", self.meta_header)
            + "\n"
            + str.join("\t", self.col_header)
            + "\n"
        )

    @property
    def csq_label(self):
        """
        String labelling the INFO field label of VEP consequence
        annotations. Will raise a HeaderError if access is attempted
        but no VEP CSQ or ANN field is present in the header.
        """
        if self.__csq_label is None:
            self.csq_fields
        return self.__csq_label

    @csq_label.setter
    def csq_label(self, c):
        self.__csq_label = c

    @property
    def csq_fields(self):
        """
        A list of CSQ field names in the order they are represented
        in CSQ INFO field entries. Set to None on initialization.
        Will raise a HeaderError if access is attempted but no VEP
        CSQ or ANN field is present in the header.
        """

        if self.__csq_fields is None:
            if self.__csq_label is None:
                try:
                    csq_header = self.metadata["INFO"]["CSQ"][-1]
                    csq = "CSQ"
                except KeyError:
                    try:
                        csq_header = self.metadata["INFO"]["ANN"][-1]
                        csq = "ANN"
                    except KeyError:
                        raise HeaderError(
                            "No CSQ or ANN field in INFO header - "
                            + "unable to retrieve consequence fields."
                        )
                self.csq_label = csq
            else:
                csq = self.__csq_label
                csq_header = self.metadata["INFO"][csq][-1]
            try:
                match = csq_header["Description"].split("'")[1].split("'")[0]
                self.__csq_fields = [label.strip() for label in match.split("|")]
            except:
                raise HeaderError(
                    "Could not parse {} Format in ".format(csq)
                    + "header. Unable to retrieve consequence "
                    + "annotations."
                )
        return self.__csq_fields

    @csq_fields.setter
    def csq_fields(self, csq):
        self.__csq_fields = csq

    def _parse_metadata(self):
        """
        Extract INFO, FORMAT, FILTER and contig information from VCF
        meta header and store in dicts
        """

        # check first line is essential fileformat line
        ff_match = self._meta_re.match(self.meta_header[0])
        if not ff_match or ff_match.group(1) != "fileformat":
            raise ParseError(
                "Error: First line of VCF should be fileformat"
                + "metaheader (e.g. ##fileformat=VCFv4.2)"
            )
        else:
            self.fileformat = ff_match.group(2)
        for h in self.meta_header:
            self._parse_header_line(h)
        for field_type in ["FORMAT", "INFO"]:
            try:
                for field in self.metadata[field_type]:
                    self._set_field_translation(field_type, field)
            except KeyError:
                if field_type == "INFO":
                    # no FORMAT field in header is common - e.g. sites only VCFs
                    warnings.warn(
                        "No '{}' field in header!".format(field_type), stacklevel=5
                    )

    def _parse_header_line(self, h):
        """
        Parse a metaheader line and assign to self.metadata dict where
        keys are the type of metadata line and values are dicts of IDs to
        lists of either dicts of key-value pairs or string values.
        """

        match_d = self._dict_re.match(h)
        match_m = self._meta_re.match(h)
        field = None
        fid = None
        if match_d:
            # line is an e.g. INFO/FORMAT/FILTER/ALT/contig with multiple keys
            field = match_d.group(1)
            fid = match_d.group(2)
            rest = match_d.group(3) or ""
            d = dict([(x, y) for (x, y) in self._subd_re.findall(rest)])
            if not field in self.metadata:
                self.metadata[field] = {fid: [d]}
            else:
                if fid in self.metadata[field]:
                    # multiple values - extend list
                    self.metadata[field][fid].append(d)
                else:
                    self.metadata[field][fid] = [d]
        elif match_m:
            field = match_m.group(1)
            fid = match_m.group(2)
            if field in self.metadata:
                self.metadata[field].append(fid)
            else:
                self.metadata[field] = [fid]
        else:
            raise HeaderError("Invalid metaheader line {}".format(h))
        if field in self._required_keys:
            # ensure ALT/FORMAT/FILTER/INFO contain required keys
            last = self.metadata[field][fid][-1]  # check entry we've just added
            for k in self._required_keys[field]:
                if not k in last:
                    raise HeaderError(
                        "Missing required key '{}' in metaheader line: {}".format(k, h)
                    )

    def _set_field_translation(self, field_type, field):
        """
        returns a tuple of variable class type (int, float or str)
        and whether the value requires splitting
        """

        f = self.metadata[field_type][field][0]
        ctype = None
        if f["Type"] == "String" or f["Type"] == "Character":
            ctype = str
        elif f["Type"] == "Float":
            ctype = float
        elif f["Type"] == "Integer":
            ctype = int
        elif f["Type"] != "Flag":
            raise ParseError(
                "Unrecognised FORMAT Type '{}' in header".format(f["Type"])
            )
        split = False
        if f["Number"].isdigit():
            if int(f["Number"]) > 1:
                split = True
        else:  # if not digit should be 'A', 'G', 'R' or '.' - split
            split = True
        if field_type == "INFO":
            setter = self._info_field_translater
        elif field_type == "FORMAT":
            setter = self._format_field_translater
        else:
            raise ParseError(
                "'{}' not recognised as a ".format(field_type)
                + "field type for translation"
            )
        setter[field] = (ctype, split)

    def add_header_field(self, name, string=None, field_type=None, dictionary=None):
        """
        Add a header field with given name and optional field type,
        and dictionary of properties.
        Args:
            name:   name of field to add
            string: string to add to field. Ignored if 'dictionary'
                    is provided.
            field_type:
                    type of field - e.g. if INFO/FILTER/FORMAT
                    field. Required if providing a dictionary.
            dictionary:
                    a dict of keys to values for the given field.
                    If 'field_type' is specified, this arg must be
                    provided and must contain all the essential keys
                    for that field type. For example, an 'INFO'
                    field must have 'Number', 'Type', and
                    'Description' keys.
        """

        if dictionary is None and string is None:
            raise Exception("Either dictionary or string argument is required")
        if field_type is not None and field_type in self._required_keys:
            if dictionary is None:
                raise Exception(
                    "Header type {} requires a dictionary.".format(field_type)
                )
        if dictionary:
            if not field_type:
                raise Exception("field_type is required for use with " + "dictionary")
            if name in self.metadata[field_type]:
                self.metadata[field_type][name].append(dictionary)
            else:
                self.metadata[field_type][name] = [dictionary]
                self._set_field_translation(field_type, name)
            h_vals = []
            if field_type in self._required_keys:
                # first append required keys in order
                for k in self._required_keys[field_type]:
                    try:
                        h_vals.append(k + "=" + dictionary[k])
                    except KeyError:
                        raise Exception(
                            "Header type '{}'".format(field_type)
                            + " requires '{}' field.".format(k)
                        )
                # then append any additional keys
                for k in dictionary:
                    if k in self._required_keys[field_type]:
                        continue
                    h_vals.append(k + "=" + dictionary[k])
            else:
                for k in dictionary:
                    h_vals.append(k + "=" + dictionary[k])
            h_string = (
                str.join(",", ["##" + field_type + "=<ID=" + name] + h_vals) + ">"
            )
        else:
            h_string = "##" + name + "=" + string
            if name in self.metadata:
                self.metadata[name].append(string)
            else:
                self.metadata[name] = [string]
        self.meta_header.append(h_string)


class VcfRecord(object):
    """
    A single record from a Vcf created by parsing a non-header line
    from a VCF file. May contain multiple alternate alleles.
    """

    _gt_splitter = re.compile(r"[\/\|]")

    __slots__ = [
        "header",
        "cols",
        "CHROM",
        "POS",
        "ID",
        "REF",
        "ALT",
        "QUAL",
        "FILTER",
        "INFO",
        "FORMAT",
        "__SPAN",
        "__CSQ",
        "samples",
        "_sample_idx",
        "__CALLS",
        "__ALLELES",
        "__INFO_FIELDS",
        "GT_FORMAT",
        "_SAMPLE_GTS",
        "_got_gts",
        "_parsed_info",
        "_parsed_gts",
        "__is_sv",
        "__var_type",
    ]

    def __init__(self, line, caller):
        """
        VcfRecord objects require a line and a related VcfReader
        object for initialization.
        Args:
            line:   a non-header line from a VCF file without
                    newline character
            caller: a VcfReader object (normally the same VcfReader
                    object that read the input line). Metadata and
                    sample information will be read from this object
                    in order to initialize the VcfRecord object.
        """
        self.cols = line.split("\t", 9)  # only collect first 9 columns initially
        # splitting whole line on a VCF with
        # lots of columns/samples is costly
        try:
            (
                self.CHROM,
                pos,
                self.ID,
                self.REF,
                self.ALT,
                qual,
                self.FILTER,
                self.INFO,
            ) = self.cols[:8]
        except ValueError as err:
            if len(self.cols) < 8:
                raise ParseError(
                    "Not enough columns for following line:\n{}".format(line)
                )
            else:
                raise err

        self.POS = int(pos)
        self.ALT = self.ALT.split(",")

        try:
            self.QUAL = float(qual)
        except ValueError:
            self.QUAL = qual
        self.SPAN = None
        self.INFO_FIELDS = None
        self.FORMAT = None
        self.GT_FORMAT = None
        self.CALLS = None
        self.ALLELES = None
        self.header = caller.header
        self.CSQ = None
        self.IS_SV = None
        self.VAR_TYPE = None
        self._SAMPLE_GTS = {}
        self._parsed_info = {}
        self._parsed_gts = {}
        self._got_gts = False  # flag indicating whether we've already
        # retrieved GT dicts for every sample
        if len(self.cols) > 8:
            self.FORMAT = self.cols[8]
            self.GT_FORMAT = self.FORMAT.split(":")

    def __str__(self):
        """
        Represent the VCF line as it should appear as a VCF record.
        This uses the values stored in self.cols to avoid
        potentially needless splitting of sample calls.
        """
        return str.join("\t", self.cols)

    def add_ids(self, ids, replace=False):
        """
        Adds given IDs to the ID field of the VCF record. If the
        record already has an ID (i.e. is not '.') these IDs are
        added to the existing value(s) unless the replace
        argument is True.
        Args:
            ids:     A list of IDs to add.
            replace: If True, existing ID values are replaced,
                     otherwise the given IDs are added to.
                     Default = False.
        """
        if replace or self.ID == ".":
            self.ID = str.join(";", ids)
        else:
            uids = set(ids + self.ID.split(";"))
            self.ID = str.join(";", uids)
        self.cols[2] = self.ID  # also change cols so is reflected in __str__

    def add_filter(self, filters, replace=False):
        """
        Add provided filters to FILTER field of VCF. If replace is
        False, filters will be added to existing FILTER annotations
        except for those annotated just 'PASS'.
        Args:
            filters:    A list of filter strings to add.
            replace:    If True replace any existing filter annotations.
        """
        if not replace and self.FILTER != "PASS":
            filters.append(self.FILTER)
        filters = set(filters)
        self.FILTER = ";".join(sorted(filters))
        self.cols[6] = self.FILTER

    @property
    def IS_SV(self):
        """True if record represents a structural variant"""
        if self.__is_sv is None:
            self.__is_sv = "SVTYPE" in self.INFO_FIELDS
        return self.__is_sv

    @IS_SV.setter
    def IS_SV(self, is_sv):
        self.__is_sv = is_sv

    @property
    def VAR_TYPE(self):
        """
        String indicating whether the ALT is a snp, mnp, ins, del or sv
        """
        if self.__var_type is None:
            ref_len = len(self.REF)
            if self.IS_SV:
                sv_type = self.INFO_FIELDS["SVTYPE"].split(":")[0].lower()
                self.__var_type = sv_type
            elif all(ref_len == len(ALT) for ALT in self.ALT):
                if ref_len == 1:
                    self.__var_type = "snp"
                else:
                    self.__var_type = "mnp"

            elif all(ref_len < len(ALT) for ALT in self.ALT):
                self.__var_type = "ins"
            elif all(ref_len > len(ALT) for ALT in self.ALT):
                self.__var_type = "del"
            else:
                self.__var_type = "indel"

        return self.__var_type

    @VAR_TYPE.setter
    def VAR_TYPE(self, t):
        self.__var_type = t

    @property
    def ALLELES(self):
        """list of REF and ALT alleles in order"""

        if self.__ALLELES is None:
            self.__ALLELES = [self.REF] + self.ALT
        return self.__ALLELES

    @ALLELES.setter
    def ALLELES(self, alleles):
        self.__ALLELES = alleles

    @property
    def INFO_FIELDS(self):
        """
        A dict of INFO field names to values. All returned values
        are Strings except for Flags which are assigned True.
        To obtain values parsed into the appropriate Type as defined
        by the VCF header, use the 'parsed_info_fields()' method.
        """
        if self.__INFO_FIELDS is None:
            self.__INFO_FIELDS = OrderedDict()
            for i in self.INFO.split(";"):
                try:
                    (f, v) = i.split("=", 1)
                except ValueError:
                    (f, v) = (i, True)
                self.__INFO_FIELDS[f] = v
        return self.__INFO_FIELDS

    @INFO_FIELDS.setter
    def INFO_FIELDS(self, i):
        self.__INFO_FIELDS = i

    def remove_info_fields(self, fields):
        """Remove given INFO fields from the record if they exists."""
        rewrite = False
        for f in fields:
            val = self.INFO_FIELDS.pop(f, None)
            rewrite = True if val is not None else rewrite
        if rewrite:
            self._rewrite_info_string()

    def add_info_fields(self, info, append_existing=False):
        """
        Requires a dict of INFO field names to a list of values.
        Adds or replaces existing INFO fields in the record with
        the items in given dict.
        Args:
            info: A dict of INFO field names to add with values
                  being list of values for the given field.
            append_existing:
                  Add values to existing INFO fields in a record.
                  If the field being added already exists and this
                  argument is True, the values provided will be
                  added to the existing values. If the Number
                  property is a fixed value, multiple values at the
                  same index will be separated by '|' characters.
                  Otherwise they will be separated by commas.
                  Default = False.
        """
        for k, v in sorted(info.items()):
            if append_existing and k in self.INFO_FIELDS:
                self._append_to_existing_info(k, v)
            else:
                self.INFO_FIELDS[k] = v
        self._rewrite_info_string()

    def _append_to_existing_info(self, field, values):
        if field in self.header.metadata["INFO"]:
            if self.header.metadata["INFO"][field][-1]["Type"] == "Flag":
                self.INFO_FIELDS[field] = True
                return
            elif self.header.metadata["INFO"][field][-1]["Number"] == ".":
                self.INFO_FIELDS[field] += "," + values
                return
            elif self.header.metadata["INFO"][field][-1]["Number"] == "1":
                self.INFO_FIELDS[field] += "|" + values
                return
        old = self.INFO_FIELDS[field].split(",")
        new = str(values).split(",")
        if len(old) != len(new):
            raise ParseError(
                "New {} INFO field '{}'".format(field, values)
                + "has differing number of values to existing "
                + "field '{}'".format(self.INFO_FIELDS[field])
            )
        self.INFO_FIELDS[field] = str.join(
            ",", (str.join("|", x) for x in zip(old, new))
        )

    def _rewrite_info_string(self):
        info = []
        for f, v in self.INFO_FIELDS.items():
            if isinstance(v, bool):  # is Flag
                info.append(f)
            if isinstance(v, list):  # join list values with commas
                info.append(f + "=" + str.join(",", [str(x) for x in v]))
            else:
                info.append(f + "=" + str(v))
        self.INFO = str.join(";", info)
        self.cols[7] = self.INFO  # also change cols so is reflected in __str__

    def parsed_info_fields(self, fields=None):
        if fields is not None:
            f_list = [x for x in fields if x in self.INFO_FIELDS]
        else:
            f_list = list(self.INFO_FIELDS)
        d = dict(
            (f, self._get_parsed_info_value(f, self.INFO_FIELDS[f])) for f in f_list
        )
        return d

    def _get_parsed_info_value(self, field, value):
        try:
            return self._parsed_info[field]
        except KeyError:
            pass
        try:
            f = self.header._info_field_translater[field]
        except KeyError:
            try:
                f = (COMMON_INFO[field]["Class"], COMMON_INFO[field]["Split"])
                self.header._info_field_translater[field] = f
            except KeyError:
                raise ParseError(
                    "Unrecognised INFO field '{}'".format(field)
                    + "at {}:{}. ".format(self.CHROM, self.POS)
                    + "Non-standard  INFO fields should be "
                    + "represented in VCF header."
                )
        # f[0] is the class type of field, f[1] = True if values should be split
        if f[0] is None:  # is a Flag
            pv = True
        else:
            try:
                conv = lambda x: None if x == "." else f[0](x)
                if f[1]:
                    if isinstance(value, str):
                        pv = list(map(conv, value.split(",")))
                    else:
                        # if newly added INFO field it may already be a list
                        pv = list(map(conv, value))
                else:
                    pv = conv(value)
            except (ValueError, TypeError, AttributeError):
                pv = ""
                # raise ParseError("Unexpected value type in value '{}' for {} "
                #                 .format(field, value) + "INFO field at " +
                #                 " {}:{}" .format(self.CHROM, self.POS))
        self._parsed_info[field] = pv
        return pv

    @property
    def SPAN(self):
        """
        Returns end position of a record according to END INFO
        field if available, or otherwise POS + len(REF) -1.
        """
        if self.__SPAN == None:
            if "END" in self.INFO_FIELDS:
                self.__SPAN = int(self.INFO_FIELDS["END"])
            else:
                self.__SPAN = self.POS + len(self.REF) - 1
        return self.__SPAN

    @SPAN.setter
    def SPAN(self, end):
        self.__SPAN = end

    def add_format_field(self, field, sample_values, default_value="."):
        """
        Add FORMAT field and values to record.
        Args:
            field:  Name of the field to add. If it already exists
                    for this record it will be overwritten.
            sample_values:
                    dict of sample IDs to values for this field.
                    Missing samples will be given the value
                    specified by 'default_value' argument.
            default_value:
                    Value to add to genotype calls for samples not
                    provided in the sample_values dict. Defaults to
                    '.'.
        """
        self._got_gts = False  # will need to reparse calls
        if field in self.GT_FORMAT:
            i = self.GT_FORMAT.index(field)
        else:
            i = len(self.GT_FORMAT)
            self.GT_FORMAT.append(field)
            self.FORMAT = ":".join(self.GT_FORMAT)
            self.cols[8] = self.FORMAT
        for s in self.header.samples:
            value = sample_values.get(s, default_value)
            call = self.CALLS[s].split(":")
            if len(call) < i:  # pad truncated calls
                call = call + ["."] * (i - len(call))
            new_call = ":".join(call[0:i] + [value] + call[i + 1 :])
            self.cols[self.header.sample_cols[s]] = new_call
            self.__CALLS[s] = new_call

    @property
    def CALLS(self):
        """
        split sample call fields and assign to self.CALLS dict of
        sample id to call string.
        self.CALLS does not get created in __init__ to save on
        overhead in VCF with many samples where we might not be
        interested in parsing sample calls
        As of python 3.6 the CALLS dict will maintain sample order,
        but to safely get a list of calls in the same order as the
        input VCF the following syntax should be used:
        >>> v = VcfReader(my_vcf)
        >>> for record in v:
        ...     calls = [record.CALLS[x] for x in record.header.samples]
        """
        if self.__CALLS is None:
            if self.header.samples:
                calls = self.cols.pop()
                self.cols.extend(calls.split("\t"))
                self.__CALLS = dict(zip(self.header.samples, self.cols[9:]))
            else:
                self.__CALLS = {}
        return self.__CALLS

    @CALLS.setter
    def CALLS(self, calls):
        self.__CALLS = calls

    def sample_calls(self):
        """
        Retrieve a dict of sample names to a dict of genotype field
        names to values. All returned values are strings. For
        values cast to the appropriate types (int/float/string/list)
        use the 'parsed_gts' function.
        >>> record.sample_calls()
        {'Sample_1': {'GT': '0/0', 'AD': '10,0', 'DP': '10', 'GQ': '30'},
        'Sample_2': {'GT': '0/1', 'AD': '6,6', 'DP': '12', 'GQ': '33'}}
        >>> d = record.sample_calls()
        >>> s1 = d['Sample_1']
        >>> s1['GQ']
        '30'
        ...or more concisely:
        >>> record.sample_calls()['Sample_1']['GQ']
        '30'
        """

        if self._got_gts:
            return self._SAMPLE_GTS
        else:
            self._got_gts = True
            # get_sample_call() sets self._SAMPLE_GTS[s] for future use
            return dict([(s, self.get_sample_call(s)) for s in self.header.samples])

    def get_sample_call(self, sample):
        """
        Retrieve a dict of genotype field names to values for a
        single sample.
        This method creates dicts as needed for given samples, so
        may be more efficient than using the 'sample_calls()'
        method when parsing a VCF with many samples and you are
        only interested in a information from a small number of
        these samples.
        All values returned are strings. For values cast to an
        appropriate type (int/float/string/list) use the
        'parsed_gts(sample=[sample])' function.
        Args:
            sample: name of the sample to retrieve (as it appears
                    in the VCF header
        Example:
            >>> s1 = record.get_sample_call('Sample_1')
            >>> s1
            {'GT': '0/0', 'AD': '10,0', 'DP': '10', 'GQ': '30'}
            >>> s1['GQ']
            '30'
        """
        try:
            return self._SAMPLE_GTS[sample]
        except KeyError:
            if sample in self.CALLS:
                d = dict(
                    [
                        (f, v)
                        for (f, v) in zip(
                            self.GT_FORMAT, (self.CALLS[sample].split(":"))
                        )
                    ]
                )
                self._SAMPLE_GTS[sample] = d
                if len(self._SAMPLE_GTS) == len(self.header.samples):
                    # if we've now set self._SAMPLE_GTS for all samples set
                    # self._got_gts to True to prevent unnecessary looping in
                    # sample_calls() method
                    self._got_gts = True
                return d
            else:
                raise ParseError("Sample {} is not in VCF".format(sample))

    def parsed_gts(self, samples=None, fields=None):
        """
        Returns a dict of GT field names to dicts of sample names
        to values. By default, values for all samples and fields
        will be retrieved, but a list of sample IDs and a list of
        FORMAT fields to retrieve can be given.
        Missing values will be assigned to None.
        Because this is a relatively costly function, you are
        advised to avoid calling this repeatedly for a single
        record - you may speed things up by only calling for a
        subset of samples and fields but in any case you probably
        want to call this function once only per record, storing
        the results in a variable.
        Args:
            samples: Optional list of sample names to retrieve
                     values for. Default = None (values retrieved
                     for all samples)
            fields:  Optional list of field names (as they appear
                     in the FORMAT field of the record) to retrieve
                     values for. Default = None (values retrieved
                     for all fields)
        >>> record.parsed_gts()
        {'GT': {'Sample_1': (0, 0), 'Sample_2': (0, 1)},
        'AD': {'Sample_1': (10, 0), 'Sample_2': (6, 6)},
        'DP': {'Sample_1': 10, 'Sample_2': 12},
        'GQ': {'Sample_1': 30, 'Sample_2': 33}}
        >>> record.parsed_gts(samples=['Sample_2'])
        {'GT': {'Sample_2': (0, 1)}, 'AD': {'Sample_2': (6, 6)},
        'DP': {'Sample_2': 12}, 'GQ': {'Sample_2': 33}}
        >>>  record.parsed_gts(fields=['GT', 'GQ'])
        {'GT': {'Sample_1': (0, 0), 'Sample_2': (0, 1)},
        'GQ': {'Sample_1': 30, 'Sample_2': 33}}
        """

        d = dict()
        if fields is not None:
            f_list = fields
        else:
            f_list = self.GT_FORMAT
        if samples is not None:
            s_list = samples
        else:
            s_list = self.header.samples
        updated = False
        for f in f_list:
            if f in self._parsed_gts:
                missing_samps = [s for s in s_list if s not in self._parsed_gts[f]]
                if not missing_samps:
                    d[f] = self._parsed_gts[f]
                    continue
                else:
                    updated = True
                    if missing_samps != s_list:  # some missing samps,but not all
                        d[f].update(
                            dict(
                                (s, self._parsed_gts[f][s])
                                for s in s_list
                                if s in self._parsed_gts[f]
                            )
                        )
                    d[f].update(
                        dict(
                            zip(
                                missing_samps,
                                self._get_parsed_gt_fields(
                                    f,
                                    (
                                        (
                                            self.sample_calls()[s][f]
                                            if f in self.sample_calls()[s]
                                            else None
                                        )
                                        for s in missing_samps
                                    ),
                                ),
                            )
                        )
                    )

            else:
                updated = True
                d[f] = dict(
                    zip(
                        s_list,
                        self._get_parsed_gt_fields(
                            f,
                            (
                                (
                                    self.sample_calls()[s][f]
                                    if f in self.sample_calls()[s]
                                    else None
                                )
                                for s in s_list
                            ),
                        ),
                    )
                )
        if updated:
            for f in f_list:
                self._parsed_gts[f].update(d[f])
        return d

    def _get_parsed_gt_fields(self, field, values=[], none_to_zeros=False):
        # TODO - make this more efficient with a cython extension
        """
        Retrieves values of genotype field parsed so that the
        returned values are of the expected type (str, int or float)
        and are split into list format if appropriate. Fields are
        handled according to the information present in the VCF
        header metadata
        Args:
            field:  genotype field to retrieve as it appears in the
                    FORMAT field of the VCF record and in the VCF
                    header.
            values: list of values from a genotype field
            none_to_zeros:
                    for float or integer FORMAT fields convert '.'
                    values to zero instead of None.
        """

        try:
            f = self.header._format_field_translater[field]
        except KeyError:
            try:
                f = (COMMON_FORMAT[field]["Class"], COMMON_FORMAT[field]["Split"])
                self.header._format_field_translater[field] = f
            except KeyError:
                raise ParseError(
                    "Unrecognised FORMAT field '{}'".format(field)
                    + "at {}:{}. ".format(self.CHROM, self.POS)
                    + "Non-standard  FORMAT fields should be "
                    + "represented in VCF header."
                )
        # f[0] is the class type of field, f[1] = True if values should be split
        pv = []
        split_default = None
        if f[1] and none_to_zeros:
            if f[0] is not None and f[0] is not str:
                split_default = 0
        for val in values:
            if field == "GT":  # GT is a special case, make tuples of alleles
                alleles = self._gt_splitter.split(val)
                try:
                    pv.append(tuple(None if x == "." else int(x) for x in alleles))
                except ValueError:
                    raise ParseErrror(
                        "Could not parse GT {} ".format(val)
                        + "at {}:{}".format(self.CHROM, self.POS)
                    )
            else:
                try:
                    if f[1]:
                        pv.append(
                            tuple(
                                f[0](x) if x != "." else split_default
                                for x in val.split(",")
                            )
                        )
                    else:
                        pv.append(f[0](val))
                except (ValueError, TypeError, AttributeError) as err:
                    if val is None or val == ".":
                        v = None
                        if none_to_zeros:
                            v = split_default
                        if f[1]:
                            pv.append((v,))
                        else:
                            pv.append(v)
                    else:
                        # raise err
                        raise ParseError(
                            "Unexpected value ('{}')".format(val)
                            + " for {} ".format(field)
                            + "FORMAT "
                            + "field at {}:{}".format(self.CHROM, self.POS)
                        )
        return pv

    @property
    def CSQ(self):
        """
        A list of dicts of CSQ/ANN annotations from VEP to values.
        Empty values are represented by empty Strings. Will raise
        a HeaderError if the associated VCF header does not contain
        CSQ/ANN information and a ParseError if the record being
        parsed does not contain a CSQ/ANN annotation in the INFO
        field.
        """
        if self.__CSQ is None:
            lbl = self.header.csq_label
            try:
                csqs = self.INFO_FIELDS[lbl].split(",")
            except KeyError:
                raise ParseError(
                    "Could not find '{}' label in ".format(lbl)
                    + "INFO field of record at {}:{}".format(self.CHROM, self.POS)
                )
            self.__CSQ = []
            alleleToNum = {}
            for c in csqs:
                d = OrderedDict(
                    [(k, v) for (k, v) in zip(self.header.csq_fields, c.split("|"))]
                )

                self.__CSQ.append(d)
        return self.__CSQ

    @CSQ.setter
    def CSQ(self, c):
        self.__CSQ = c

    def in_cis_with(self, sample, allele, other, other_allele):
        """
        Returns True if the two alleles are physically phased
        according to the PID and PGT fields of both records.
        Args:
            sample: Sample ID to check phasing data for.
            allele: Allele number of this record.
            other:  Other record to compare with this record.
            other_allele:
                    Allele number for other record.
        """
        if "PID" not in self.GT_FORMAT or "PID" not in other.GT_FORMAT:
            return False
        if "PGT" not in self.GT_FORMAT or "PGT" not in other.GT_FORMAT:
            return False
        try:
            pid1 = self.sample_calls()[sample]["PID"]
            pid2 = other.sample_calls()[sample]["PID"]
            pgt1 = self.sample_calls()[sample]["PGT"]
            pgt2 = other.sample_calls()[sample]["PGT"]
        except KeyError:
            # when joining VCFs together only some samples may have PID/PGT
            return False
        if pid1 != pid2:
            return False
        if pgt1 == "." or pgt2 == ".":
            return False
        try:
            phase1 = pgt1.split("|").index(str(allele))
            phase2 = pgt2.split("|").index(str(other_allele))
            return phase1 == phase2
        except ValueError:  # allele might not be in phase group
            # TODO: for some multiallelic variants I've spotted that the GT is
            # '0/2' while the PGT is '0|1' - is this a bug in HaplotypeCaller
            # or will the 'ALT' always be 1 in the PGT?
            return False


class HeaderError(ValueError):
    pass


class ParseError(ValueError):
    pass
