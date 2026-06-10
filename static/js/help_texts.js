var helpTexts = {
 "setup":
   "Refer to the manual (top right) for more information.",

 "cmp-trt":
  "<table><tr><th>Zoom in/out of subplots</th><td>Right-click subplot</td></tr>" +
  "<tr><th>Zoom in/out of correlation</th><td>Right-click plot column</td></tr></table>",
 "core":
 "Selects a sample set of specified <i>Size</i> that represents the greatest diversity within the samples.<br>" +
 "Based on genetic distance data from chosen <i>Marker set</i>.<br>" +
 "Select samples to <i>Include</i> or <i>Exclude</i> in the core set using dropdowns below.",
 "transform":
  "<table><tr><th>Input traits</th><td>All <i>Traits</i> and <i>Trait groups</i> selected under <i>Compare traits</i></td></tr>" +
  "<tr><th>PCA transform</th><td>PCs are named <i>&lt;trait name prefix&gt;_PC&lt;number&gt;</i></td></tr>" +
  "<tr><th>All other transforms</th><td>Traits are named <i>&lt;original name&gt;_&lt;Trait name suffix&gt;</i></td></tr>" +
  "<tr><th>More info</th><td>In the manual (top right)</td></tr></table>",

// "trt-corr":
//  "Correct phenotype values for sample relatedness<br>Values regressed with kinship (random effect) and principal compponents (fixed effect),<br> residuals from the model are saved as corrected values",

// "fa-input":
//  "<table><tr><th>Coding regions</th><td>Start to end of each gene in query</td></tr>" +
//  "<tr><th>Promoter regions</th><td><i>n</i> bp upstream of each gene in query</td></tr>" +
//  "<tr><th>Entire region</th><td>Entire span of the query</td></tr></table>",
// "fa-impact":
//  "<table><tr><th>High</th><td>Likely disruptive to the protein</td></tr>" +
//  "<tr><th>Moderate</th><td>Might change protein effectiveness</td></tr>" +
//  "<tr><th>Low</th><td>Likely to have no effect on the protein</td></tr>" +
//  "<tr><th>Modifier</th><td>Usually non-coding, difficult to predict impact</td></tr></table>",
// "fa-filter":
//  "<table><tr><th>Sample groups</th><td>Created in Phenotype/sample view</td></tr>" +
//  "<tr><th>Types</th><td>Structural variant types may be annotated by the variant caller used<br>SNP, MNP and indel variants are annotated based on length</td></tr>" +
//  "<tr><th>Consequence</th><td>H: High, M: Moderate, L: Low, -: Modifier<br>Sequence Ontology nomenclature, see <a target='_blank' href='http://www.sequenceontology.org/'>here</a></td></tr></table>",
 "fa-skip":
  "<table><tr><th>Reference-only samples</th><td>Skip samples identical to the reference across the region specified</td></tr>" +
  "<tr><th>Reference-only variants</th><td>Skip variants where each sample is identical to the reference</td></tr>" +
  "<tr><th>Absent-only samples</th><td>Skip samples with no genotypes called across the region specified</td></tr>" +
  "<tr><th>Absent-only variants</th><td>Skip variants where each sample has no genotype called</td></tr>" +
  "<tr><th>Variants without variation</th><td>Skip variants where each sample is genotyped the same<br>Useful for finding differences between a pair of accessions</td></tr></table>",

 "vars-table":
  "<table><tr><th>Genotype encoding</th><td>" +
    "<div style='padding-right: 15px; display: inline-block;'>" +
      "<div style='width: 28px; display: inline-block; text-align: center; padding-right: 7px;'>⸺</div>Reference homozygous<br>" +
      "<div style='width: 28px; display: inline-block; text-align: center; padding-right: 7px;'>⟩</div>Part of preceding deletion</div>" +
    "<div style='padding-right: 5px; display: inline-block;'>" +
      "<div style='width: 28px; display: inline-block; text-align: center; padding-right: 7px;'>—/X</div>REF/X heterozygous<br>" +
      "<div style='width: 28px; display: inline-block; text-align: center; padding-right: 7px;'>·</div>Missing call</div></td></tr>" +
  "<tr><th>Additional information</th><td>Right-click <i>Reference, Annotation, Severity score</i> or genotype<br>" +
  "More info about interpreting genotype qualities <a target='_blank' href='https://samtools.github.io/hts-specs/VCFv4.1.pdf#subsubsection.1.4.2'>here</a></td></tr>" +
  "<tr><th>Color-coding</th><td>" +
  "Genotypes by error probability: <br>" +
   "<table><tr class='center-text'><td style='width:45px; background-color: rgba(204, 0, 0, 0.65);'>100%</td>" +
   "<td style='width:45px; background-color: rgba(204, 80, 0, 0.65);'>10%</td>" +
   "<td style='width:45px; background-color: rgba(204, 160, 0, 0.65);'>1%</td>" +
   "<td style='width:45px; background-color: rgba(168, 204, 0, 0.65);'>0.1%</td>" +
   "<td style='width:45px; background-color: rgba(88, 204, 0, 0.65);'>0.01%</td>" +
   "<td style='width:45px; background-color: rgba(8, 204, 0, 0.65);'><0.001%</td></tr></table>" +
  "Misingness and heterozygosity by value: <br>" +
   "<table><tr class='center-text'><td style='width:45px; background-color: rgba(0, 204, 0, 0.65);'>0</td>" +
   "<td style='width:45px; background-color: rgba(102, 204, 0, 0.65);'>0.05</td>" +
   "<td style='width:45px; background-color: rgba(204, 204, 0, 0.65);'>0.1</td>" +
   "<td style='width:45px; background-color: rgba(204, 102, 0, 0.65);'>0.15</td>" +
   "<td style='width:45px; background-color: rgba(204, 0, 0, 0.65);'>>=0.2</td></tr></table>" +
  "Minor allele frequency by value: <br>" +
   "<table><tr class='center-text'><td style='width:45px; background-color: rgba(204, 0, 0, 0.65);'>0</td>" +
   "<td style='width:45px; background-color: rgba(204, 102, 0, 0.65);'>0.125</td>" +
   "<td style='width:45px; background-color: rgba(204, 204, 0, 0.65);'>0.25</td>" +
   "<td style='width:45px; background-color: rgba(102, 204, 0, 0.65);'>0.375</td>" +
   "<td style='width:45px; background-color: rgba(0, 204, 0, 0.65);'>0.5</td></tr></table>" +
  "</td></tr></table>",

// "filt-ext":
//  "Use an existing marker set as the input for filtering<br>Use for reapplying filters or constructing haplotype block from existing marker sets",
 "filt-indel":
  "Specify maximum size of 1 to filter for only SNPs.",
// "ld-prune":
//  "For each marker, removes all markers within <i>Window</i> with r-squared above <i>Threshold</i>",
// "impute":
//  "Imputation and phasing carried out by <a target='_blank' href='http://faculty.washington.edu/browning/beagle/beagle.html'>BEAGLE</a>",
// "haps":
//  "<table><tr><th>LD-based</th><td>Constructed using the <a target='_blank' href='https://doi.org/10.1126/science.1069424'>Gabriel 2002</a> method<br>Blocks larger than <i>Max size</i> will be split.</td></tr>" +
//  "<tr><th>SNP number window</th><td>Blocks of <i>span</i> markers, moving by <i>step</i> each time</td></tr>" +
//  "<tr><th>Base-pair window</th><td>Blocks of <i>span</i> base pairs, moving by <i>step</i> each time</td></tr></table>",

 "assoc-qc":
  "<table><tr><th>View <i>Marker set</i> filters</th><td>Right-click <i>Marker set</i> checkbox label</td></tr>" +
  "<tr><th>View members of<br><i>Sample</i> or <i>Region</i> group</th><td>Right-click checkbox label</td></tr></table>",

 "assoc-opts":
  "<table><tr><th>View <i>Marker set</i> filters</th><td>Right-click <i>Marker set</i> checkbox label</td></tr>" +
  "<tr><th>View members of<br><i>Sample</i> or <i>Trait</i> group</th><td>Right-click checkbox label</td></tr></table>",

// "model":
//  "<table><tr><th>GLM</th><td>Generalized linear model via <a target='_blank' href='https://doi.org/10.1038/ng.2310'>GEMMA</a><br>Does not correct for kinship</td></tr>" +
//  "<tr><th>MLM</th><td>Mixed linear model via GEMMA<br>Must correct for kinship</td></tr>" +
//  "<tr><th>Multi-trait MLM</th><td>Models association with multiple traits at once, via GEMMA</td></tr>" +
//  "<tr><th>BLINK</th><td>Bayesian multi-locus model, via <a target='_blank' href='https://zzlab.net/blink/'>BLINK</a><br>Must correct for kinship, calculates kinship internally</td></tr>" +
//  "<tr><th>BayesB</th><td>Bayesian mixed linear model, via <a target='_blank' href=https://cnsgenomics.com/software/gctb//'>GCTB</a><br>Assumes most markers have zero effects</td></tr>" +
//  "<tr><th>BayesC/th><td>Bayesian mixed linear model, via GCTB<br>Assumes most markers have non-zero effects</td></tr>" +
//  "<tr><th>Random forest</th><td>Calculates marker importance based on random forest regression<br>Averages the importances of <i>n</i> iterations</td></tr></table>", //<br>Corrects for spurious associations using a permutation approach</td></tr></table>",
// "cov":
//  ".txt or .tsv file, tab-separated<br>First column is for sample IDs, all samples must be present<br>No header of column names must be present",
// "rf-options":
//  "<table><tr><th>Optimization</th><td>Input <i>mtry</i> and <i>ntree</i> values, comma and space (, ) separated<br>Look for minimum value that gives PVE equal to that of higher values</td></tr>" +
//  "<tr><th>Test markers</th><td><i>mtry</i>, number of markers to consider in each split</td></tr>" +
//  "<tr><th>Trees per forest</th><td><i>ntree</i>, number of trees in the forest</td></tr>" +
//  "<tr><th>Iterations</th><td>Number of forests to aggregate for importance estimation</td></tr>" +
//  "<tr><th>Permutations</th><td>Number of forests to aggregate for correcting spurious associations</td></tr></table>",

 "assoc-res":
  "<table><tr><th>View <i>Marker set</i> filters</th><td>Right-click <i>Marker set</i> checkbox label</td></tr>" +
  "<tr><th>View <i>Run</i> parameters</th><td>Right-click plot legend dataset</td></tr>" +
  "<tr><th>Add <i>Selected regions</i></th><td>Right-click marker or gene, or<br>zoom into a single chromosome region</td></tr>" +
  "<tr><th>Highlight gene</th><td>Left-click gene</td></tr>" +
  "<tr><th>View per-marker LD</th><td>Left-click marker</td></tr></table>",
 "as-group":
  "<table><tr><th>Plot as groups</th><td>Groups results per <i>Window</i><br>Plots as a single marker; size represents the number of runs<br>with significant association within the window (at or above <i>Threshold</i>)<br>Save traits as a <i>Trait group</i> to plot</td></tr>" +
  "<tr><th>Plot averages</th><td>Plot the average significance and average position of peak markers<br>(by default, plots the position and significance of the most significant peak marker)</td></tr>" +
  "<tr><th>Combinations separately</th><td>Plot all combinations of <i>Trait group, Marker set</i> and <i>Model</i> separately<br>(by default, all results for a <i> Trait group</i> are grouped into one)</td></tr>" +
  "<tr><th>Number of traits</th><td>Minimum number of traits with significant association within <i>Window</i> to plot</td></tr></table>",
// "finemap":
//  "<table><tr><th>PiMASS</th><td>Determines markers with highest posterior probability of being causal, via <a target='_blank' href='https://doi.org/10.1214/11-AOAS455'>PiMASS</a></td></tr>" +
//  "<th>Random forest</th><td>Calculates marker importance based on random forest regression<br>Averages the importances of <i>n</i> iterations<br>Corrects for spurious associations using a permutation approach</td></tr></table>",
// "pimass-options":
//  "Options respectively correspond to:<br><i> -w, -s, -num, -hmin, -hmax, -pmin, -pmax, -nstart, -smin, -smax</i><br>Check <a target='_blank' href='https://www.haplotype.org/download/pimass-manual.pdf#page=7'>PiMASS manual</a> for more info.",

 "make-group":
  "Select at least one item to group.",
 "delete-group":
  "Select at least one group to delete.",
 "select-column":
  "Make a selection from the field above.",
 "select-gene":
  "Choose one or more genes from the Genes table.",
 "select-var":
  "Choose one or more variants from the Variants table.",
 "select-phe":
  "Choose a phenotype to plot from the field above.",
};
