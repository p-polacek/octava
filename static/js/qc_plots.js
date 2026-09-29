function percentile(values, p, scaled=false) {
    const sorted = [...values].sort((a, b) => a - b);

    const pos = (p / 100) * (sorted.length - 1);
    const lower = Math.floor(pos);
    const upper = Math.ceil(pos);
    var perc;

    if (lower === upper) {
        perc = sorted[lower];
    } else {
      const weight = pos - lower;
      perc = sorted[lower] * (1 - weight) + sorted[upper] * weight;
    }

    if (scaled) {
      return (perc - Math.min(...values)) / Math.max(...values)
    }
    return perc
}

function fetchChr() {
  var filtFrom = document.querySelector("input[name='filt_from']:checked");
  if (filtFrom) {
    markers = filtFrom.value;
  }
  if (typeof(haps) !== "undefined" && haps[markers]) {
    var plotChoice = document.getElementById("hap-choice-field").value;
  } else {
    var plotChoice = document.getElementById("choice-field").value;
  }
  var chrDiv = document.getElementById("chr-div");
  var pcaPcDiv = document.getElementById("pca-pc-div");
  var pcaSampDiv = document.getElementById("pca-sample-div");
  if (plotChoice === "ld") {
    chrDiv.classList.remove("collapsed");
    pcaPcDiv.classList.add("collapsed");
    pcaSampDiv.classList.add("collapsed");
    var chr = document.getElementById("chr-field").value;
    if (markers == null) {
      markers = "";
    }
    var requests = [];
    var request = $.ajax({
      url: "/chr_ld",
      method: "POST",
      data: {
        genome: genome,
        markers: markers,
        chrom: chr,
      },
      success: function (data) {
        chrData = data;
      }
    });
    requests.push(request);

    $.when.apply($, requests).done(function() {
      showPlot();
    })
  } else if (plotChoice == "bdist") {
    chrDiv.classList.remove("collapsed");
    pcaPcDiv.classList.add("collapsed");
    pcaSampDiv.classList.add("collapsed");
    var chr = document.getElementById("chr-field").value;
    var requests = [];
    var request = $.ajax({
      url: "/chr_bdist",
      method: "POST",
      data: {
        genome: genome,
        markers: markers,
        chrom: chr,
      },
      success: function (data) {
        chrData = data;
      }
    });
    requests.push(request);

    $.when.apply($, requests).done(function() {
      showPlot();
    })
  } else {
    chrDiv.classList.add("collapsed");
    pcaPcDiv.classList.add("collapsed");
    pcaSampDiv.classList.add("collapsed");
    if (plotChoice === "pca") {
      pcaPcDiv.classList.remove("collapsed");
      pcaSampDiv.classList.remove("collapsed");
    } else if (plotChoice === "kin") {
      pcaSampDiv.classList.remove("collapsed");
    }
    showPlot();
  }
}

// Plot
function showPlot() {
  if (typeof(haps) !== "undefined" && haps[markers]) {
    var plotChoice = document.getElementById("hap-choice-field").value;
  } else {
    var plotChoice = document.getElementById("choice-field").value;
  }
  if (histData == null) {
    return;
  }
  var data = histData[plotChoice];
  var traces = [];

  if (plotChoice === "dens") {
    var origSizes = Object.fromEntries(
      names.map((ch, i) => [ch, sizes[i]]));
    var chrSizes = Object.fromEntries(
      names.map((ch, i) => [ch, Math.floor(sizes[i] / mdbin) * mdbin + mdbin]));
    var binOffsets = {};
    chrNames = names.filter(ch => !concatNames.includes(ch));

    var otherSize = concatNames.reduce((s, ch) => {
      s += origSizes[ch];
      delete chrSizes[ch];
      return s;
    }, 0);

    chrNames.forEach((ch, idx) => {
      binOffsets[ch] = chrNames.slice(0, idx).reduce((s, c) => {
        s += chrSizes[c];
        return s;
      }, 0);
    });

    if (concatNames.length > 0) {
      binOffsets["Other"] = chrNames.reduce((s, c) => {
        s += chrSizes[c];
        return s;
      }, 0);
      chrSizes["Other"] = otherSize;
      chrNames.push("Other");
    }

    chrSizes = chrNames.map( c => chrSizes[c]);
    binOffsets = chrNames.map( c => binOffsets[c]);

    var chrTicks = [];
    for (let i = 0; i < chrNames.length; i++) {
      chrTicks.push(binOffsets[i] + Math.floor(chrSizes[i] / 2));
    }

    var trace = {
      type: "bar",
      x: data[1],
      y: data[0],
      text: data[2],
      textposition: "none",
      hoverinfo: "text+y",
      hovertemplate: "%{text}<br>Count: %{y}",
      hoverlabel: {
        namelength: 0,
        bgcolor: "white",
      },
      marker: { color: "#60cdcd" },
      visible: true,
    };
    traces.push(trace);

    var lastChrIdx = chrNames.length - 1;
    var xMax = sizes.reduce((s, a) => s + a, 0);
    var yMax = Math.max(...data[0]) * 1.05;
    var layout = {
      //plot_bgcolor: "rgba(0, 0, 0, 0)",
      //paper_bgcolor: "rgba(0, 0, 0, 0)",
      xaxis: {
        tickvals: chrTicks,
        ticktext: chrNames,
        showgrid: false,
        range: [0, xMax],
        minallowed: 0,
        maxallowed: xMax,
      },
      yaxis: {
        title: "Count",
        range: [0, yMax],
        minallowed: 0,
        maxallowed: yMax,
      },
      margin: {
        l: 80,
        r: 40,
        b: 50,
        t: 70,
        pad: 4
      },
      bargap: 0,
      hovermode: "x",
    };
    var spanShapes = [];
    for (i = 0; i < chrNames.length; i+= 2) {
      var shape = {
        type: "rect",
        xref: "x",
        yref: "paper",
        x0: binOffsets[i],
        y0: 0,
        x1: binOffsets[i] + chrSizes[i],
        y1: 1,
        fillcolor: "#000",
        opacity: 0.05,
        line: {
          width: 0,
        },
      };
      spanShapes.push(shape);
    }
    layout.shapes = spanShapes;
  } else if (plotChoice === "pca") {
    var[plottedSamples, PCAdata, PVEs] = data;
    var sampleGroupChoice = [];
    var showLabels = document.getElementById("show-labels").checked;
    var pcNumChoice = document.querySelectorAll("input[name='pc_nums']:checked");
    var pcNums = [...pcNumChoice].map(v => parseInt(v.value));
    document.querySelectorAll(".sample_group-checkbox").forEach( function(cbox) {
      if (cbox.checked) {
        // Keep both group contents and name
        var samplesList = JSON.parse(cbox.value);
        sampleGroupChoice.push(
          [
            samplesList,
            cbox.nextElementSibling.textContent]);
      }
    });

    var trait = document.querySelector("input[name='by_phe']:checked");
    var sampleChoice;
    if (sampleGroupChoice.length > 0 && trait != null) {
      sampleChoice = sampleGroupChoice.map(v => v[0]).flat();
      sampleChoice = plottedSamples.filter(s => sampleChoice.includes(s));
    } else {
      sampleChoice = plottedSamples;
    }

    var legend = (sampleGroupChoice.length > 0 && trait == null);
    var opacity = legend ? 0.6 : 1;

    var trace = {
      x: PCAdata[pcNums[0]].filter((v, i) => sampleChoice.includes(plottedSamples[i])),
      y: PCAdata[pcNums[1]].filter((v, i) => sampleChoice.includes(plottedSamples[i])),
      name: "<i>All samples</i>",
      type: "scatter",
      traceType: "PCA",
      mode: "markers",
      text: sampleChoice,
      opacity: opacity,
      showlegend: legend,
      hovertemplate: "%{text}",
      hoverlabel: {
        namelength: 0,
        bgcolor: "white",
      },
      marker: {
         size: 6,
         color: "#60cdcd",
      },
      visible: true,
    }

    if (pcNums.length > 2) {
      trace.z = PCAdata[pcNums[2]].filter((v, i) => sampleChoice.includes(plottedSamples[i]));
      trace.type = "scatter3d";
      trace.marker.size = 4;
    }
    traces.push(trace);

    // Plot of sample groups
    if (sampleGroupChoice.length > 0 && trait == null) {
      for (var [idx, group] of sampleGroupChoice.entries()) {
        var [members, groupName] = group;
        var trace = {
          x: PCAdata[pcNums[0]].filter((_, idx) => members.includes(plottedSamples[idx])),
          y: PCAdata[pcNums[1]].filter((_, idx) => members.includes(plottedSamples[idx])),
          name: groupName,
          type: "scatter",
          traceType: "PCA",
          mode: "markers",
          text: plottedSamples.filter(v => members.includes(v)),
          opacity: opacity,
          showlegend: legend,
          hovertemplate: "%{text}",
          hoverlabel: {
            namelength: 0,
            bgcolor: "white",
          },
          marker: {
             color: colorCycle[idx],
             size: 6,
          },
          visible: true,
        }

        if (pcNums.length > 2) {
          trace.z = PCAdata[pcNums[2]].filter((_, idx) => members.includes(plottedSamples[idx]));
          trace.type = "scatter3d";
          trace.marker.size = 4;
        }
        traces.push(trace);
      }
    }

    if (showLabels) {
      // Traces of labels
      var trace = {
        x: PCAdata[pcNums[0]].filter((v, i) => sampleChoice.includes(plottedSamples[i])),
        y: PCAdata[pcNums[1]].filter((v, i) => sampleChoice.includes(plottedSamples[i])),
        name: "<i>All samples</i>",
        type: "scatter",
        mode: "text",
        text: sampleChoice,
        showlegend: legend,
        textposition: "top right",
        textfont: {
          size: 11
        },
        hoverinfo: "skip",
        visible: true,
      };
      if (pcNums.length > 2) {
        trace.z = PCAdata[pcNums[2]].filter((v, i) => sampleChoice.includes(plottedSamples[i]));
        trace.type = "scatter3d";
      }
      traces.push(trace);

      if (sampleGroupChoice.length > 0 && trait == null) {
        for (var [idx, group] of sampleGroupChoice.entries()) {
          var [members, groupName] = group;
          var trace = {
            x: PCAdata[pcNums[0]].filter((_, idx) => members.includes(plottedSamples[idx])),
            y: PCAdata[pcNums[1]].filter((_, idx) => members.includes(plottedSamples[idx])),
            name: groupName,
            type: "scatter",
            mode: "text",
            text: plottedSamples.filter(v => members.includes(v)),
            showlegend: legend,
            textposition: "top right",
            textfont: {
              size: 11
            },
            hoverinfo: "skip",
            visible: true,
          }

          if (pcNums.length > 2) {
            trace.z = PCAdata[pcNums[2]].filter((_, idx) => members.includes(plottedSamples[idx]));
            trace.type = "scatter3d";
          }
          traces.push(trace);
        }
      }
    }

    var layout;
    if (pcNums.length > 2) {
      layout = {
        legend: {
          x: 1.1,
          y: 0.95,
          xanchor: "left",
          yanchor: "top",
        },
        showlegend: true,
        scene: {
          camera: {
            eye: {
              x: 0,
              y: 0,
              z: 2,
            }
          },
          xaxis: {
            title: `PC ${pcNums[0]+1}`,
            showgrid: false,
          },
          yaxis: {
            title: `PC ${pcNums[1]+1}`,
            showgrid: false,
           },
          zaxis: {
            title: `PC ${pcNums[2]+1}`,
            showgrid: false,
          },
        },
        margin: {
          t: 30,
          b: 10,
          l: 10,
          r: 100,
          pad: 4
        },
      };
    } else {
      layout = {
        legend: {
          x: 1.1,
          y: 0.95,
          xanchor: "left",
          yanchor: "top",
        },
        showlegend: true,
        xaxis: {
          title: `PC ${pcNums[0]+1}`,
        },
        yaxis: {
          title: `PC ${pcNums[1]+1}`,
        },
        margin: {
          t: 80,
          b: 40,
          l: 50,
          r: 100,
          pad: 4
        },
      };
    }

    var annotations = [];
    annotations.push({
      text: "PVE:",
      x: 1,
      y: 0.95,
      showarrow: false,
      yanchor: "middle",
      xanchor: "left",
      valign: "top",
      xref: "paper",
      yref: "paper",
      font: { family: "monospace", size: 14 },
      visible: true,
    });
    PVEs.forEach(function(PVE, idx) {
      idx ++;
      annText = `PC${idx}: `;
      if (idx < 10 && PVE < 10) {
        annText += " ";
      }
      annText += `${PVE.toFixed(2)}%`;
      annotations.push({
        text: annText,
        x: 1,
        y: 0.95 - (idx * Math.min(0.05, (1 / PVEs.length))),
        showarrow: false,
        yanchor: "middle",
        xanchor: "left",
        valign: "top",
        xref: "paper",
        yref: "paper",
        font: { family: "monospace", size: 14 },
        visible: true,
      });
    });
    layout.annotations = annotations;

    if (trait != null) {
      traitIndex = pheCols[trait.value];
      if (traitIndex !== -1) {
        traceText = sampleChoice.map(s => {
          if (pheData[s][traitIndex] !== "") {
            return `${s}<br>Value: ${parseFloat(pheData[s][traitIndex]).toFixed(3)}`;
          } else {
            return s;
          }
        });
        trace.marker.color = sampleChoice.map(s => pheData[s][traitIndex]);
        trace.text = traceText;
        trace.marker.colorscale = [
          [0, "rgb(242,242,187)"],
          [0.5, "rgb(245,184,78)"],
          [1, "rgb(238,62,10)"]
        ];
        trace.marker.colorbar = {
          title: {
            text: "Phenotype",
            font: {
              size: 12,
            },
          },
          x: 1.1,
          xref: "paper",
        };
        trace.marker.showscale = true;
      }
    }

  } else if (plotChoice === "kin") {
    var zVals = data[0].flat();

    var sampleGroupChoice = [];
    document.querySelectorAll(".sample_group-checkbox").forEach( function(cbox) {
      if (cbox.checked) {
        var samplesList = JSON.parse(cbox.value);
        sampleGroupChoice.push(
          [
            samplesList,
            cbox.nextElementSibling.textContent]);
      }
    });
    var groupNum = sampleGroupChoice.length;
    var trait = document.querySelector("input[name='by_phe']:checked");

    var trace = {
      type: "heatmap",
      z: data[0],
      x: data[1],
      y: data[1],
      zmin: percentile(zVals, 1),
      zmax: percentile(zVals, 99),
      showscale: true,
      //colorscale: "Jet",
      colorbar: {
        title: {
          text: "Kinship",
          font: {
            size: 12,
          },
        },
      },
    };
    traces.push(trace);

    if (groupNum > 0) {
      const groupSets = sampleGroupChoice.map(g => new Set(g[0]));

      var groupMatrix = data[1].map(samp =>
          groupSets.map(group => group.has(samp) ? 1 : 0)
      );

      var groupTrace = {
        type: "heatmap",
        z: groupMatrix,
        x: sampleGroupChoice.map(g => g[1]),
        y: data[1],
        showscale: false,
        zmin: 0,
        zmax: 1,
        colorscale: [
          [0, "white"],
          [1, "red"]
        ],
        xaxis: "x2",
        toggleable: false,
      };
      traces.push(groupTrace);
    }
    if (trait != null) {

      traitIndex = pheCols[trait.value];
      if (traitIndex !== -1) {
        var traceMatrix = data[1].map(s => [pheData[s][traitIndex]]);

        var pheTrace = {
          type: "heatmap",
          z: traceMatrix,
          x: ["Phenotype"],
          y: data[1],
          showscale: true,
          colorscale: [
            [0, "rgb(242,242,187)"],
            [0.5, "rgb(245,184,78)"],
            [1, "rgb(238,62,10)"]
          ],
          xaxis: "x3",
          colorbar: {
            title: {
              text: "Phenotype",
              font: {
                size: 12,
              },
            },
            x: 1.15,
            xref: "paper",
          },
        };
        traces.push(pheTrace);
      }
    }

    var layout = {
      xaxis: {
        tickangle: -45,
        minallowed: -0.5,
        maxallowed: data[1].length - 0.5,
      },
      yaxis: {
        tickangle: -45,
        minallowed: -0.5,
        maxallowed: data[1].length - 0.5,
      },
      margin: {
        t: 30,
        pad: 4
      },
    }

    var numCols = groupNum + (trait == null ? 0 : 1);
    var offset = numCols / 50;
    if (groupNum > 0) {
      offset += 0.02;
    }
    if (trait != null) {
      offset += 0.02;
    }

    if (offset) {
      var xBound = Math.max(1 - offset, 0.5);
      layout.xaxis.domain = [0, xBound];

      if (groupNum > 0) {
        var x2Bound = trait == null ? 1 : 0.98 - (0.96 - xBound) / numCols;
        layout.xaxis2 = {
          domain: [xBound + 0.02, x2Bound],
          tickangle: -45,
          type: "category",
          tickmode: "array",
          tickvals: sampleGroupChoice.map(g => g[1]),
          ticktext: sampleGroupChoice.map(g => g[1]),
        }
      }
      if (trait != null) {
        var x3Bound = groupNum === 0 ? xBound : x2Bound;
        layout.xaxis3 = {
          domain: [x3Bound + 0.02, 1],
          tickangle: -45,
        }
      }
    }
  } else if (plotChoice === "bdist") {
    var chr = document.getElementById("chr-field").value;
    var xMax = sizes[names.indexOf(chr)];
    var yMax = Math.max(...chrData.map(v => v[2])) * 1.05;
    var trace = {
      type: "bar",
      x: chrData.map(v => (v[0] + v[1]) / 2),
      y: chrData.map(v => v[2]),
      textposition: "none",
      customdata: chrData.map(v => [v[0], v[1]].join("‒").replace(/\B(?=(\d{3})+(?!\d))/g, ",")),
      hovertemplate: "%{customdata}",
      hoverlabel: {
        namelength: 0,
        bgcolor: "white",
      },
      width: chrData.map(v => v[2]),
      marker: { color: "#60cdcd" },
      visible: true,
    };
    traces.push(trace);
    var layout = {
      xaxis: {
        title: "Position",
        showgrid: false,
        range: [0, xMax],
        minallowed: 0,
        maxallowed: xMax,
      },
      yaxis: {
        title: "Span",
        range: [0, yMax],
        minallowed: 0,
        maxallowed: yMax,
      },
      margin: {
        t: 30,
        pad: 4
      },
      hovermode: "x",
    }
  } else if (plotChoice === "ld") {
    var trace = {
      type: "heatmap",
      z: chrData[0],
      x: chrData[1],
      y: chrData[1],
      colorscale: [
       [0, "rgb(255,255,255)"],
       [1, "rgb(255,0,0)"]],
      zmin: 0,
      zmax: 1,
      autocolorscale: false,
      showscale: true,
      colorbar: {
        title: {
          text: "LD",
          font: {
            size: 12,
          },
        },
      },
    };
    traces.push(trace);
    var layout = {
      xaxis: {
        tickangle: -45,
        automargin: true,
      },
      yaxis: {
        tickangle: -45,
        automargin: true,
      },
      margin: {
        t: 30,
        pad: 0
      },
      plot_bgcolor: "rgb(150,150,150)",
      //title: "Random selection of 1,000 markers within chromosome",
      //plot_bgcolor: "rgb(0,0,0)",
    }
  } else if (plotChoice === "baf") {
    var [hgts, bins] = data;
    for (let i = 0; i < hgts.length; i++) {
      var trace = {
        type: "bar",
        name: `${ordinal(i+1)} haplotype`,
        x: bins,
        y: hgts[i],
        opacity: 0.4,
        textposition: "none",
        visible: true,
      }
      traces.push(trace);
    };
    var yMax = Math.max(...hgts.reduce((acc, hgt) => acc.concat(hgt), [])) + 1;
    var layout = {
      barmode: "overlay",
      bargap: 0,
      xaxis: {
        title: "Value",
        showgrid: false,
        range: [0, 1],
        minallowed: 0,
        maxallowed: 1,
      },
      yaxis: {
        title: "Count",
        //range: [0, yMax],
        minallowed: 0,
        //maxallowed: yMax,
      },
      showlegend: true,
    }
  } else {
    var txt = data[0].map((val, idx) => { return `${data[1][idx]}-${data[1][idx+1]}` });
    var xMax = Math.max(...data[1]);
    var yMax = Math.max(...data[0]) * 1.05;
    var trace = {
      type: "bar",
      x: data[1],
      y: data[0],
      text: txt,
      textposition: "none",
      hovertemplate: "%{text}<br>Count: %{y}",
      hoverlabel: {
        namelength: 0,
        bgcolor: "white",
//        font: {
//          size: 13,
//        },
      },
      marker: { color: "#60cdcd" },
      offset: 0,
      visible: true,
    };

    if (data[0].length === 1) {
      xMax = 0.1;
      trace.width = 0.01;
    }

    traces.push(trace);
    var layout = {
      xaxis: {
        title: "Value",
        range: [0, xMax],
        minallowed: 0,
        maxallowed: xMax,
      },
      yaxis: {
        title: "Count",
        range: [0, yMax],
        minallowed: 0,
        maxallowed: yMax,
      },
      hovermode: "x",
      margin: {
        l: 80,
        r: 40,
        b: 50,
        t: 70,
        pad: 4
      },
      bargap: 0,
    };

    var thresh = document.getElementById(plotChoice + "-thresh");
    if (thresh != null) {
      thresh = thresh.value;
      layout.shapes = [{
        type: "line",
        x0: thresh,
        x1: thresh,
        y0: 0,
        y1: 1,
        yref: "paper",
        line: {
          color: "red",
          width: 2
        }
      }];
    }
  }

  var config = {
    editable: false,
    modeBarButtons: [
      ["toImage", buttonsToAdd["save_svg"], buttonsToAdd["save_html"]],
      ["zoom2d", "pan2d", "select2d", "lasso2d", buttonsToAdd["clear_sel"]],
      ["zoomIn2d", "zoomOut2d", "autoScale2d", "resetScale2d"]],
    modeBarButtonsToRemove: [],
  };
  if (plotChoice === "kin") {
    config.modeBarButtons.push(
      ["edit", "legend"].map(key => buttonsToAdd[key])
    );
  } else {
    config.modeBarButtons.push(
      ["edit", "recolor", "legend"].map(key => buttonsToAdd[key])
    );
  }
  if (plotChoice === "dens") {
    //config.modeBarButtons[3].push(buttonsToAdd["reshade"]);
  } else if (plotChoice === "pca") {
    config.modeBarButtons[3].push(buttonsToAdd["opacity"]);
    config.modeBarButtons[3].push(buttonsToAdd["pve"]);
    config.modeBarButtons[3].push(buttonsToAdd["m_size"]);
    config.modeBarButtons[3].push(buttonsToAdd["m_symbol"]);
    config.modeBarButtons[3].push(buttonsToAdd["move_down"]);
    config.modeBarButtons[3].push(buttonsToAdd["_delete"]);
  } else if (plotChoice === "baf") {
    config.modeBarButtons[3].push(buttonsToAdd["opacity"]);
    config.modeBarButtons[3].push(buttonsToAdd["group"]);
    config.modeBarButtons[3].push(buttonsToAdd["move_down"]);
    config.modeBarButtons[3].push(buttonsToAdd["_delete"]);
  }

  document.getElementById("plot-div").innerHTML = "";
  //Plotly.purge("plot-div");
  Plotly.newPlot("plot-div", traces, layout, config);
};

function resizePlot() {
  plotDiv = document.getElementById("plot-div");
  if (plotDiv.innerHTML) {
    Plotly.Plots.resize("plot-div");
  }
}

function fetchPlots(markers, isHap, callback) {
  $.ajax({
    url: "/hist",
    method: "POST",
    data: {
      "genome": genome,
      "markers": markers,
      "hap": isHap,
    },
    success: function (results) {
      callback(results);
    }
  });
}

function ordinal(i) {
    let j = i % 10, k = i % 100;
    if (j == 1 && k != 11) {
        return i + "st";
    }
    if (j == 2 && k != 12) {
        return i + "nd";
    }
    if (j == 3 && k != 13) {
        return i + "rd";
    }
    return i + "th";
}



//////////////////////////////////////////////////////////////////////////////
// Unused table functions
//////////////////////////////////////////////////////////////////////////////
//
//var xLabels  = {
// "smiss": "Missingness",
// "vmiss": "Missingness",
// "shet":  "Heterozygosity",
// "vhet":  "Heterozygosity",
// "maf":   "Minor AF",
// "hwe":   "-log10(P) HWE"}
//
//function getMaxRows(plotChoice) {
//  if ((plotChoice == "smiss") | (plotChoice == "shet")) {
//    return sampNum;
//  }
//  else {
//    return varNum;
//  }
//};
//
// In/decrement table
//function increment() {
//  var plotChoice = document.getElementById("choice-field").value;
//  var maxRows    = getMaxRows(plotChoice);
//  var tableFrom  = document.getElementById("table-from");
//  var tableTo    = document.getElementById("table-to");
//  var fromValue  = parseInt(tableFrom.innerText, 10);
//  var toValue    = parseInt(tableTo.innerText, 10);
//
//  if (toValue + 10 > maxRows) {
//    tableFrom.innerText = Math.max(1, maxRows - 10);
//    tableTo.innerText = maxRows;
//  }
//  else {
//    tableFrom.innerText = fromValue + 10;
//    tableTo.innerText = toValue + 10;
//  }
//  if (fromValue != parseInt(tableFrom.innerText, 10)) {
//    fillTable(plotChoice);
//  }
//}
//
//function decrement() {
//  var plotChoice = document.getElementById("choice-field").value;
//  var maxRows    = getMaxRows(plotChoice);
//  var tableFrom  = document.getElementById("table-from");
//  var tableTo    = document.getElementById("table-to");
//  var fromValue  = parseInt(tableFrom.innerText, 10);
//  var toValue    = parseInt(tableTo.innerText, 10);
//
//  if (fromValue - 10 < 1) {
//    tableFrom.innerText = 1;
//    tableTo.innerText = Math.min(10, maxRows);
//  }
//  else {
//    tableFrom.innerHTML = fromValue - 10;
//    tableTo.innerHTML = toValue - 10;
//  }
//  if (fromValue != parseInt(tableFrom.innerText, 10)) {
//    fillTable(plotChoice);
//  }
//}
//
// Fill table
//function fillTable(plotChoice) {
//  var fromVal = document.getElementById("table-from").innerText;
//  var toVal = document.getElementById("table-to").innerText;
//
//  $.ajax({
//    url: "/filt",
//    method: "POST",
//    data: {
//     genome: genome,
//     markers: markers,
//     metric: plotChoice,
//     start_row: fromVal,
//     end_row: toVal },
//    success: function (data) {
//      var table = document.getElementById("plot-table");
//      table.innerHTML = "";
//      if ((plotChoice == "smiss") | (plotChoice == "shet")) {
//        var idColumn = "Sample";
//      }
//      else {
//        var idColumn = "Variant";
//      }
//      var cols = [idColumn, xLabels[plotChoice]];
//
//      var thead = document.createElement("thead");
//      var tr = document.createElement("tr");
//      cols.forEach(function (col) {
//        var th = document.createElement("th");
//        var tdDiv = document.createElement("div");
//        tdDiv.classList.add("disp-flex");
//        var textDiv = document.createElement("div");
//        textDiv.innerHTML = col;
//        textDiv.classList.add("flex");
//        var resizDiv = document.createElement("div");
//        resizDiv.classList.add("resizable", "resiz-render");
//        tdDiv.appendChild(textDiv);
//        tdDiv.appendChild(resizDiv);
//        th.appendChild(tdDiv);
//        tr.appendChild(th);
//      });
//      thead.appendChild(tr);
//
//      var tbody = document.createElement("tbody");
//      data.forEach(function (row) {
//        var tr = document.createElement("tr");
//        row.forEach(function (cell) {
//          var td = document.createElement("td");
//          td.textContent = cell;
//          td.classList.add("qc-col");
//          tr.appendChild(td);
//        });
//        tbody.appendChild(tr);
//      });
//      table.appendChild(thead);
//      table.appendChild(tbody);
//    },
//    error: function () {
//      alert("Failed to fetch table.");
//    }
//  });
//};
//
//  //document.getElementById("table-from").innerHTML = 1;
//  //document.getElementById("table-to").innerHTML = 10;
//
//  var tableDiv = document.getElementById("table-div");
//  if (["dens", "pca", "kin", "ld"].includes(plotChoice)) {
//    tableDiv.classList.add("collapsed");
//  } else {
//    tableDiv.classList.remove("collapsed");
//    fillTable(plotChoice);
//  }
//
//////////////////////////////////////////////////////////////////////////////
// Unused PCA scree plot functions
//////////////////////////////////////////////////////////////////////////////
//
//    var pcs = Array.from({ length: data.length }, (_, i) => i + 1);
//    var trace = {
//      type: "scatter",
//      x: pcs,
//      y: data,
//      mode: "lines+markers",
//      line: {
//        color: "#60cdcd",
//      },
//      marker: {
//        color: "#014e60",
//        size: 10
//      },
//    };
//    traces.push(trace);
//    var layout = {
//      title: "Scree plot",
//      xaxis: {
//        title: "Principal component",
//        tickvals: pcs,
//        ticktext: pcs,
//      },
//      yaxis: { 
//        title: "Proportion of variance explained",
//        //range: [0, 100],
//      },
//      margin: {
//        t: 70,
//        pad: 4
//      },
//    }
//
//////////////////////////////////////////////////////////////////////////////
// Unused plot titles
//////////////////////////////////////////////////////////////////////////////
//
//var pTitles = {
// "smiss": "Sample missingness",
// "vmiss": "Marker missingness",
// "shet":  "Sample heterozygosity",
// "vhet":  "Marker heterozygosity",
// "maf":   "Minor allele frequency",
// "hwe":   "Hardy-Weinberg equilibrium",
// "dens":  "Marker density"}
//
//      title: {
//        text: pTitles[plotChoice],
//        x: 0.5,
//      },
//