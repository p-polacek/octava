function fetchChr() {
  if (typeof(haps) !== "undefined" && haps[markers]) {
    var plotChoice = document.getElementById("hap-choice-field").value;
  } else {
    var plotChoice = document.getElementById("choice-field").value;
  }
  var chrDiv = document.getElementById("chr-div");
  var pcaSampDiv = document.getElementById("pca-sample-div");
  if (plotChoice === "ld") {
    chrDiv.classList.remove("collapsed");
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
    pcaSampDiv.classList.add("collapsed");
    if (plotChoice === "pca") {
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
      marker: { color: "#47b023" },
      visible: true,
    };
    traces.push(trace);

    var lastChrIdx = chrNames.length - 1;
    var yMax = Math.max(...data[0]) * 1.05;
    var layout = {
      //plot_bgcolor: "rgba(0, 0, 0, 0)",
      //paper_bgcolor: "rgba(0, 0, 0, 0)",
      xaxis: {
        tickvals: chrTicks,
        ticktext: chrNames,
        showgrid: false,
        range: [0, xMax],
        minallowed: -10000000,
        maxallowed: data[1][data[1].length - 1] + mdbin / 2,
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
    document.querySelectorAll(".samples_group-checkbox").forEach( function(cbox) {
      if (cbox.checked) {
        // Keep both group contents and name
        var samplesList = JSON.parse(cbox.value);
        sampleGroupChoice.push(
          [
            samplesList,
            cbox.nextElementSibling.textContent]);
      }
    });
    var legend = (sampleGroupChoice.length > 0);
    var opacity = legend ? 0.5 : 1;

    var trace = {
      x: PCAdata[0],
      y: PCAdata[1],
      z: PCAdata[2],
      name: "<i>All samples</i>",
      type: "scatter3d",
      mode: "markers",
      text: plottedSamples,
      opacity: opacity,
      showlegend: legend,
      hovertemplate: "%{text}",
      hoverlabel: {
        namelength: 0,
        bgcolor: "white",
      },
      marker: {
         size: 4,
      },
      visible: true,
    }
    traces.push(trace);

    // Plot of sample groups
    if (sampleGroupChoice.length > 0) {
      for (var [idx, group] of sampleGroupChoice.entries()) {
        var [members, groupName] = group;
        var trace = {
          x: PCAdata[0].filter((_, idx) => members.includes(plottedSamples[idx])),
          y: PCAdata[1].filter((_, idx) => members.includes(plottedSamples[idx])),
          z: PCAdata[2].filter((_, idx) => members.includes(plottedSamples[idx])),
          name: groupName,
          type: "scatter3d",
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
             color: plotColors[idx],
             size: 4,
          },
          visible: true,
        }
        traces.push(trace);
      }
    }

    var annotations = [];
    annotations.push({
      text: "PVE:",
      x: 1,
      y: 1,
      showarrow: false,
      yanchor: "middle",
      valign: "top",
      xref: "paper",
      yref: "paper",
      font: { family: "monospace", size: 14 },
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
        y: 1 - (idx * Math.min(0.05, (1 / PVEs.length))),
        showarrow: false,
        yanchor: "middle",
        valign: "top",
        xref: "paper",
        yref: "paper",
        font: { family: "monospace", size: 14 },
      });
    });
    var layout = {
      annotations: annotations,
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
          title: "PC1",
          showgrid: false,
        },
        yaxis: {
          title: "PC2",
          showgrid: false,
         },
        zaxis: {
          title: "PC3",
          showgrid: false,
        },
      },
      margin: {
        t: 30,
        b: 10,
        l: 10,
        r: 10,
        pad: 4
      },
    }
  } else if (plotChoice === "kin") {
    var trace = {
      type: "heatmap",
      z: data[0],
      x: data[1],
      y: data[1],
      colorbar: {
        title: {
          text: "Hamming<br>distance",
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
      },
      yaxis: {
        tickangle: -45,
      },
      margin: {
        t: 30,
        pad: 4
      },
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
      width: chrData.map(v => v[2]),
      marker: { color: "#47b023" },
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
        range: [0, yMax],
        minallowed: 0,
        maxallowed: yMax,
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
      marker: { color: "#47b023" },
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
      hovermode: "closest",
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
    modeBarButtonsToAdd: 
      ["edit", "recolor"]
      .map(key => buttonsToAdd[key]),
    modeBarButtonsToRemove: [],
  };
  if (["pca", "baf"].includes(plotChoice)) {
    config.modeBarButtonsToAdd.push(buttonsToAdd["legend"]);
    config.modeBarButtonsToAdd.push(buttonsToAdd["move_down"]);
    config.modeBarButtonsToAdd.push(buttonsToAdd["_delete"]);
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
//        color: "#47b023",
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