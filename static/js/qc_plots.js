// Static data
var pTitles = {
 "smiss": "Sample missingness",
 "vmiss": "Marker missingness",
 "shet":  "Sample heterozygosity",
 "vhet":  "Marker heterozygosity",
 "maf":   "Minor allele frequency",
 "hwe":   "Hardy-Weinberg equilibrium",
 "dens":  "Marker density"}
var xLabels  = {
 "smiss": "Missingness",
 "vmiss": "Missingness",
 "shet":  "Heterozygosity",
 "vhet":  "Heterozygosity",
 "maf":   "Minor AF",
 "hwe":   "-log10(P) HWE"}

// Helper functions
function checkName(event) {
  var setName = document.getElementById("name").value.trim();
  if (setName == "") {
    alert("Set name is required.");
    event.preventDefault();
  } else if (savedSets.includes(setName)) {
    alert("Set with this name already exists.");
    event.preventDefault();
  }
}

function getMaxRows(plotChoice) {
  if ((plotChoice == "smiss") | (plotChoice == "shet")) {
    return sampNum;
  }
  else {
    return varNum;
  }
};

function fetchChrLd() {
var plotChoice = document.getElementById("choice-field")
  var ldChrDiv = document.getElementById("ld-chr-div");
  if (plotChoice.value === "ld") {
    ldChrDiv.classList.remove("collapsed");
    var ldChr = document.getElementById("ld-chr-field").value;
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
        chrom: ldChr,
      },
      success: function (data) {
        chrLdData = data;
      }
    });
    requests.push(request);

    $.when.apply($, requests).done(function() {
      showPlot();
    })
  } else {
    ldChrDiv.classList.add("collapsed");
    showPlot();
  }
}

// In/decrement table
function increment() {
  var plotChoice = document.getElementById("choice-field").value;
  var maxRows    = getMaxRows(plotChoice);
  var tableFrom  = document.getElementById("table-from");
  var tableTo    = document.getElementById("table-to");
  var fromValue  = parseInt(tableFrom.innerText, 10);
  var toValue    = parseInt(tableTo.innerText, 10);

  if (toValue + 10 > maxRows) {
    tableFrom.innerText = Math.max(1, maxRows - 10);
    tableTo.innerText = maxRows;
  }
  else {
    tableFrom.innerText = fromValue + 10;
    tableTo.innerText = toValue + 10;
  }
  if (fromValue != parseInt(tableFrom.innerText, 10)) {
    fillTable(plotChoice);
  }
}

function decrement() {
  var plotChoice = document.getElementById("choice-field").value;
  var maxRows    = getMaxRows(plotChoice);
  var tableFrom  = document.getElementById("table-from");
  var tableTo    = document.getElementById("table-to");
  var fromValue  = parseInt(tableFrom.innerText, 10);
  var toValue    = parseInt(tableTo.innerText, 10);

  if (fromValue - 10 < 1) {
    tableFrom.innerText = 1;
    tableTo.innerText = Math.min(10, maxRows);
  }
  else {
    tableFrom.innerHTML = fromValue - 10;
    tableTo.innerHTML = toValue - 10;
  }
  if (fromValue != parseInt(tableFrom.innerText, 10)) {
    fillTable(plotChoice);
  }
}

// Fill table
function fillTable(plotChoice) {
  var fromVal = document.getElementById("table-from").innerText;
  var toVal = document.getElementById("table-to").innerText;

  $.ajax({
    url: "/filt",
    method: "POST",
    data: {
     genome: genome,
     metric: plotChoice,
     start_row: fromVal,
     end_row: toVal },
    success: function (data) {
      var table = document.getElementById("plot-table");
      table.innerHTML = "";
      if ((plotChoice == "smiss") | (plotChoice == "shet")) {
    var idColumn = "Sample";
      }
      else {
    var idColumn = "Variant";
      }
      var cols = [idColumn, xLabels[plotChoice]];

      var thead = document.createElement("thead");
      var tr = document.createElement("tr");
      cols.forEach(function (col) {
        var th = document.createElement("th");
        var tdDiv = document.createElement("div");
        tdDiv.classList.add("disp-flex");
        var textDiv = document.createElement("div");
        textDiv.innerHTML = col;
        textDiv.classList.add("flex");
        var resizDiv = document.createElement("div");
        resizDiv.classList.add("resizable", "resiz-render");
        tdDiv.appendChild(textDiv);
        tdDiv.appendChild(resizDiv);
        th.appendChild(tdDiv);
        tr.appendChild(th);
      });
      thead.appendChild(tr);

      var tbody = document.createElement("tbody");
      data.forEach(function (row) {
        var tr = document.createElement("tr");
        row.forEach(function (cell) {
          var td = document.createElement("td");
          td.textContent = cell;
          td.classList.add("qc-col");
          tr.appendChild(td);
        });
        tbody.appendChild(tr);
      });
      table.appendChild(thead);
      table.appendChild(tbody);
    },
    error: function () {
      alert("Failed to fetch table.");
    }
  });
};

// Plot
function showPlot() {
  var plotChoice = document.getElementById("choice-field").value;
  if (histData == null) {
    return;
  }
  var data = histData[plotChoice];

  if (plotChoice === "dens") {
    if (visNames.length > 0) {
      var visIdxs = chrNames.filter( n => {
        return visNames.includes(n);
      }).map( n => {
        return chrNames.indexOf(n);
      });
      var visSizes = visIdxs.map( i => {
        return Math.floor(chrSizes[i] / 1000000) * 1000000 + 1000000;
      });
      var visOffsets = visIdxs.map( i => {
        return Math.floor(chrOffsets[i] / 1000000) * 1000000 + 1000000;
      });
      var visTicks = visOffsets.map( (o, i) => {
        return visOffsets[i] + Math.floor(visSizes[i] / 2);
      });
      visNames.push("Other");
      visOffsets.push(
       (visOffsets[
        visOffsets.length - 1]
        + visSizes[
         visSizes.length - 1]));
      visSizes.push(
       Math.max(
       Math.floor(
        chrSizes[chrSizes.length - 1] / 1000000)
        * 1000000
        + 1000000,
       visOffsets[visOffsets.length - 1] + 1000000));

      visTicks.push(
       visOffsets[visOffsets.length - 1]
       + Math.floor(visSizes[visSizes.length - 1] / 2));
    } else {
      visNames       = chrNames;
      var visTicks   = chrTicks;
      var visOffsets = chrOffsets;
      var visSizes   = chrSizes;
    }

    var trace = {
      type: "bar",
      x: data[1],
      y: data[0],
      text: data[2],
      hoverinfo: "text+y",
      hovertemplate: "%{text}",
      marker: { color: "#47b023" },
    };
  
    var layout = {
      title: {
        text: pTitles[plotChoice],
        x: 0.5,
      },
      plot_bgcolor: "rgba(0, 0, 0, 0)",
      paper_bgcolor: "rgba(0, 0, 0, 0)",
      xaxis: {
        tickvals: visTicks,
        ticktext: visNames,
        gridcolor: "rgba(0,0,0,0.1)",
      },
      yaxis: {
        title: "Count",
        gridcolor: "rgba(0,0,0,0.1)",
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
    for (i = 0; i < visNames.length; i+= 2) {
      var shape = {
        type: "rect",
        xref: "x",
        yref: "paper",
        x0: visOffsets[i],
        y0: 0,
        x1: visOffsets[i] + visSizes[i],
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
    var pcs = Array.from({ length: data.length }, (_, i) => i + 1);
    var trace = {
      type: "scatter",
      x: pcs,
      y: data,
      mode: "lines+markers",
      line: {
        color: "#47b023",
      },
      marker: {
        color: "#014e60",
        size: 10
      },
    };
    var layout = {
      title: "Scree plot",
      xaxis: {
        title: "Principal component",
        tickvals: pcs,
        ticktext: pcs,
      },
      yaxis: { 
        title: "Proportion of variance explained",
        //range: [0, 100],
      },
      margin: {
        t: 70,
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
        title: "Distance",
      },
    };
    var layout = {
      xaxis: {
        tickangle: -45,
      },
      yaxis: {
        tickangle: -45,
      },
      margin: {
        t: 10,
        pad: 4
      },
    }
  } else if (plotChoice === "ld") {
    var ldChr = document.getElementById("ld-chr-field").value;
    var trace = {
      type: "heatmap",
      z: chrLdData[0],
      x: chrLdData[1],
      y: chrLdData[1],
      colorscale: [
       [0, "rgb(255,255,255)"],
       [1, "rgb(255,0,0)"]],
      colorbar: {
        title: "LD",
      },
    };
    var layout = {
      xaxis: {
        tickangle: -45,
      },
      yaxis: {
        tickangle: -45,
      },
      margin: {
        t: 10,
        pad: 0
      },
      plot_bgcolor: "rgb(150,150,150)",
      //plot_bgcolor: "rgb(0,0,0)",
    }
  } else {
    var trace = {
      type: "bar",
      x: data[1],
      y: data[0],
      marker: { color: "#47b023" },
    };

    var layout = {
      title: {
        text: pTitles[plotChoice],
        x: 0.5,
      },
      plot_bgcolor: "rgba(0, 0, 0, 0)",
      paper_bgcolor: "rgba(0, 0, 0, 0)",
      xaxis: {
        title: xLabels[plotChoice],
        gridcolor: "rgba(0,0,0,0.1)",
      },
      yaxis: {
        title: "Count",
        gridcolor: "rgba(0,0,0,0.1)",
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

  document.getElementById("table-from").innerHTML = 1;
  document.getElementById("table-to").innerHTML = 10;

  var tableDiv = document.getElementById("table-div");
  if (["dens", "pca", "kin", "ld"].includes(plotChoice)) {
    tableDiv.classList.add("collapsed");
  } else {
    tableDiv.classList.remove("collapsed");
    fillTable(plotChoice);
  }

  document.getElementById("plot-div").innerHTML = "";
  //Plotly.purge("plot-div");
  Plotly.newPlot("plot-div", [trace], layout);
};

function resizePlot() {
  plotDiv = document.getElementById("plot-div");
  if (plotDiv.innerHTML) {
    Plotly.Plots.resize("plot-div");
  }
}

function fetchPlots(markers, callback) {
  payload = {};
  payload["genome"] = genome;
  payload["markers"] = markers

  $.ajax({
    url: "/hist",
    method: "POST",
    data:   payload,
    success: function (results) {
      callback(results);
    }
  });
}