var customOrder = [];
var customColors = {};

// Buttons:
// clear_sel    # Clear selection
// legend       # Toggle legend
// pve          # Toggle PVE
// edit         # Edit plot
// recolor      # Recolor visible traces
// reshade      # Swap chromosome shading
// opacity      # Change transparency of visible traces
// reset_color  # Reset colors of visible traces
// move_down    # Move visible traces down
// group        # Combine visible traces
// group_mp     # Combine visible traces
// _delete      # Delete visible traces


var buttonsToAdd = {
  clear_sel: {
    name: "Clear selection",
    icon: Plotly.Icons.eraseshape,
    click: function(gd) {
      Plotly.update(gd, { selectedpoints: null }, { selections: [] });
    }
  },
  legend: {
    name: "Toggle legend",
    icon: {
      "width": 100,
      "height": 100,
      "path": "M 12 14 A 10 10 0 0 0 2 25 A 10 10 0 0 0 12 35 A 10 10 0 0 0 23 25 A 10 10 0 0 0 12 14 z M 39 20 C 36 20 34 22 34 25 C 34 28 36 30 39 30 L 93 30 C 96 30 98 28 98 25 C 98 22 96 20 93 20 L 39 20 z M 12 45 A 10 10 0 0 0 2 56 A 10 10 0 0 0 12 66 A 10 10 0 0 0 23 56 A 10 10 0 0 0 12 45 z M 39 51 C 36 51 34 53 34 56 C 34 59 36 61 39 61 L 93 61 C 96 61 98 59 98 56 C 98 53 96 51 93 51 L 39 51 z M 12 76 A 10 10 0 0 0 2 87 A 10 10 0 0 0 12 97 A 10 10 0 0 0 23 87 A 10 10 0 0 0 12 76 z M 39 82 C 36 82 34 84 34 87 C 34 90 36 92 39 92 L 93 92 C 96 92 98 90 98 87 C 98 84 96 82 93 82 L 39 82 z ",
    },
    click: function(gd) {
      Plotly.relayout(gd, { showlegend: !gd.layout.showlegend });
      var colorBar = gd.data.filter(t => t.showscale != null);
      var colorMarkerBar = gd.data.filter(t => t.marker != null && t.marker.showscale != null);
      if (colorBar.length > 0) {
        Plotly.restyle(gd, { showscale: [!colorBar[0].showscale] });
      }
      if (colorMarkerBar.length > 0) {
        Plotly.restyle(gd, { "marker.showscale": [!colorMarkerBar[0].marker.showscale] });
      }
    }
  },
  pve: {
    name: "Toggle PVE",
    icon: {
      "width": 100,
      "height": 100,
      "path": "M 28 16 A 11 11 0 0 0 16 28 A 11 11 0 0 0 28 40 A 11 11 0 0 0 40 28 A 11 11 0 0 0 28 16 z M 75 16 A 7 7 0 0 0 70 19 L 19 70 A 7 7 0 0 0 19 80 A 7 7 0 0 0 29 80 L 80 29 A 7 7 0 0 0 80 19 A 7 7 0 0 0 75 16 z M 71 59 A 11 11 0 0 0 59 71 A 11 11 0 0 0 71 83 A 11 11 0 0 0 83 71 A 11 11 0 0 0 71 59 z",
    },
    click: function(gd) {
      Plotly.relayout(gd, { annotations: gd.layout.annotations.map(a => ({ ...a, visible: !a.visible })) });
    }
  },
  edit: {
    name: "Edit plot",
    icon: Plotly.Icons.pencil,
    click: function(gd) {
      Plotly.newPlot(
        gd.id,
        gd.data,
        gd.layout,
        { ...gd._context, editable: !gd._context.editable, edits: {shapePosition: false} });
      if (typeof setPlotEventListeners === "function") {
        setPlotEventListeners();
      }
    }
  },
  recolor: {
    name: "Recolor visible traces",
    icon: {
      "width": 100,
      "height": 100,
      "path": "m 11,11 c -3,0 -5,2 -5,5 v 19 c 0,3 2,5 5,5 H 68 c 3,0 5,-2 5,-5 v -5 h 7 V 48 L 41,62 c -1,0 -2,2 -2,4 v 25 c 3e-4,2 1,4 4,4 2,-3e-4 4,-1 4,-4 V 70 l 40,-13 c 1,-0 2,-2 2,-4 V 26 C 91,24 89,22 86,22 H 74 v -5 c 0,-3 -2,-5 -5,-5 z",
    },
    click: function(gd) {
      if (gd.data.length === 0) {
        return;
      }
      var clickEvent = event;
      var traceColors;
      var traceIdxs = [];
      var traceIds = [];
      gd.data.forEach(function(trace, idx) {
        if (trace.visible === true) {
          traceIdxs.push(idx);
          traceIds.push(trace.Id || "");
        }
      });

      if (traceIdxs.length === 0) {
        return;
      }

      var currColor;
      if (gd.data[traceIdxs[0]].marker != null && gd.data[traceIdxs[0]].marker.color != null) {
        currColor = gd.data[traceIdxs[0]].marker.color;
      } else if (gd.data[traceIdxs[0]].line != null && gd.data[traceIdxs[0]].line.color != null) {
        currColor = gd.data[traceIdxs[0]].line.color;
      } else {
        currColor = "#ff0000";
      }
      // If hex value includes alpha
      var isAlpha = false;
      if (currColor.length === 9) {
        isAlpha = true;
        currColor = currColor.slice(0, 7);
      }
      showHelp(clickEvent, null, `<input type="color" id="colorPicker" value="${currColor}" />`);

      document.getElementById("colorPicker").click();
      document.getElementById("help-popup").classList.add("collapsed");
      document.getElementById("colorPicker").addEventListener("input", function() {
        var colorChoice = this.value;

        for (let i = 0; i < traceIdxs.length; i++) {
          var recolor = {
            "marker.color": [colorChoice],
            "line.color": [colorChoice],
          }
          if (gd.data[traceIdxs[i]].type !== "box") {
            recolor["fillcolor"] = [colorChoice];
          }
          if (isAlpha) {
            recolor["marker.color"] = [colorChoice + "80"];
            recolor["marker.line.color"] = [colorChoice];
          }

          Plotly.restyle(
            gd,
            recolor,
            [traceIdxs[i]]);
        }
      });
    },
  },
  reshade: {
    name: "Swap chromosome shading",
    icon: {
      "width": 100,
      "height": 100,
      "path": "M 10,9 A 6,6 0 0 0 3,15 v 68 a 6,6 0 0 0 6,6 h 79 a 6,6 0 0 0 6,-6 V 15 A 6,6 0 0 0 89,9 Z M 16,21 H 83 V 78 H 16 Z m 11,3 a 6,6 0 0 0 -6,6 v 36 a 6,6 0 0 0 6,6 6,6 0 0 0 6,-6 V 31 a 6,6 0 0 0 -6,-6 z m 17,0 a 6,6 0 0 0 -6,6 v 36 a 6,6 0 0 0 6,6 6,6 0 0 0 6,-6 V 31 a 6,6 0 0 0 -6,-6 z",
    },
    click: function(gd) {
      var shapes = gd.layout.shapes;
      for (shape of shapes) {
        if (shape.shapeType === "chrom") {
          shape.visible = !shape.visible;
        }
      }
      Plotly.relayout("live-plot", {shapes: shapes});
    }
  },
  opacity: {
    name: "Change transparency of visible traces",
    icon: {
      "width": 100,
      "height": 100,
      "path": "M 50 6 C 26 6 6 26 6 50 C 6 73 26 93 50 93 C 73 93 93 73 93 50 C 93 26 73 6 50 6 z M 50 19 C 66 19 80 33 80 50 C 80 66 66 80 50 80 C 33 80 19 66 19 50 C 19 33 33 19 50 19 z M 60 28 A 5 5 0 0 0 56 30 L 30 56 A 5 5 0 0 0 30 64 A 5 5 0 0 0 38 64 L 64 37 A 5 5 0 0 0 64 30 A 5 5 0 0 0 60 28 z M 66 49 A 5 5 0 0 0 62 51 L 51 63 A 5 5 0 0 0 51 70 A 5 5 0 0 0 59 70 L 70 59 A 5 5 0 0 0 70 51 A 5 5 0 0 0 66 49 z",
    },
    click: function(gd) {
      var clickEvent = event;
      var traceOpacities = [];
      var traceIdxs = [];
      var traceIds = [];
      gd.data.forEach(function(trace, idx) {
        if (trace.visible === true) {
          traceIdxs.push(idx);
          traceIds.push(trace.Id || "");
          traceOpacities.push(trace.opacity || 1);
        }
      });
      var opacityValue = 1;
      if ([...new Set(traceOpacities)].length === 1) {
        opacityValue = traceOpacities[0];
      }

      showHelp(clickEvent, null, `<input type="range" id="opacitySlider" min="0" max="1" step="0.01" value="${opacityValue}" />`);
      document.getElementById("opacitySlider").addEventListener("input", function() {
        var opacityChoice = this.value;
        Plotly.restyle(
          gd,
          {opacity: opacityChoice,
           line: {
             opacity: opacityChoice},
          },
          traceIdxs);
      });
    },
  },
  reset_color: {
    name: "Reset colors of visible traces",
    icon: {
      "width": 100,
      "height": 100,
      "path": "m 51,8 c -10,0 -20,4 -27,12 a 6,6 0 0 0 0,9 6,6 0 0 0 9,-0 c 7,-8 20,-10 30,-3 9,6 12,19 7,29 -4,9 -15,13 -25,11 L 54,55 11,54 31,91 39,79 C 56,85 75,78 84,62 93,46 88,25 72,14 67,10 60,8 53,8 c -0,-0 -1,-0 -2,-0 z",
    },
    click: function(gd) {
      var traceIds = [];
      var traceIdxs = [];
      var traceColors = [];

      gd.data.forEach(function(trace, idx) {
        if (trace.visible === true) {
          traceIdxs.push(idx);
          delete customColors[trace.Id || ""];
          traceColors.push(colorCycle[idx % 10]);
        }
      });

      Plotly.restyle(
        gd, {
          "marker.color": traceColors,
          "line.color": traceColors},
        traceIdxs);
    },
  },
  move_down: {
    name: "Move visible traces down",
    icon: {
      "width": 100,
      "height": 100,
      "path": "m 48,12 c -2,10e-7 -5,2 -5,5 V 61 L 33,50 c -1,-1 -5,-1 -7,0 l -1,1 c -1,1 -1,5 -10e-7,7 L 44,80 c 0,0 2,1 3,1 0,0 0,-0 0,-0 0,0 0,0 0,0 H 47 13 c -2,0 -5,2 -5,5 l 0,2 c -7e-7,2 2,5 5,5 l 72,2e-6 c 2,10e-7 5,-2 5,-5 l -1e-6,-2 c 0,-2 -2,-5 -5,-5 l -34,2e-6 -0,-1e-6 c 0,1e-6 0,-0 0,-0 0,0 0,0 0,0 1,0 2,-0 3,-1 L 75,59 c 1,-1 1,-5 10e-7,-7 l -1,-1 c -1,-1 -5,-1 -7,10e-7 l -10,10 0,-43 c -10e-7,-2 -2,-5 -5,-5 z",
    },
    click: function(gd) {
      var visibleIdxs = [];
      gd.data.forEach( function(trace, idx) {
        if (trace.visible === true) {
          visibleIdxs.push(idx);
        }
      });
      Plotly.moveTraces(gd, visibleIdxs);

      // Save custom order
      var traceIds = gd.data.map( t => t.Id);
      customOrder = traceIds.concat(customOrder.filter( i => !traceIds.includes(i) ));
    },
  },
  group: {
    name: "Combine visible traces",
    icon: {
      "width": 100,
      "height": 100,
      "path": "M 8 9 A 6 6 0 0 0 1 15 L 1 84 A 6 6 0 0 0 8 90 L 28 90 A 6 6 0 0 0 35 84 A 6 6 0 0 0 28 78 L 14 78 L 14 21 L 28 21 A 6 6 0 0 0 35 15 A 6 6 0 0 0 28 9 L 8 9 z M 71 9 A 6 6 0 0 0 64 15 A 6 6 0 0 0 71 21 L 85 21 L 85 78 L 71 78 A 6 6 0 0 0 64 84 A 6 6 0 0 0 71 90 L 91 90 A 6 6 0 0 0 98 84 L 98 15 A 6 6 0 0 0 91 9 L 71 9 z M 50 27 A 6 6 0 0 0 43 33 L 43 43 L 33 43 A 6 6 0 0 0 27 50 A 6 6 0 0 0 33 56 L 43 56 L 43 66 A 6 6 0 0 0 50 72 A 6 6 0 0 0 56 66 L 56 56 L 66 56 A 6 6 0 0 0 72 50 A 6 6 0 0 0 66 43 L 56 43 L 56 33 A 6 6 0 0 0 50 27 z",
    },
    click: function(gd) {
      var visibleIdxs = [];
      var visibleData = [];
      gd.data.forEach( function(trace, idx) {
        if (trace.visible === true) {
          visibleIdxs.push(idx);
          visibleData.push(trace);
        }
      });
      if (visibleData.length < 2) {
        return;
      }
      // Create a combined trace
//      var yAxes = Object.keys(gd.layout)
//        .filter(v => v.startsWith("yaxis"))
//        .map(v => v.replace("axis", ""));
      var tracesByY = {};
      visibleData.forEach((trace, index) => {
        let yAxis = trace.yaxis || "y";
        if (!tracesByY[yAxis]) {
          tracesByY[yAxis] = [];
        }
        tracesByY[yAxis].push(trace);
      });

      var newTraces = [];
      for (const [y, traces] of Object.entries(tracesByY)) {
        var cmb = traces[0];
        for (let i = 1; i < traces.length; i++) {
          cmb.x = cmb.x.concat(traces[i].x);
          cmb.y = cmb.y.concat(traces[i].y);
          if (cmb.hovertext != null) {
            cmb.hovertext = cmb.hovertext.concat(traces[i].hovertext);
          }
        }
        cmb.name = "Combined";
        newTraces.push(cmb);
      }

      Plotly.deleteTraces(gd, visibleIdxs);
      Plotly.addTraces(gd, newTraces);
    }
  },
  group_mp: {
    name: "Combine visible traces",
    icon: {
      "width": 100,
      "height": 100,
      "path": "M 8 9 A 6 6 0 0 0 1 15 L 1 84 A 6 6 0 0 0 8 90 L 28 90 A 6 6 0 0 0 35 84 A 6 6 0 0 0 28 78 L 14 78 L 14 21 L 28 21 A 6 6 0 0 0 35 15 A 6 6 0 0 0 28 9 L 8 9 z M 71 9 A 6 6 0 0 0 64 15 A 6 6 0 0 0 71 21 L 85 21 L 85 78 L 71 78 A 6 6 0 0 0 64 84 A 6 6 0 0 0 71 90 L 91 90 A 6 6 0 0 0 98 84 L 98 15 A 6 6 0 0 0 91 9 L 71 9 z M 50 27 A 6 6 0 0 0 43 33 L 43 43 L 33 43 A 6 6 0 0 0 27 50 A 6 6 0 0 0 33 56 L 43 56 L 43 66 A 6 6 0 0 0 50 72 A 6 6 0 0 0 56 66 L 56 56 L 66 56 A 6 6 0 0 0 72 50 A 6 6 0 0 0 66 43 L 56 43 L 56 33 A 6 6 0 0 0 50 27 z",
    },
    click: function(gd) {
    
      // Combine single traces
      var visibleIdxs = [];
      var traces = [];
      gd.data.forEach( function(trace, idx) {
        if (trace.visible === true && trace.traceType === "mp_sgl") {
          visibleIdxs.push(idx);
          traces.push(trace);
        }
      });
      if (traces.length >= 2) {
        var cmb = traces[0];
        for (let i = 1; i < traces.length; i++) {
          cmb.x = cmb.x.concat(traces[i].x);
          cmb.y = cmb.y.concat(traces[i].y);
          if (cmb.text != null) {
            cmb.text = cmb.text.concat(traces[i].text);
          }
        }
        cmb.name = "Combined";
        cmb.Id = null;
        var legendGroup = crypto.randomUUID();
        cmb.legendgroup = legendGroup;

        traces = [cmb];

        var hapTraces = [];
        gd.data.forEach( function(trace, idx) {
          if (trace.visible === true && trace.traceType === "hapline") {
            visibleIdxs.push(idx);
            hapTraces.push(trace);
          }
        });

        if (hapTraces.length > 0) {
          var cmbHap = hapTraces[0];
          for (let i = 1; i < hapTraces.length; i++) {
            cmbHap.x = cmbHap.x.concat(hapTraces[i].x);
            cmbHap.y = cmbHap.y.concat(hapTraces[i].y);
          }
          cmbHap.legendgroup = legendGroup;
          cmbHap.line.color = cmb.marker.color;
          traces.push(cmbHap);
        }

        Plotly.deleteTraces(gd, visibleIdxs);
        Plotly.addTraces(gd, traces);
      }

      // Combine bubble plot traces
      var visibleIdxs = [];
      var traces = [];
      gd.data.forEach( function(trace, idx) {
        if (trace.visible === true && trace.traceType === "mp_bbl") {
          visibleIdxs.push(idx);
          traces.push(trace);
        }
      });
      if (traces.length >= 2) {
        var cmb = traces[0];
        for (let i = 1; i < traces.length; i++) {
          var tr = traces[i];
          var cmbX = [...new Set([...cmb.x, ...tr.x])];

          var cmbPos = [...cmb.positions, ...tr.positions.filter(v => !cmb.positions.includes(v))];

          var cmbRuns = cmbX
            .map( v => [
              cmb.runs[cmb.x.indexOf(v)],
              tr.runs[tr.x.indexOf(v)]])
            .map( arr => arr.filter(item => item != undefined).flat());

          var cmbSize = cmbRuns.map(v => v.length * 4 + 2);

          var cmbY = [];
          for (v of cmbX) {
            var cmbVal = [];
            if (cmb.x.indexOf(v) !== -1) {
              cmbVal.push(cmb.y[cmb.x.indexOf(v)] * cmb.runs[cmb.x.indexOf(v)].length);
            }
            if (tr.x.indexOf(v) !== -1) {
              cmbVal.push(tr.y[tr.x.indexOf(v)] * tr.runs[tr.x.indexOf(v)].length);
            }
            cmbY.push(cmbVal.reduce((a, b) => a + b) / cmbRuns[cmbX.indexOf(v)].length);
          }

          var cmbText = cmbX.map( (_, i) => {
            return `${cmbPos[i]}<br><br>${cmbRuns[i].join("<br>")}`;
          });

//        cmbText = cmbX
//          .map( v => [
//            cmb.text[cmb.x.indexOf(v)],
//            tr.text[tr.x.indexOf(v)]])
//          .map( arr => arr.filter(item => item != undefined));

          cmb.x = cmbX;
          cmb.y = cmbY;
          cmb.runs = cmbRuns;
          cmb.marker.size = cmbSize;
          cmb.positions = cmbPos;
          cmb.text = cmbText;
        }
        cmb.name = "Combined";
 
        Plotly.deleteTraces(gd, visibleIdxs);
        Plotly.addTraces(gd, [cmb]);
      }
    }
  },
  _delete: {
    name: "Delete visible traces",
    icon: {
      "width": 100,
      "height": 100,
      "path": "m 24,17 a 7,7 0 0 0 -5,2 7,7 0 0 0 0,10 L 39,49 19,70 a 7,7 0 0 0 0,10 7,7 0 0 0 10,0 L 50,60 70,80 a 7,7 0 0 0 10,0 7,7 0 0 0 0,-10 L 60,49 80,29 a 7,7 0 0 0 0,-10 7,7 0 0 0 -5,-2 7,7 0 0 0 -5,2 L 50,39 29,19 a 7,7 0 0 0 -5,-2 z",
    },
    click: function(gd) {
      var visibleIdxs = [];
      gd.data.forEach( function(trace, idx) {
        if ((trace.visible === true) && (!trace.permanent)) {
          visibleIdxs.push(idx);
        }
      });
      Plotly.deleteTraces(gd, visibleIdxs);
    }
  },
}

