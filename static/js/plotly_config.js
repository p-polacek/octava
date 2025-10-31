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
// m_size       # Change marker size of visible traces
// log          # Switch to log y scale
// reset_color  # Reset colors of visible traces
// move_down    # Move visible traces down
// group        # Combine visible traces
// group_mp     # Combine visible traces
// _delete      # Delete visible traces


var buttonsToAdd = {
  save_html: {
    name: "Save plot as interactive HTML",
    //icon: Plotly.Icons.disk,
    icon: {
      "width": 100,
      "height": 100,
      "path": "M 12.433594,93.725087 C 8.9142385,93.704507 6.1017715,90.226466 6.3242125,86.80868 c 0.0099,-24.882663 -0.02896,-49.766614 0.02761,-74.648469 0.219507,-3.3953736 3.541788,-6.1949366 6.9243505,-5.8726396 17.186411,0.0086 34.374176,-0.02664 51.559734,0.02678 3.187934,0.277819 6.218174,1.68241 8.719094,3.635051 2.173919,1.9561836 4.124919,4.1422636 6.225551,6.1771306 3.549089,3.596333 7.16511,7.132583 10.649291,10.78931 2.128662,2.745219 3.335993,6.231174 3.241736,9.712617 -0.002,17.08326 0.03261,34.168292 -0.02425,51.25092 -0.165738,3.261728 -3.229243,6.103009 -6.519591,5.855509 -24.89805,0.0059 -49.796092,0.0025 -74.69414,-0.0098 z m 8.113281,-9.265625 c 0.01851,-7.611364 -0.05529,-15.226672 0.05108,-22.835529 0.373893,-3.312998 3.660096,-5.965796 6.980947,-5.609801 15.366891,0.01012 30.735142,-0.0293 46.101179,0.02803 3.307451,0.207826 6.12703,3.405215 5.77304,6.730579 V 84.45946 h 5.044922 c 0.0036,-16.255861 7.45e-4,-32.51172 -0.0078,-48.767578 C 84.134416,33.947846 83.199578,32.33186 81.835834,31.185122 76.920296,26.322837 72.049319,21.4076 67.093839,16.5906 c -0.575698,-0.375954 -1.203316,-0.677099 -1.863361,-0.871375 -0.01978,7.568882 0.058,15.142345 -0.0533,22.708321 -0.309243,3.191628 -3.368741,5.856528 -6.592796,5.580838 -10.790367,-0.0112 -18.936619,0.03839 -29.725897,-0.03426 -3.418694,-0.28557 -5.945543,-3.773512 -5.66577,-7.099435 V 15.564923 H 15.50196 v 68.894532 h 5.044922 z m 49.728516,0 V 65.293446 H 29.724609 V 84.459462 Z M 52.835938,34.730946 c 0.968528,-0.06855 1.984996,0.109173 2.923828,-0.115234 0.59064,-0.41344 0.165835,-1.220309 0.292092,-1.814453 -0.01775,-5.604403 0.04424,-11.212327 -0.03819,-16.814454 -0.411442,-0.745365 -1.433828,-0.310945 -2.114453,-0.42218 -3.173344,0.02409 -6.358066,-0.0582 -9.524219,0.05499 -0.751035,0.286511 -0.31726,1.188301 -0.426945,1.776953 0.02698,5.646878 -0.06218,11.302309 0.05195,16.94375 0.410712,0.693252 1.363864,0.275961 2.014844,0.391094 2.273698,9.35e-4 4.547395,0.0023 6.821094,-4.69e-4 z",
    },
    click: function(gd) {
      async function exportPlotlyHTML(gd) {
        const fig = { data: gd.data, layout: gd.layout };

        const html = `
          <!DOCTYPE html>
          <html>
          <head>
            <meta charset="utf-8">
            <script
             src="https://cdn.plot.ly/plotly-2.35.2.min.js">
            </script>
            <style>
              html, body {
                margin: 0;
                padding: 0;
                width: 100%;
                height: 100%;
                overflow: hidden;
              }
              #plot {
                width: 100%;
                height: 100%;
              }
            </style>
          </head>
          <body>
            <div id="plot"></div>
            <script>
              const fig = ${JSON.stringify(fig, null, 2)};
              Plotly.newPlot('plot', fig.data, fig.layout);
              window.addEventListener('resize', () => {
                Plotly.Plots.resize('plot');
              });
            </script>
          </body>
          </html>
        `;

        const blob = new Blob([html], { type: 'text/html' });
        const url = URL.createObjectURL(blob);
        const a = document.createElement('a');
        a.href = url;
        a.download = 'octava_plot.html';
        a.click();
        URL.revokeObjectURL(url);
      }
      exportPlotlyHTML(gd);
    } 
  }, 
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
      var editable = !gd._context.editable;
      // Enable rescaling axes
      var axes = Object.keys(gd.layout).filter(k => k.startsWith("xaxis") || k.startsWith("yaxis"));
      for (var ax of axes) {
        ax = gd.layout[ax];
        ax.minallowed = editable ? false : ax.range[0];
        ax.maxallowed = editable ? false : ax.range[1];
      }
      Plotly.newPlot(
        gd.id,
        gd.data,
        gd.layout,
        { ...gd._context, editable: editable, edits: {shapePosition: false} });

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
        currColor = "#1f77b4ff";
      }

      showHelp(clickEvent, null, `<div id="colorPickerContainer"></div>`);
      const popup = document.getElementById("help-popup");
      const rect = popup.getBoundingClientRect();
      popup.style.left = Math.max(window.innerWidth, rect.left) - 280 + "px";

      const picker = new Picker({
        parent: document.getElementById("colorPickerContainer"),
        color: currColor,
        popup: false,
        editor: true,
        alpha: true,
        onChange: (color) => {
          const colorChoice = color.hex;

          for (let i = 0; i < traceIdxs.length; i++) {
            var recolor = {
              "marker.color": [colorChoice],
              "line.color": [colorChoice],
            }
            if (gd.data[traceIdxs[i]].type !== "box") {
              recolor["fillcolor"] = [colorChoice];
            }
            //if (isAlpha) {
            //  recolor["marker.color"] = [colorChoice + "80"];
            //  recolor["marker.line.color"] = [colorChoice];
            //}

            Plotly.restyle(
              gd,
              recolor,
              [traceIdxs[i]]);
          }
        },
        onDone: (color) => {
          picker.destroy();
          document.getElementById("help-popup").classList.add("collapsed");
        },
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
//  opacity: {
//    name: "Change transparency of visible traces",
//    icon: {
//      "width": 100,
//      "height": 100,
//      "path": "M 50 6 C 26 6 6 26 6 50 C 6 73 26 93 50 93 C 73 93 93 73 93 50 C 93 26 73 6 50 6 z M 50 19 C 66 19 80 33 80 50 C 80 66 66 80 50 80 C 33 80 19 66 19 50 C 19 33 33 19 50 19 z M 60 28 A 5 5 0 0 0 56 30 L 30 56 A 5 5 0 0 0 30 64 A 5 5 0 0 0 38 64 L 64 37 A 5 5 0 0 0 64 30 A 5 5 0 0 0 60 28 z M 66 49 A 5 5 0 0 0 62 51 L 51 63 A 5 5 0 0 0 51 70 A 5 5 0 0 0 59 70 L 70 59 A 5 5 0 0 0 70 51 A 5 5 0 0 0 66 49 z",
//    },
//    click: function(gd) {
//      var clickEvent = event;
//      var traceOpacities = [];
//      var traceIdxs = [];
//      var traceIds = [];
//      gd.data.forEach(function(trace, idx) {
//        if (trace.visible === true) {
//          traceIdxs.push(idx);
//          traceIds.push(trace.Id || "");
//          traceOpacities.push(trace.opacity || 1);
//        }
//      });
//      var opacityValue = 1;
//      if ([...new Set(traceOpacities)].length === 1) {
//        opacityValue = traceOpacities[0];
//      }
//
//      showHelp(clickEvent, null, `<input type="range" id="opacitySlider" min="0" max="1" step="0.01" value="${opacityValue}" />`);
//      document.getElementById("opacitySlider").addEventListener("input", function() {
//        var opacityChoice = this.value;
//        Plotly.restyle(
//          gd,
//          {opacity: opacityChoice,
//          },
//          traceIdxs);
//      });
//    },
//  },
  m_size: {
    name: "Change marker size of visible traces",
    icon: {
      "width": 100,
      "height": 100,
      "path": "M 65 6 C 54 7 43 14 39 25 C 35 39 42 55 57 60 C 71 64 86 57 91 42 L 91 42 C 96 28 88 12 74 8 C 71 7 68 6 65 6 z M 22 61 C 15 61 9 65 7 72 C 4 80 9 89 17 92 C 25 95 34 90 37 82 L 37 82 C 40 74 35 65 27 62 C 25 62 24 61 22 61 z",
    },
    click: function(gd) {
      var clickEvent = event;
      var traceSizes = [];
      var traceIdxs = [];
      var traceColors = [];
      var editLine = false;
      gd.data.forEach(function(trace, idx) {
        if (trace.visible === true && (trace.traceType == null || ["mp_sgl", "hapline"].includes(trace.traceType))) {
          if (trace.type.includes("scatter") || trace.type === "box") {
            if (trace.traceType != null) {
                editLine = true;
            }
            traceIdxs.push(idx);
            if (trace.marker != null) {
              traceSizes.push(trace.marker.size || 6);
              traceColors.push(trace.marker.color);
            } else if (trace.line != null) {
              traceSizes.push(trace.line.size || 6)
              traceColors.push(trace.line.color);
            }
          }
        }
      });
      var sizeValue = 6;
      if ([...new Set(traceSizes)].length === 1) {
        sizeValue = traceSizes[0];
      }

      showHelp(clickEvent, null, `<input type="range" id="sizeSlider" min="1" max="20" step="1" value="${sizeValue}" />`);
      document.getElementById("sizeSlider").addEventListener("input", function() {
        var sizeChoice = this.value;
        var updateDict = {
            "marker.size": Array(traceIdxs.length).fill(sizeChoice),
            "marker.color": traceColors,
        }
        if (editLine) {
            updateDict["line.width"] = Array(traceIdxs.length).fill(sizeChoice);
            updateDict["line.color"] = traceColors;
        }

        Plotly.restyle(
          gd,
          updateDict,
          traceIdxs
        );
      });
    },
  },
  log: {
    name: "Toggle log-scale y-axis",
    icon: {
      "width": 100,
      "path": "M 62 6 A 6 6 0 0 0 55 13 A 6 6 0 0 0 62 19 L 80 19 L 80 86 A 6 6 0 0 0 86 93 A 6 6 0 0 0 93 86 L 93 13 A 6 6 0 0 0 86 6 L 62 6 z M 51 35 C 44 35 40 37 36 41 C 32 45 31 51 31 58 L 31 70 C 31 77 32 83 36 87 C 40 91 44 93 51 93 C 57 93 62 91 65 87 C 69 83 71 77 71 70 L 71 58 C 71 51 69 45 65 41 C 62 37 57 35 51 35 z M 6 37 L 6 91 L 22 91 L 22 37 L 6 37 z M 51 50 C 53 50 55 51 55 55 L 55 73 C 55 77 53 78 51 78 C 48 78 47 77 47 73 L 47 55 C 47 51 48 50 51 50 z",
      "height": 100,
    },
    click: function(gd) {
      Plotly.relayout("genexp-plot", {"yaxis.type": gd.layout.yaxis.type === "linear" ? "log" : "linear"});
    }
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
      var maxesByY = {};
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
          if (cmb.type === "bar") {
            cmb.y = cmb.y.map((_,idx) => { return cmb.y[idx] + traces[i].y[idx] });
          } else {
            cmb.y = cmb.y.concat(traces[i].y);
          }
          if (cmb.hovertext != null) {
            cmb.hovertext = cmb.hovertext.concat(traces[i].hovertext);
          }
        }
        cmb.name = "Combined";
        newTraces.push(cmb);
      }

      Plotly.deleteTraces(gd, visibleIdxs);
      Plotly.addTraces(gd, newTraces);
      //Plotly.relayout(
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

