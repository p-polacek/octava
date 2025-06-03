var customOrder = [];
var customColors = {};

var plotColors = [
  //"#1f77b4",  // muted blue
  "#ff7f0e",  // safety orange
  //"#2ca02c",  // cooked asparagus green
  "#d62728",  // brick red
  "#9467bd",  // muted purple
  "#8c564b",  // chestnut brown
  "#e377c2",  // raspberry yogurt pink
  "#7f7f7f",  // middle gray
  "#bcbd22",  // curry yellow-green
  "#17becf"   // blue-teal
];

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
        { ...gd._context, editable: !gd._context.editable});
      setPlotEventListeners();
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
      var colorGradients;
      var traceIdxs = [];
      var traceIds = [];
      gd.data.forEach(function(trace, idx) {
        if (trace.visible === true) {
          traceIdxs.push(idx);
          traceIds.push(trace.Id || "");
        }
      });

      showHelp(clickEvent, null, '<input type="color" id="colorPicker" value="#ff0000" />');
      document.getElementById("colorPicker").click();
      document.getElementById("help-popup").classList.add("collapsed");
      document.getElementById("colorPicker").addEventListener("input", function() {
        var colorChoice = this.value;
        traceColors = [];
        colorGradients = [];
        gd.data.forEach(function(trace, idx) {
          if (trace.visible === true) {
            traceColors.push(colorChoice);
            customColors[trace.Id || ""] = colorChoice;
          }
        });

        var recolor = {
          "marker.color": traceColors,
          "line.color": traceColors,
        }
        if (gd.data[0].type !== "box") {
          recolor["fillcolor"] = traceColors;
        }

        Plotly.restyle(
          gd,
          recolor,
          traceIdxs);
      });
    },
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
      var traceIdxs = [];
      var traceIds = [];
      gd.data.forEach(function(trace, idx) {
        if (trace.visible === true) {
          traceIdxs.push(idx);
          traceIds.push(trace.Id || "");
        }
      });

      showHelp(clickEvent, null, '<input type="range" id="opacitySlider" min="0" max="1" step="0.01" value="1" />');
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
        cmb.legendgroup = "Combined";
        newTraces.push(cmb);
      }

      Plotly.deleteTraces(gd, visibleIdxs);
      Plotly.addTraces(gd, newTraces);
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
        if (trace.visible === true) {
          visibleIdxs.push(idx);
        }
      });
      Plotly.deleteTraces(gd, visibleIdxs);
    }
  },
}

