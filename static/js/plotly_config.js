const scriptUrl = "https://cdn.plot.ly/plotly-2.35.2.min.js";
const script = document.createElement("script");
script.src = scriptUrl;
script.type = "text/javascript";
script.async = true;
document.head.appendChild(script);

var config = {
  editable: false,
  modeBarButtonsToAdd: [
    {
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
    {
      name: "Recolor traces (double-click to change mode)",
      icon: {
        "width": 100,
        "height": 100,
        "path": "m 11,11 c -3,0 -5,2 -5,5 v 19 c 0,3 2,5 5,5 H 68 c 3,0 5,-2 5,-5 v -5 h 7 V 48 L 41,62 c -1,0 -2,2 -2,4 v 25 c 3e-4,2 1,4 4,4 2,-3e-4 4,-1 4,-4 V 70 l 40,-13 c 1,-0 2,-2 2,-4 V 26 C 91,24 89,22 86,22 H 74 v -5 c 0,-3 -2,-5 -5,-5 z",
      },
      click: function(gd) {
        var clickEvent = event;
        if (clickTimeout) {
          clearTimeout(clickTimeout);
          recolorMode = recolorMode === "cycle" ? "custom" : "cycle";
          clickTimeout = null;
        } else {
          clickTimeout = setTimeout(function() {
            if (recolorMode === "cycle") {
              var t;
              var color;
              var colors = colorCycle.slice();
              var traceIdxs = [];
              var traceColors = [];
              var legendGroups = [];

              for (let i = 0; i < gd.data.length; i++) {
                t = gd.data[i];
                if (t.traceType === "__fm__" || t.traceType == null) {
                  if (t.legendgroup == null || !legendGroups.includes(t.legendgroup)) {
                    if (colors.length === 0) {
                      colors = colorCycle.slice();
                    }
                    color = colors.splice(Math.floor(Math.random() * colors.length), 1)[0];
                    //var [r,g,b] = Array(3).fill().map(() => Math.floor(Math.random() * 208));
                    //color = `rgb(${r}, ${g}, ${b})`;
                    traceIdxs.push(i);
                    traceColors.push(color);
                    if (t.legendgroup != null) {
                      legendGroups.push(t.legendgroup);
                    }
                  } else {
                    traceIdxs.push(i);
                    traceColors.push(traceColors[legendGroups.indexOf(t.legendgroup)]);
                  }
                }
              }
            } else {
              var visibleTraces = [];
              var visibleIdx = 0;
              gd.data.forEach( function(trace, idx) {
                if (
                 trace.showlegend !== false &&
                 (!trace.hasOwnProperty("visible") ||
                  trace.visible !== "legendonly")) {
                  visibleTraces.push(trace.Id);
                  visibleIdx = idx;
                }
              });
              if (visibleTraces.length !== 1) {
                showHelp(clickEvent, null, "Select exacly one trace to recolor.");
                clickTimeout = null;
                return false;
              } else {
                var traceColor = "";
                while (traceColor === "") {
                  var traceColor = prompt("Trace color:");
                }
              }
              if (traceColor == null) {
                clickTimeout = null;
                return false;
              }
              traceColors = [traceColor];
              traceIdxs = [visibleIdx];
            }
            Plotly.restyle(
              "live-plot",
              {"line.color": traceColors,
               "marker.color": traceColors},
              traceIdxs);

            clickTimeout = null;
          }, delay);
        }
      }
    },
    {
      name: "Edit plot",
      icon: Plotly.Icons.pencil,
      click: function(gd) {
        Plotly.newPlot(
          gd.id,
          gd.data,
          gd.layout,
          { ...gd._context, editable: !gd._context.editable});
      }
    },
  ],
  modeBarButtonsToRemove: [],
};
