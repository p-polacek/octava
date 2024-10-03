function contextMenu(event, payload) {
  event.preventDefault();
  payload["genome"] = genome;
  $.ajax({
    url: "/contextmenu",
    method: "POST",
    data:   payload,
    success: function (results) {
      var varData = results[0];
      var keys = results[1];
      console.log(results);
      var fieldText = "<table>";
      for (let i = 0; i < varData[0].length; i++) {
        if (keys.length > 0) {
           fieldText += "<tr><th>" + keys[i] + "</th>";
        }
        varData.forEach( function(elem) {
          fieldText += "<td>" + elem[i] + "</td>";
        });
        fieldText += "</tr>";
      }
      fieldText += "</table>";
      showHelp(event, null, fieldText);
      return;
    }
  });
}