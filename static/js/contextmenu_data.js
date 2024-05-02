function contextmenuData(event, payload) {
	event.preventDefault();
	payload["genome"] = genome;

	$.ajax({
		url: "/contextmenu_data",
		method: "POST",
		data:   payload,
		success: function (results) {
			var varData = results[0];
			var keys = results[1];
			//alert(JSON.stringify(varData));

			if (typeof varData === "object") {
				var fieldText = "<table>";
				for (let i = 0; i < keys.length; i++) {
					fieldText += "<tr><th>" + keys[i] + "</th>";
					varData.forEach( function(elem) {
						fieldText += "<td>" + elem[i] + "</td>";
					});
					fieldText += "</tr>";
				}
				fieldText += "</table>";
			} else {
				var fieldText = "<table><th>" + varData + "</th></table>";
			}
			showHelp(event, null, fieldText);
			return;
		}
	});
}