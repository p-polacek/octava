function getVarData(event, chrom, pos, field) {
	event.preventDefault();

	$.ajax({
		url: "/get_var_data",
		method: "POST",
		data: { genome: genome,
				chrom:  chrom,
				pos:    pos,
				field:  field },
		success: function (data) {
			var varData = data[0];
			var keys = data[1];

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