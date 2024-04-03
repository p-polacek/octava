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
				keys.forEach( function(row) {
					fieldText += "<tr><th>" + row + "</th>";
					varData.forEach( function(elem) {
						fieldText += "<td>" + elem[row] + "</td>";
					});
					fieldText += "</tr>";
				});
				fieldText += "</table>";
			} else {
				var fieldText = "<table><th>" + varData + "</th></table>";
			}
			showHelp(event, null, fieldText);
			return;
		}
	});
}