function validateForm() {
	if (document.getElementById("fmt-field").value == "gene") {
		if (0 && document.getElementById("gene-field").value == "") {
			alert("Specify a gene name.");
			return false;
		}
	}
	else if (document.getElementById("fmt-field").value == "region") {
		var start = document.getElementById("start-field").value;
		var end = document.getElementById("end-field").value;
		if (start == "") {
			alert("Specify region start.");
			return false;
		}
		else if (end == "") {
			alert("Specify region end.");
			return false;
		}
		else if (parseInt(end) < parseInt(start)) {
			alert("Region start cannot be smaller than end.");
			return false;
		}
	}
	var impacts = document.querySelectorAll(".impact-checkbox");
	var atLeastOneChecked = Array.from(impacts).some(checkbox => checkbox.checked);
	if (!atLeastOneChecked) {
		alert("Specify target impact.");
		return false;
	}
return true;
}
