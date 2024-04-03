function makeGroup(event, setToGroup) {
	if (setToGroup === "genes") {
		var selectedItems = document.getElementById("gene-field").value.split(",");
		for (let i = 0; i < selectedItems.length; i++) {
			selectedItems[i] = selectedItems[i].trim();
		}
	} else if (setToGroup === "expruns") {
		var runsTable = document.getElementById("runs-table")
		var runRows = runsTable.querySelector("tbody").querySelectorAll("tr");
		var selectedItems = [];
		for (let i = 0; i < runRows.length; i++) {
			if (!runRows[i].classList.contains("collapsed")) {
				selectedItems.push(runRows[i].querySelector("td").innerHTML);
			}
		}
	} else {
		var allItems = document.querySelectorAll("." + setToGroup + "-checkbox");
		var selectedItems = [];
		for (var item of allItems) {
			if (item.checked) {
				selectedItems.push(item.value);
			}
		}
	}

	if (selectedItems.length !== 0) {
		var groupName = "";
		while (groupName === "") {
			var groupName = window.prompt("Group name:");
		}
		if (groupName !== null) {
			var data = JSON.stringify({user_id: user_id,
									   set: setToGroup,
									   name: groupName,
									   items: selectedItems});
			var request = new XMLHttpRequest();
			request.open("POST", "/makeDeleteGroup", true);
			request.send(data);

			request.onreadystatechange = function() {
				if (request.readyState === XMLHttpRequest.DONE) {
					if (request.status === 204) {
						var dropdown = document.getElementById(setToGroup + "_group-dropdown");
						var checkboxRow = document.createElement("div");
						checkboxRow.classList.add("checkbox-row");
						checkboxRow.classList.add("group-option");
						var checkbox = document.createElement("input")
						checkbox.type = "checkbox";
						checkbox.classList.add(setToGroup + "_group-checkbox");
						checkbox.id = setToGroup + "_group-dropdown" + groupName;
						checkbox.value = selectedItems.join(",");

						checkbox.addEventListener("change", function(event) {
							updateCount(this.parentElement.parentElement);
							if (setToGroup === "expruns") {
								filterRunsTable();
							}
						});

						var label = document.createElement("label");
						label.classList.add("checkbox-label")
						label.htmlFor = checkbox.id;
						label.appendChild(document.createTextNode(groupName));

						showGroupContents(checkboxRow);

						checkboxRow.appendChild(checkbox);
						checkboxRow.appendChild(document.createTextNode(" "));
						checkboxRow.appendChild(label);
						dropdown.appendChild(checkboxRow);

					} else {
						alert(request.responseText);
					}
				}
			};
		}
	} else {
		showHelp(event, "make-group");
	}
	event.preventDefault();
}

function deleteGroup(event, setToDelete) {
	var allItems = document.querySelectorAll("." + setToDelete + "_group-checkbox");
	var selectedItems = [];
	for (var item of allItems) {
		if (item.checked) {
			selectedItems.push(item);
		}
	}
	if (selectedItems.length !== 0) {
		var deletePrompt = window.confirm("Delete selected groups?");
		if (deletePrompt == true) {
			var names = [];
			for (var item of selectedItems) {
				names.push(item.nextElementSibling.textContent);
			}
			var data = JSON.stringify({user_id: user_id,
									   set: setToDelete,
									   items: names});

			var request = new XMLHttpRequest();
			request.open("POST", "/makeDeleteGroup", true);
			request.send(data);

			// Remove checkboxes from the dropdown
			for (let i = 0; i < selectedItems.length; i++) {
				selectedItems[i].parentElement.remove();
			}

			var dropdown = document.getElementById(setToDelete + "_group-dropdown");
			updateCount(dropdown);
		}
	} else {
		showHelp(event, "delete-group");
	}
	event.preventDefault();
}