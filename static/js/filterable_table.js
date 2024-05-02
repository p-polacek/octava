function parseTSV(tsvStream, sep) {
	// Parses a string into a data object based on a separator (tab by default)
	var lines = tsvStream.trim().split(/\r\n|\r|\n/);
	if (sep == null) {
		for (var i = 0; i < lines.length; i++) {
			lines[i] = lines[i].split(/\t/);
		}
	} else {
		for (var i = 0; i < lines.length; i++) {
			lines[i] = lines[i].split(sep);
		}
	}
	return lines;
}

function createTable(tableData, table) {
	// Creates a filterable table given a data object and a table element
	parentElem = table.parentElement;
	parentElem.removeChild(table);

	table.innerHTML = "";
	var tableId = table.id;

	var cols = tableData[0];
	var thead = document.createElement("thead");
	var tbody = document.createElement("tbody");
	var tr = document.createElement("tr");

	for (var i = 0; i < cols.length; i++) {
		var th = document.createElement("th");
		th.classList.add("header-top");

		var tdDiv = document.createElement("div");
		tdDiv.style.display = "flex";

		var textDiv = document.createElement("div");
		textDiv.style.flex = "1";
		textDiv.innerHTML = cols[i] + "&nbsp;&#9013;";
		textDiv.id = tableId + "_" + i;
		textDiv.classList.add("clickable");
		textDiv.addEventListener("click", function() {
			dropDown(event, this);
		});

		var resizDiv = document.createElement("div");
		resizDiv.classList.add("resizable");
		resizDiv.style.flex = "0 0 auto";
		resizDiv.style.width = "17px";

		tdDiv.appendChild(textDiv);
		tdDiv.appendChild(resizDiv);

		th.appendChild(tdDiv);
		tr.appendChild(th);
	}
	thead.appendChild(tr);
	table.appendChild(thead);

	tableData.slice(1).forEach(function(rowData) {
		var tr = document.createElement("tr");

		rowData.forEach(function(cellData) {
			var td = document.createElement("td");
			td.classList.add("default-col");
			if (typeof cellData === "object") {
				const subcellWidth = 100 / cellData.length;

				for (var i = 0; i < cellData.length; i++) {
					const subcell = document.createElement("div");
					subcell.className = "subcell";
					subcell.style.width = subcellWidth + "%";
					subcell.textContent = cellData[i];
					td.appendChild(subcell);
				}				
			} else {
				td.textContent = cellData;
			}
			tr.appendChild(td);
		});
		tbody.appendChild(tr);
	});
	table.appendChild(tbody);

	table.style.opacity = "0";
	table.style.position = "absolute";
	document.body.appendChild(table);

	// Set correct padding for header cells
	table.querySelectorAll("th").forEach( function(header) {
		var pad = window.getComputedStyle(header).getPropertyValue("padding");
		header.querySelector(".clickable").style.padding = [pad, "0px", pad, pad].join(" ");
		header.querySelector(".resizable").style.margin = [pad, pad, pad, "0px"].join(" ");
		header.style.padding = "0px";
	});

	// Set title for shortened cells
	table.querySelectorAll("td").forEach( function(cell) {
		if (cell.offsetWidth < cell.scrollWidth) {
			cell.title = cell.innerText;
		}
	});

	document.body.removeChild(table);
	table.style.position = "";
	table.style.opacity = "";
	parentElem.appendChild(table);
}

function getUniqueValues(table, colIndex) {
	// Returns all unique values in a table column
	var uniqueValues = [];
	var rows = table.querySelector("tbody").querySelectorAll("tr");
	for (var i = 0; i < rows.length; i++) {
		var cell = rows[i].children[colIndex].textContent;
		if (!uniqueValues.includes(cell)) {
			uniqueValues.push(cell);
		}
	}
	return uniqueValues;
}

function findParentTable(elem) {
	let par = elem.parentNode;
    while (par !== null && par.tagName !== "TABLE") {
        par = par.parentNode;
    }
    return par;
}

function makeDropdown(tableId, colIndex, uniqueValues) {
	// Create new dropdown if it doesn't exist
	dropdown = document.createElement("div");
	dropdown.id = tableId + "dropdown_" + colIndex;
	dropdown.classList.add("dropdown-menu");
	dropdown.classList.add("collapsed");

	// Create search bar
	var searchBar = document.createElement("input");
	searchBar.type = "text";
	searchBar.classList.add("search-bar");
	searchBar.placeholder = "Search...";

	(function(dropdown) {
		searchBar.addEventListener("keyup", function() {
			filterDropdown([dropdown, tableId]);
			sortGray(dropdown);
		});
	})(dropdown);

	// Create (un)select all button
	var toggleAllButton = document.createElement("button");
	toggleAllButton.type = "button";
	toggleAllButton.textContent = "Select all";

	(function(dropdown) {
		toggleAllButton.addEventListener("click", function() {
			toggleAll([dropdown, tableId], true);
		});
	})(dropdown);

	var header = document.createElement("div");
	header.classList.add("sticky");
	header.classList.add("non-transparent");
	header.classList.add("small-padding");
	header.appendChild(searchBar);
	header.insertAdjacentHTML("beforeend", "&nbsp;");
	header.appendChild(toggleAllButton);
	dropdown.appendChild(header);

	// Get unique values for the column
	// Create checkbox for each unique value
	uniqueValues.forEach(function(value) {
		var checkboxRow = document.createElement("div");
		checkboxRow.classList.add("checkbox-row");
		var checkbox = document.createElement("input");
		checkbox.type = "checkbox";
		checkbox.id = tableId + colIndex + value;
		checkbox.value = value;
		checkbox.checked = true;

		(function(dropdown, val) {
			checkbox.addEventListener("change", function(event) {
				applyFilters([dropdown, tableId], [val], event.target.checked);
				updateCount(dropdown);
			});
		})(dropdown, value);

		var label = document.createElement("label");
		label.htmlFor = checkbox.id;
		label.classList.add("checkbox-label");
		if (!value) {
			label.appendChild(document.createTextNode("(blank)"));
			label.classList.add("gray-text");
		} else {
			label.appendChild(document.createTextNode(value));
		}

		checkboxRow.appendChild(checkbox);
		checkboxRow.appendChild(label);
		dropdown.appendChild(checkboxRow);
	});

	// Prevent hiding the popup when clicked within
	popupClickedRetain(dropdown);

	document.getElementById(tableId).appendChild(dropdown);
	dropdowns[tableId][colIndex] = dropdown;

	return dropdown;
}

function applyFilters(dropdownData, checkedVals, checked) {
	var dropdown = dropdownData[0];
	var tableId = dropdownData[1];

	var colIndex = dropdown.id.split(/[_]+/).pop();
	var table = document.getElementById(tableId);
	var rows = table.getElementsByTagName("tbody")[0].getElementsByTagName("tr");

	// Show rows for checked checkbox
	// Not applicable if the changed checkbox was unchecked
	if (checked) {
		// Loop through each row and show if filtered in
		for (var i = 0; i < rows.length; i++) {
			var cells = rows[i].children;
			var showRow = false;

			// If row is to be filtered back in
			if (checkedVals.includes(cells[colIndex].textContent)) {
				showRow = true;

				// Re-apply remaining filters to the row
				// Loop through each filter checkbox for each dropdown menu
				for (const [index, dropdown] of Object.entries(dropdowns[tableId])) {
					if (index === colIndex) { continue }

					var checkboxes = dropdown.querySelectorAll('input[type="checkbox"]');
					for (checkbox of checkboxes) {
						if (!checkbox.checked && cells[index].textContent === checkbox.value) {
							showRow = false;
							break;
						}
					}
				}
			}
			// Apply hide/show to row
			if (showRow) {
				rows[i].classList.remove("collapsed");
			} else {
				rows[i].classList.add("collapsed");
			}
		}
	}

	// Apply filters from projects-table as well, if filtering run-table
	if (tableId === "runs-table") {
		var projectsTable = document.getElementById("projects-table");
		var projectRows = projectsTable.querySelector("tbody").querySelectorAll("tr");
		var visibleBPs = [];

		for (let i = 0; i < projectRows.length; i++) {
			var row = projectRows[i];
			if (!row.classList.contains("collapsed")) {
				visibleBPs.push(row.querySelector("td").textContent);
			}
		}
	} else {
		var visibleBPs = null;
	}

	// Loop through each row and hide based on filters
	for (var i = 0; i < rows.length; i++) {
		var cells = rows[i].children;
		var showRow = true;

		if (visibleBPs !== null) {
			if (!visibleBPs.includes(cells[1].innerHTML)) {
				showRow = false;
			}
		}

		// Loop through each filter checkbox for each dropdown menu
		Object.keys(dropdowns[tableId]).forEach(function(colIndex) {
			dropdowns[tableId][colIndex].querySelectorAll('input[type="checkbox"]').forEach(function(checkbox) {
				if (!checkbox.checked && cells[colIndex].textContent === checkbox.value) {
					showRow = false;
				}
			});
		});

		// Apply hide/show to row
		if (showRow) {
			rows[i].classList.remove("collapsed");
		} else {
			rows[i].classList.add("collapsed");
		}
	}

	// Highlight columns with filters applied
	var header = table.getElementsByTagName("th");
	Object.keys(dropdowns[tableId]).forEach(function(colIndex) {
		var toHighlight = false;
		checkboxes = dropdowns[tableId][colIndex].querySelectorAll('input[type="checkbox"]')
		for (let i = 0; i < checkboxes.length; i++) {
			var checkbox = checkboxes[i];
			var checkboxRow = checkbox.parentElement;

			if (!checkbox.checked && !checkboxRow.classList.contains("grayed-out")) {
				toHighlight = true;
				break;
			}
		}
		header[colIndex].style.color = toHighlight ? "#ff9c78" : "";
	});
}

function sortGray(dropdown) {
	// Sorts checkboxes to plase greyed-out values to the bottom
	var checkboxRows = Array.from(dropdown.querySelectorAll(".checkbox-row"));
	checkboxRows.sort(function(a, b) {
		var classA = a.classList.contains("grayed-out") ? 1 : -1;
		var classB = b.classList.contains("grayed-out") ? 1 : -1;
		return classA - classB;
	});

	while (dropdown.childNodes.length > 1) {
		dropdown.removeChild(dropdown.lastChild);
	}
	// Append the sorted divs back to the container
	checkboxRows.forEach(function(checkboxRow) {
		dropdown.appendChild(checkboxRow);
	});
}
