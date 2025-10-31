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
function createTable(tableData, table, undraggable) {
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
    if (!undraggable.includes(i)) {
      th.classList.add("draggable");
    }
    var tdDiv = document.createElement("div");
    tdDiv.classList.add("disp-flex");
    var textDiv = document.createElement("div");
    textDiv.innerHTML = cols[i] + "&nbsp;&#9013;";
    textDiv.id = tableId + "_" + i;
    textDiv.classList.add("clickable", "flex");
    textDiv.addEventListener("click", function() {
      dropDown(event, this);
    });
    var resizDiv = document.createElement("div");
    resizDiv.classList.add("resizable", "resiz-render");
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
          subcell.style.paddingRight = "2px";
          subcell.style.paddingLeft = "2px";
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
    header.querySelector(".clickable").style.padding = [
     pad,
     "0px",
     pad,
     pad].join(" ");
    header.querySelector(".resizable").style.margin = [
     pad,
     pad,
     pad,
     "0px"].join(" ");
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
var dropdowns = {}; // Array to store dropdown menus for each column
function dropDown(event, clickedElem) {
  var table = findParentTable(clickedElem);
  var tableId = table.id

  var colIdNum = clickedElem.id.split(/[_]+/).pop();
  var colIndex = [...table.querySelector("thead").querySelectorAll("th")]
    .map(h => h.querySelector("div").querySelector("div"))
    .indexOf(clickedElem);

  if (dropdowns[tableId]) {
    var dropdown = dropdowns[tableId][colIdNum];
  } else {
    dropdowns[tableId] = {};
  }

  // Check if dropdown already exists
  if (!dropdown) {
    var uniqueValues = getUniqueValues(table, colIndex).sort();
    var dropdown = makeDropdown(tableId, colIdNum, uniqueValues);

    // Remove starting rows from "Filter samples" dropdown
    if (tableId === "vars-table" && colIndex == 0) {
      deleteDropdownRows(dropdown, [0,1,2,3,4,5,6]);
    }

    // Update runs-table and projects-table counters on filtering
    if (["runs-table", "projects-table", "samples-table"].includes(tableId)) {
      var checkboxes = dropdown.querySelectorAll("input[type='checkbox']");
      checkboxes.forEach(function(checkbox) {
        checkbox.addEventListener("change", function(event) {
          updateTableCount(tableId);
        });
      });
    }

    // Update runs-table on filtering project-table
    //if (tableId === "projects-table") {
    //  var checkboxes = dropdown.querySelectorAll("input[type='checkbox']");
    //  checkboxes.forEach(function(checkbox) {
    //    checkbox.addEventListener("change", function(event) {
    //      filterRunsTable();
    //    });
    //  });
    //}
  }

  // Position the dropdown menu
  var headerRect = clickedElem.getBoundingClientRect();
  dropdown.style.top = (headerRect.bottom + window.scrollY) + "px";
  dropdown.style.left = (headerRect.left + window.scrollX) + "px";

  var clickedIsVisible = !dropdown.classList.contains("collapsed");
  if (visibleDropdown) {
    hideDropdown(visibleDropdown);
  }

  if (!clickedIsVisible) {
    dropdown.classList.remove("collapsed");
    visibleDropdown = [dropdown, tableId];
    event.stopPropagation();

    // Gray-out checkbox values only available for hidden rows
    // Get visible unique values for the column
    var visibleValues = [];
    var rows = table.querySelector("tbody").getElementsByTagName("tr");
    for (let i = 0; i < rows.length; i++) {
      if (!rows[i].classList.contains("collapsed")) {
        var cell = rows[i].children[colIndex].textContent;
        if (!visibleValues.includes(cell)) {
          visibleValues.push(cell);
        }
      }
    }

    // Loop through checkboxes and gray-out
    var checkboxes = dropdown.querySelectorAll('input[type="checkbox"]');
    checkboxes.forEach(function(checkbox) {
      var checkboxRow = checkbox.parentElement;
      if (!(visibleValues.includes(checkbox.value)) && checkbox.checked) {
        checkboxRow.classList.add("grayed-out");
      } else {
        checkboxRow.classList.remove("grayed-out");
      }
    });

    // Focus on the search bar
    dropdown.querySelector("input[type='text']").focus();

    sortGray(dropdown);
    event.stopPropagation();
  }
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
function makeDropdown(tableId, colIdNum, uniqueValues) {
  // Create new dropdown if it doesn't exist
  dropdown = document.createElement("div");
  dropdown.id = tableId + "dropdown_" + colIdNum;
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
    checkbox.id = tableId + colIdNum + value;
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
  dropdowns[tableId][colIdNum] = dropdown;
  return dropdown;
}
function applyFilters(dropdownData, checkedVals, checked) {
  var dropdown = dropdownData[0];
  var tableId = dropdownData[1];
  var table = document.getElementById(tableId);
  var rows = table.getElementsByTagName("tbody")[0].getElementsByTagName("tr");
  //console.log(dropdownData);
  var colIdNum = dropdown.id.split(/[_]+/).pop();
  var colId = tableId + "_" + colIdNum;
  var colIndex = [...table.querySelector("thead").querySelectorAll("th")]
    .map(h => h.querySelector("div").querySelector("div"))
    .indexOf(document.getElementById(colId));

  // Show rows for checked checkbox
  // Not applicable if the changed checkbox was unchecked
  if (checked) {
    // Loop through each row and show if filtered in
    for (var i = 0; i < rows.length; i++) {
      var cells = rows[i].children;
      var showRow = false;
      // If row is to be filtered back in
      if (checkedVals.includes(cells[colIdNum].textContent)) {
        showRow = true;
        // Re-apply remaining filters to the row
        // Loop through each filter checkbox for each dropdown menu
        for (const [index, dropdown] of Object.entries(dropdowns[tableId])) {
          if (index === colIdNum) { continue }
          var checkboxes = dropdown.querySelectorAll('input[type="checkbox"]');
          for (checkbox of checkboxes) {
            if (
             !checkbox.checked &&
             cells[index].textContent === checkbox.value) {
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
  //if (tableId === "runs-table") {
  //  var projTable = document.getElementById("projects-table");
  //  var projRows = projTable.querySelector("tbody").querySelectorAll("tr");
  //  var visibleBPs = [];
  //  for (let i = 0; i < projRows.length; i++) {
  //    var row = projRows[i];
  //    if (!row.classList.contains("collapsed")) {
  //      visibleBPs.push(row.querySelector("td").textContent);
  //    }
  //  }
  //} else {
  //  var visibleBPs = null;
  //}
  // Loop through each row and hide based on filters
  for (var i = 0; i < rows.length; i++) {
    var cells = rows[i].children;
    var showRow = true;
    //if (visibleBPs !== null) {
    //  if (!visibleBPs.includes(cells[1].innerHTML)) {
    //    showRow = false;
    //  }
    //}
    // Loop through each filter checkbox for each dropdown menu
    Object.keys(dropdowns[tableId]).forEach(function(colIdNum) {
      var colId = tableId + "_" + colIdNum;
      var colIndex = [...table.querySelector("thead").querySelectorAll("th")]
        .map(h => h.querySelector("div").querySelector("div"))
        .indexOf(document.getElementById(colId));

      var colDropdown = dropdowns[tableId][colIdNum];
      var colBoxes = colDropdown.querySelectorAll('input[type="checkbox"]');
      colBoxes.forEach(function(checkbox) {
        if (
         !checkbox.checked &&
         cells[colIndex].textContent === checkbox.value) {
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
  Object.keys(dropdowns[tableId]).forEach(function(colIdNum) {
    var colId = tableId + "_" + colIdNum;
    var colIndex = [...table.querySelector("thead").querySelectorAll("th")]
      .map(h => h.querySelector("div").querySelector("div"))
      .indexOf(document.getElementById(colId));

    var toHighlight = false;
    var colDropdown = dropdowns[tableId][colIdNum];
    colBoxes = colDropdown.querySelectorAll('input[type="checkbox"]')
    for (let i = 0; i < colBoxes.length; i++) {
      var checkbox = colBoxes[i];
      var checkboxRow = checkbox.parentElement;
      if (
       !checkbox.checked &&
       !checkboxRow.classList.contains("grayed-out")) {
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
function moveColumn(table, from, to) {
  for (const row of table.rows) {
    const cells = Array.from(row.cells);
    const moved = cells.splice(from, 1)[0];
    cells.splice(to, 0, moved);
    row.innerHTML = "";
    for (const cell of cells) row.appendChild(cell);
  }
}


function makeDraggable(table) {
  const header = table.tHead?.rows[0];
  if (!header) return; // skip tables without header
  let draggingIndex = null;

  for (const th of header.cells) {
    if (th.classList.contains("draggable")) {
      th.draggable = true;
      th.addEventListener("dragstart", e => {
        draggingIndex = th.cellIndex;
        th.classList.add("dragging");
        e.dataTransfer.effectAllowed = "move";
        });

      th.addEventListener("dragend", () => {
        th.classList.remove("dragging");
      });

      th.addEventListener("dragover", e => e.preventDefault());

      th.addEventListener("drop", e => {
        e.preventDefault();
        const targetIndex = th.cellIndex;
        moveColumn(table, draggingIndex, targetIndex);
      });
    }
  }
}
