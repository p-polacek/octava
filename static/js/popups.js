function showHelp(event, field, helpText) {
  // Show help popup when clickig help button
  const helpPopup = document.getElementById("help-popup");
  if (helpText == null) {
    helpText = helpTexts[field];
  }
  helpPopup.innerHTML = helpText;

  helpPopup.classList.remove("collapsed");
  const popupWidth = helpPopup.offsetWidth + 20;
  const popupHeight = helpPopup.offsetHeight + 20;
  const windowWidth = window.innerWidth;
  const windowHeight = window.innerHeight;
  let topPos = event.clientY + window.scrollY;
  let leftPos = event.clientX + window.scrollX;

  // Ensure the popup doesn't extend beyond the bottom or left of the window
  if (topPos + popupHeight > windowHeight + window.scrollY) {
    topPos = Math.max(0, windowHeight + window.scrollY - popupHeight);
  }
  if (leftPos + popupWidth > windowWidth + window.scrollX) {
    leftPos = Math.max(0, windowWidth + window.scrollX - popupWidth);
  }
  helpPopup.style.top = topPos + "px";
  helpPopup.style.left = leftPos + "px";
  event.stopPropagation();
  event.preventDefault();
}
function showDropdown(event, clickedElem) {
  // Shows dropdown below the clicked element
  var dropdown = clickedElem.parentElement.querySelector(".dropdown-menu");
  var isHidden = dropdown.classList.contains("collapsed");
  hideDropdown(visibleDropdown);
  if (isHidden) {
    dropdown.classList.remove("collapsed");
    visibleDropdown = [dropdown, null];
    event.stopPropagation();
    // Style filter dropdown to popup directly below the clicked element
    var clickedRect = clickedElem.getBoundingClientRect();
    dropdown.style.top = (clickedRect.bottom + window.scrollY) + "px";
    dropdown.style.left = (clickedRect.left + window.scrollX) + "px";
    // Focus on the search bar
    var searchBar = dropdown.querySelector("input[type='text']");
    if (searchBar != null) {
      searchBar.focus();
    }
  }
}
function filterDropdown(dropdownData) {
  // Filter checkboxes in dropdown menu
  if (dropdownData != null) {
    var dropdown = dropdownData[0];
    var term = dropdown.querySelector("input[type='text']").value;
    if (term[0] === "\\") {
      term = term.slice(1);
      var options = dropdown.querySelectorAll(".checkbox-row");

      if (term.length === 0) {
        for (let i = 0; i < options.length; i++) {
          options[i].classList.remove("collapsed");
        }
      } else {
        try {
          re = new RegExp(term);
        }  catch (e) {
          return;
        }

        for (let i = 0; i < options.length; i++) {
          var label = options[i].textContent || options[i].innerText;
          options[i].classList.toggle("collapsed", !re.test(label));
        }
      }
    } else {
      term = term.toUpperCase();
      var term = term.split(" ").map(i => {
        return i.trim();
      }).filter(i => {
        return i !== "";
      });
      var options = dropdown.querySelectorAll(".checkbox-row");  
      if (term.length === 0) {
        for (let i = 0; i < options.length; i++) {
          options[i].classList.remove("collapsed");
        } 
      } else {
        for (let i = 0; i < options.length; i++) {
          var label = options[i].textContent || options[i].innerText;
          options[i].classList.add("collapsed");
          for (let j = 0; j < term.length; j++) {
            if (label.toUpperCase().indexOf(term[j]) > -1) {
              options[i].classList.remove("collapsed");
              break;
            }
          }
        }
      }
    }
  }
}
function hideDropdown(dropdownData) {
  // Hide dropdown menu, reset search bar and filtering, hide help popup
  if (dropdownData != null) {
    var dropdown = dropdownData[0];
    dropdown.classList.add("collapsed");
    var searchBar = dropdown.querySelector("input[type='text']");
    if (searchBar != null) {
      searchBar.value = "";
      filterDropdown(dropdownData);
    }
    var counter = dropdown.parentElement.querySelector("span");
    if ((counter != null) &&
        (counter.innerHTML.includes("selected"))) {
      updateCount(dropdown);
    }
    visibleDropdown = null;
  }
  document.getElementById("help-popup").classList.add("collapsed");
}
function updateCount(dropdown, counter) {
  if (counter == null) {
    var counter = dropdown.parentElement.querySelector("span.counter");
  }
  if (counter != null) {
    var checkboxes = dropdown.querySelectorAll("input[type='checkbox']");
    var count = 0;
    checkboxes.forEach(function(checkbox) {
      if (checkbox.checked) { count++ }
    });
    counter.innerHTML = count + " selected";
  }
}
function updateTracker(dropdown, choice) {
  var tracker = dropdown.parentElement.querySelector("span.tracker");
  if (tracker != null) {
    var choices = dropdown.querySelectorAll("input[type='radio']:checked");
    if (choices.length === 0) {
      choice = "";
    } else {
      label = dropdown.querySelector('label[for="' + choices[0].id + '"]');
      choice = label.textContent;
    }
    tracker.innerHTML = "&nbsp;&nbsp;" + choice;
    hideDropdown(visibleDropdown);
  }
}
function toggleAll(dropdownData, filter) {
  var dropdown = dropdownData[0];
  var tableId = dropdownData[1];
  var checkboxes = dropdown.querySelectorAll('input[type="checkbox"]');
  var visibleCheckboxes = [];
  var visibleLabels = [];
  for (var i = 0; i < checkboxes.length; i++) {
    var checkboxRow = checkboxes[i].parentElement;
    if (!checkboxRow.classList.contains("collapsed")) {
      visibleCheckboxes.push(checkboxes[i]);
      visibleLabels.push(checkboxes[i].value);
    }
  }
  if (visibleCheckboxes[0] != null) {
    // If first checkbox is checked, all are unchecked and vice versa
    var selectAll = !visibleCheckboxes[0].checked;
    visibleCheckboxes.forEach(function(checkbox) {
      checkbox.checked = selectAll;
    });
    if (filter === true) {
      var checkedVals = [];
      visibleCheckboxes.forEach(function(checkbox) {
          if (checkbox.checked) {
            checkedVals.push(checkbox.value);
          }
      });
      applyFilters(dropdownData, checkedVals, selectAll);
    }
  }
  updateCount(dropdown);
  if (["runs-table", "projects-table", "samples-table"].includes(tableId)) {
    updateTableCount(tableId);
  }
  if (tableId === "projects-table") {
    filterRunsTable();
  }
}
function deleteDropdownRows(dropdown, rowsToDelete) {
  var checkboxRows = dropdown.querySelectorAll(".checkbox-row");
  rowsToDelete.forEach( function(row) {
    dropdown.removeChild(checkboxRows[row]);
  });
}
function popupClickedRetain(popup) {
  // Prevent hiding the popup when clicked within
  popup.addEventListener("click", function(event) {
    event.stopPropagation();
  });
}
function showGroupContents(checkboxRow) {
  checkboxRow.addEventListener("contextmenu", function(event) {
    var items = JSON.parse(checkboxRow.children[0].value);
    showHelp(event, null, "<b>" + checkboxRow.innerText + ":</b><br>" + items.join("<br>"));
  });
}


document.addEventListener("DOMContentLoaded", function() {
  document.querySelectorAll(".group-option").forEach(opt => {
    showGroupContents(opt);
  });
  document.querySelectorAll(".scrollable").forEach( container => {
    container.addEventListener("scroll", () => hideDropdown(visibleDropdown));
  });
  document.querySelectorAll(".dropdown-menu").forEach(dropdown => {
    // Count checked boxes per dropdown
    var checkboxes = dropdown.querySelectorAll("input[type='checkbox']");
    if (checkboxes.length != 0) {
      updateCount(dropdown);
    }
    // Prevent hiding the popup when clicked within
    popupClickedRetain(dropdown);
  });
  // Style auxilliary buttons to round shape
  document.querySelectorAll(".aux-button").forEach(button => {
    const clone = button.cloneNode(true);
    clone.style.visibility = "hidden";
    clone.style.position = "absolute";
    document.body.appendChild(clone);
    var shapeLim = Math.max(clone.offsetWidth, clone.offsetHeight);
    shapeLim = Math.floor(shapeLim / 2) * 2;
    button.style.height = shapeLim + "px";
    button.style.width = shapeLim + "px";
    document.body.removeChild(clone);
  });
});
document.addEventListener("keydown", function(event) {
  if (event.which === 27) {
    hideDropdown(visibleDropdown);
  }
});