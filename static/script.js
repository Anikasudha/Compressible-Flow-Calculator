document.addEventListener("DOMContentLoaded", function () {
    console.log("🚀 Script Loaded");

    const btnIsen = document.getElementById("isen_btn");
    const btnNor = document.getElementById("nor_btn");
    const btnObl = document.getElementById("obl_btn");
    const btnReset = document.getElementById("reset_btn"); 
    const outputSection = document.getElementById("output"); 
    const isenForm = document.getElementById("isentropic");
    const norForm = document.getElementById("normal");
    const oblForm = document.getElementById("oblique");

    function hideAllForms() {
        isenForm.style.display = "none";
        norForm.style.display = "none";
        oblForm.style.display = "none";
    }

    function showForm(form, button) {
        hideAllForms();
        form.style.display = "flex";
        setActiveButton(button);
        localStorage.setItem("selectedFlow", button.id); // Store last selection
    }

    function setActiveButton(clickedButton) {
        [btnIsen, btnNor, btnObl].forEach(btn => btn.classList.remove("active"));
        clickedButton.classList.add("active");
    }

    btnIsen.addEventListener("click", () => showForm(isenForm, btnIsen));
    btnNor.addEventListener("click", () => showForm(norForm, btnNor));
    btnObl.addEventListener("click", () => showForm(oblForm, btnObl));

    // Loading last form
    const savedFlow = localStorage.getItem("selectedFlow");
    if (savedFlow === "nor_btn") {
        showForm(norForm, btnNor);
    } else if (savedFlow === "obl_btn") {
        showForm(oblForm, btnObl);
    } else {
        showForm(isenForm, btnIsen);
    }

    // Reset button
    btnReset.addEventListener("click", function () {
        const currentForm = document.querySelector(".form:not([style*='display: none'])");
        if (currentForm) {
            currentForm.querySelectorAll("input").forEach(input => input.value = "");
            currentForm.querySelectorAll("select").forEach(select => select.value = "");
        }
        outputSection.style.display = "none";
    });

    const limits = {
        "M": { min: 0.1, max: 25, placeholder: "Mach number (0.1 - 25)" },
        "Mstar": { min: 0.1, max: 25, placeholder: "Characteristic Mach Number (0.1 - 25)" },
        "P0_P": { min: 1, max: 100, placeholder: "Pressure Ratio (>1)" },
        "rho0_rho": { min: 1, max: 10, placeholder: "Density Ratio (>1)" },
        "T0_T": { min: 1, max: 10, placeholder: "Temperature Ratio (>1)" },
        "M1": { min: 0.1, max: 25, placeholder: "Upstream Mach Number (0.1 - 25)" },
        "M2": { min: 0.1, max: 25, placeholder: "Downstream Mach Number (0.1 - 25)" },
        "P2_P1": { min: 1, max: 100, placeholder: "Pressure Ratio (>1)" },
        "rho2_rho1": { min: 1, max: 10, placeholder: "Density Ratio (>1)" },
        "T2_T1": { min: 1, max: 10, placeholder: "Temperature Ratio (>1)" },
        "theta": { min: 0, max: 90, placeholder: "Deflection Angle (0 - 90°)" },
        "beta": { min: 0, max: 90, placeholder: "Shock Angle (0 - 90°)" },
        "gamma": { min: 1, max: 1.67, placeholder: "Ratio of Specific Heats (1 - 1.67)" },
        "R": { min: 100, max: 400, placeholder: "Gas Constant (100 - 400 J/kg·K)" }
    };

    function updateInput(selectElement, inputElement, errorElement) {
        let selectedValue = selectElement.value;
        console.log(`Dropdown changed: ${selectedValue}`);
        if (limits[selectedValue]) {
            let { min, max, placeholder } = limits[selectedValue];
            inputElement.min = min;
            inputElement.max = max;
            inputElement.placeholder = placeholder;
            inputElement.value = ""; 
            errorElement.textContent = ""; 
            console.log(`Updated input: min=${min}, max=${max}, placeholder=${placeholder}`);
        } else {
            console.log("No matching limit found!");
        }
    }

    function validateInput(inputElement, errorElement) {
        let min = parseFloat(inputElement.min);
        let max = parseFloat(inputElement.max);
        let value = parseFloat(inputElement.value);
        if (isNaN(value) || value < min || value > max) {
            errorElement.innerHTML = `Value must be between ${min} - ${max}`;
            inputElement.style.border = "2px solid red";
        } else {
            errorElement.textContent = "";
            inputElement.style.border = "";
        }
    }

    document.querySelectorAll("select").forEach(select => {
        let input = select.parentElement.querySelector("input");
        if (input) {
            let errorMessage = document.createElement("span");
            errorMessage.style.color = "red";
            errorMessage.style.fontSize = "12px";
            errorMessage.classList.add("error-message");
            input.parentElement.appendChild(errorMessage);
            updateInput(select, input, errorMessage);
            select.addEventListener("change", function () {
                updateInput(select, input, errorMessage);
            });
            input.addEventListener("blur", function () {
                validateInput(input, errorMessage);
            });
        }
    });
});