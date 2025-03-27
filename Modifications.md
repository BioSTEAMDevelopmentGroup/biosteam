# **Installation Guide**  

1. **Clone this repository**  
   ```sh
   git clone https://github.com/omarbay/biosteam.git
   cd biosteam
   ```

2. **Ensure Python Version < 3.13**  
   Check your Python version:  
   ```sh
   python --version
   ```  
   If you have Python 3.13 or later, install Python 3.12.9 and ensure it's used in the virtual environment.  

3. **Create and Activate a Virtual Environment**  
   ```sh
   python -m venv venv  
   source venv/bin/activate  # On macOS/Linux  
   venv\Scripts\activate     # On Windows  
   ```  

4. **Install a Compatible Thermosteam Version for This Biosteam Version**  
   Run the following commands:  
   ```sh
   git submodule init  
   git submodule update  # (This may take a long time. To speed it up, consider removing 'Bioindustrial-park' and 'How2STEAM' from the .gitmodules file after initialization.)
   cd thermosteam  
   pip install .
   cd ..
   ```

5. **Install Dependencies from `requirements.txt`**  
   ```sh
   pip install -r requirements.txt
   ```

6. **Run the Example Script**  
   ```sh
   python binary-distillation-example.py
   ```

# Issues Encountered  

### **Biosteam Module Compatibility**  
**Issue:** Some requirements in the Biosteam module require a Python version lower than 3.13 (latest: 3.13.2).  
**Fix:** Install Python 3.12.9.  

### **Thermosteam Version Compatibility**  
**Issue:** The Thermosteam version available on pip (0.51.5) is not compatible with the Biosteam version on GitHub, which requires Thermosteam 0.51.6.  
**Fix:** Clone the Thermosteam submodule from Git and uninstall the pip version.  
**Commands:**  
```sh
git submodule init  
git submodule update  # (This may take too long. To speed up future updates, consider removing 'Bioindustrial-park' and 'How2STEAM' from the .gitmodules file after initialization.)
cd thermosteam  
pip install .
```

---

# **Modifications**  

### **Modified Functions (All in `distillation.py`)**  
- `Distillation.__init__`  
- `Distillation._set_distillation_product_specifications`  
- `BinaryDistillation._run`  
- `Distillation._run_binary_distillation_mass_balance`  
- `Distillation._update_distillate_and_bottoms_temperature`  
- `BinaryDistillation._design`  
- `BinaryDistillation._run_McCabeThiele`  

### **New Functions**  
- **`BinaryDistillation._calculate_stages`**  
  - Extracts the code responsible for calculating the number of stages from `BinaryDistillation._run_McCabeThiele`.  
  - Returns the computed number of stages.  

- **`BinaryDistillation._mccabe_thiele_find_X_bot`**  
  - Performs a binary search on `x_bot`, validating it based on the provided number of stages.  

- **`BinaryDistillation._mccabe_thiele_find_Y_top`**  
  - Performs a binary search on `y_top`, validating it based on the provided number of stages.  

### **Modifications** 

### **1. `Distillation.__init__`**  
- Added `N_stages` as a parameter in the class constructor (default: `None`).  
- Passes `N_stages` to `Distillation._set_distillation_product_specifications`.  

### **2. `x_bot` Setter**  
- Added an assertion to allow `x_bot == None`.  

### **3. `y_top` Setter**  
- Added an assertion to allow `y_top == None`.  

### **4. `Distillation._set_distillation_product_specifications`**  
- Added `N_stages` as a parameter.  
- Introduced cases:  
  - If `N_stages` and `y_top` are provided but `x_bot` is not.  
  - If `N_stages` and `x_bot` are provided but `y_top` is not.  

### **5. `BinaryDistillation._run`**  
- Added a condition to prevent calling `Distillation._update_distillate_and_bottoms_temperature` if `x_bot` or `y_top` is `None`.  

### **6. `Distillation._run_binary_distillation_mass_balance`**  
- Added a condition to check if `x_bot` and `y_top` are not `None` before performing mass balance calculations.  
- Replaced `else` with `elif spec != 'Composition'`.  
- Moved the code from the previous `else` block into the above two cases.  

### **7. `BinaryDistillation._run_McCabeThiele`**  
- Now handles three cases:  
  - **Case 1:** If `N_stages` is `None` but `x_bot` and `y_top` are defined → Calls `BinaryDistillation._calculate_stages`.  
  - **Case 2:** If `x_bot` is `None` but `N_stages` and `y_top` are defined → Calls `BinaryDistillation._mccabe_thiele_find_X_bot`.  
  - **Case 3:** If `y_top` is `None` but `N_stages` and `x_bot` are defined → Calls `BinaryDistillation._mccabe_thiele_find_Y_top`.  
