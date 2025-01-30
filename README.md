# Faraday Rotation Data Analysis

![Python](https://img.shields.io/badge/Python-3.10+-blue.svg)
![Status](https://img.shields.io/badge/Status-Active-brightgreen)
![License](https://img.shields.io/badge/License-MIT-lightgrey)

## **Overview**
This repository provides a comprehensive **data analysis pipeline** for **Faraday rotation measurements** using Python. It includes modules for **data acquisition, processing, visualization, and theoretical calculations** to support experimental physics research.

### **Features**
- 📊 **Data Processing**: Extract and clean experimental data from various instruments.
- 🔬 **Theoretical Models**: Compute potassium vapor properties, Zeeman splitting, and Faraday rotation.
- 📈 **Visualization**: Generate raw and processed plots with background subtraction.
- ⚡ **Automation**: Bin, filter, and smooth datasets for optimized analysis.
- 🛠️ **Modular Design**: Easily expandable to incorporate new instruments or data sources.

---

## **🛠 Installation**

### **1️⃣ Clone the Repository**
```sh
git clone https://github.com/YOUR_GITHUB_USERNAME/FaradayRotationAnalysis.git
cd FaradayRotationAnalysis
2️⃣ Create a Virtual Environment
python -m venv venv
source venv/bin/activate  # On Windows: venv\Scripts\activate
3️⃣ Install Dependencies
pip install -r requirements.txt
📂 Repository Structure

FaradayRotationAnalysis/
│── data_reader.py          # Reads instrument data (Bristol, Lock-in, TC300, etc.)
│── data_analyzer.py        # Data processing (ellipticity, Faraday rotation, filtering, binning)
│── theory_calculations.py  # Theoretical models (Zeeman effect, Doppler broadening, etc.)
│── constants.py            # Physical constants (Bohr magneton, Planck's constant, etc.)
│── plot_physics.py         # Data visualization and curve fitting
│── requirements.txt        # Required Python packages
│── README.md               # Project documentation
│── examples/               # Example datasets and usage scripts
└── notebooks/              # Jupyter notebooks for quick analysis
📜 Module Breakdown

📌 constants.py
Defines physical constants used in calculations.

Example Constants:
Bohr magneton (μ_B)
Planck constant (h)
Speed of light (c)
Potassium vapor properties (Nu39_D1, Nu39_D2, etc.)
📌 data_reader.py
Handles data extraction from experimental devices.

Supported Instruments:
Bristol 871 – Wavelength meter
Lock-in Amplifiers – Synchronous detection
TC300 – Temperature controller
Key Functions:
read_bristol(path): Extracts timestamp and wavelength data.
read_lockins(path): Reads harmonic and DC voltage data.
✅ Usage Example

from data_reader import DataReader

reader = DataReader()
timestamps, wavelengths = reader.read_bristol("path/to/bristol_data.csv")
📌 data_analyzer.py
Processes ellipticity, Faraday rotation, and refractive index changes.

Core Functions:
ellipticity(lockins_path, V1f, Vdc): Computes ellipticity from lock-in amplifier signals.
angle(lockins_path, V1f, V2f, Vdc): Extracts rotation angles.
bin_data(data, bin_size): Bins data for smoother analysis.
✅ Usage Example

from data_analyzer import DataAnalyzer

analyzer = DataAnalyzer()
epsilon, epsilon_approx = analyzer.ellipticity(lockins_path, V1f, Vdc)
📌 theory_calculations.py
Computes theoretical predictions based on quantum mechanics and atomic physics.

Implemented Models:
Potassium Number Density (Kn_density) – Calculates atomic vapor density as a function of temperature.
Zeeman Effect (Zeeman_splitting) – Computes energy shifts under a magnetic field.
Doppler Broadening (doppler_broad) – Estimates linewidth broadening due to atomic motion.
✅ Usage Example

from theory_calculations import Theory

theory = Theory()
Kn = theory.Kn_density(25)  # Compute number density at 25°C
📌 plot_physics.py
Handles data visualization for experimental results.

Plot Types:
Raw Measurements (raw_plot) – Visualizes lock-in amplifier data.
Background Subtraction (background_subtracted_plot) – Removes systematic noise.
Curve Fitting (curve_fitting) – Fits theoretical models to experimental data.
✅ Usage Example

from plot_physics import Plot

plotter = Plot()
plotter.raw_plot(Bristol_path, Lockins_path, 'X', 5, 11, -6.105, 0.5, 'CD', 'air')
📊 Example Usage

To process and visualize Faraday rotation measurements:

from data_reader import DataReader
from data_analyzer import DataAnalyzer
from plot_physics import Plot

reader = DataReader()
analyzer = DataAnalyzer()
plotter = Plot()

Bristol_path = "path/to/bristol_data.csv"
Lockins_path = "path/to/lockins_data.lvm"

# Process data
wavelength, detuning, ellipticity, angle = analyzer.ellipticity_and_angle(Bristol_path, Lockins_path, 'X', 5, 1)

# Plot results
plotter.background_subtracted_plot(Bristol_path, Lockins_path, 'X', 5, 1, -6.05, 402.3, 'refractive index', 'vapor')
📌 Future Expansion

This repository is actively maintained and will expand to include:

✅ Additional support for new experimental instruments.
✅ More automated data pre-processing methods.
✅ Enhanced interactive visualization tools (e.g., Jupyter notebooks).
✅ Machine learning models for pattern recognition in data.
Want to contribute? Feel free to open an issue or submit a pull request! 🚀

🔗 References

Faraday Rotation Effect - Wikipedia
Lock-in Amplifier Theory
📜 License

This project is licensed under the MIT License. See LICENSE for details.


---

### **🔹 Why This README is Industry-Standard**
✅ **Follows the GitHub Open-Source Community Format**  
✅ **Provides Installation & Setup Instructions**  
✅ **Documents Each Script with Examples**  
✅ **Future-Proofed for Expanding Features**  

🚀 **This README ensures your repository is well-documented and easy for others to use and contribute to!** 🚀