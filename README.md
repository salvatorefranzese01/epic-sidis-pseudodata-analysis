# ePIC SIDIS pseudodata analysis

Code used to process ePIC simulated events and extract four-dimensional SIDIS cross-section pseudodata for charged pions and kaons in the 9 × 275 GeV and 9 × 130 GeV beam configurations. 
The analysis considers the 26.07.1 SIDIS/Pythia6-eic samples. The associated analysis note describes the event selections, binning, efficiency correction, uncertainties, and limitations.

## Files

| File | Role |
| --- | --- |
| `epic_studies.cpp` | Reads simulated events and writes generated and reconstructed electron and hadron ROOT trees, including truth-matched quantities needed for efficiency studies. |
| `epic_extraction_cleaned.cpp` | Extracts cross-section tables for 9 × 275 GeV, using the corresponding MC sample cross sections and a target luminosity of 2.5 fb⁻¹. |
| `epic_extraction_9x130.cpp` | Equivalent extraction for 9 × 130 GeV, with its own MC sample cross sections and a target luminosity of 1.0 fb⁻¹. The function inside this file is also named `epic_extraction_cleaned`. |
| `merge_all.py` | Filters and merges the four separate Q² productions into CSV and YAML tables for each hadron species. |

