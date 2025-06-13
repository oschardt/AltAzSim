# Alt Az Simulator

This project simulates multi-session sky scans using the Green Bank Telescope (GBT), specifically tailored for mapping the entire night sky. It is **not** intended to be a universal telescope simulation tool—deliberate simplifications and assumptions were made to reflect the capabilities and constraints of the GBT. These design choices were carefully reviewed to ensure realistic behavior within the GBT's operational envelope. That said, for telescopes that function similarly to the GBT, only slight—or even no—modifications may be needed.

Significant effort went into making the user-facing interface clean, intuitive, and well-organized. For a demonstration of how to use this package, see `GBT_Project_Observation.ipynb`, where we simulate a six-night scanning observation.

## Dependencies

This project uses the following Python packages:

- **astropy** – time handling, coordinate transformations, observatory location  
- **healpy** – HEALPix map creation, pixel queries, and sky visualization  
- **numpy** – numerical operations and array handling  
- **matplotlib** – plotting and PDF export  
- **jupyter** – for interactive demo and experimentation (not required for backend code)

Install all dependencies with:

```bash
pip install -r requirements.txt
```
