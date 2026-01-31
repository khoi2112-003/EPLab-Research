# Argon ICP Global Plasma Model

A physics-based global model for Argon Inductively Coupled Plasma (ICP) thrusters, with an integrated neural network surrogate for real-time predictions.

## Overview

This project implements and validates the global plasma model from:

> **"A Global Model Study of Plasma Chemistry and Propulsion Parameters of a Gridded Ion Thruster Using Argon as Propellant"**  
> Magaldi et al., 2022

The model solves a system of stiff ODEs representing mass and energy conservation within the discharge chamber. Additionally, I have extended this work to develop a machine learning surrogate that can replace computationally expensive ODE integration for real-time control applications.

## Validation Status

| Output | Status | Notes |
|--------|--------|-------|
| Species Densities (Ar, Ar⁺, Arᵐ, Arʳ, Arᵖ) | ✅ Validated | < 5% mean error at 2 mTorr |
| Electron Temperature (Tₑ) | ⚠️ In Progress | ~3-5% overestimation in 0D model |
| Gas Temperature (Tg) | ⚠️ In Progress | Investigating wall-loss sensitivity |

## Project Structure

```
Argon_ICP_Model/
├── Code/
│   ├── Ar_script.m              # Main driver script & training data generation
│   ├── rhs_global.m             # Core ODE system (species + energy balance)
│   ├── RateCoefficients_Ar.m    # Argon reaction rate coefficients
│   ├── adapt_rates.m            # Rate coefficient scaling factors
│   ├── electron_energy_rhs.m    # Electron energy balance equation
│   ├── neutral_energy_rhs.m     # Gas/neutral energy balance equation
│   ├── icp_circuit_power.m      # ICP circuit power absorption model
│   ├── calc_rind.m              # Induced resistance calculation
│   ├── kn_wall_paper.m          # Neutral wall-loss rate (diffusion model)
│   ├── do_sweep.m               # Power sweep utility with warm starts
│   ├── validate_model.m         # Model validation against reference data
│   ├── inverse_model.m          # Parameter optimization via inverse modeling
│   └── NN_Surrogate.ipynb       # Neural network surrogate development
```

## Physics-Based Global Model

The core solver implements a 7-state ODE system:

**State Vector:** `[nAr, nArm, nArr, nArp, ne, Te, Tg]`

- `nAr` — Ground state Argon density (m⁻³)
- `nArm` — Metastable Argon density (m⁻³)
- `nArr` — Resonant Argon density (m⁻³)
- `nArp` — 4p state Argon density (m⁻³)
- `ne` — Electron/ion density (m⁻³)
- `Te` — Electron temperature (eV)
- `Tg` — Gas temperature (K)

### Key Physics

- **20 reaction channels** including ionization, excitation, de-excitation, and radiative transitions
- **ICP circuit model** with induced resistance and coil losses
- **Wall losses** via Bohm flux (ions) and diffusion (neutrals/excited states)
- **Energy balance** coupling absorbed RF power to electron heating and wall losses

## Neural Network Surrogate

To enable rapid parametric studies and model-predictive control, I am developing a data-driven surrogate trained on the physics-based solver output.

### Pipeline

1. **Data Generation:** ODE solver runs across parameter space (Power: 200–1600 W, Pressure: 1–6 mTorr)
2. **Architecture:** Multi-Layer Perceptron (MLP) mapping inputs (P, p) → steady-state outputs
3. **Outputs:** `nAr, nArm, nArr, nArp, ne, Tg` (6 targets; Tₑ excluded due to model limitations)

### Goal

Reduce prediction time from **minutes (ODE integration)** to **milliseconds (NN inference)** for real-time thruster control.

## Quick Start

```bash
# Clone the repository
git clone https://github.com/khoi2112-003/EPLab-Research.git

# Navigate to the directory
cd Argon_ICP_Model/Code

# Run the main solver script (MATLAB)
Ar_script
```

### Requirements

**MATLAB:**
- R2020a or later
- No additional toolboxes required

**Python (for NN surrogate):**
- Python 3.8+
- NumPy, Pandas, Scikit-learn, TensorFlow/PyTorch
- Jupyter Notebook

## File Descriptions

| File | Description |
|------|-------------|
| `Ar_script.m` | Main entry point. Runs power×pressure sweep and generates NN training data |
| `rhs_global.m` | Assembles the full ODE right-hand side for `ode15s` |
| `RateCoefficients_Ar.m` | Temperature-dependent reaction rate coefficients (k₁–k₂₀) |
| `adapt_rates.m` | Applies scaling factors to rate coefficients for sensitivity studies |
| `electron_energy_rhs.m` | Computes dTₑ/dt from power absorption and collisional losses |
| `neutral_energy_rhs.m` | Computes dTg/dt from elastic heating and wall conduction |
| `icp_circuit_power.m` | Models RF power coupling with impedance matching effects |
| `calc_rind.m` | Calculates plasma-induced resistance from Bessel functions |
| `kn_wall_paper.m` | Diffusive wall-loss rate for neutrals (Eq. 8 from paper) |
| `do_sweep.m` | Utility for running power sweeps with solution continuation |
| `validate_model.m` | Compares model outputs against reference data with error metrics |
| `inverse_model.m` | Optimizes model parameters (σ, γ, Rcoil) to minimize validation error |
| `NN_Surrogate.ipynb` | Jupyter notebook for training and evaluating the neural network surrogate |

## Optimized Parameters

These parameters were obtained via inverse optimization and achieve < 5% mean error on all species densities at 2 mTorr:

```matlab
params.Rcoil = 0.543183;           % Coil resistance (Ω)
params.sigma_heating = 4.852e-18;  % Ion-neutral heating cross-section (m²)
params.sigma_nn = 8.237e-19;       % Neutral-neutral collision cross-section (m²)
params.sigma_i = 7.447e-19;        % Ion-neutral collision cross-section (m²)
params.gamma_g = 1.000e-04;        % Ground state wall-loss probability
params.gamma_m = 0.001110;         % Metastable wall-loss probability
params.gamma_r = 0.001945;         % Resonant wall-loss probability
params.gamma_p = 0.005982;         % 4p state wall-loss probability
```

## References

1. Magaldi, B. V., et al. (2022). "A Global Model Study of Plasma Chemistry and Propulsion Parameters of a Gridded Ion Thruster Using Argon as Propellant." *Applied Sciences*, 12(19), 9972.

## Author

**Khoi Nguyen**  
Undergraduate Research Assistant  
Electric Propulsion Laboratory (EPLab)  
University of Illinois at Urbana-Champaign

## License

This project is for academic and research purposes.
