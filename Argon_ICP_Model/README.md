In this directory is my code for the replication and validation of the global plasma model for argon inductively coupled plasma (ICP) thruster from the paper "A Global Model Study of Plasma Chemistry and Propulsion Parameters of a Gridded Ion Thruster Using Argon as Propellant" by (Magaldi et al., 2022). So far, I have managed to produce the species densities matching those from the paper, while I am still debugging the electron temperature and gas temperature plots. 
Furthermore, I have independently extended this project to combine computational physics with machine learning. Currently, I am using the physics-based model to generate high-fidelity training data for a neural network (NN) surrogate. THe goal of this personal initiative is to demonstrate how data-driven surrogates can replace computationally expensive ODE integration, reducing prediction time from minutes to milliseconds for real-time control applications.

1. Physics-Based Global Model

The core solver implements a set of stiff Ordinary Differential Equations (ODEs) representing the conservation of mass and energy within the discharge chamber.
Key Capabilities:
- Species Balance: Tracks the time-evolution of Ground state Argon, Metastable Argon ($Ar_m$), Resonant Argon ($Ar_r$), and Ions ($Ar^+$).
- Energy Balance: Solves the electron temperature ($T_e$) and gas temperature ($T_g$) equations coupled with power deposition logic.
- Validation: The model successfully reproduces the species density profiles reported in the reference literature.

Current Validation Status:
- ✅ Species Densities: Successfully validated against Magaldi et al. (2022). The reaction rates and particle balance are consistent with published results.
- ⚠️ Energy Balance: Currently investigating discrepancies in the electron and gas temperature evolution ($T_e$, $T_g$) relative to the COMSOL outputs in the reference paper. This investigation focuses on the sensitivity of wall-loss probabilities and power coupling efficiency.

2. Neural Network Surrogate (In Development)

To enable rapid parametric studies and model-predictive control, I am developing a data-driven surrogate model trained on the physics-based solver's output.
The Pipeline:
  - Data Generation: The ODE solver is run across a sweep of input parameters (Pressure: 1–50 mTorr, Power: 100–1000 W).
  - Architecture: A Multi-Layer Perceptron (MLP) is trained to map inputs $(P, \mathcal{P}_{abs})$ directly to steady-state plasma properties $(n_{species}, T_e)$.
  - Objective: Minimize the Mean Squared Error (MSE) between the ODE solution and the NN prediction.

To run,

# Clone the repository
git clone https://github.com/khoi2112-003/EPLab-Research.git

# Navigate to the directory
cd Argon_ICP_Model/Code

# Run the main solver script
Ar_script
