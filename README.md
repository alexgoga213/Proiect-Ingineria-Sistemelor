# ⚙️ System Identification of Dynamic Systems (MATLAB / Simulink)

> Building mathematical models of **second-order dynamic systems** (a hydraulic circuit and a nonlinear pendulum) from noisy, quantized input–output data, using **classical non-parametric methods** and **modern parametric methods**.

![MATLAB](https://img.shields.io/badge/MATLAB-R2022b-0076A8?logo=mathworks&logoColor=white)
![Simulink](https://img.shields.io/badge/Simulink-models-E16737)
![System Identification Toolbox](https://img.shields.io/badge/System%20Identification-Toolbox-0076A8)

---

## ✨ Highlights

- **2 physical processes**, simulated in Simulink in realistic conditions: **white noise** on input and output, **DAC/ADC quantization** (9–13 bit) and **saturation**
- **Time-domain identification** from the step response, with time constants estimated by least-squares regression
- **Frequency-domain identification** from a chirp input: resonance, phase lag and damping
- **Parametric identification** with a PRBS input and 6 model structures: ARX, ARMAX, OE, IV4, N4SID, SSEST
- **Rigorous validation** on separate data: residual whiteness/correlation tests, fit comparison and the normalized error **ε<sub>MPN</sub>**

## 🧪 What's inside

### 1. Hydraulic circuit: step response (time domain)
[`setup_circuit_hidraulic.m`](setup_circuit_hidraulic.m) · [`circuit_hidraulic_R2022b.slx`](circuit_hidraulic_R2022b.slx)

- Steady-state gain from the step response: **K = (y<sub>st</sub> − y<sub>0</sub>) / (u<sub>st</sub> − u<sub>0</sub>)**
- **Dominant time constant T₁** by linear regression on `ln(y_st − y(t))` (normal equations / least squares)
- **Non-dominant time constant T₂** by solving the inflection-point equation graphically
- Model: **H(s) = K / ((T₁s + 1)(T₂s + 1))**, simulated in state space and compared to the measured data

### 2. Hydraulic circuit: chirp input (frequency domain)
[`setup_circuit_hidraulic_chirp.m`](setup_circuit_hidraulic_chirp.m)

- Swept-sine (chirp) excitation over two decades around the system bandwidth
- Read the **amplitude ratio M** and the **phase lag φ** at resonance directly from the signals
- Estimate the **natural frequency ω<sub>n</sub>** and the **damping ratio ζ**, giving
  **H(s) = K·ω<sub>n</sub>² / (s² + 2ζω<sub>n</sub>s + ω<sub>n</sub>²)**

### 3. Nonlinear pendulum: parametric identification
[`setup_pendul.m`](setup_pendul.m) · [`pendul_R2022b.slx`](pendul_R2022b.slx)

- Excitation with a **PRBS (SPAB) signal** around an operating point, tuned to the system's dynamics
- Data split into **identification** and **validation** sets, with means removed and decimation applied
- Models compared:

| Structure | Type | MATLAB |
|---|---|---|
| ARX | equation error, least squares | `arx` |
| ARMAX | ARX + colored noise model | `armax` |
| OE | output error | `oe` |
| IV | instrumental variables | `iv4` |
| N4SID | subspace state-space | `n4sid` |
| SSEST | prediction-error state-space | `ssest` |

- Every model is validated with **`resid`** (whiteness of the residuals and their correlation with the input) and **`compare`** (fit on unseen data)

## 📐 Performance metrics

- **Root mean squared error:** J = (1/√N) · ‖y − ŷ‖
- **Normalized error:** ε<sub>MPN</sub> = ‖y − ŷ‖ / ‖y − ȳ‖ · 100%

The full results, plots and discussion are in the report: **[Goga_Alex_Andrei_IS.pdf](Goga_Alex_Andrei_IS.pdf)**.

## ▶️ How to run

Requirements: **MATLAB R2022b+**, **Simulink**, **System Identification Toolbox**

1. Clone the repository and open the folder in MATLAB.
2. Run one of the scripts, for example:
   ```matlab
   setup_pendul
   ```
   The script sets the process parameters, simulates the Simulink model, extracts `t, u, y` and runs the identification step by step (each `%%` section can be run on its own).

> `setup_pendul.m` uses `generate_input_signal` to build the PRBS input. This helper comes with the course materials and is not included here.

## 🧠 Skills demonstrated

`system identification` · `signal processing` · `least squares` · `frequency response` · `state-space models` · `model validation` · `MATLAB/Simulink`

---

<sub>University project · Systems Engineering · author: Goga Alex-Andrei</sub>
