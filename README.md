# Agent-Based Modeling of Financial-Macroeconomic Systems  
**Empirical Estimation of Wealth, Cost, and Misperception Channels**  

## 📘 Overview

This repository contains the R code developed for the seminar paper _"Agent-Based Modeling of Financial-Macroeconomic Systems: Empirical Estimation of Wealth, Cost, and Misperception Channels". The project implements an integrated agent-based model combining a heuristic-switching macroeconomic framework with a financial market agent-based model. Parameters are estimated via a genetic algorithm using simulated data.

The model is primarily based on:

- **Macroeconomic model**: Kukacka & Zila (2024), *Wealth, Cost, and Misperception: Empirical Estimation of Three Interaction Channels in a Financial-Macroeconomic Agent-Based Model*.
- **Financial model**: Franke & Westerhoff (2011, 2012)

---

## 📂 Repository Structure

```
project-root/
├── main.R                 # Entry point: runs the full pipeline (simulation + GA optimization)
├── setup_test.yaml        # YAML config file specifying model and GA parameters
└── functions/             # All core model components and supporting scripts
    ├── f_FinModel.R       # Financial Agent-Based Model
    ├── f_MacroModel.R     # Macroeconomic Model with heuristic switching
    ├── f_GenData.R        # Pseudo-data generation functions
    ├── f_LoadData.R       # Loads config and initializes simulation parameters
    ├── f_lol.R            # Objective function (negative log-likelihood)
    └── f_GA.R             # Genetic Algorithm implementation
```

---

## 🧠 Model Description

### 🔹 Integrated Framework

The simulation integrates two key components:

1. **Macroeconomic Model (New Keynesian)**  
   - Based on Kukacka & Zila (2024) and De Grauwe's heuristic-switching framework.
   - Three-equation NKM with bounded rationality and adaptive expectations.
   - Agents use three heuristics:  
     - Adaptive (ADA)  
     - Trend-Following (TR)  
     - Learning Anchoring & Adjustment (LAA)  
   - Heuristics are selected via a multinomial discrete choice mechanism based on past forecast accuracy.
   - Incorporates two macro-financial interaction channels:  
     - **Wealth effect** (`c₁`)  
     - **Cost effect** (`c₂`)

2. **Financial Agent-Based Model**  
   - Based on Franke & Westerhoff (DCA-HPM framework).
   - Two trader types: Fundamentalists and Chartists.
   - Dynamic switching based on price misalignment, herding, and predisposition.
   - Introduces a third channel:  
     - **Misperception effect** (`h`), linking macro output gap to perceived asset value.

---

## ⚙️ How It Works

### 1. Configuration  
Edit `setup_test.yaml` to define:
- Initial values for parameters
- Parameter bounds for optimization
- Genetic Algorythm settings (population size, iterations, etc.)

### 2. Simulation Workflow

- `main.R`: Entry point that runs the complete pipeline:
  - Loads configuration (`f_LoadData.R`)
  - Generates pseudo-data (`f_GenData.R`)
  - Runs the genetic algorithm to estimate parameters (`f_GA.R`)
  - Evaluates fitness via negative log-likelihood (`f_lol.R`)
  - Uses models from `f_MacroModel.R` and `f_FinModel.R`

### 3. Optimization Objective

- Implemented in `f_lol.R`
- Computes **negative log-likelihood** by:
  - Simulating both macroeconomic and financial data
  - Smoothing simulated vs. observed distributions with kernel density estimation
- Goal: Find parameters `c₁`, `c₂`, and `h` that minimize the negative log-likelihood

---

## 🧪 Parameter Estimation

- Parameters estimated:
  - `c₁`: Wealth effect (∈ [–1, 1.2])
  - `c₂`: Cost effect (∈ [0, 30])
  - `h`: Misperception effect (∈ [0, 1.5])

- Estimation via **Genetic Algorithm**:
  - Binary-encoded chromosomes (16 bits per parameter)
  - Population: 50 individuals
  - Iterations: 100 generations
  - Crossover: 90%
  - Mutation: ~2%
  - Selection pressure: 1

- Best result achieved:  
  `f(c₁ = 0.3054, c₂ = 16.3202, h = 0.4053) = 0.5911`

---

## 📚 References

- Kukacka, J., & Zila, E. (2024). *Wealth, Cost, and Misperception: Empirical Estimation of Three Interaction Channels in a Financial-Macroeconomic Agent-Based Model*. IES Working Paper 22/2024.
- Bask, M. (2012). Asset price misalignments and monetary policy. International Journal of Finance & Economics 17 (3), 221–241.
- Bernanke, B. S., & Gertler, M. (2000). Monetary policy and asset price volatility. Working Paper 7559, National Bureau of Economic Research.
- Brock , W. A., & Hommes, C. H. (1998). Heterogeneous beliefs and routes to chaos in a simple asset pricing model. Journal of Economic Dynamics & Control 22, 1235–1274.
- Cho , D., & Jang, T. S. (2019). Asset market volatility and new keynesian macroeconomics: A game-theoretic approach. Computational Economics 54 (1), 245–266.
- De Grauwe, P. (2010). Top-down versus bottom-up macroeconomics. CESifo Economic Studies 56 (4), 465–497.
- De Grauwe, P. (2011). Animal spirits and monetary policy. Economic Theory 47 (2), 423–457.
- De Grauwe, P. (2012a). Booms and busts in economic activity: A behavioral explanation. Journal of Economic Behavior & Organization 83 (3), 484–501.
- De Grauwe, P. (2012b). Lectures on Behavioral Macroeconomics, Princeton University Press.
- De Grauwe, P., & Ji, Y. (2019). Behavioural Macroeconomics: Theory and Policy. Oxford University Press.
- De Grauwe, P., & Ji, Y. (2020). Structural reforms, animal spirits, and monetary policies. European Economic Review 124, 103395.
- De Grauwe, P., & Kaltwasser, P. R. (2012). Animal spirits in the foreign exchange market. Journal of Economic Dynamics and Control 36 (8), 1176–1192.
- Franke, R., & Westerhoff, F. (2011). Estimation of a structural stochastic volatility model of asset pricing. Computational Economics 38 (1), 53–83.
- Franke, R., & Westerhoff, F. (2012). Structural stochastic volatility in asset pricing dynamics: Estimation and model contest. Journal of Economic Dynamics and Control 36 (8), 1193–1211.
- Gali, J. (2015). Monetary Policy, Inflation, and the Business Cycle: An Introduction to the New Keynesian Framework and Its Applications (2nd ed.). Princeton: Princeton University Press.
- Kontonikas , A., & Montagnoli, A. (2006). Optimal monetary policy and asset price misalignments. Scottish Journal of Political Economy 53 (5), 636–654.
- Kukacka, J., & Zila, E. (2024). Wealth, Cost, and Misperception: Empirical Estimation of Three Interaction Channels in a Financial-Macroeconomic Agent-Based Model. IES Working Papers, 22/2024
- Lengnick , M., & Wohltmann, H.-W. (2013). Agent-based financial markets and New Keynesian macroeconomics: a synthesis. Journal of Economic Interaction and Coordination 8 (1), 1–32.
- Lengnick , M., & Wohltmann, H.-W. (2013). Optimal monetary policy in a new Keynesian model with animal spirits and financial markets. Journal of Economic Dynamics and Control 64, 148–165.
- Naimzada , A., & Pireddu, M. (2013). Dynamic behavior of real and stock markets with a varying degree of interaction. Milan: University of Milan Bicocca Department of Economics, Management and Statistics.
- Simon, H. A. (1955). A behavioral model of rational choice. The Quarterly Journal of Economics 69 (1), 99–118.
- Tversky , A., & Kahneman, D. (1974). Judgment under uncertainty: Heuristics and biases. Science 185 (4157), 1124–1131.
- Wand, M., & Jones, M. C. (1995). Kernel Smoothing. London: Chapman and Hall Ltd.
- Westerhoff, F. (2012). Interactions between the real economy and the stock market: A simple agent-based approach. Discrete Dynamics in Nature and Society 2012, 504840.



---

