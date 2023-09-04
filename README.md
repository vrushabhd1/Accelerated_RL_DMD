# Accelerated_RL_DMD
This repository provides a Matlab code for the research paper Accelerated Reinforcement Learning via Dynamic Mode Decomposition.

[![Open in MATLAB Online](https://www.mathworks.com/images/responsive/global/open-in-matlab-online.svg)](https://matlab.mathworks.com/open/github/v1?repo=vrushabhd1/Accelerated_RL_DMD)
# Discrete-Time Reinforcement Learning with Dynamic Mode Decomposition

Welcome to the repository for our work on applying the Dynamic Mode Decomposition (DMD) principle to discrete-time reinforcement learning for solving optimal control problems in a network of subsystems. In this project, we tackle control design as a linear quadratic regulator graphical problem, where the performance function couples the dynamics of subsystems.

## Problem Overview

Optimal control problems in networked systems can be complex, especially for larger networks. To address this challenge, we present two distinct approaches:

### 1. Model-Free Reinforcement Learning

We start with a model-free discrete-time reinforcement learning algorithm. This approach is based on online behaviors without relying on explicit knowledge of system dynamics. While effective, it can be time-consuming for larger networks.

### 2. Model-Free Reinforcement Learning with Dynamic Mode Decomposition

To overcome the computational burden of the first approach, we introduce an efficient model-free reinforcement learning algorithm based on dynamic mode decomposition (DMD). DMD reduces the dimensionality of the measured data while preserving the essential dynamic information of the original network. This algorithm is designed for online implementation, making it suitable for real-time control.

## Validation

We validate our proposed methodology using examples from two different network types:

1. **Consensus Network**: We demonstrate the effectiveness of our approach in achieving consensus within a network of subsystems.

2. **Power System Network**: We apply our methodology to a power system network to illustrate its utility in solving control problems in complex, real-world scenarios.

## MATLAB Code

You can find the MATLAB code for  DMD-based reinforcement learning algorithms in this repository. The code is organized to facilitate easy implementation and experimentation.

Feel free to explore the code and adapt it to your specific control problems. We hope that our work and the provided code will be valuable resources for your research and practical applications.

Please refer to the documentation (https://ieeexplore.ieee.org/abstract/document/10076803) within the repository for detailed instructions on using the code and replicating the results presented in our work. If you have any questions or encounter issues, don't hesitate to open an issue or reach out to us for assistance.

Thank you for your interest in our research, and we look forward to your contributions and feedback!
