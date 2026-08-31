# iLQRSolver

A C++ implementation of the **iterative Linear Quadratic Regulator (iLQR)** trajectory optimization algorithm for rigid-body systems, built on top of [Pinocchio](https://github.com/stack-of-tasks/pinocchio) for dynamics and kinematics computations.

This project was developed as part of an undergraduate research (Iniciação Científica) project on trajectory optimization for robotic manipulators and underactuated systems.

## Overview

iLQR is a trajectory optimization method for finite-horizon, discrete-time nonlinear systems. Given an initial state and a goal state, it iteratively linearizes the system dynamics around a nominal trajectory and solves a local LQR problem (via a Riccati-like backward pass) to compute a feedback correction. A forward rollout with line search then updates the trajectory, and the process repeats until convergence.

This implementation:

- Uses **Pinocchio** to compute the rigid-body dynamics (via the Articulated Body Algorithm, `aba`) and their analytical derivatives (`computeABADerivatives`) for any robot model loaded from a URDF file (or built programmatically).
- Discretizes the continuous-time dynamics with a **4th-order Runge-Kutta (RK4)** integrator, and propagates the corresponding state-transition Jacobians (A, B matrices) through the same RK4 stages.
- Solves the backward pass with Riccati-style recursions, including regularization of the control-cost Hessian (`Guu`) to keep it positive definite.
- Performs a backtracking **line search** on the forward rollout to guarantee cost decrease at each iteration.
- Supports an **actuator/underactuation selection matrix** `S` (useful for systems such as the Acrobot, where only a subset of joints is actuated) and optional **viscous joint damping**.

## Repository Structure

```
.
├── src/
│   ├── iLQRSolver.cpp          # Solver implementation
│   └── include/iLQRSolver.hpp  # Solver class declaration
├── test/
│   ├── acrobot.cpp             # Swing-up example (underactuated 2-link Acrobot)
│   └── iiwa_14.cpp             # Point-to-point reaching example (KUKA iiwa 14 arm)
├── utils/
│   ├── iiwa_14.urdf            # KUKA iiwa 14 URDF model
│   ├── iiwa_14_sim.urdf        # Variant used for simulation
│   └── iiwa_14_sim_adulterado.urdf  # Modified/perturbed variant (model-mismatch experiments)
├── plot_acrobot.py             # Plots the Acrobot trajectory/torque CSV output
├── plot_iiwa.py                # Plots the iiwa 14 trajectory/torque CSV output
└── Makefile
```

## Dependencies

- A C++11-capable compiler (`g++`)
- [Pinocchio](https://github.com/stack-of-tasks/pinocchio) (rigid-body dynamics library), discoverable via `pkg-config`
- Eigen (pulled in as a Pinocchio dependency)
- Boost (`boost_system`), `pthread`, and `urdfdom` (for URDF parsing)
- Python 3 with `pandas`, `numpy`, and `matplotlib` (only needed to run the plotting scripts)

## Building

The provided `Makefile` builds both example executables:

```bash
make
```

This produces `bin/acrobot` and `bin/iiwa_14`, compiled against `src/iLQRSolver.cpp` with the flags Pinocchio reports through `pkg-config --cflags/--libs pinocchio`.

> **Note:** The `bin/` and `obj/` directories are expected to exist beforehand (only `.gitkeep` placeholders are tracked in version control).

To build and run a single example manually:

```bash
mkdir -p obj bin
make
./bin/acrobot
./bin/iiwa_14
```

## Usage

The `iLQRSolver` class wraps a Pinocchio `Model` (built either directly or from a URDF path) and exposes a simple setup/solve API:

```cpp
iLQRSolver ilqr(urdf_path, h, Nt);        // or: iLQRSolver ilqr(model, h, Nt);

ilqr.setCosts(Q, R, Qn);                  // stage-cost state/control weights and terminal weight
ilqr.setBoundaryConditions(x0, xgoal);    // initial state and target state
ilqr.setViscousDamping(c);                // optional
ilqr.setGravity(g_z);                     // optional

bool converged = ilqr.solve(max_iters, tol);

auto xtraj = ilqr.getXTraj();  // state trajectory
auto utraj = ilqr.getUTraj();  // control trajectory
auto Ktraj = ilqr.getKTraj();  // time-varying feedback gains
```

- `h` is the integration time step and `Nt` is the number of time steps in the horizon (`Tfinal = (Nt - 1) * h`).
- State vectors follow Pinocchio's convention: `x = [q; v]`, where `q` are generalized positions and `v` are generalized velocities, so `Nx = model.nq + model.nv`.
- The optional selection matrix `S` maps the control vector `u` to joint torques (`tau = S * u`), which is what allows modeling underactuated systems like the Acrobot (see `test/acrobot.cpp`).

### Examples

- **`test/acrobot.cpp`** — swing-up of a 2-link underactuated Acrobot (only the shoulder joint is actuated), built directly in code (no URDF). Exports the resulting trajectory to `trajetoria.csv`, which can be visualized with `plot_acrobot.py`.
- **`test/iiwa_14.cpp`** — point-to-point joint-space reaching task for a 7-DOF KUKA iiwa 14 arm, loaded from `utils/iiwa_14.urdf`. Exports states, controls, and feedback gains to `trajetoria_iiwa.csv`, visualizable with `plot_iiwa.py`.

To visualize a result after running an example:

```bash
python3 plot_acrobot.py     # after running ./bin/acrobot
python3 plot_iiwa.py        # after running ./bin/iiwa_14
```

## Known Limitations / Work in Progress

- The include list in `iLQRSolver.hpp` is broader than strictly necessary and may be trimmed in the future.
- Some cost configurations for the iiwa 14 example (commented out in `test/iiwa_14.cpp`) converge but can produce torques beyond realistic actuator limits — actuator/torque constraints are not currently enforced by the solver.
- Some in-code comments and console output are in Portuguese, reflecting the project's origin as an academic research exercise.

## License

This project is licensed under the MIT License — see [LICENSE](LICENSE) for details.
