# iLQRSolver

A C++ implementation of the **iterative Linear Quadratic Regulator (iLQR)** for trajectory optimization of rigid-body robots.

The solver is organized around interchangeable dynamics backends: the iLQR algorithm is written once, while the computation of the robot dynamics and of their gradients, typically a large share of the computation time [1], can be swapped. This makes it possible to compare a general-purpose baseline built on [Pinocchio](https://github.com/stack-of-tasks/pinocchio) with specialized, hardware-optimized implementations.

> This repository is part of an undergraduate research project (*Iniciação Científica*).

## Motivation and goals

The gradient of the rigid-body dynamics is a central operation in iLQR and in many other planning and control algorithms. Plancher et al. [1] report that it accounts for roughly 30% to 90% of the total computation time of typical nonlinear MPC implementations, and propose optimized implementations of it for a CPU, a GPU and an FPGA.

The goals of this project are to:

1. integrate those optimized dynamics-gradient implementations (CPU and GPU) into an iLQR solver;
2. validate them in simulation, using a Pinocchio-based solver as the baseline;
3. use them to control a real KUKA LBR iiwa 14 manipulator.

See the [Roadmap](#roadmap) for the current status.

## Solver backends

All backends share the same algorithm and the same public API. They differ only in how they evaluate the continuous-time dynamics and its Jacobians.

| Backend | Class | Dynamics and gradients | Supported robots | Parallelism | Status |
|---|---|---|---|---|---|
| **Pinocchio** (baseline) | `iLQRSolverPinocchio` | Pinocchio's Articulated Body Algorithm and its analytical derivatives [3] | Any model loaded from a URDF (or built in code) with `nq = nv`; optional actuation matrix, viscous damping and gravity | Serial | Working |
| **Optimized** (CPU) | `iLQRSolverOptimized` | Rigid-body dynamics and analytical gradient code from Plancher et al. [1] | KUKA LBR iiwa 14 only (7 DoF, hard-coded) | OpenMP across knot points | Working |
| **GPU** (CUDA) | `iLQRSolverGPU` | n/a | n/a | n/a | Skeleton only, not built yet |

## How the solver works

The solver implements iLQR [4, 5]: it finds a control sequence that drives the system from an initial state `x0` towards a goal state `xgoal` over a horizon of `Nt` knot points spaced `h` seconds apart. The state is `x = [q; v]` (joint positions and velocities) and the control `u` is a vector of joint torques. The cost is quadratic around the goal: a running cost penalizing state error (weight `Q`) and control effort (weight `R`), plus a terminal cost (weight `Qn`).

Each iteration of the algorithm has three steps:

1. **Linearization.** The continuous dynamics are discretized with a 4th-order Runge-Kutta (RK4) integrator and linearized along the current trajectory. The Jacobians are propagated through the RK4 stages, so they are the exact derivatives of the discretized dynamics. As in any iLQR, second derivatives of the dynamics are ignored (unlike full DDP).
2. **Backward pass.** A Riccati-like recursion, from the last knot point to the first, computes a feedforward correction `d_k` and a feedback gain `K_k` for every knot point. If the control-cost Hessian is not positive definite, a regularization term is added to its diagonal.
3. **Forward rollout with line search.** A new trajectory is simulated with the nonlinear dynamics using `u_k - alpha * d_k - K_k * (x_new_k - x_k)`. The step size `alpha` starts at 1 and is halved until the cost decreases sufficiently.

The iteration stops when the largest feedforward correction, `max |d_k|`, falls below `tol` (default `1e-3`) or when `max_iters` (default 1000) is reached.

### Code organization

`iLQRSolverBase` implements everything above. A backend only has to implement two virtual functions:

```cpp
// Continuous-time dynamics: x' = f(x, u)
Eigen::VectorXd continuousDynamics(const Eigen::VectorXd& x, const Eigen::VectorXd& u);

// Jacobians of the continuous dynamics: Jx = df/dx and Ju = df/du
void getContinuousJacs(const Eigen::VectorXd& x, const Eigen::VectorXd& u,
                       Eigen::MatrixXd& Jx, Eigen::MatrixXd& Ju);
```

and set `Nx` and `Nu` before calling `iLQRSolverBase::initializeMemory()`. Adding a new backend (such as the GPU one) means filling in these two functions.

## Parallelism and MPC

**Parallel derivatives.** The linearization of the dynamics at each knot point is independent of the others. When `use_parallel_derivatives` is enabled, the backward pass computes all of them in an OpenMP parallel loop. This is the coarse-grained parallelism (many independent gradient evaluations) described in [1]. The Optimized backend enables it; the Pinocchio backend does not, because it reuses a single `pinocchio::Data` workspace and is not safe to call concurrently as written. The number of threads is currently the compile-time constant `N_THREADS = 4` in `src/iLQRSolverBase.cpp`.

**MPC-like simulation.** The `*_mpc` examples run a receding-horizon loop in simulation, in a real-time-iteration style: at each control step the solver is warm-started with the previous solution (`advanceMPCStep`), runs a **single** iLQR iteration (`solve(1)`), and the first control is applied. The examples report the average solve time per control step.

This test is preliminary and is not a sufficient validation of the solver for MPC. It only asks the arm to reach a single fixed goal, so the receding horizon never has to follow a time-varying reference. The simulated state is the solver's own predicted next state, so the plant is exactly the model (no model mismatch, noise or disturbances), and the simulation neither enforces a time budget per control step nor accounts for the computation delay. The planned test follows Plancher and Kuindersma [2] (Section 5.3): the manipulator has to track a moving reference trajectory (a figure-eight path for the end effector in [2]), with the solver re-planning at every control step within a fixed time budget (10 ms in [2]) in a simulation that runs in real time. See the [Roadmap](#roadmap).

## Repository structure

```
.
├── include/
│   ├── iLQRSolverBase.hpp           # Backend-independent iLQR algorithm
│   ├── iLQRSolverPinocchio.hpp      # Pinocchio backend
│   ├── iLQRSolverOptimized.hpp      # Optimized CPU backend (iiwa 14)
│   └── iLQRSolverGPU.hpp            # GPU backend (skeleton)
├── src/
│   ├── iLQRSolverBase.cpp
│   ├── iLQRSolverPinocchio.cpp
│   ├── iLQRSolverOptimized.cpp
│   └── iLQRSolverGPU.cu             # Skeleton, not part of the Makefile yet
├── test/                            # Example programs, one binary each (see below)
├── utils/
│   ├── dynamicsGradient.h           # Third-party code from [1] (MIT license);
│   ├── math.h                       #   see "Third-party code" below
│   ├── minv.h                       #
│   └── iiwa_14.urdf                 # KUKA LBR iiwa 14 model, also from [1] (used by the Pinocchio backend)
├── scripts/
│   ├── plot_acrobot.py              # Plots the CSV written by the acrobot example
│   └── plot_iiwa.py                 # Plots the CSV written by the iiwa examples
├── Makefile
├── LICENSE
└── THIRD_PARTY_NOTICES.md
```

## Dependencies

- Linux with `g++` (C++11, OpenMP) and `pkg-config`
- [Eigen 3](https://eigen.tuxfamily.org)
- [Pinocchio](https://github.com/stack-of-tasks/pinocchio) with URDF support (`urdfdom`) and Boost (`boost_system`). Required by `make`, since every example is linked against the Pinocchio backend. The Optimized backend itself only needs Eigen and OpenMP (see below).
- Python 3 with `pandas`, `numpy` and `matplotlib` (only for the plotting scripts)
- CUDA: not required yet

Developed and tested with Pinocchio 3.9.0 on Linux Mint 21.3.

## Building

From the repository root:

```bash
make
```

This builds all examples into `bin/`. Use `make clean` to remove the build products.

### Building only the Optimized backend (without Pinocchio)

```bash
g++ -std=c++11 -O3 -fopenmp -I include -c src/iLQRSolverBase.cpp      -o obj/iLQRSolverBase.o
g++ -std=c++11 -O3 -fopenmp -I include -c src/iLQRSolverOptimized.cpp -o obj/iLQRSolverOptimized.o
g++ -std=c++11 -O3 -fopenmp -I include test/opt_iiwa_14.cpp \
    obj/iLQRSolverBase.o obj/iLQRSolverOptimized.o -o bin/opt_iiwa_14
```

If the compiler cannot find Eigen, add `-I /usr/include/eigen3` to each command.

## Running the examples

Run the binaries **from the repository root**: the examples load `utils/iiwa_14.urdf` and write their CSV output using paths relative to the current directory.

| Binary | Backend | Description |
|---|---|---|
| `acrobot` | Pinocchio | Swing-up of an underactuated 2-link Acrobot (only the first joint is actuated, through the actuation matrix `S`). The model is built in code, with no URDF. Writes `trajetoria.csv`. |
| `iiwa_14` | Pinocchio | Point-to-point reaching on the iiwa 14, from the zero configuration at rest to a fixed goal configuration (5 s horizon, `h = 0.01 s`, gravity disabled). Prints the final state, the final error and the solve time, and writes `trajetoria_iiwa.csv`. |
| `opt_iiwa_14` | Optimized | Same problem as `iiwa_14`, solved with the Optimized backend. |
| `iiwa_14_mpc`, `opt_iiwa_14_mpc` | Pinocchio, Optimized | Preliminary MPC-like simulation (see above; not a sufficient test) of the same reaching task: 64 knot points, `h = 0.01 s`, 500 control steps (5 s), one iLQR iteration per step. Prints the final error and the average time per step, and writes `trajetoria_iiwa.csv`. |
| `iiwa_14_printless`, `opt_iiwa_14_printless` | Pinocchio, Optimized | The reaching problem without console output or CSV export (e.g., for profiling with `perf`). |

Plotting:

```bash
./bin/opt_iiwa_14
python3 scripts/plot_iiwa.py        # reads trajetoria_iiwa.csv

./bin/acrobot
python3 scripts/plot_acrobot.py     # reads trajetoria.csv
```

The iiwa CSV has one row per time step with the columns `x1..x14` (state), `u1..u7` (torques) and `K1..K98` (the 7x14 feedback gain matrix, flattened row by row). The last row holds the final state with zero control and gains.

## Usage

```cpp
#include "iLQRSolverOptimized.hpp"   // or "iLQRSolverPinocchio.hpp"

iLQRSolverOptimized ilqr(h, Nt);                      // Optimized backend (iiwa 14)
// iLQRSolverPinocchio ilqr("utils/iiwa_14.urdf", h, Nt);   // Pinocchio backend

ilqr.setCosts(Q, R, Qn);                  // state, control and terminal weights
ilqr.setBoundaryConditions(x0, xgoal);    // initial and goal states
ilqr.setGravity(0.0);                     // optional; the iiwa examples disable gravity

ilqr.solve(max_iters, tol);

auto xtraj = ilqr.getXTraj();             // Nt states
auto utraj = ilqr.getUTraj();             // Nt-1 controls
auto Ktraj = ilqr.getKTraj();             // Nt-1 feedback gains (Nu x Nx)
```

MPC-like loop:

```cpp
ilqr.setBoundaryConditions(x0, xgoal);
for (int i = 0; i < steps; ++i)
{
    if (i > 0) ilqr.advanceMPCStep(x_current);   // new initial state; shifts the control sequence (warm start)
    ilqr.solve(1);                               // a single iLQR iteration per control step
    u = ilqr.getUTraj()[0];                      // control to apply
    x_current = ilqr.getXTraj()[1];              // here: predicted next state (ideal model)
}
```

| Method | Available in | Description |
|---|---|---|
| `setCosts(Q, R, Qn)` | all | Cost weights (default: zero matrices) |
| `setBoundaryConditions(x0, xgoal)` | all | Initial and goal states |
| `solve(max_iters = 1000, tol = 1e-3)` | all | Runs iLQR |
| `advanceMPCStep(x_current)` | all | Sets a new initial state and shifts the control sequence by one step |
| `getXTraj()`, `getUTraj()`, `getKTraj()` | all | Trajectories and feedback gains |
| `getNx()`, `getNu()` | all | State and control dimensions |
| `setGravity(g)` | Pinocchio, Optimized | Sets the gravity acceleration (the iiwa examples disable it with `setGravity(0.0)`) |
| `setViscousDamping(c)` | Pinocchio | Per-joint viscous damping coefficients |
| `setSelectionMatrix(S)` | Pinocchio | Actuation matrix: joint torques are `S * u`, so `Nu = S.cols()`. If `Nu` changes, call `setCosts` again. |

Notes:

- `solve()` returns `false` only when the backward pass fails because the control-cost Hessian could not be regularized. Otherwise it returns `true`, **including when `max_iters` was reached without convergence**.

## Roadmap

- [x] Pinocchio-based baseline solver
- [x] CPU-optimized backend using the dynamics-gradient code from [1], with multi-threaded derivative evaluation
- [x] Preliminary MPC-like simulation test (point-to-point reaching only; not a sufficient test)
- [ ] Detect the number of physical CPU cores at runtime instead of hard-coding 4 threads
- [ ] Replace the current MPC-like test with a path-tracking test in the style of [2]: the manipulator follows a moving reference trajectory (a figure-eight end-effector path in [2]) while the solver re-plans at every control step under a fixed time budget
- [ ] GPU (CUDA) backend, possibly reusing code from the repositories of [1] and [2]
- [ ] Real-time control of a real KUKA LBR iiwa 14

## Known limitations

- **The Optimized backend is specific to the iiwa 14.** The number of joints (7), the inertial parameters, the kinematics and the viscous damping (0.5 per joint) are hard-coded in `utils/dynamicsGradient.h`, and the actuation matrix is the identity. The inertial parameters and the damping are the same as in `utils/iiwa_14.urdf`.
- **The Pinocchio backend assumes `nq = nv`**, so it does not support floating-base or quaternion joints.
- **There are no constraints.** Torque and joint limits are not enforced, and some aggressive cost settings (commented out in the examples) converge to solutions that exceed realistic actuator limits.
- **Line search.** If no step size gives a sufficient cost decrease, the smallest step tried is still accepted.
- **Initial guess.** There is no API for an initial control guess: the first `solve()` starts from zero controls. Later calls are warm-started from the previous solution.
- **The MPC-like test is preliminary and insufficient**: a single fixed goal, a plant identical to the model, and no time budget per control step (see [Parallelism and MPC](#parallelism-and-mpc)).
- **Fixed goal only.** The cost is quadratic around a single goal state `xgoal`; time-varying references, needed for path tracking, are not supported yet.
- **Thread count** is hard-coded to 4 (`N_THREADS`).
- Some code comments and console messages are in Portuguese.
- In general, the code still lacks a proper amount of comments

## Third-party code

`utils/dynamicsGradient.h`, `utils/math.h` and `utils/minv.h` come from [plancherb1/fast-rbd-gradients](https://github.com/plancherb1/fast-rbd-gradients), the code accompanying [1], released under the MIT license (Copyright (c) 2021 Brian Plancher). `utils/iiwa_14.urdf` is the iiwa 14 model used in that same repository (inline comments in the file refer to the original Drake URDF). The full notice is in [THIRD_PARTY_NOTICES.md](THIRD_PARTY_NOTICES.md).

Pinocchio and Eigen are external dependencies and are not distributed with this repository.

## References

1. B. Plancher, S. M. Neuman, T. Bourgeat, S. Kuindersma, S. Devadas and V. J. Reddi. *Accelerating Robot Dynamics Gradients on a CPU, GPU, and FPGA.* IEEE Robotics and Automation Letters, 2021. Code: <https://github.com/plancherb1/fast-rbd-gradients>
2. B. Plancher and S. Kuindersma. *A Performance Analysis of Parallel Differential Dynamic Programming on a GPU.* International Workshop on the Algorithmic Foundations of Robotics (WAFR). Code: <http://bit.ly/ParallelDDP>
3. J. Carpentier and N. Mansard. *Analytical Derivatives of Rigid Body Dynamics Algorithms.* Robotics: Science and Systems (RSS), 2018.
4. W. Li and E. Todorov. *Iterative Linear Quadratic Regulator Design for Nonlinear Biological Movement Systems.* International Conference on Informatics in Control, Automation and Robotics (ICINCO), 2004.
5. Y. Tassa, T. Erez and E. Todorov. *Synthesis and Stabilization of Complex Behaviors through Online Trajectory Optimization.* IEEE/RSJ International Conference on Intelligent Robots and Systems (IROS), 2012.

## License

This project is licensed under the MIT License. See [LICENSE](LICENSE) for details. Third-party code is covered by its own license (see above).
