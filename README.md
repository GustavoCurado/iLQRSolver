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
| **Pinocchio** (baseline) | `iLQRSolverPinocchio` | Pinocchio's Articulated Body Algorithm and its analytical derivatives [3] | Any model loaded from a URDF (or built in code) with `nq = nv`; optional actuation matrix, viscous damping and gravity | OpenMP across knot points (off by default) | Working |
| **Optimized** (CPU) | `iLQRSolverOptimized` | Rigid-body dynamics and analytical gradient code from Plancher et al. [1] | KUKA LBR iiwa 14 only (7 DoF, hard-coded) | OpenMP across knot points (on by default) | Working |
| **GPU** (CUDA) | `iLQRSolverGPU` | n/a | n/a | n/a | Skeleton only, not built yet |

## How the solver works

The solver implements iLQR [4, 5]: it finds a control sequence that drives the system from an initial state `x0` towards a goal state `xgoal` over a horizon of `Nt` knot points spaced `h` seconds apart. The state is `x = [q; v]` (joint positions and velocities) and the control `u` is a vector of joint torques. The cost is quadratic around the goal: a running cost penalizing state error (weight `Q`) and control effort (weight `R`), plus a terminal cost (weight `Qn`).

Each iteration of the algorithm has three steps:

1. **Linearization.** The continuous dynamics are discretized with a 4th-order Runge-Kutta (RK4) integrator and linearized along the current trajectory. The Jacobians are propagated through the RK4 stages, so they are the exact derivatives of the discretized dynamics. As in any iLQR, second derivatives of the dynamics are ignored (unlike full DDP).
2. **Backward pass.** A Riccati-like recursion, from the last knot point to the first, computes a feedforward correction `d_k` and a feedback gain `K_k` for every knot point. If the control-cost Hessian is not positive definite, a regularization term is added to its diagonal.
3. **Forward rollout with line search.** A new trajectory is simulated with the nonlinear dynamics using `u_k - alpha * d_k - K_k * (x_new_k - x_k)`. The step size `alpha` starts at 1 and is halved until the cost decreases sufficiently.

The iteration stops when the largest feedforward correction, `max |d_k|`, falls below `tol` (default `1e-3`), when `max_iters` (default 1000) is reached or, if a time limit `time_limit_ms` is given, when the elapsed time exceeds it (the limit is checked before each iteration after the first, so at least one iteration always runs).

### Code organization

`iLQRSolverBase` implements everything above. A backend only has to implement two virtual functions:

```cpp
// Continuous-time dynamics: x' = f(x, u)
Eigen::VectorXd continuousDynamics(const Eigen::VectorXd& x, const Eigen::VectorXd& u);

// Evaluates x' = f(x, u) and its Jacobians Jx = df/dx and Ju = df/du
void getContinuousJacs(const Eigen::VectorXd& x, const Eigen::VectorXd& u,
                       Eigen::VectorXd& xdot, Eigen::MatrixXd& Jx, Eigen::MatrixXd& Ju);
```

and set `Nx` and `Nu` before calling `iLQRSolverBase::initializeMemory()`. Adding a new backend (such as the GPU one) means filling in these two functions. Since `getContinuousJacs` also returns `x'`, the RK4 linearization does not need to evaluate the dynamics a second time. When parallel derivatives are enabled, both functions are called concurrently from several threads, so a backend must not share mutable workspaces between calls (both current backends keep one workspace per thread, indexed by `omp_get_thread_num()`).

## Parallelism and MPC

**Parallel derivatives.** The linearization of the dynamics at each knot point is independent of the others. When parallel derivatives are enabled, the backward pass computes all of them in an OpenMP parallel loop. This is the coarse-grained parallelism (many independent gradient evaluations) described in [1]. Both backends support it: the Optimized backend keeps one copy of the transformation and inertia arrays per thread, and the Pinocchio backend keeps one `pinocchio::Data` per thread. It is enabled by default in the Optimized backend and disabled by default in the Pinocchio backend; `setParallelDerivatives(true/false)` changes it. The number of threads is currently the compile-time constant `N_THREADS = 4` in `include/iLQRSolverBase.hpp`.

**MPC-like simulation.** The `*_mpc` examples run a receding-horizon loop in simulation: at each control step the solver is warm-started with the previous solution (`advanceMPCStep`), re-optimizes, and the first control is applied. The examples report the average solve time per control step.

The MPC-like tests are not adequate yet and should not be used as a validation of the solver. They will be revised in the future. The planned test follows Plancher and Kuindersma [2] (Section 5.3): the manipulator has to track a moving reference trajectory (a figure-eight path for the end effector in [2]), with the solver re-planning at every control step within a fixed time budget (10 ms in [2]) in a simulation that runs in real time. See the [Roadmap](#roadmap).

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
├── test/                            # Example and benchmark programs (see below)
├── utils/
│   ├── dynamicsGradient.h           # Third-party code from [1] (MIT license);
│   ├── math.h                       #   see "Third-party code" below
│   ├── minv.h                       #
│   └── iiwa_14.urdf                 # KUKA LBR iiwa 14 model, also from [1] (used by the Pinocchio backend)
├── scripts/
│   ├── plot_acrobot.py              # Plots the CSV written by the acrobot example
│   ├── plot_iiwa.py                 # Plots the CSV written by the iiwa examples
│   └── run_benchmarks.sh            # Runs the benchmarks with the CPU configured for stable timings
├── Makefile
├── LICENSE
└── THIRD_PARTY_NOTICES.md
```

## Dependencies

- Linux with `g++` (C++11, OpenMP) and `pkg-config`
- [Eigen 3](https://eigen.tuxfamily.org) (located through `pkg-config`)
- [Pinocchio](https://github.com/stack-of-tasks/pinocchio) with URDF support (`urdfdom`) and Boost (`boost_system`). Required by the Pinocchio executables, and therefore by a plain `make`. The Optimized executables only need Eigen and OpenMP (see below).
- An x86-64 CPU with AVX, since the Optimized executables are compiled with `-march=native -mavx`
- Python 3 with `pandas`, `numpy` and `matplotlib` (only for the plotting scripts)
- `sudo` and `cpupower` (only for `scripts/run_benchmarks.sh`)
- CUDA: not required yet

Developed and tested with Pinocchio 3.9.0 on Linux Mint 21.3.

## Building

From the repository root:

```bash
make
```

This builds all examples and benchmarks into `bin/`. Use `make clean` to remove the build products.

### Building only the Optimized backend (without Pinocchio)

```bash
make bin/opt_iiwa_14 bin/opt_iiwa_14_mpc bin/opt_iiwa_14_printless bin/opt_benchmark_1 bin/opt_benchmark_2
```

These `opt_*` targets only need Eigen and OpenMP, and any subset of them can be built.

## Running the examples

Run the binaries **from the repository root**: the examples load `utils/iiwa_14.urdf` and write their CSV output using paths relative to the current directory.

| Binary | Backend | Description |
|---|---|---|
| `acrobot` | Pinocchio | Swing-up of an underactuated 2-link Acrobot (only the first joint is actuated, through the actuation matrix `S`). The model is built in code, with no URDF. Writes `trajetoria.csv`. |
| `iiwa_14` | Pinocchio | Point-to-point reaching on the iiwa 14, from the zero configuration at rest to a fixed goal configuration (5 s horizon, `h = 0.01 s`, gravity disabled). Prints the final state, the final error and the solve time, and writes `trajetoria_iiwa.csv`. |
| `opt_iiwa_14` | Optimized | Same problem as `iiwa_14`, solved with the Optimized backend. |
| `iiwa_14_mpc` | Pinocchio | MPC-like simulation (not adequate yet, see above) of the reaching task: 64 knot points, `h = 0.01 s`, 500 control steps (5 s), one iLQR iteration per step. Prints the final error and the average time per step, and writes `trajetoria_iiwa.csv`. |
| `opt_iiwa_14_mpc` | Optimized | MPC-like simulation (not adequate yet, see above) in which the joint-space goal moves sinusoidally in time and is updated at every step with `setGoal`: 64 knot points, `h = 0.01 s`, 1000 control steps (10 s), up to 100 iterations per step with a 2 ms time limit, and the plant simulated with `simulatePlant`. Prints the final error, the average number of iterations and the average time per step, and writes `trajetoria_iiwa.csv`. |
| `iiwa_14_printless`, `opt_iiwa_14_printless` | Pinocchio, Optimized | The reaching problem without console output or CSV export (e.g., for profiling with `perf`). |

Plotting:

```bash
./bin/opt_iiwa_14
python3 scripts/plot_iiwa.py        # reads trajetoria_iiwa.csv

./bin/acrobot
python3 scripts/plot_acrobot.py     # reads trajetoria.csv
```

The iiwa CSV has one row per time step with the columns `x1..x14` (state), `u1..u7` (torques) and `K1..K98` (the 7x14 feedback gain matrix, flattened row by row). The last row holds the final state with zero control and gains.

## Benchmarks

There are two benchmark programs, each built for both backends (`benchmark_1` and `benchmark_2` for Pinocchio, `opt_benchmark_1` and `opt_benchmark_2` for Optimized):

- **`benchmark_1`** measures the dynamics derivatives. It reports the latency of a single call to `continuousDynamics`, `getContinuousJacs` and `computeDerivatives` (average of 100000 calls), and the time to run `computeDerivatives` on batches of N = 16, 32, 64, 128 and 500 random points (the first four are the values of N used in [1]), with one thread and with `N_THREADS` OpenMP threads. For each batch it prints the median, mean, standard deviation, minimum and 99th percentile, and the multi-thread speedup.
- **`benchmark_2`** measures the full solver. It solves the iiwa 14 reaching problem (`Nt = 501`, `h = 0.01 s`, until convergence or 10000 iterations) 30 times, single-threaded and multi-threaded, and reports the time spent in the derivatives, in the backward (Riccati) recursion and in the forward rollouts, as a share of the total `solve()` time.

`scripts/run_benchmarks.sh`, run from the repository root after `make`, executes them:

```bash
./scripts/run_benchmarks.sh            # all four executables
./scripts/run_benchmarks.sh 1          # benchmark_1 and opt_benchmark_1
./scripts/run_benchmarks.sh 2          # benchmark_2 and opt_benchmark_2
./scripts/run_benchmarks.sh no_turbo   # benchmark_2 and opt_benchmark_2 with Intel Turbo Boost disabled
```

The script uses `sudo` to set the CPU frequency governor to `performance` (through `cpupower`) and sets `OMP_WAIT_POLICY=PASSIVE`. The `no_turbo` mode also disables Intel Turbo Boost through `intel_pstate` while the benchmarks run and re-enables it afterwards. The governor is not restored to its previous value.

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
    ilqr.setGoal(xgoal_i);                          // optional: new goal for this step
    if (i > 0) ilqr.advanceMPCStep(x_current);      // new initial state; shifts the control sequence (warm start)
    ilqr.solve(max_iters, tol, time_limit_ms);      // optionally limited in time
    u = ilqr.getUTraj()[0];                         // control to apply
    x_current = ilqr.simulatePlant(x_current, h);   // simulates one control period
}
```

| Method | Available in | Description |
|---|---|---|
| `setCosts(Q, R, Qn)` | all | Cost weights (default: zero matrices) |
| `setBoundaryConditions(x0, xgoal)` | all | Initial and goal states |
| `solve(max_iters = 1000, tol = 1e-3, time_limit_ms = 0)` | all | Runs iLQR. If `time_limit_ms > 0`, it also stops when the elapsed time exceeds the limit. |
| `advanceMPCStep(x_current)` | all | Sets a new initial state and shifts the control sequence by one step |
| `setGoal(xgoal)` | all | Changes only the goal state (e.g., between MPC steps) |
| `resetTrajectories()` | all | Resets the trajectories (zero controls), so that the next `solve()` starts from scratch |
| `simulatePlant(x, dt_total, substeps = 100)` | all | Integrates the solver's own dynamics model (RK4, `substeps` sub-steps) for `dt_total` seconds from `x`, applying the current controls with feedback, `u_k - K_k (x - x_k)`, and returns the final state |
| `computeDerivatives(x, u, A, B)` | all | Discrete-time (RK4) Jacobians `A = df/dx` and `B = df/du` at one point |
| `setParallelDerivatives(enable)` | all | Enables or disables the OpenMP evaluation of the derivatives (default: on in Optimized, off in Pinocchio) |
| `getTimeDerivativesUs()`, `getTimeRiccatiUs()`, `getTimeRolloutUs()`, `getTimeTotalUs()`, `getLastIters()` | all | Time (in microseconds) spent in the derivatives, in the backward (Riccati) recursion, in the rollouts and in total during the last `solve()`, and the number of iterations of that call |
| `getXTraj()`, `getUTraj()`, `getKTraj()` | all | Trajectories and feedback gains |
| `getNx()`, `getNu()` | all | State and control dimensions |
| `setGravity(g)` | Pinocchio, Optimized | Sets the gravity acceleration (the iiwa examples disable it with `setGravity(0.0)`) |
| `setViscousDamping(c)` | Pinocchio | Per-joint viscous damping coefficients |
| `setSelectionMatrix(S)` | Pinocchio | Actuation matrix: joint torques are `S * u`, so `Nu = S.cols()`. If `Nu` changes, call `setCosts` again. |

Notes:

- `solve()` returns `false` only when the backward pass fails because the control-cost Hessian could not be regularized. Otherwise it returns `true`, **including when `max_iters` was reached or the time limit expired without convergence**.

## Roadmap

- [x] Pinocchio-based baseline solver
- [x] CPU-optimized backend using the dynamics-gradient code from [1], with multi-threaded derivative evaluation
- [x] Multi-threaded derivative evaluation in the Pinocchio backend
- [x] Benchmarks of the dynamics derivatives and of the full solver
- [x] Preliminary MPC-like simulation tests (not adequate yet)
- [ ] Detect the number of physical CPU cores at runtime instead of hard-coding 4 threads
- [ ] Replace the current MPC-like tests with a path-tracking test in the style of [2]: the manipulator follows a moving reference trajectory (a figure-eight end-effector path in [2]) while the solver re-plans at every control step under a fixed time budget
- [ ] GPU (CUDA) backend, possibly reusing code from the repositories of [1] and [2]
- [ ] Real-time control of a real KUKA LBR iiwa 14

## Known limitations

- **The Optimized backend is specific to the iiwa 14.** The number of joints (7), the inertial parameters, the kinematics and the viscous damping (0.5 per joint) are hard-coded in `utils/dynamicsGradient.h`, and the actuation matrix is the identity. The inertial parameters and the damping are the same as in `utils/iiwa_14.urdf`.
- **The Pinocchio backend assumes `nq = nv`**, so it does not support floating-base or quaternion joints.
- **There are no constraints.** Torque and joint limits are not enforced, and some aggressive cost settings (commented out in the examples) converge to solutions that exceed realistic actuator limits.
- **Line search.** If no step size gives a sufficient cost decrease, the smallest step tried is still accepted.
- **Initial guess.** There is no API for an initial control guess: the first `solve()` starts from zero controls. Later calls are warm-started from the previous solution.
- **The MPC-like tests are not adequate yet** and should not be used as a validation of the solver (see [Parallelism and MPC](#parallelism-and-mpc)).
- **Constant goal along the horizon.** The cost is quadratic around a single goal state `xgoal`. It can be changed between solves with `setGoal`, but a reference that varies along the horizon is not supported yet.
- **Thread count** is hard-coded to 4 (`N_THREADS`, in `include/iLQRSolverBase.hpp`).
- Some code comments and console messages are in Portuguese.
- In general, the code still lacks a proper amount of comments

## Third-party code

`utils/dynamicsGradient.h`, `utils/math.h` and `utils/minv.h` come from [plancherb1/fast-rbd-gradients](https://github.com/plancherb1/fast-rbd-gradients), the code accompanying [1], released under the MIT license (Copyright (c) 2021 Brian Plancher). `utils/iiwa_14.urdf` is the iiwa 14 model used in that same repository (inline comments in the file refer to the original Drake URDF). The `time_delta_us` timing macro in `include/iLQRSolverBase.hpp` is also taken from that repository. The full notice is in [THIRD_PARTY_NOTICES.md](THIRD_PARTY_NOTICES.md).

Pinocchio and Eigen are external dependencies and are not distributed with this repository.

## References

1. B. Plancher, S. M. Neuman, T. Bourgeat, S. Kuindersma, S. Devadas and V. J. Reddi. *Accelerating Robot Dynamics Gradients on a CPU, GPU, and FPGA.* IEEE Robotics and Automation Letters, 2021. Code: <https://github.com/plancherb1/fast-rbd-gradients>
2. B. Plancher and S. Kuindersma. *A Performance Analysis of Parallel Differential Dynamic Programming on a GPU.* International Workshop on the Algorithmic Foundations of Robotics (WAFR). Code: <http://bit.ly/ParallelDDP>
3. J. Carpentier and N. Mansard. *Analytical Derivatives of Rigid Body Dynamics Algorithms.* Robotics: Science and Systems (RSS), 2018.
4. W. Li and E. Todorov. *Iterative Linear Quadratic Regulator Design for Nonlinear Biological Movement Systems.* International Conference on Informatics in Control, Automation and Robotics (ICINCO), 2004.
5. Y. Tassa, T. Erez and E. Todorov. *Synthesis and Stabilization of Complex Behaviors through Online Trajectory Optimization.* IEEE/RSJ International Conference on Intelligent Robots and Systems (IROS), 2012.

## License

This project is licensed under the MIT License. See [LICENSE](LICENSE) for details. Third-party code is covered by its own license (see above).
