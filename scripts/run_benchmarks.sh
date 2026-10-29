#!/bin/bash

set -e

sudo cpupower frequency-set -g performance > /dev/null
export OMP_WAIT_POLICY=PASSIVE

if [ "$1" == "1" ]; then
	./bin/benchmark_1
	./bin/opt_benchmark_1
elif [ "$1" == "2" ]; then
	./bin/benchmark_2
	./bin/opt_benchmark_2
elif [ "$1" == "no_turbo" ]; then
	echo 1 | sudo tee /sys/devices/system/cpu/intel_pstate/no_turbo > /dev/null
	./bin/benchmark_2
	./bin/opt_benchmark_2
	echo 0 | sudo tee /sys/devices/system/cpu/intel_pstate/no_turbo > /dev/null
else
	./bin/benchmark_1
	./bin/opt_benchmark_1
	./bin/benchmark_2
	./bin/opt_benchmark_2
fi
