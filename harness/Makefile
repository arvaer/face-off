# Makefile for virial artifact harness
CXX     = c++
CXXFLAGS = -O2 -std=c++17 -Wall

.PHONY: all clean run

all: harness

harness: virial_artifact_harness.cpp
	$(CXX) $(CXXFLAGS) $< -o $@

# Quick smoke test: nstlist=20, rlist=1.2, 1000 steps, 200 particles
run: harness
	./harness 20 1.2 1000 200 | head -20

# Full artifact scan (requires matplotlib)
analyze: harness
	python3 analyze_artifact.py --harness ./harness --nsteps 5000 --N 500

clean:
	rm -f harness *.png *.csv
