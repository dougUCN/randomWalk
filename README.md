# random walk

This simulation is a quick calculation used to determine the impact on stored UCN counts caused by a thin aluminum window in the region of the polarizing magnet for the Los Alamos neutron electric dipole moment experiment.

A description of the methodology used is found in the "Monte Carlo simulations of UCN transport and spin transport" chapter in this [thesis](https://github.com/dougUCN/thesis).

## Getting started

```
cmake .
make
```

Dependencies: [Boost C++ libraries](https://www.boost.org/), HDF5

Once the executable is built, run

```
./randomWalk_t.x --help
```

to see available parameter options

## Utility scripts

**parallel.py**: [Uses GNU parallel](https://www.gnu.org/software/parallel/) to run multiple random walks together

**parse.py**: Parse and plot the results of a `*.h5` file
