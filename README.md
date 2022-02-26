# Binary Star Simulator

## About

This repository contains a tool for simulating and visualizing the orbital 
mechanics of p-type binary star systems. Additionally, orbiting planets and 
satellites can be simulated as well as circumbinary habitable zones. All 
parameters are freely configurable in conformity with the laws of physics.

## Getting Started

Use these instructions to get a copy of the project running in a virtual 
environment on your local machine.

### Prerequisites

The simulator is written in Python. It is highly recommended using a virtual 
environment (e.g. [virtualenv](https://docs.python.org/3/tutorial/venv.html)) 
to install additional packages while keeping them separated from your other 
Python projects.

### Installation

To get started, clone this repository, set up a virtual environment (at 
least Python 3.8) and install the dependencies by running `pip install -r 
requirements.txt` from the root directory of this project.

## Usage

The root directory contains several scripts that simulate different aspects 
of p-type binary star systems. At the beginning of each script you will find 
a small section where various simulation parameters can be set. Here you can 
also choose from a set of pre-parameterized star systems by specifying the 
corresponding name. The available systems are listed in `star_catalogue.yml`. 
Define your own parameterization by adding new star systems to the catalogue. 
Note that certain combinations of parameters lead to physical unstable 
solutions and therefore cannot be simulated.

## License

This project is published under the [MIT license](LICENSE.txt).
