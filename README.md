# Single-Point-Positioning: determining the position of the receiver based on the given observation and navigation RINEX files
# Satellite Position Based on Almanac Ephemeris

This repository contains the Python code for calculating the position of GNSS satellites based on almanac ephemeris data. The satellite coordinates are presented in the Earth-fixed coordinate system. The task is part of the Satellite Navigation Systems course at the Faculty of Geodesy and Cartography, Warsaw University of Technology.

## Task Overview

The aim of the task is to calculate the position of a GNSS satellite based on the data from the almanac. The satellite coordinates should be presented in the Earth-fixed coordinate system.

## Repository Structure

This repository consists of a single Python script:

1. **main.py**: This script contains all the steps from loading the almanac data, to processing the data and calculating the satellite's position.

## Instructions

1. Clone this repository to your local machine.
2. Ensure that you have Python installed. If not, you can install it from [here](https://www.python.org/downloads/).
3. Navigate to the cloned repository and open `main.py` to view the code.
4. To run the code, ensure that you have all the necessary Python packages installed. You can install the necessary packages using pip (e.g., `pip install numpy matplotlib`).
5. Before running the script, please ensure that you have the correct path to the almanac data. Modify the file paths in the script as needed.
6. Run the script using a Python interpreter (e.g., `python main.py`).

## Dependencies

The code in this repository requires the following Python packages:

- Numpy
- Matplotlib

## Contact

If you encounter any issues or have questions, please open an issue in this GitHub repository.
