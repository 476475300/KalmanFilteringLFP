# Using Kalman filtering techniques to obtain a better State of Charge-estimation of the battery pack of a Solar Car
Demo repository for the thesis "Using Kalman filtering techniques to obtain a better State of Charge-estimation of the battery pack of a Solar Car".

Author: Frederik Vanmaele

Promotor: Prof. Dr. Ir. Bart De Moor

Supervisor: Ir. Mauricio Agudelo


# Data
Input-output sequences that are used in the Kalman Filtering or Modelling demos.

# Parameters
Averaged parameters for the best set of each model, as specified in the thesis. 
OCV_SOC is the relationship based on the OCV-tests conducted on the three sample cells.

# Modelling
The folder "Models" contains the different Simulink models and can be inspected separately.
The output of the different models on the race data can be viewed by running the respective files. The input-output sequence of each sample cell can be used by selecting the right input in the file.

# Kalman Filtering
The files that do not start with UKF or EKF are helper functions.
## EKF
EKF_main contains the framework of the Extended Kalman Filtering algorithm.
EKF_script runs the EKF_main file for the BWSC race data obtained in the lab on the sample cells.
EKF_Marokko runs the EKF_main file with the real data of the Solar Challenge Morocco. 
## UKF
Completely analogue.


