
#Orbit_Estimation_tool
This workpackage consists of two primary tools which includes Orbit Estimation tool and State Vector to TLE conversion tool.
Orbit Estimation tool, uses GPS data in terms of position and velocity vectors from satellites, and then Batch Least Square Estimator algorithm is engineered to generate one single value for all the parameters in the dataset which can be further used to generate Keplerian parameters and find the exact position of the satellite. The tool is validated using satellite simulation tool with an accuracy of 99.8% 
The libraries used in this program are astropy libraries from which we can retrieve the required constant values to be used in the calculation of Keplerian parameters. In order for proper execution of the file ensure that the file path and file name of the dataset is changed according to the users' computer. 


#State_to_TLE
The State_to_TLE tool is an important tool for satellite operators to determine the initial position of the satellite immediately after launch and point the antennas in the right direction to establish a down link. This tool uses State vectors which are poition and velocity vectors along with Unix date and time to convert the vectors into TLEs(Two Line Elements) where the parameters are represented in a specific order representing the Keplerian parameters required to point the antenna in the right direction. The tool is currently in operation for all upcoming missions, with an accuracy of 99.9%


Execution of the file
Both these tools have very simplified execution process involved i.e executing the files on a python enviornment, provided the user changes the name of the file and file path based on their dataset. 
