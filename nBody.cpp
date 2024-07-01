#include <iostream>
#include <vector>
#include <algorithm>
#include <random> 
#include <thread>
#include <condition_variable>
#include <mutex>
#include <chrono>
#include <iostream>
#include <fstream>


//#define DEBUG_ACC
#define FILE

using matrix = std::vector<std::vector<double>>;

// N-body simulation

// Simulation parameters

struct Config {
	size_t numberOfThreads;
	size_t numberOfBodies;
	std::string ouputFileName;
	double startTime;
	double endTime;
	double dt;
	static double softening;
	static double G;	

	double numberOfTimeSteps(){ return (endTime - startTime) / dt; }
	double bodiesPerThread(){ return numberOfBodies/numberOfThreads; }	
	double totalMass(){ return 10*numberOfBodies; }	
};

double Config::softening = 11.6;
double Config::G = 9.8;
Config config {1, 3, "nBody.txt", 0, 20, 0.1};

constexpr size_t numberOfDimensions{ 3 };

std::thread threads[65];
std::vector<double> masses[65];
matrix initialPositions[65];
matrix initialVelocities[65];

// shared between threads
double culumtativeBodiesMasses[65];
std::vector<double> culumtativeBodiesPositions[65];

// synchronization
size_t centerOfMassCalculatedIteration[65];
std::mutex centerOfMassCalculatedIterationMutex[65];
std::condition_variable centerOfMassCalculatedIterationCondVariable[65];
std::condition_variable test;
std::mutex printMutex;
std::mutex fileWriteMutex;
std::ofstream outputFile;


void initSharedVariables()
{
	for (size_t i = 0; i < config.numberOfThreads; i++)
	{
		culumtativeBodiesMasses[i] = 0.0;
		std::vector<double> v{ 0.0, 0.0, 0.0 };
		culumtativeBodiesPositions[i] = v;
	}
}

matrix getAcceleration(const matrix& planetPositions, const std::vector<double>& masses, const double gravityConstant, const double softening, const size_t currentThreadIdx, const int currentIteration)
{
	matrix accelerations;
	accelerations.reserve(config.bodiesPerThread());

	for (size_t i = 0; i <config.bodiesPerThread(); i++)
	{
		accelerations.push_back({ 0.0, 0.0, 0.0 } );
	}

	double dx, dy, dz, inv_r3;

	for (size_t i = 0; i < config.bodiesPerThread(); i++)
	{
		// Effect of bodies in the same thread region
		for (size_t j = 0; j < config.bodiesPerThread(); j++)
		{
			if (i == j)
			{
				continue;
			}

			dx = planetPositions[j][0] - planetPositions[i][0];
			dy = planetPositions[j][1] - planetPositions[i][1];
			dz = planetPositions[j][2] - planetPositions[i][2];

			inv_r3 = std::pow((dx * dx + dy * dy + dz * dz + softening * softening), -1.5);

			accelerations[i][0] += gravityConstant * (dx * inv_r3) * masses[j];
			accelerations[i][1] += gravityConstant * (dy * inv_r3) * masses[j];
			accelerations[i][2] += gravityConstant * (dz * inv_r3) * masses[j];

#ifdef DEBUG_ACC
			std::cout << "dx: " << dx << std::endl;
			std::cout << "dy: " << dy << std::endl;
			std::cout << "dz: " << dz << std::endl;

			std::cout << "inv_r3: " << inv_r3 << std::endl;
			std::cout << "G: " << G << std::endl;
			std::cout << "masses[j]: " << masses[j] << std::endl;

			std::cout << "accelerations[i][0]: " << accelerations[i][0] << std::endl;
			std::cout << "accelerations[i][1]: " << accelerations[i][1] << std::endl;
			std::cout << "accelerations[i][2]: " << accelerations[i][2] << std::endl;
#endif // DEBUG_ACC
		}

		// Effect of bodies outside the thread region
		for (size_t threadIdx = 0; threadIdx < config.numberOfThreads; threadIdx++)
		{
			if (currentThreadIdx == threadIdx)
			{
				continue;
			}
			
			{
				std::unique_lock<std::mutex> lock(centerOfMassCalculatedIterationMutex[threadIdx]);


				centerOfMassCalculatedIterationCondVariable[threadIdx].wait(lock, [threadIdx, currentIteration] { return centerOfMassCalculatedIteration[threadIdx] >= currentIteration; });
			}

			dx = culumtativeBodiesPositions[threadIdx][0] - planetPositions[i][0];
			dy = culumtativeBodiesPositions[threadIdx][1] - planetPositions[i][1];
			dz = culumtativeBodiesPositions[threadIdx][2] - planetPositions[i][2];


			inv_r3 = std::pow((dx * dx + dy * dy + dz * dz + softening * softening), -1.5);

			accelerations[i][0] += gravityConstant * (dx * inv_r3) * culumtativeBodiesMasses[threadIdx];
			accelerations[i][1] += gravityConstant * (dy * inv_r3) * culumtativeBodiesMasses[threadIdx];
			accelerations[i][2] += gravityConstant * (dz * inv_r3) * culumtativeBodiesMasses[threadIdx];
		}

	}

	return accelerations;
}

struct Energy {
	double kinetic;
	double potential;
};

Energy getEnergy(matrix planetPositions, matrix velocities, std::vector<double> masses, const double gravityConstant)
{
	Energy result{0,0};

	// Kinetic Energy :
	for (size_t i = 0; i < config.numberOfBodies; i++)
	{
		for (size_t j = 0; j < numberOfDimensions; j++)
		{
			result.kinetic += 0.5 * masses[i] * velocities[i][j] * velocities[i][j];
		}
	}

	double dx, dy, dz, inv_r3;

	// Potential Energy :
	for (size_t i = 0; i < config.numberOfBodies; i++)
	{
		for (size_t j = i+1; j < config.numberOfBodies; j++)
		{
			dx = planetPositions[j][0] - planetPositions[i][0];
			dy = planetPositions[j][1] - planetPositions[i][1];
			dz = planetPositions[j][2] - planetPositions[i][2];


			inv_r3 = std::pow((dx * dx + dy * dy + dz * dz), 0.5);
			result.potential += (gravityConstant * masses[i] * masses[j]) / inv_r3;
		}
	}

	result.potential *= -1;
	return result;
}

void printMatrix(const matrix& m, const std::string& message)
{
	{
		std::lock_guard<std::mutex> lock(printMutex);
		std::cout << message << std::endl;

		std::cout << "[" << std::endl;
		for (size_t i = 0; i < m.size(); i++)
		{
			std::cout << "  [";
			for (size_t j = 0; j < m[i].size(); j++)
			{
				std::cout << m[i][j];
				if (j != m[i].size() - 1)
				{
					std::cout << ", ";
				}
			}
			std::cout << "]" << std::endl;
		}
		std::cout << "]" << std::endl << std::endl;
	}
}

void initialVelocityValues(matrix mat[65])
{
	for (size_t i = 0; i < config.numberOfThreads; i++)
	{
		mat[i] = matrix();
	}
	mat[0].push_back({ -0.1,-0.1,0.1 });
	mat[0].push_back({ -0.1,0.3,0.3 });
	mat[0].push_back({ -0.2,0.2,0.2 });
}

void initialPlanetPositionValues(matrix mat[65])
{
	for (size_t i = 0; i < config.numberOfThreads; i++)
	{
		mat[i] = matrix();
	}
	mat[0].push_back({0,0,0});
	mat[0].push_back({ 40,40,10 });
	mat[0].push_back({ -50,40,10 });
}

void nBodyWorkerThread(const size_t threadIdx);


int main(int argc, char *argv[])
{
	if(argc >= 2)
	{
		config.numberOfThreads = std::atoi(argv[1]);
	}
	if(argc >= 3)
	{
		config.numberOfBodies = std::atoi(argv[2]);
	}
	if(argc >= 4)
	{
		config.ouputFileName = argv[3];
	}
	if(argc >= 5)
	{
		config.startTime = std::atoi(argv[4]);
	}
	if(argc >= 6)
	{
		config.endTime = std::atoi(argv[5]);
	}
	if(argc >= 7)
	{
		config.dt = std::stod(argv[6]);
	}
	if(argc >= 8)
	{
		Config::softening = std::stod(argv[7]);
	}
	if(argc >= 9)
	{
		Config::G = std::stod(argv[8]);
	}

	outputFile.open(config.ouputFileName);
	std::uniform_real_distribution<double> unifs[65];
	std::uniform_real_distribution<double> unifsVel[65];
	outputFile << config.numberOfBodies << "," << config.numberOfBodies << "," << config.dt << "\n";

	
	int start = -10000, end = 10000;
	// int interval = (end - start) / (config.numberOfThreads * 2);
	// for (int i = 0; i < config.numberOfThreads; i++)
	// {
	// 	unifs[i] = std::uniform_real_distribution<double>(start + 2*i*interval, start + interval + 2*i*interval);
	// }

	unifs[0] = std::uniform_real_distribution<double>(start, end);
	unifsVel[0] = std::uniform_real_distribution<double>(-100, 100);
	std::default_random_engine re;


	for (size_t i = 0; i < config.numberOfThreads; i++)
	{
		masses[i] = std::vector<double>();
		initialPositions[i] = std::vector<std::vector<double>>();
		initialVelocities[i] = std::vector<std::vector<double>>();

		masses[i].reserve(config.bodiesPerThread());
		initialPositions[i].reserve(config.bodiesPerThread());
		initialVelocities[i].reserve(config.bodiesPerThread());
	}

	for (size_t i = 0; i < config.numberOfThreads; i++)
	{
		for (size_t j = 0; j < config.bodiesPerThread(); j++)
		{
			masses[i].push_back(config.totalMass() / config.numberOfBodies);

			std::vector<double> v1, v2;
			for (size_t k = 0; k < numberOfDimensions; k++)
			{
				v1.push_back(unifs[0](re));
				v2.push_back(unifsVel[0](re));
			}

			initialPositions[i].push_back(v1);
			initialVelocities[i].push_back(v2);
		}
	}

	//initialPlanetPositionValues(initialPositions);
	//initialVelocityValues(initialVelocities);

	initSharedVariables();
	auto t1 = std::chrono::high_resolution_clock::now();
	for (int i = 0; i < config.numberOfThreads; ++i)
	{
		threads[i] = std::thread(nBodyWorkerThread, i);
	}

	for (int i = 0; i < config.numberOfThreads; ++i)
	{
		threads[i].join();
	}
	auto t2 = std::chrono::high_resolution_clock::now();
	std::chrono::duration<double, std::milli > ms_double = t2 - t1;
	std::cout << ms_double.count() / 1000 << "ms\n";
	outputFile.close();
}

double calculateTotalLocalMass(const size_t threadIdx)
{
	double totalLocalMass{ 0.0 };

	for (size_t i = 0; i < config.bodiesPerThread(); i++)
	{
		totalLocalMass += masses[threadIdx][i];
	}

	return totalLocalMass;
}

std::vector<double> calculateCenterOfMassPosition(const size_t threadIdx, const double totalMass)
{
	std::vector<double> centerOfMassPosition{0,0,0};

	for (size_t i = 0; i < config.bodiesPerThread(); i++)
	{
		for (size_t j = 0; j < numberOfDimensions; j++)
		{
			centerOfMassPosition[j] += initialPositions[threadIdx][i][j] * masses[threadIdx][i];
		}
	}

	for (size_t i = 0; i < numberOfDimensions; i++)
	{
		centerOfMassPosition[i] /= totalMass;
	}

	return centerOfMassPosition;
}

void nBodyWorkerThread(const size_t threadIdx)
{
	//printMatrix(initialPositions[threadIdx], "Initial position");
	
	matrix accelerations;
	double totalLocalMass;
	std::vector<double> centerOfMassPosition;

	// step -1
	totalLocalMass = calculateTotalLocalMass(threadIdx);
	centerOfMassPosition = std::move(calculateCenterOfMassPosition(threadIdx, totalLocalMass));

	{
		std::lock_guard<std::mutex> lock(centerOfMassCalculatedIterationMutex[threadIdx]);
		centerOfMassCalculatedIteration[threadIdx] = 0;
		culumtativeBodiesMasses[threadIdx] = totalLocalMass;
		culumtativeBodiesPositions[threadIdx] = centerOfMassPosition;
		centerOfMassCalculatedIterationCondVariable[threadIdx].notify_all();
	}

	accelerations = std::move(getAcceleration(initialPositions[threadIdx], masses[threadIdx], Config::G, Config::softening, threadIdx, 0));
	

	for (size_t i = 1; i <= config.numberOfTimeSteps(); i++)
	{
#ifdef FILE
		{
			std::lock_guard<std::mutex> l(fileWriteMutex);

			if (threadIdx == 0)
			{
				outputFile << config.startTime + config.dt*i << ":";
				for (size_t m = 0; m < config.numberOfThreads; m++)
				{
					for (size_t j = 0; j < config.bodiesPerThread(); j++)
					{
						outputFile << "(";
						for (size_t k = 0; k < numberOfDimensions; k++)
						{
							outputFile << initialPositions[m][j][k];
							if (k < numberOfDimensions - 1)
							{
								outputFile << ",";
							}
							else
							{
								outputFile << "),";
							}
						}
					}
				}
				outputFile << "\n";
			}
		}
#endif
		//std::cout << "Center of mass position: " << centerOfMassPosition[0] << ", " << centerOfMassPosition[1] << ", " << centerOfMassPosition[2] << std::endl;

		for (size_t j = 0; j < config.bodiesPerThread(); j++)
		{
			for (size_t k = 0; k < numberOfDimensions; k++)
			{
				initialVelocities[threadIdx][j][k] += accelerations[j][k] * config.dt / 2.0;
				initialPositions[threadIdx][j][k] += initialVelocities[threadIdx][j][k] * config.dt;
			}
		}
		//if (threadIdx == 0)
			//printMatrix(initialVelocities[0], "Vel");
		
		totalLocalMass = calculateTotalLocalMass(threadIdx);
		centerOfMassPosition = std::move(calculateCenterOfMassPosition(threadIdx, totalLocalMass));

		{
			std::lock_guard<std::mutex> lock(centerOfMassCalculatedIterationMutex[threadIdx]);
			centerOfMassCalculatedIteration[threadIdx] = i;
			culumtativeBodiesMasses[threadIdx] = totalLocalMass;
			culumtativeBodiesPositions[threadIdx] = centerOfMassPosition;
			centerOfMassCalculatedIterationCondVariable[threadIdx].notify_all();
		}

		//printMatrix(initialPositions[threadIdx], "Planet positions");
		accelerations = std::move(getAcceleration(initialPositions[threadIdx], masses[threadIdx], Config::G, Config::softening, threadIdx, i));
		
		for (size_t j = 0; j < config.bodiesPerThread(); j++)
		{
			for (size_t k = 0; k < numberOfDimensions; k++)
			{
				initialVelocities[threadIdx][j][k] += accelerations[j][k] * config.dt / 2.0;
			}
		}
	}

	//printMatrix(initialPositions[threadIdx], "Planet positions");
}