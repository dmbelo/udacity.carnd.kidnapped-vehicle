/*
 * particle_filter.cpp
 *
 *  Created on: Dec 12, 2016
 *      Author: Tiffany Huang
 */

#include <random>
#include <algorithm>
#include <iostream>
#include <numeric>

#include "particle_filter.h"

using namespace std;



void ParticleFilter::Initialize(double x, double y, double theta, double std[]) {
	// TODO: Set the number of particles. Initialize all particles to first position (based on estimates of 
	//   x, y, theta and their uncertainties from GPS) and all weights to 1. 
	// Add random Gaussian noise to each particle.
	// NOTE: Consult particle_filter.h for more information about this method (and others in this file).

	default_random_engine gen;
  	normal_distribution<double> N_x(x, std[0]);
	normal_distribution<double> N_y(y, std[1]);
	normal_distribution<double> N_theta(theta, std[2]);

	num_particles = 500;
	Particle p;

	for (int i = 0; i < num_particles; i++) {

		p.x = N_x(gen);
		p.y = N_y(gen);
		p.theta = N_theta(gen);
		p.weight = 1;
		particles.push_back(p);

	}

}

void ParticleFilter::Prediction(double dt, double std[], double v, double r) {
	// TODO: Add measurements to each particle and add random Gaussian noise.
	// NOTE: When adding noise you may find std::normal_distribution and std::default_random_engine useful.
	//  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
	//  http://www.cplusplus.com/reference/random/default_random_engine/

	// Add noise to account for uncertainty in control inputs
	default_random_engine gen;

	for (int i = 0; i < num_particles; i++) {

		double x = particles[i].x;
		double y = particles[i].y;
		double theta = particles[i].theta;
		normal_distribution<double> N_x(x, std[0]);
		normal_distribution<double> N_y(y, std[1]);
		normal_distribution<double> N_theta(theta, std[2]);

		// Avoid divide by zero
		if (fabs(r) > 1e-3) {

			double tmp1 = v / r;
			double tmp2 = theta + r * dt;
			particles[i].x += tmp1 * (sin(tmp2) - sin(theta)) + N_x(gen);
			particles[i].y += tmp1 * (cos(theta) - cos(tmp2)) + N_y(gen);
			particles[i].theta += r * dt + N_theta(gen);

		}
		else {

			particles[i].x += v * dt * cos(theta) + N_x(gen);
			particles[i].y += v * dt * sin(theta) + N_y(gen);
			particles[i].theta += theta + N_theta(gen);

		}
		

	}

}

void ParticleFilter::DataAssociation(vector<LandmarkObs> predicted, vector<LandmarkObs>& observations) {
	// TODO: Find the predicted measurement that is closest to each observed measurement and assign the 
	//   observed measurement to this particular landmark.
	// NOTE: this method will NOT be called by the grading code. But you will probably find it useful to 
	//   implement this method and use it as a helper during the updateWeights phase.

	double distance;
	double min_distance = 1e6;

	for (int i = 0; i < observations.size(); i++) {

		for (int j = 0; j < predicted.size(); j++) {

			distance = dist(observations[i].x, observations[i].y,
			                predicted[i].x, predicted[i].y);

			if (distance < min_distance) {
				min_distance = distance;
			}

		}

	}


}

void ParticleFilter::UpdateWeights(double sensor_range, double std_landmark[], 
		vector<LandmarkObs> observations, Map map_landmarks) {
	// TODO: Update the weights of each particle using a mult-variate Gaussian distribution. You can read
	//   more about this distribution here: https://en.wikipedia.org/wiki/Multivariate_normal_distribution
	// NOTE: The observations are given in the VEHICLE'S coordinate system. Your particles are located
	//   according to the MAP'S coordinate system. You will need to transform between the two systems.
	//   Keep in mind that this transformation requires both rotation AND translation (but no scaling).
	//   The following is a good resource for the theory:
	//   https://www.willamette.edu/~gorr/classes/GeneralGraphics/Transforms/transforms2d.htm
	//   and the following is a good resource for the actual equation to implement (look at equation 
	//   3.33. Note that you'll need to switch the minus sign in that equation to a plus to account 
	//   for the fact that the map's y-axis actually points downwards.)
	//   http://planning.cs.uiuc.edu/node99.html


}

void ParticleFilter::Resample() {
	// TODO: Resample particles with replacement with probability proportional to their weight. 
	// NOTE: You may find std::discrete_distribution helpful here.
	//   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution

}

void ParticleFilter::Write(string filename) {

	// You don't need to modify this file.
	ofstream dataFile;
	dataFile.open(filename, ios::app);
	
	for (int i = 0; i < num_particles; ++i) {

		dataFile << particles[i].x << " " << particles[i].y << " " << particles[i].theta << "\n";

	}

	dataFile.close();

}
