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
#include <map>

#include "particle_filter.h"

void ParticleFilter::init(double x, double y, double theta, double std[]) {
	// TODO: Set the number of particles. Initialize all particles to first position (based on estimates of 
	//   x, y, theta and their uncertainties from GPS) and all weights to 1. 
	// Add random Gaussian noise to each particle.
	// NOTE: Consult particle_filter.h for more information about this method (and others in this file).

	num_particles = 10;

	// Create normal distributions
	std::default_random_engine generator;
	std::normal_distribution<double> dist_x(x, std[0]);
	std::normal_distribution<double> dist_y(y, std[1]);
	std::normal_distribution<double> dist_theta(theta, std[2]);

	// Pre-allocate vectors
	particles.resize(num_particles);
	weights.resize(num_particles);

	for (unsigned int i = 0; i < num_particles; i++) {

		particles[i].id = i;
		particles[i].x = dist_x(generator);
		particles[i].y = dist_y(generator);
		particles[i].theta = dist_theta(generator);
		particles[i].weight = 1.0;
		weights[i] = 1.0;

	}

	is_initialized = true;

}

void ParticleFilter::prediction(double delta_t, double std_pos[], double velocity, double yaw_rate) {
	// TODO: Add measurements to each particle and add random Gaussian noise.
	// NOTE: When adding noise you may find std::normal_distribution and std::default_random_engine useful.
	//  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
	//  http://www.cplusplus.com/reference/random/default_random_engine/

	double r = yaw_rate;
	double v = velocity;
	double d_t = delta_t;
	double x_0, y_0, theta_0;
	double x_f, y_f, theta_f;

	std::default_random_engine generator;

	// Address divide by zero when yaw rate is 0
	for (unsigned int i = 0; i < num_particles; i++) {

		x_0 = particles[i].x;
		y_0 = particles[i].y;
		theta_0 = particles[i].theta;

		if (fabs(r) > 1e-6) {

			x_f = x_0 + v/r * (sin(theta_0 + r * d_t) - sin(theta_0));
			y_f = y_0 + v/r * (cos(theta_0) - cos(theta_0 + r * d_t));
			theta_f = theta_0 + r * d_t;

		} else {

			x_f = x_0 + v * d_t * cos(theta_0);
			y_f = y_0 + v * d_t * sin(theta_0);
			theta_f = theta_0;

		}

		// Update the particle position with the prediction and add gaussian noise
		std::normal_distribution<double> dist_x(x_f, std_pos[0]);
		std::normal_distribution<double> dist_y(y_f, std_pos[1]);
		std::normal_distribution<double> dist_theta(theta_f, std_pos[2]);

		particles[i].x = dist_x(generator);
		particles[i].y = dist_y(generator);
		particles[i].theta = dist_theta(generator);

	}
	
}

void ParticleFilter::dataAssociation(std::vector<LandmarkObs> predicted, std::vector<LandmarkObs>& observations) {
	// TODO: Find the predicted measurement that is closest to each observed measurement and assign the 
	//   observed measurement to this particular landmark.
	// NOTE: this method will NOT be called by the grading code. But you will probably find it useful to 
	//   implement this method and use it as a helper during the updateWeights phase.

	double distance, min_distance;

	for (unsigned int i = 0; i < observations.size(); i++) {

		min_distance = 1e10;

		for (unsigned int j = 0; j < predicted.size(); j++) {

			distance = dist(observations[i].x, observations[i].y,
			                predicted[j].x, predicted[j].y);

			if (distance < min_distance) {

				observations[i].id = predicted[j].id;
				min_distance = distance;

			}

		}


	}

}

void ParticleFilter::updateWeights(double sensor_range, double std_landmark[], 
		std::vector<LandmarkObs> observations, Map map_landmarks) {
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

	// Create a map to be able to access landmarks by id
	std::map<int, Map::single_landmark_s> landmarks_id_map;

	for (unsigned int i = 0; i < map_landmarks.landmark_list.size(); i++) {

		landmarks_id_map.insert(std::make_pair(map_landmarks.landmark_list[i].id_i,
		                                       map_landmarks.landmark_list[i]));

	}

	// Start the main particle loop
	for (unsigned int i = 0; i < num_particles; i++) {

		// Extract particle position
		double x_p = particles[i].x;
		double y_p = particles[i].y;
		double theta_p = particles[i].theta;
		double cos_t = cos(theta_p);
		double sin_t = sin(theta_p);

		// Convert observations from vehicle to map coordinate system
		double x_m, y_m;
		std::vector<LandmarkObs> observations_transformed;
		observations_transformed.resize(observations.size());

		for (unsigned int j = 0; j < observations.size(); j++) {

			x_m = x_p + observations[j].x * cos_t - observations[j].y * sin_t; 
			y_m = y_p + observations[j].x * sin_t + observations[j].y * cos_t;
			observations_transformed[j] = LandmarkObs{-1, x_m, y_m};	

		}

		// Find landmarks in range
		std::vector<LandmarkObs> landmarks_in_range;
		for (unsigned int j = 0; j < map_landmarks.landmark_list.size(); j++) {

			double distance  = dist(x_p, y_p,
			                        map_landmarks.landmark_list[j].x_f,
									map_landmarks.landmark_list[j].y_f);

			if (distance < sensor_range) {

				landmarks_in_range.push_back(LandmarkObs{map_landmarks.landmark_list[j].id_i,
				                                         map_landmarks.landmark_list[j].x_f,
														 map_landmarks.landmark_list[j].y_f });

			}

		}

		// Calculate weights
		if (landmarks_in_range.size() > 0) {

			// Associate each measurement with the respective nearest neighbor landmark in range
			dataAssociation(landmarks_in_range, observations_transformed);

			particles[i].weight = 1.0; //reset the weight of the particle

			for (const auto obs:observations_transformed) {

				double x = landmarks_id_map[obs.id].x_f;
				double y = landmarks_id_map[obs.id].y_f;
				double x_mu = obs.x;
				double y_mu = obs.y;

				double dx = x - x_mu;
				double dy = y - y_mu;
				double std_x = std_landmark[0];
				double std_y = std_landmark[1];
				double num = exp(-0.5 * (dx * dx / (std_x * std_x) + dy * dy / (std_y * std_y)));
				double den = 2 * M_PI * std_x * std_y;
				double prob = num/den;

				particles[i].weight *= prob;

			}

			weights[i] = particles[i].weight; // Assign to weight of ith particle

		} else {

			weights[i] = 0.0;

		}

	}

}

void ParticleFilter::resample() {
	// TODO: Resample particles with replacement with probability proportional to their weight. 
	// NOTE: You may find std::discrete_distribution helpful here.
	//   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution
	
	// std::random_device rd;
    // std::mt19937 gen(rd());
	std::default_random_engine generator;
    std::discrete_distribution<> dist(weights.begin(), weights.end());
    std::vector<Particle> particles_new;

    particles_new.resize(num_particles);
    for (unsigned int i = 0; i < num_particles; i++) {
        particles_new[i] = particles[dist(generator)];
    }

    particles = particles_new;

}

void ParticleFilter::write(std::string filename) {
	// You don't need to modify this file.
	std::ofstream dataFile;
	dataFile.open(filename, std::ios::app);
	for (int i = 0; i < num_particles; ++i) {
		dataFile << particles[i].x << " " << particles[i].y << " " << particles[i].theta << "\n";
	}
	dataFile.close();
}
