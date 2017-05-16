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

using namespace std;



void ParticleFilter::init(double x, double y, double theta, double std[]) {
	// TODO: Set the number of particles. Initialize all particles to first position (based on estimates of 
	//   x, y, theta and their uncertainties from GPS) and all weights to 1. 
	// Add random Gaussian noise to each particle.
	// NOTE: Consult particle_filter.h for more information about this method (and others in this file).

	default_random_engine gen;
  	normal_distribution<double> N_x(x, std[0]);
	normal_distribution<double> N_y(y, std[1]);
	normal_distribution<double> N_theta(theta, std[2]);

	num_particles = 10;
	Particle p;

	for (int i = 0; i < num_particles; i++) {

		p.x = N_x(gen);
		p.y = N_y(gen);
		p.theta = N_theta(gen);
		p.weight = 1;
		particles.push_back(p);
		weights.push_back(0.0);

	}

	is_initialized = true;

}

void ParticleFilter::prediction(double dt, double std[], double v, double r) {
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

void ParticleFilter::dataAssociation(std::vector<LandmarkObs> predicted, vector<LandmarkObs>& observations) {
	// TODO: Find the predicted measurement that is closest to each observed measurement and assign the 
	//   observed measurement to this particular landmark.
	// NOTE: this method will NOT be called by the grading code. But you will probably find it useful to 
	//   implement this method and use it as a helper during the updateWeights phase.

	double distance;
	double min_distance;

	for (int i = 0; i < predicted.size(); i++) {

		min_distance = 1e6;

		for (int j = 0; j < observations.size(); j++) {

			distance = dist(observations[j].x, observations[j].y,
			                predicted[i].x, predicted[i].y);

			if (distance < min_distance) {

				predicted[i].id = observations[j].id;
				min_distance = distance;

			}

		}

	}


}

void ParticleFilter::updateWeights(double sensor_range, double std_landmark[], 
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

	// LandmarkObs landmark;

	double x_p, y_p, theta_p, costheta, sintheta;

	for (int i = 0; i < particles.size(); i++) {

		x_p = particles[i].x;
		y_p = particles[i].y;
		theta_p = particles[i].theta;
		costheta = cos(theta_p);
		sintheta = sin(theta_p);

		// Convert observations from vehicle to global coordinate system
		std::map<int, LandmarkObs> transformed_observations_map;
		std::vector<LandmarkObs> transformed_observations;
		transformed_observations.resize(observations.size());

		for (int j = 0; j < observations.size(); j++) {

			LandmarkObs tmp;
			tmp.id = observations[j].id;
			tmp.x = x_p + observations[j].x * costheta - observations[j].y * sintheta; 
			tmp.y = y_p + observations[j].x * sintheta + observations[j].y * costheta;
			transformed_observations[j] = tmp;
			transformed_observations_map.insert(std::make_pair(tmp.id, tmp));

		}

		// Find landmarks in sensor range
		std::vector<LandmarkObs> landmarks_in_range;
		for (int j = 0; j < map_landmarks.landmark_list.size(); j++) {

			if (dist(x_p, map_landmarks.landmark_list[j].x_f,
					 y_p, map_landmarks.landmark_list[j].y_f) < sensor_range) {

				landmarks_in_range.push_back(LandmarkObs{map_landmarks.landmark_list[j].id_i,
				                                         map_landmarks.landmark_list[j].x_f,
														 map_landmarks.landmark_list[j].y_f });

			}

		}

		dataAssociation(landmarks_in_range, transformed_observations);

		particles[i].weight = 1.0; //reset the weight of the particle

		for (const auto landmark:landmarks_in_range) {

			double x = landmark.x;
			double y = landmark.y;
			double x_mu = transformed_observations_map[landmark.id].x;
			double y_mu = transformed_observations_map[landmark.id].y;

			double dx = x - x_mu;
			double dy = y - y_mu;
			double std_x = std_landmark[0];
			double std_y = std_landmark[1];
			double num = exp(-0.5 * (dx * dx / (std_x * std_x) + dy * dy / (std_y * std_y)));
			double den = 2 * M_PI * std_x * std_y;

			particles[i].weight *= num / den;

		}

		weights[i] = particles[i].weight; // Assign to weight of ith particle

	}	

}

void ParticleFilter::resample() {
	// TODO: Resample particles with replacement with probability proportional to their weight. 
	// NOTE: You may find std::discrete_distribution helpful here.
	//   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution

	std::random_device rd;
    std::mt19937 gen(rd());
    std::discrete_distribution<> d(weights.begin(), weights.end());
    vector<Particle> particles_new;

    particles_new.resize(num_particles);
    for (int n=0; n < num_particles; ++n) {
        particles_new[n] = particles[d(gen)];
    }

    particles = particles_new;

}

void ParticleFilter::write(string filename) {

	// You don't need to modify this file.
	ofstream dataFile;
	dataFile.open(filename, ios::app);
	
	for (int i = 0; i < num_particles; ++i) {

		dataFile << particles[i].x << " " << particles[i].y << " " << particles[i].theta << "\n";

	}

	dataFile.close();

}
