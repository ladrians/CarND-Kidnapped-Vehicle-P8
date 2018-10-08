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
#include <math.h> 
#include <iostream>
#include <sstream>
#include <string>
#include <iterator>

#include "particle_filter.h"

using namespace std;

void ParticleFilter::init(double x, double y, double theta, double std[]) {
	// TODO: Set the number of particles. Initialize all particles to first position (based on estimates of 
	//   x, y, theta and their uncertainties from GPS) and all weights to 1. 
	// Add random Gaussian noise to each particle.
	// NOTE: Consult particle_filter.h for more information about this method (and others in this file).

	if (initialized() == true)
		return;

	cout << "Initialization" << endl;
	num_particles = 75; // tested with 10, 50, 75 (speed vs efficiency)

	default_random_engine gen;

	normal_distribution<double> dist_x(x, std[0]);
	normal_distribution<double> dist_y(y, std[1]);
	normal_distribution<double> dist_theta(theta, std[2]);

	for (int i=0; i < num_particles; i++)
	{
		Particle p;
		p.id = i;
		p.x = dist_x(gen);
		p.y = dist_y(gen);
		p.theta = dist_theta(gen);
		p.weight = 1.0;

		particles.push_back(p);
		weights.push_back(p.weight);
	}

	is_initialized = true;
}

void ParticleFilter::prediction(double delta_t, double std_pos[], double velocity, double yaw_rate) {
	// TODO: Add measurements to each particle and add random Gaussian noise.
	// NOTE: When adding noise you may find std::normal_distribution and std::default_random_engine useful.
	//  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
	//  http://www.cplusplus.com/reference/random/default_random_engine/

	default_random_engine gen;

	for (int i=0; i < num_particles; i++)
	{
		double new_x;
		double new_y;
		double new_theta;

		if (fabs(yaw_rate) < 1e-5)
		{
			new_x = particles[i].x + velocity*delta_t*cos(particles[i].theta);
			new_y = particles[i].y + velocity*delta_t*sin(particles[i].theta);
			new_theta = particles[i].theta;
		}
		else
		{
			new_x = particles[i].x + (velocity/yaw_rate)*(sin(particles[i].theta+(yaw_rate*delta_t))-sin(particles[i].theta));
			new_y = particles[i].y + (velocity/yaw_rate)*(cos(particles[i].theta)-cos(particles[i].theta+(yaw_rate*delta_t)));
			new_theta = particles[i].theta + yaw_rate*delta_t;
		}
		
		particles[i].x = new_x;
		particles[i].y = new_y;
		particles[i].theta = new_theta;
		//cout << "** Prediction(" << particles[i].x << ", " << particles[i].y << ", " << particles[i].theta << ")" << endl;
	}
}

void ParticleFilter::dataAssociation(std::vector<LandmarkObs> predicted, std::vector<LandmarkObs>& observations) {
	// TODO: Find the predicted measurement that is closest to each observed measurement and assign the 
	//   observed measurement to this particular landmark.
	// NOTE: this method will NOT be called by the grading code. But you will probably find it useful to 
	//   implement this method and use it as a helper during the updateWeights phase.

	// Data association on the updateWeights method
}

void ParticleFilter::updateWeights(double sensor_range, double std_landmark[], 
		const std::vector<LandmarkObs> &observations, const Map &map_landmarks) {
	// TODO: Update the weights of each particle using a mult-variate Gaussian distribution. You can read
	//   more about this distribution here: https://en.wikipedia.org/wiki/Multivariate_normal_distribution
	// NOTE: The observations are given in the VEHICLE'S coordinate system. Your particles are located
	//   according to the MAP'S coordinate system. You will need to transform between the two systems.
	//   Keep in mind that this transformation requires both rotation AND translation (but no scaling).
	//   The following is a good resource for the theory:
	//   https://www.willamette.edu/~gorr/classes/GeneralGraphics/Transforms/transforms2d.htm
	//   and the following is a good resource for the actual equation to implement (look at equation 
	//   3.33
	//   http://planning.cs.uiuc.edu/node99.html

/* Steps to follow from tip https://discussions.udacity.com/t/dont-get-the-point-of-transformation/248991/5
  for each particle
    for each observation  
        observation transform to map coordinates  
        for each landmarks   
             calculate the euclidiean distance from the observation to the landmark
        assign landmark to the observation with the shortest distance
        calculate weight*/

	double normalizer = (1 / (2 * M_PI*std_landmark[0] * std_landmark[1]));
	double x2 = 2 * std_landmark[0] * std_landmark[0];
	double y2 = 2 * std_landmark[1] * std_landmark[1];

	int particles_count = particles.size();
	for (int i=0; i < particles_count; i++)
	{	
		std::vector<int> associations;
		std::vector<double> sense_x;
		std::vector<double> sense_y;

		particles[i].weight = 1.0;

		// observation transform to map coordinates
		LandmarkObs landmark_t;
		int observations_count = observations.size();
		for (int j=0; j < observations_count; j++)
		{	
			landmark_t = observations[j];
			LandmarkObs obs_t;
			
			obs_t.id = landmark_t.id;
			obs_t.x = particles[i].x + (cos(particles[i].theta) * landmark_t.x) - (sin(particles[i].theta) * landmark_t.y);
			obs_t.y = particles[i].y + (sin(particles[i].theta) * landmark_t.x) + (cos(particles[i].theta) * landmark_t.y);

			double min_distance = sensor_range;
			int min_distance_idx = -1;

			for (int k = 0; k < map_landmarks.landmark_list.size(); k++)
			{
				// calculate the euclidiean distance from the observation to the landmark
				double dists = dist(map_landmarks.landmark_list[k].x_f, map_landmarks.landmark_list[k].y_f, obs_t.x, obs_t.y);

				if (dists < min_distance)
				{
					min_distance = dists;
					min_distance_idx = k;
				}
			}

			if (min_distance_idx != -1)
			{ // assign landmark to the observation with the shortest distance
				
				double mu_x = map_landmarks.landmark_list[min_distance_idx].x_f;
				double mu_y = map_landmarks.landmark_list[min_distance_idx].y_f;
				
				double exponent_term = (obs_t.x - mu_x)*(obs_t.x - mu_x) / x2 + (obs_t.y - mu_y)*(obs_t.y - mu_y) / y2;

				long double new_weight = normalizer*exp(-exponent_term);
				
				particles[i].weight *= new_weight;
			}
			else
				particles[i].weight = 0;

			associations.push_back(min_distance_idx + 1);
			sense_x.push_back(obs_t.x);
			sense_y.push_back(obs_t.y);
		}
		particles[i] = SetAssociations(particles[i], associations, sense_x, sense_y);
	}
}

void ParticleFilter::resample() {
	// TODO: Resample particles with replacement with probability proportional to their weight. 
	// NOTE: You may find std::discrete_distribution helpful here.
	//   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution

	default_random_engine gen;
	discrete_distribution<int> dist_w(weights.begin(), weights.end());

	std::vector<Particle> resample_particles;

	for (int i=0; i < num_particles; i++)
	{
		resample_particles.push_back(particles[dist_w(gen)]);
	}
	particles = resample_particles;
}

Particle ParticleFilter::SetAssociations(Particle& particle, const std::vector<int>& associations, 
                                     const std::vector<double>& sense_x, const std::vector<double>& sense_y)
{
    //particle: the particle to assign each listed association, and association's (x,y) world coordinates mapping to
    // associations: The landmark id that goes along with each listed association
    // sense_x: the associations x mapping already converted to world coordinates
    // sense_y: the associations y mapping already converted to world coordinates

	// Clear the previous associations
	particle.associations.clear();
	particle.sense_x.clear();
	particle.sense_y.clear();

	// Associate
    particle.associations = associations;
    particle.sense_x = sense_x;
    particle.sense_y = sense_y;

    return particle;
}

string ParticleFilter::getAssociations(Particle best)
{
	vector<int> v = best.associations;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<int>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
string ParticleFilter::getSenseX(Particle best)
{
	vector<double> v = best.sense_x;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<float>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
string ParticleFilter::getSenseY(Particle best)
{
	vector<double> v = best.sense_y;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<float>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
