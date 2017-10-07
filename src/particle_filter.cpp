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

static default_random_engine gen;

void ParticleFilter::init(double x, double y, double theta, double std[]) {
	// TODO: Set the number of particles. Initialize all particles to first position (based on estimates of 
	//   x, y, theta and their uncertainties from GPS) and all weights to 1. 
	// Add random Gaussian noise to each particle.
	// NOTE: Consult particle_filter.h for more information about this method (and others in this file).

	// set particle number
	num_particles = 100;
	//Set standard deviations for x, y, and theta.
	double std_x = std[0];
	double std_y = std[1];
	double std_theta = std[2];
	// creates a normal (Gaussian) distribution for x,y,theta
	normal_distribution<double> dist_x(x, std_x);
	normal_distribution<double> dist_y(y, std_y);
	normal_distribution<double> dist_theta(theta, std_theta);
	
	for (int i = 0; i < num_particles; ++i) {
		double sample_x, sample_y, sample_theta;
	// construct particle struct
		Particle p;
		p.id = i;
		p.x = dist_x(gen);
		p.y = dist_y(gen);
		p.theta = dist_theta(gen);
		p.weight = 1.0;

		particles.push_back(p);
	}

	is_initialized = true;

}


void ParticleFilter::prediction(double delta_t, double std_pos[], double velocity, double yaw_rate) {
	// TODO: Add measurements to each particle and add random Gaussian noise.
	// NOTE: When adding noise you may find std::normal_distribution and std::default_random_engine useful.
	//  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
	//  http://www.cplusplus.com/reference/random/default_random_engine/
	
	//Update particle postions
	for (int i = 0; i < num_particles; ++i) {
			//adding noise to sensor information

		if (fabs(yaw_rate) < 0.00001) {
			particles[i].x += velocity * delta_t * cos(particles[i].theta);
			particles[i].y += velocity * delta_t * sin(particles[i].theta);
		} 
		else {
			particles[i].x += velocity / yaw_rate * (sin(particles[i].theta + yaw_rate*delta_t) - sin(particles[i].theta));
			particles[i].y += velocity / yaw_rate * (cos(particles[i].theta) - cos(particles[i].theta + yaw_rate*delta_t));
			particles[i].theta += yaw_rate * delta_t;
		}
		normal_distribution<double> dist_x(particles[i].x, std_pos[0]);
		normal_distribution<double> dist_y(particles[i].y, std_pos[1]);
		normal_distribution<double> dist_theta(particles[i].theta, std_pos[2]);

		particles[i].x = dist_x(gen);
		particles[i].y = dist_y(gen);
		particles[i].theta = dist_theta(gen);
	}

}

void ParticleFilter::dataAssociation(std::vector<LandmarkObs> predicted, std::vector<LandmarkObs>& observations) {
	// TODO: Find the predicted measurement that is closest to each observed measurement and assign the 
	//   observed measurement to this particular landmark.
	// NOTE: this method will NOT be called by the grading code. But you will probably find it useful to 
	//   implement this method and use it as a helper during the updateWeights phase.
	
	for (int i=0; i<observations.size(); ++i){

		// observations Vector of landmark observations from particle
		LandmarkObs obs=observations[i];

		int near_id = 0;
		double distance_thr = 10000000;

		// predicted Vector of predicted landmark observations from map
		for (int j=0; j<predicted.size(); j++){
			LandmarkObs pred=predicted[j];

			double distance = dist(obs.x,obs.y,pred.x,pred.y);
			
			if (distance < distance_thr){
				distance_thr = distance;
				near_id = pred.id;
			}
		}

		observations[i].id = near_id;
	}
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
	//   and the following is a good resouusrce for the actual equation to implement (look at equation 
	//   3.33
	//   http://planning.cs.uiuc.edu/node99.html


	for (int i=0; i<num_particles; ++i){
		// read in particle position
		double p_x = particles[i].x;
		double p_y = particles[i].y;
		double p_theta = particles[i].theta;

		
		vector<LandmarkObs> predictions;

		// calculte landmarks in the range of sensor (predictions)
		for (int j=0; j<map_landmarks.landmark_list.size(); j++){

			float x_mark = map_landmarks.landmark_list[j].x_f;
			float y_mark = map_landmarks.landmark_list[j].y_f;
			int id_mark = map_landmarks.landmark_list[j].id_i;

		// consider sensor range
			if (fabs(x_mark - p_x) < sensor_range && fabs(y_mark - p_y) < sensor_range){
				LandmarkObs predict;
				predict.x = x_mark;
				predict.y = y_mark;
				predict.id = id_mark;
				predictions.push_back(predict);
			}
		}

		// Transformation
		vector<LandmarkObs> observation_trans;    
		for (int j=0; j<observations.size(); j++){
			double x_transformed = p_x+ cos(p_theta)*observations[j].x - sin(p_theta)*observations[j].y;
			double y_transformed = p_y+ sin(p_theta)*observations[j].x + cos(p_theta)*observations[j].y;
			LandmarkObs o_trans;
			o_trans.x = x_transformed;
			o_trans.y = y_transformed;
			o_trans.id = observations[j].id;
			observation_trans.push_back(o_trans);
		}
				
		// data ossociation, to link observed objects with landmarks, output is make observatrion id equeal to landmark id
		dataAssociation(predictions, observation_trans);

		// Update weight of each particle
		double sig_x = std_landmark[0];
		double sig_y = std_landmark[1];
		//double updated_weight;
		particles[i].weight = 1.0;

		for (int j =0; j<observation_trans.size();j++){

			double x_obs,y_obs,mu_x,mu_y;
			
			x_obs = observation_trans[j].x;
			y_obs = observation_trans[j].y;
			for (int k=0; k<predictions.size(); k++){
				if (predictions[k].id == observation_trans[j].id){
					mu_x = predictions[k].x;
					mu_y = predictions[k].y;
				}
			}
			double gauss_norm = (1/(2 * M_PI * sig_x * sig_y));
			double exponent = ((x_obs - mu_x)*(x_obs - mu_x))/(2 * sig_x*sig_x) + ((y_obs - mu_y)*(y_obs - mu_y))/(2 * sig_y*sig_y);
			double updated_weight =  gauss_norm * exp(-exponent);
			
			//total observations weight
			particles[i].weight *= updated_weight;
		}
	}
}


void ParticleFilter::resample() {
	vector<Particle> new_particles;
	vector<double> weights;

	for (int i=0; i<num_particles; ++i){
	weights.push_back(particles[i].weight);
	}
	 // max weight
	double wmax=*max_element(weights.begin(), weights.end());
	 // Generate random  
	uniform_real_distribution<double> unirealdist(0.0, 1.0);
	double rand_p = unirealdist(gen);

	 // max weight used for spin wheel
	double wmax_temp = wmax*2*rand_p;

	 // create random index
	uniform_int_distribution<int> uniintdist(0, num_particles-1);
	int index = uniintdist(gen);

	 //spin wheel
	for (int i=0; i<num_particles; ++i){
		
		while (wmax_temp>weights[index]){
			wmax_temp=wmax_temp-weights[index];
			index=(index+1)%num_particles;
		}
		new_particles.push_back(particles[index]);
	}
	 //copy to particles
	particles = new_particles;
}

Particle ParticleFilter::SetAssociations(Particle particle, std::vector<int> associations, std::vector<double> sense_x, std::vector<double> sense_y)
{
	//particle: the particle to assign each listed association, and association's (x,y) world coordinates mapping to
	// associations: The landmark id that goes along with each listed association
	// sense_x: the associations x mapping already converted to world coordinates
	// sense_y: the associations y mapping already converted to world coordinates

	//Clear the previous associations
	particle.associations.clear();
	particle.sense_x.clear();
	particle.sense_y.clear();

	particle.associations= associations;
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
