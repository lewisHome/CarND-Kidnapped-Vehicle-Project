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
#include <map>

#include "particle_filter.h"

using namespace std;

void ParticleFilter::init(double x, double y, double theta, double std[]) {
	// TODO: Set the number of particles. Initialize all particles to first position (based on estimates of 
	//   x, y, theta and their uncertainties from GPS) and all weights to 1. 
	// Add random Gaussian noise to each particle.
	// NOTE: Consult particle_filter.h for more information about this method (and others in this file).
	
	//could probably put this in header file
	default_random_engine gen;
	
	//set standard deviations for x, y and theta
	float std_x;
	float std_y;
	float std_theta;
	
	std_x = std[0];
	std_y = std[1];
	std_theta = std[2];
	
	//generater normal gaussian distribution about GPS points
  normal_distribution<double> dist_x(x, std_x);
  normal_distribution<double> dist_y(y, std_y);
  normal_distribution<double> dist_theta(theta, std_theta);

  num_particles = 100;
  double weight_=1.0/(float)num_particles;
  cout<<"Initial weight "<<weight_<<endl;
  
  for(int i =0; i<num_particles; ++i)
  {
    Particle particle_;
    
    particle_.id     = i;
    particle_.x      = dist_x(gen);
    particle_.y      = dist_y(gen);
    particle_.theta  = dist_theta(gen);
    particle_.weight = weight_;
    particles.push_back(particle_);

    weights.push_back(weight_);
  }
  cout<<"Init complete.\nNum particles: "<<num_particles<<endl; 
	is_initialized=true;
	particles_new.resize(num_particles);
}

void ParticleFilter::prediction(double delta_t, double std_pos[], double velocity, double yaw_rate) {
	// TODO: Add measurements to each particle and add random Gaussian noise.
	// NOTE: When adding noise you may find std::normal_distribution and std::default_random_engine useful.
	//  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
	//  http://www.cplusplus.com/reference/random/default_random_engine/

  //generater normal gaussian distribution for adding noise to predicted position
  default_random_engine gen;	
  normal_distribution<double> dist_x(0.0, std_pos[0]);
  normal_distribution<double> dist_y(0.0, std_pos[1]);
  normal_distribution<double> dist_theta(0.0, std_pos[2]);
	
	//predict new postion for each particle
	for (int i=0; i<num_particles; ++i)
	{
	  //avoid division by 0
	  if (fabs(yaw_rate) < 0.00001)
	  {
  	  particles[i].x += velocity*delta_t*cos(particles[i].theta);
  	  particles[i].y += velocity*delta_t*sin(particles[i].theta);
	  }
	  else
	  {
  	  particles[i].x += (velocity/yaw_rate)*(sin(particles[i].theta + yaw_rate*delta_t)-sin(particles[i].theta));
  	  particles[i].y += (velocity/yaw_rate)*(cos(particles[i].theta) -cos(particles[i].theta + yaw_rate * delta_t));
      particles[i].theta += yaw_rate * delta_t;  
    }
    
    //add noise    
    particles[i].x      += dist_x(gen);
    particles[i].y      += dist_y(gen);
    particles[i].theta  += dist_theta(gen);
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
	//   and the following is a good resource for the actual equation to implement (look at equation 
	//   3.33
	//   http://planning.cs.uiuc.edu/node99.html
	//transfer observations to map co-ordinates
  //calc normalisation term
  double gauss_norm = 1.0/(2.0* M_PI * std_landmark[0] * std_landmark[1]);


	for (int i = 0; i < num_particles; ++i)
	{
	  //clear previous time step data
  	long double weight_ = 1.0;
    particles[i].sense_x.clear();
    particles[i].sense_y.clear();
    particles[i].associations.clear();
    
    //find landmarks in range
    vector<LandmarkObs> predictions;
    for(unsigned int j=0; j<map_landmarks.landmark_list.size();++j)
    {
      double lm_x = map_landmarks.landmark_list[j].x_f;
      double lm_y = map_landmarks.landmark_list[j].y_f;
      int lm_id = map_landmarks.landmark_list[j].id_i;
      
      double landmark_range = pow(pow((lm_x - particles[i].x),2.0)+pow((lm_y - particles[i].y),2.0),0.5);
      
      if (landmark_range<=sensor_range)
      {
        predictions.push_back(LandmarkObs{lm_id, lm_x, lm_y});
      }
    }
    
	  //transform observations to map co-ordinates
	  for (unsigned int j=0; j < observations.size(); ++j)
	  {
	    double sense_x_ = particles[i].x + (cos(particles[i].theta)*observations[j].x-sin(particles[i].theta)*observations[j].y);
	    
	    double sense_y_ = particles[i].y + (sin(particles[i].theta)*observations[j].x+cos(particles[i].theta)*observations[j].y);

	  //Associate particles with landmarks by nearest neighbour search
      double best_distance = sensor_range;
      double x_min, y_min;
      int id;
	    for(unsigned int k=0; k<predictions.size(); ++k)
	    {
	      double distance = dist(predictions[k].x, predictions[k].y,sense_x_,sense_y_);
	      if (distance< best_distance)
	      {
	        best_distance = distance;
	        id = predictions[k].id;
	        x_min=predictions[k].x;
	        y_min=predictions[k].y;
	      }
	    }	    
      
      //store particle observation data for debug
	    particles[i].sense_x.push_back(sense_x_);
	    particles[i].sense_y.push_back(sense_y_);
	    particles[i].associations.push_back(id);

      //calc particle weight
      double Xexponent = 0.5 * pow((sense_x_ - x_min)/std_landmark[0],2.0);
      double Yexponent = 0.5 * pow((sense_y_ - y_min)/std_landmark[1],2.0);
      weight_*=gauss_norm*exp(-(Xexponent+Yexponent));
    }
    particles[i].weight = weight_;
    weights[i]= weight_;
//    cout<<"particle Number: "<<i<<" Weight val: "<<weights[i]<<endl;
  }
}


void ParticleFilter::resample() {
	// TODO: Resample particles with replacement with probability proportional to their weight. 
	// NOTE: You may find std::discrete_distribution helpful here.
	//   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution
//std::random_device rd;
  default_random_engine gen;
  discrete_distribution<int> index(weights.begin(),weights.end()); 

  for (int n=0; n<num_particles; ++n)
  {
  
    int genID=index(gen);
    particles_new[n] =particles[genID];
//    cout<<"particle "<<n<<" taken from old particle "<<genID<<endl;        
  }
  particles = particles_new;
}

Particle ParticleFilter::SetAssociations(Particle& particle, const std::vector<int>& associations, 
                                     const std::vector<double>& sense_x, const std::vector<double>& sense_y)
{
    //particle: the particle to assign each listed association, and association's (x,y) world coordinates mapping to
    // associations: The landmark id that goes along with each listed association
    // sense_x: the associations x mapping already converted to world coordinates
    // sense_y: the associations y mapping already converted to world coordinates

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
