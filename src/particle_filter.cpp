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
 
   num_particles = 10;
   double std_x = std[0];
   double std_y = std[1];
   double std_theta = std[2];
   std::random_device rd;
   std::mt19937 gen(rd());
   normal_distribution<double> dist_x(x, std_x);
   normal_distribution<double> dist_y(y, std_y);
   normal_distribution<double> dist_theta(theta, std_theta);

   for (int i = 0; i < num_particles; ++i) {
      double sample_x, sample_y, sample_theta;
      sample_x = dist_x(gen);
      sample_y = dist_y(gen);
      sample_theta = dist_theta(gen);
      Particle p = { i, sample_x ,sample_y ,sample_theta, 1.0 };
      particles.push_back(p);
      weights.push_back(1.0);
   }
   is_initialized = true;
}

void ParticleFilter::prediction(double delta_t, double std_pos[], double velocity, double yaw_rate) {
   // TODO: Add measurements to each particle and add random Gaussian noise.
   // NOTE: When adding noise you may find std::normal_distribution and std::default_random_engine useful.
   //  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
   //  http://www.cplusplus.com/reference/random/default_random_engine/

   //cout << __FUNCTION__ << endl;
   std::random_device rd;
   std::mt19937 gen(rd());
   for (int i = 0; i < num_particles; ++i) {
      Particle& p = particles[i];
      double x_new, y_new, theta_new;
      x_new = y_new = theta_new = 0.0;
      if (fabs(yaw_rate) < 0.00001) 
      {
         x_new = p.x + velocity * delta_t * cos(p.theta);
         y_new = p.y + velocity * delta_t * sin(p.theta);
         theta_new = p.theta;
      }
      else
      {
         x_new = p.x + velocity / yaw_rate*(sin(p.theta + yaw_rate*delta_t) - sin(p.theta));
         y_new = p.y + velocity / yaw_rate*(cos(p.theta) - cos(p.theta + yaw_rate*delta_t));
         theta_new = p.theta + yaw_rate*delta_t;
      }
      double std_x = std_pos[0];
      double std_y = std_pos[1];
      double std_theta = std_pos[2];
      normal_distribution<double> dist_x(x_new, std_x);
      normal_distribution<double> dist_y(y_new, std_y);
      normal_distribution<double> dist_theta(theta_new, std_theta);
      p.x = dist_x(gen);
      p.y = dist_y(gen);
      p.theta = dist_theta(gen);
   //   cout << "particle= " << i << ", delta_t=" << delta_t << ", vel=" << velocity << ", yaw_rate=" << yaw_rate << endl;
   //   cout << "- orig(" << p.x << "," << p.x << "," << p.theta << ")" << endl;
   //   cout << "- new (" << x_new << "," << y_new << "," << theta_new << ")" << endl;
   //   cout << "- pred(" << p.x << "," << p.y << "," << p.theta << ")" << endl;
   }
}

void ParticleFilter::dataAssociation(std::vector<LandmarkObs> predicted, std::vector<LandmarkObs>& observations) {
   // TODO: Find the predicted measurement that is closest to each observed measurement and assign the 
   //   observed measurement to this particular landmark.
   // NOTE: this method will NOT be called by the grading code. But you will probably find it useful to 
   //   implement this method and use it as a helper during the updateWeights phase.
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

   // cout << __FUNCTION__ << endl;
   // factor needed for calculate of multi-variate Gaussian distribution  
   const double a = 1 / (2 * M_PI * std_landmark[0] * std_landmark[1]);
   // update weights of all particles
   for (int i = 0; i < num_particles; ++i) {
      Particle& p = particles[i];
      p.weight = 1.0; // initialize weight to 1
      // create lists for associating observations to landmarks
      vector<int> associations;
      vector<double> sense_x;
      vector<double> sense_y;  
      // transform observations from car coordinates to map coordinates
      std::vector<LandmarkObs> trans_observations;
      for (unsigned j = 0; j < observations.size(); ++j) {
         LandmarkObs obs = observations[j];
         LandmarkObs transformed;
         transformed.x = p.x + cos(p.theta)*obs.x - sin(p.theta)*obs.y;
         transformed.y = p.y + sin(p.theta)*obs.x + cos(p.theta)*obs.y;
         //cout << "-- obs: x=" << transformed.x << ", y=" << transformed.y << endl;
         trans_observations.push_back(transformed);
      }
      // associate observations with map landmarks
      for (unsigned j = 0; j < map_landmarks.landmark_list.size(); ++j) {
         Map::single_landmark_s lm = map_landmarks.landmark_list[j];
         // is landmark in range ?
         double dist_landmark = dist(p.x, p.y, lm.x_f, lm.y_f);
         if (dist_landmark <= sensor_range)
         {
            // cout << "-- lm: x=" << lm.x_f << ", y=" << lm.y_f << endl;
            // get nearest observation
            double dist_min = 1000.0;
            LandmarkObs nearest;
            for (unsigned k = 0; k < trans_observations.size(); ++k) {
               LandmarkObs transformed = trans_observations[k];
               double d = dist(lm.x_f, lm.y_f, transformed.x, transformed.y);
               if (d < dist_min)
               {
                  dist_min = d;
                  nearest = transformed;
               }
            }
            //cout << "-- nn: x=" << nearest.x << ", y=" << nearest.y << ", id=" << lm.id_i << endl;

            // calculate multi-variate Gaussian distribution
            double x_diff = nearest.x - lm.x_f;
            double y_diff = nearest.y - lm.y_f;
            double b = ((x_diff * x_diff) / (2 * std_landmark[0] * std_landmark[0])) + ((y_diff * y_diff) / (2 * std_landmark[1] * std_landmark[1]));
            p.weight *= a * exp(-b);
            // remember nearest observation for landmark
            associations.push_back(lm.id_i);
            sense_x.push_back(nearest.x);
            sense_y.push_back(nearest.y);
         }
      }
      // update particle weights with combined multi-variate Gaussian distribution
      weights[i] = p.weight;
      // set assosiations between nearest observations and landmarks for current particle
      p = SetAssociations(p, associations, sense_x, sense_y);
   }
}

void ParticleFilter::resample() {
   // TODO: Resample particles with replacement with probability proportional to their weight. 
   // NOTE: You may find std::discrete_distribution helpful here.
   //   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution

   //std::default_random_engine gen;
   std::random_device rd;
   std::mt19937 gen(rd());
   std::vector<Particle> new_particles;

    // generate random starting index for resampling wheel
   std::uniform_int_distribution<int> uniintdist(0, num_particles - 1);
   auto index = uniintdist(gen);

   // get max weight
   double max_weight = *std::max_element(weights.begin(), weights.end());

   // uniform random distribution [0.0, max_weight)
   std::uniform_real_distribution<double> unirealdist(0.0, max_weight);

   // spin the resample wheel!
   double beta = 0.0;
   for (int i = 0; i < num_particles; i++) {
      beta += unirealdist(gen) * 2.0;
      while (beta > weights[index]) {
         beta -= weights[index];
         index = (index + 1) % num_particles;
      }
      new_particles.push_back(particles[index]);
   }
   particles.swap(new_particles);

   //// Use discrete distribution to return particles by weight
   //for (int i = 0; i < num_particles; ++i) {
   //   discrete_distribution<int> index(weights.begin(), weights.end());
   //   int index_sample = index(gen);
   //   new_particles.push_back(particles[index_sample]);
   //}
   //// Replace particles
   //particles = new_particles;
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

   particle.associations = associations;
   particle.sense_x = sense_x;
   particle.sense_y = sense_y;

   return particle;
}

string ParticleFilter::getAssociations(Particle best)
{
   vector<int> v = best.associations;
   stringstream ss;
   copy(v.begin(), v.end(), ostream_iterator<int>(ss, " "));
   string s = ss.str();
   s = s.substr(0, s.length() - 1);  // get rid of the trailing space
   return s;
}
string ParticleFilter::getSenseX(Particle best)
{
   vector<double> v = best.sense_x;
   stringstream ss;
   copy(v.begin(), v.end(), ostream_iterator<double>(ss, " "));
   string s = ss.str();
   s = s.substr(0, s.length() - 1);  // get rid of the trailing space
   return s;
}
string ParticleFilter::getSenseY(Particle best)
{
   vector<double> v = best.sense_y;
   stringstream ss;
   copy(v.begin(), v.end(), ostream_iterator<double>(ss, " "));
   string s = ss.str();
   s = s.substr(0, s.length() - 1);  // get rid of the trailing space
   return s;
}
