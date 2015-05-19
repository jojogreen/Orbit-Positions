// Orbit Positions.cpp : Defines the entry point for the console application.
//
// This code is designed to calculate the postion of a planet at a certain date.
// The user inputs the body and date they wish to find the position of the planet
// and the output will be a radius and an angle. It utilizes a bisection method to
// solve Kepler's equation. Additionally, it uses a 2-body approximation.
//
// Coded by Evan Kaiser

#include "stdafx.h"
#include <iostream>
#include <iomanip>
#include <math.h>
#include <string>

using namespace std;

const double pi = 3.14159265358979323846;
// Semi-Major Axis (meters)
double A[9] = { 57909050000, 108939000000, 149598261000, 227939100000, 778547200000, 1433449370000, 2870671400, 4498542600000, 5874000000000 };
//e_e = eccentricity
double e_e[9] = { .2056, .0068, 0.0167, .0934, .0484, .0542, .0472, .0086, .24188 };
//I_i = inclination degrees from the ecliptic plane
double I_i[9] = { 7.005, 3.3947, 0, 1.851, 1.305, 2.484, .77, 1.769, 17.142 };
//m_0 = mean anomoly radians
double m_0[9] = { 9.89293 * 10 * pi / 180, 4.10610 * 10 * pi / 180, 1.2608 * 100 * pi / 180, 8.001909 * 10 * pi / 180,
1.261179 * 100 * pi / 180, 1.4399 * 100 * pi / 180, 2.0842 * 100 * pi / 180, 2.708158 * 100 * pi / 180, 3.87325 * 10 * pi / 180 };
//argument of periapse radians
double omega_mat[9] = { 0.5091, 0.9536, 4.5120, 5.0009, 4.7826, 5.9367, 1.6901, 5.1615, 1.9536 };
//Longitude of Ascending Node radians
double OMEGA_MAT[9] = { 4.8311 * 10 * pi / 180, 7.6639 * 10 * pi / 180, 2.06567 * 100 * pi / 180, 4.95128 * 10 * pi / 180, 1.005 * 100 * pi / 180,
1.13544 * 100 * pi / 180, 7.3943 * 10 * pi / 180, 1.31914 * 100 * pi / 180, 1.10357 * 100 * pi / 180, };
//True Anomaly radians
double THETA[9] = { 1.20775 * 100 * pi / 180, 1.2077500 * 100 * pi / 180, 4.15732 * 10 * pi / 180, 127.603,
9.07213 * 10 * pi / 180, 1.3047 * 100 * pi / 180, 1.4742 * 100 * pi / 180, 2.05873 * 100 * pi / 180,
2.69886 * 100 * pi / 180, };
double M = 1.7565670*pow(10, 28);

//This solves Kepler's equation using the bisection method
double Bisection(double m, double e){
	double EHi = 2 * pi, ELO = 0, E, Error = 500;
	for (int i = 0; Error > pow(10, -5) || i < 200; i++)
	{
		E = (EHi + ELO) / 2;
		Error = E - e*sin(E) - m;
		if (Error == 0){
			break;
		}
		else if (Error > 0){
			EHi = E;
		}
		else if (Error < 0){
			ELO = E;
		}
	}
	return E;
}

// This calculates the position of the planet at a time "T"
void Position(string Plan, char day){
	double i; //inclination
	double e; //eccentricity
	double a; // semi-major axis
	double OMEGA;
	double omega;
	double G = 6.67384*pow(10, -11);
	double m_o, th;

	//get orbital parameters from the array
	if (Plan.compare("Mercury") == 0 || Plan.compare("mercury") == 0){
		e = e_e[0];
		a = A[0];
		i = I_i[0];
		OMEGA = OMEGA_MAT[0];
		omega = omega_mat[0];
		m_o = m_0[0];
		th = THETA[0];

	}
	else if (Plan.compare("Venus") == 0 || Plan.compare("venus") == 0){
		e = e_e[1];
		a = A[1];
		i = I_i[1];
		OMEGA = OMEGA_MAT[1];
		omega = omega_mat[1];
		m_o = m_0[1];
		th = THETA[1];

	}
	else if (Plan.compare("Earth") == 0 || Plan.compare("earth") == 0){
		e = e_e[2];
		a = A[2];
		i = I_i[2];
		OMEGA = OMEGA_MAT[2];
		omega = omega_mat[2];
		m_o = m_0[2];
		th = THETA[2];

	}
	else if (Plan.compare("Mars") == 0 || Plan.compare("mars") == 0){
		e = e_e[3];
		a = A[3];
		i = I_i[3];
		OMEGA = OMEGA_MAT[3];
		omega = omega_mat[3];
		m_o = m_0[3];
		th = THETA[3];

	}
	else if (Plan.compare("Jupiter") == 0 || Plan.compare("jupiter") == 0){
		e = e_e[4];
		a = A[4];
		i = I_i[4];
		OMEGA = OMEGA_MAT[4];
		omega = omega_mat[4];
		m_o = m_0[4];
		th = THETA[4];

	}
	else if (Plan.compare("Saturn") == 0 || Plan.compare("saturn") == 0){
		e = e_e[5];
		a = A[5];
		i = I_i[5];
		OMEGA = OMEGA_MAT[5];
		omega = omega_mat[5];
		m_o = m_0[5];
		th = THETA[5];

	}
	else if (Plan.compare("Uranus") == 0 || Plan.compare("uranus") == 0){
		e = e_e[6];
		a = A[6];
		i = I_i[6];
		OMEGA = OMEGA_MAT[6];
		omega = omega_mat[6];
		m_o = m_0[6];
		th = THETA[6];

	}
	else if (Plan.compare("Neptune") == 0 || Plan.compare("neptune") == 0){
		e = e_e[7];
		a = A[7];
		i = I_i[7];
		OMEGA = OMEGA_MAT[7];
		omega = omega_mat[7];
		m_o = m_0[7];
		th = THETA[7];

	}
	else if (Plan.compare("Pluto") == 0 || Plan.compare("pluto") == 0){
		e = e_e[8];
		a = A[8];
		i = I_i[8];
		OMEGA = OMEGA_MAT[8];
		omega = omega_mat[8];
		m_o = m_0[8];
		th = THETA[8];

	}
	else{
		cout << endl << "WARNING!!! SPELLING ERROR IN BODY. ASSUMING EARTH" << endl << endl;
		e = e_e[2];
		a = A[2];
		i = I_i[2];
		OMEGA = OMEGA_MAT[2];
		omega = omega_mat[2];
		m_o = m_0[2];

	}

	//Output eccentricity inclination, etc
	cout << "eccentricity (e) = " << e << endl;
	cout << "inclination (i) = " << i << "degrees" << endl;
	//find "n"
	double n = sqrt(G*M / (pow(a, 3)));
	//put time in seconds
	double t = day * 3600 * 24;
	//E-eSin(E) = m is keplers equation
	double m = n*t + m_o;
	m = fmod(m, 2 * pi);
	//solve kepler's equation using the bisection method
	double E = Bisection(m, e);
	//calculate true anomaly
	double theta = 2 * atan(sqrt((1 + e) / (1 - e))*tan(E / 2)) + 2 * pi;
	theta = fmod(theta * 180 / pi, 360);
	// calculate radius from eccentricity, semi-major axis, ant true anomaly
	double radius = (a*(1 - pow(e, 2))) / (1 + e*cos(theta));
	cout << "True Anomaly = " << theta << " (deg)" << endl;
	cout << "Radius from Central body = " << radius << " (m)" << endl;

}




int _tmain(int argc, _TCHAR* argv[])
{
	char Planet[15];
	int day;
	int deltaT;
	bool cont = true;
	cout << "ORBIT POSITION DATA" << endl;
	//run the code until user exits
	while (cont == true)
	{
		//Get planet and time we wish to find true anomaly and radius of
		cout << "%-----------------------------------------------------------------%" << endl;
		cout << "What is the orbiting body? ";
		cin >> setw(15) >> Planet;
		cin.ignore();
		cout << "Orbiting body is " << Planet << endl;
		cout << "What day do you want the position for (Julian Date)? ";
		cin >> day;
		while (!cin)
		{
			cin.clear();
			cin.ignore();
			cout << "What day do you want the position (Julian Date)? ";
			cin >> day;
		}
		deltaT = 2457156 - day;
		cout << "Position Day is " << day << endl;
		string planet(Planet);
		//call position function
		Position(planet, deltaT);
		cin.ignore();
		cout << "%-----------------------------------------------------------------%" << endl;
		cout << endl;
	}
	return 0;
}

