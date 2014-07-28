// Filename:  KukaRMC-rtnetcomponent.hpp
// Copyright: 2014 ISIR-CNRS
// Author:  Vincent SAY 
// Description: Orocos component using RTNET to compute the kuka model
//              from fri data and using a restricted model controller.

#ifndef KUKA_KUKA_RMC_RTNET_COMPONENT_HPP
#define KUKA_KUKA_RMC_RTNET_COMPONENT_HPP

#include <friRTNetExampleAbstract.hpp>
#include "kukafixed.h"
#include <Eigen/Dense>
#include <fstream>
#include <qpOASES.hpp>

class KukaRMControllerRTNET : public FriRTNetExampleAbstract{
	private:
		/* modele */
		kukafixed *m_model;
		
		/* Vecteurs des coordonnees */
		Eigen::VectorXd m_q_mes;
		Eigen::VectorXd m_qp_mes;
		Eigen::VectorXd m_qpp_mes;
		Eigen::Vector3d m_X_des;
		Eigen::Vector3d m_Xp_des;
		Eigen::Vector3d m_Xpp_des;
		Eigen::Vector3d m_X_mes;
		Eigen::Vector3d m_Xp_mes;
		Eigen::Vector3d m_Xpp_mes;
		Eigen::VectorXd m_q_max;
		Eigen::VectorXd m_qp_max;
		
		/* Vecteurs des efforts */
		Eigen::VectorXd m_tau_mes;
		Eigen::VectorXd m_tau_FRI;
		Eigen::VectorXd m_tau_max;		
		
		/* Vecteurs des PID */
		Eigen::Vector3d m_pid;
		Eigen::Vector3d m_Kp;
		Eigen::Vector3d m_Kd;
		
		/* Valeurs de la trajectoire desiree */
		Eigen::Vector3d m_lemniCenter;
		double m_a;
		double m_theta;
		
		/* Variables pour QPOASES */
		qpOASES::SQProblem m_sqp;
		qpOASES::real_t *m_realH;
		qpOASES::real_t *m_realA;
		qpOASES::real_t *m_realg;
		qpOASES::real_t *m_reallb;
		qpOASES::real_t *m_realub;
		qpOASES::real_t *m_reallbA;
		qpOASES::real_t *m_realubA;
		qpOASES::real_t *m_x0;
		qpOASES::real_t *m_cputime;
		qpOASES::Options m_options;
		int m_nWSR;
		
		/* Matrice et vecteurs dynamiques */
		Eigen::MatrixXd m_H;
		Eigen::MatrixXd m_J;
		Eigen::MatrixXd m_Jp;
		Eigen::VectorXd m_c;
		Eigen::VectorXd m_g;
		
		/* Variables pour estimation de F */
		Eigen::MatrixXd m_FTab;
		Eigen::VectorXd m_bufSum;
		double m_FSat;
		unsigned int m_Fidx;
		bool m_FBoucle;
		unsigned int m_delta; //Fenetre de moyennage
		/* Variables pour AMAF */
		unsigned int m_SL;
		unsigned int m_L;
		unsigned int m_FL;
		Eigen::MatrixXd m_buf;
		Eigen::VectorXd m_prec;
		Eigen::VectorXd m_AMA;
		double m_fastLen;
		double m_slowLen;
		
		
		double m_dt; //Temps echantillonnage
		double m_time; //Temps actuel
		bool m_writeInFile;
		std::ofstream m_fluxErr;
		std::string m_fluxErrName;
		std::ofstream m_fluxTrq;
		std::string m_fluxTrqName;
		std::ofstream m_flux;
		
	public:	
		KukaRMControllerRTNET(std::string const& name);
		~KukaRMControllerRTNET();
		void updateHook();
		void nextDesVal();
		void estimateF(Eigen::Vector3d& F);
		void estimateAMAF(Eigen::Vector3d& F);
		void fdm(std::vector<double> const& qp_new);
		void pidCalc();
		void updateQpVar(Eigen::Vector3d F);
		
		void setTheta(double const Theta);
		void setA(double const A);
		void setXinit(std::vector<double> const& Xinit);
		void setName(std::string const& fluxErrName, std::string const& fluxTrqName);
		void setAllowWritingData(bool writeInFile);
		void setMAFDelta(unsigned int const& delta);
		void setAMAFLength(unsigned int const& FL, unsigned int const& L, unsigned int const& SL);
		void setJointImpedance(std::vector<double> &stiffness, std::vector<double> &damping);
};

#endif
