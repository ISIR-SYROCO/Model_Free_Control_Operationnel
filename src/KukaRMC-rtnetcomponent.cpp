// Filename:  KukaRMControllerRTNET-rtnetcomponent.cpp
// Copyright: 2014 ISIR-CNRS
// Author:  Vincent SAMY (vsamy95@gmail.com) 
// Description: This is a control law in cartesian space designs for the KUKA LWR

#include "KukaRMC-rtnetcomponent.hpp"
#include <rtt/Component.hpp>
#include <iostream>
#include <Eigen/LU>

#define PI	3.141592653589793238462643383280 /* pi */
#define CARTDOF 3 //Orientation is not taken here

KukaRMControllerRTNET::KukaRMControllerRTNET(std::string const& name) 
	: m_dt(0.001) 							//Sampling time
	, m_time(0.0) 							//Actual time
	, m_Fidx(0) 							//Index for moving average
	, m_FBoucle(false)						//Loop for moving average
	, m_delta(30)							//Window of the moving average
	, m_FSat(200.0)							//Limits of F
	, m_SL(30)								//AMAF's SlowLength
	, m_L(10)								//AMAF's Length
	, m_FL(2)								//AMAF's FastLength
	, m_fastLen(1)							//AMAF's actual fast length
	, m_slowLen(1)							//AMAF's actual slow length
	, m_a(0.0)								//Lemniscate constant (width)
	, m_theta(0.0)							//Lemniscate constant (angular frequency)
	, m_writeInFile(false)					//Enable/Disable data saving
	, m_sqp(14, 7)							//Initialization of qp solver
	, FriRTNetExampleAbstract(name)			//Initialization of inherited class
{
	/* Add functions to be used by the .ops file */
    this->addOperation("setTheta", &KukaRMControllerRTNET::setTheta, this, RTT::OwnThread);
    this->addOperation("setA", &KukaRMControllerRTNET::setA, this, RTT::OwnThread);
    this->addOperation("setXinit", &KukaRMControllerRTNET::setXinit, this, RTT::OwnThread);
    this->addOperation("setName", &KukaRMControllerRTNET::setName, this, RTT::OwnThread);
    this->addOperation("setAllowWritingData", &KukaRMControllerRTNET::setAllowWritingData, this, RTT::OwnThread);
    this->addOperation("setMAFDelta", &KukaRMControllerRTNET::setMAFDelta, this, RTT::OwnThread);
    this->addOperation("setAMAFLength", &KukaRMControllerRTNET::setAMAFLength, this, RTT::OwnThread);
    this->addOperation("setJointImpedance", &KukaRMControllerRTNET::setJointImpedance, this, RTT::OwnThread);
 	
    Eigen::MatrixXd J(6,LWRDOF);
    Eigen::MatrixXd Jp(6,LWRDOF);
    
	/* Resize all vectors */
	m_q_mes.resize(LWRDOF);
	m_qp_mes.resize(LWRDOF);
	m_qpp_mes.resize(LWRDOF);
	m_tau_max.resize(LWRDOF);
	m_tau_mes.resize(LWRDOF);
	m_tau_FRI.resize(LWRDOF);
	m_q_max.resize(LWRDOF);
	m_qp_max.resize(LWRDOF);
	m_H.resize(LWRDOF,LWRDOF);
	m_J.resize(CARTDOF,LWRDOF); 
	m_Jp.resize(CARTDOF,LWRDOF);
	m_c.resize(LWRDOF);
	m_g.resize(LWRDOF);
	m_FTab.resize(CARTDOF,m_delta);
	m_buf.resize(CARTDOF,m_L);

	/* PID's constants */
	m_Kp << 100.0, 100.0, 100.0; //Comment if Orientation is considered
	m_Kd << 0.1, 0.1, 0.1;
	//m_Kp << 100.0, 100.0, 100.0, 10.0, 10.0, 10.0; //Uncomment if Orientation is considered
	//m_Kd << 0.1, 0.1, 0.1, 0.1, 0.1, 0.1;
	
   	/* KUKA's constrains */
   	m_tau_max << 176.0, 176.0, 100.0, 100.0, 100.0, 30.0, 30.0;
   	//m_tau_max << 30.0, 30.0, 20.0, 20.0, 20.0, 15.0, 15.0;
	m_q_max << 165.0*PI/180.0, 115.0*PI/180.0, 165.0*PI/180.0, 115.0*PI/180.0, 165.0*PI/180.0, 115.0*PI/180.0, 165.0*PI/180.0;
	m_qp_max << 112.5*PI/180.0, 112.5*PI/180.0, 112.5*PI/180.0, 112.5*PI/180.0, 180.0*PI/180.0, 112.5*PI/180.0, 112.5*PI/180.0;
	//m_qp_max << 30.5*PI/180.0, 30.5*PI/180.0, 30.5*PI/180.0, 30.5*PI/180.0, 30.0*PI/180.0, 30.5*PI/180.0, 30.5*PI/180.0;

	/* Initialization to zero */
	m_pid.setZero();
	m_X_des.setZero();
	m_Xp_des.setZero();
	m_Xpp_des.setZero();
	m_X_mes.setZero();
	m_Xp_mes.setZero();
	m_Xpp_mes.setZero();
	m_qp_mes.setZero();
	m_qpp_mes.setZero();
	m_tau_mes.setZero();
	m_tau_FRI.setZero();
	m_bufSum.setZero();
	m_FTab.setZero();
	m_prec.setZero();
	m_AMA.setZero();
	m_buf.setZero();
	m_J.setZero();
	m_Jp.setZero();
	m_lemniCenter.setZero();
    
    /* Get the KUKA lWR model */
    m_model = new kukafixed("kuka");
    
    /* Model dynamics */
    m_H 	= m_model->getInertiaMatrix();
    m_c 	= m_model->getNonLinearTerms();
    m_g 	= m_model->getGravityTerms();
    J 		= m_model->getSegmentJacobian(7);
    Jp 		= m_model->getSegmentJdot(7);
    m_J 	= J.bottomRows(3);
    m_Jp 	= Jp.bottomRows(3);
    
    /* QP pointers */
    m_realH 	= new qpOASES::real_t[14*14];
    m_realA 	= new qpOASES::real_t[7*14];
    m_realg	 	= new qpOASES::real_t[14];
    m_realub 	= new qpOASES::real_t[14];
    m_reallb 	= new qpOASES::real_t[14];
    m_realubA	= new qpOASES::real_t[7];
    m_reallbA 	= new qpOASES::real_t[7];
    m_x0 		= new qpOASES::real_t[14];
    m_cputime 	= new qpOASES::real_t;
    m_nWSR 		= 15;
    *m_cputime 	= 0.0005;
    qpOASES::Bounds guessedBounds(14);
    qpOASES::Constraints guessedConstraints(7);
    guessedBounds.setupAllLower();
	guessedConstraints.setupAllInactive();
    
    m_x0[0] = 0.0;
    m_x0[0] = -5.03;
    m_x0[0] = 0.01;
    m_x0[0] = -6.02;
    m_x0[0] = -0.01;
    m_x0[0] = -0.07;
    m_x0[0] = 0.0;
    for(unsigned int i=7;i<14;i++){
    	m_x0[i] = 0.0;
    }
    
    /* Update QP matrices */
    updateQpVar(Eigen::Vector3d::Zero());
    
    /* Set QP options */
    m_options.setToReliable();
    m_options.printLevel = qpOASES::PL_DEBUG_ITER;
    m_sqp.setOptions(m_options);
    
    /* First QP solving */
    qpOASES::real_t *nullpointer = NULL; //Can take all time it needs
    int initWSR = 100; //Can do 100 iterations
    qpOASES::returnValue callReturn;
    callReturn = m_sqp.init(m_realH, m_realg, m_realA, m_reallb, m_realub, m_reallbA, m_realubA, initWSR, nullpointer, m_x0, nullpointer, &guessedBounds, &guessedConstraints);
    
    /* Check successfulness */
    if(callReturn==qpOASES::SUCCESSFUL_RETURN)
    {
    	std::cout << std::endl << "=========>QP initialization is a SUCCESS<=========" << std::endl << std::endl;
    	std::cout << "The found value is x_0 =" << std::endl;
    	m_sqp.getPrimalSolution(m_x0);
    	for(unsigned int i=0;i<14;i++){
    		std::cout << m_x0[i] << std::endl;
    	}
    	std::cout << std::endl;
    	m_options.printLevel = qpOASES::PL_NONE;
    	m_sqp.setOptions(m_options);
    }
	else if(callReturn==qpOASES::RET_MAX_NWSR_REACHED){
		std::cout << "QP could not be solved within nWSR value, it must be higher" << std::endl;
		return;
    }
	else if(callReturn==qpOASES::RET_INIT_FAILED){
		std::cout << "QP initialization fails" << std::endl;
		return;
    }
}

KukaRMControllerRTNET::~KukaRMControllerRTNET()
{
	/* Close files */
	if(m_writeInFile)
	{
		m_fluxErr.seekp(4,std::ios::end);
		m_fluxErr << "];%";
		m_fluxErr.close();
		m_fluxTrq.seekp(4,std::ios::end);
		m_fluxTrq << "];%";
		m_fluxTrq.close();
		m_flux.close();
	}
	/* delete pointers */
	delete[] m_model;
	m_model = 0;
	delete[] m_realH;
	m_realH = 0;
	delete[] m_realA;
	m_realA = 0;
	delete[] m_realg;
    m_realg = 0;
    delete[] m_realub;
    m_realub = 0;
    delete[] m_reallb;
    m_reallb = 0;
    delete[] m_realubA;
    m_realubA = 0;
    delete[] m_reallbA;
    m_reallbA = 0;
    delete[] m_x0;
    m_x0 = 0;
    delete[] m_cputime;
    m_cputime = 0;
}

/* This function is called by the FRI every 1ms */
void KukaRMControllerRTNET::updateHook()
{
	/* Get robot's state */
    std::string fri_mode("e_fri_unkown_mode");
    bool fri_cmd_mode = false;
    RTT::FlowStatus fs_event = iport_events.read(fri_mode);
    if (fri_mode == "e_fri_cmd_mode")
        fri_cmd_mode = true;
    else if (fri_mode == "e_fri_mon_mode")
        fri_cmd_mode = false;

    Eigen::VectorXd tau_dyn(LWRDOF);
    Eigen::Vector3d F;
    Eigen::MatrixXd J(6,LWRDOF);
    Eigen::MatrixXd Jp(6,LWRDOF);
    qpOASES::returnValue callReturn;
	
	std::vector<double> jntPos(LWRDOF);
	std::vector<double> jntVel(LWRDOF);
	std::vector<double> jntTrq(LWRDOF);
	std::vector<double> tau_FRI(LWRDOF);
	std::vector<double> jntDesPos(LWRDOF);
	geometry_msgs::Pose cartPos;
	
	/* Get sensors measure */
	RTT::FlowStatus joint_pos_fs = iport_msr_joint_pos.read(jntPos);
	RTT::FlowStatus joint_vel_fs = iport_msr_joint_vel.read(jntVel);
	RTT::FlowStatus joint_trq_fs = iport_msr_joint_trq.read(jntTrq);
	RTT::FlowStatus cartPos_fs   = iport_cart_pos.read(cartPos);
	
	if(fri_cmd_mode) //Check if the component has the control of the KUKA
	{
		/* Update of the mesures */
		if(joint_trq_fs==RTT::NewData)
		{
			for(unsigned int i=0;i<LWRDOF;i++){
				m_tau_mes[i] = jntTrq[i];
			}
		}
		if(joint_pos_fs==RTT::NewData)
		{
			for(unsigned int i=0;i<LWRDOF;i++){
				m_q_mes[i] = jntPos[i];
			}
			m_model->setJointPositions(m_q_mes); //Update model positions
		}
		if(joint_vel_fs==RTT::NewData)
		{
			fdm(jntVel); //calculate the measured acceleration with finite difference
			for(unsigned int i=0;i<LWRDOF;i++){
				m_qp_mes[i] = jntVel[i];
			}
			m_model->setJointVelocities(m_qp_mes); //Update model velocities
		}
		if(cartPos_fs==RTT::NewData)
		{
        	m_X_mes(0) = (double)cartPos.position.x;
        	m_X_mes(1) = (double)cartPos.position.y;
        	m_X_mes(2) = (double)cartPos.position.z;
        }		
		
		/* Update variables */
		tau_dyn = m_tau_mes-m_tau_FRI;
		m_H = m_model->getInertiaMatrix();
		m_c = m_model->getNonLinearTerms();
		m_g = m_model->getGravityTerms();
		J = m_model->getSegmentJacobian(m_model->nbSegments());
    	Jp = m_model->getSegmentJdot(m_model->nbSegments());
    	m_J = J.bottomRows(3);
    	m_Jp = Jp.bottomRows(3);
		m_Xpp_mes = m_J*m_qpp_mes+m_Jp*m_qp_mes; //Xpp = J*qpp+Jp*qp
		m_Xp_mes = m_J*m_qp_mes; //Xp = J*qp
		
		/* Estimate F */
		F = m_J*m_qpp_mes-m_J*m_H.inverse()*(m_tau_FRI-m_c-m_g-m_J.transpose()*m_pid);
		//Xpp-beta-Jv/H*u+LambdaInv*pid

		estimateF(F); //Comment if you want to use MAF method
		//estimateAMAF(F); //Uncomment if you want to use AMAF method
		
		/* Update time */
		m_time = m_time+m_dt; 
		
		/* Calculate positions, velocities and desired accelerations */
		nextDesVal();
		
		/* Calculate the PID */
		pidCalc();
		m_Xpp_des = m_pid;
		
		/* Calculate the control */
		updateQpVar(F);
		callReturn = m_sqp.hotstart(m_realH,m_realg,m_realA,m_reallb,m_realub,m_reallbA,m_realubA,m_nWSR,m_cputime);
		if(callReturn==qpOASES::SUCCESSFUL_RETURN)
		{
			m_sqp.getPrimalSolution(m_x0);
			for(unsigned int i=0;i<LWRDOF;i++){
				m_tau_FRI[i] = m_x0[i];
			}
		}
		else{
			std::cout << "QP in-loop fails error" << std::endl;
		}
		
		/* Save datas */
		if(m_writeInFile)
		{
			m_flux << std::endl;
			m_flux << "temps : " << m_time << std::endl;
			//m_flux << "X_des : " << m_q_des.transpose() << std::endl;
			m_flux << "X_mes : " << m_q_mes.transpose() << std::endl;
			//m_flux << "Xp_des : " << m_qp_des.transpose() << std::endl;
			//m_flux << "Xp_mes : " << m_qp_mes.transpose() << std::endl;
			//m_flux << "Xpp_des : " << m_qp_des.transpose() << std::endl;
			//m_flux << "Xpp_mes : " << m_qp_mes.transpose() << std::endl;
			//m_flux << "F_est_sat : " << F.transpose() << std::endl;
			//m_flux << "PID : " << m_pid.transpose() << std::endl;
			//m_flux << "Hqpp : " << (m_H*m_qpp_des).transpose() << std::endl;
			//m_flux << "coriolis : " << m_c.transpose() << std::endl;
			//m_flux << "gravity : " << m_g.transpose() << std::endl;
			//m_flux << "tau_FRI : " << m_tau_FRI.transpose() << std::endl;
			m_flux << "f_dynamics : " << tau_dyn.transpose() << std::endl;
			m_flux << "gravity : " << m_g.transpose() << std::endl;
			//m_flux << "tau_mes : " << m_tau_mes.transpose() << std::endl;
			m_fluxTrq << std::endl << m_tau_mes.transpose() << ";...";
			//m_fluxErr << std::endl << (m_X_des-m_X_mes).transpose() << ";...";
		}
		
		/* Output datas to the KUKA */
		for(unsigned int i=0;i<LWRDOF;i++){
			tau_FRI[i] = m_tau_FRI[i];
		}
		//oport_add_joint_trq.write(tau_FRI);
		oport_joint_position.write(jntPos);
	}
}

/* Function that calculates deesired position */
void KukaRMControllerRTNET::nextDesVal()
{
	Eigen::VectorXd period(1);
	Eigen::VectorXd a(1);
	
	period << m_theta*m_time;
	
	if(m_time<1.0){
		a << m_a*m_time;
	}
	else{
		a << m_a;
	}
	
	m_X_des(1) = a(0)*(period.array().sin())(0)+m_lemniCenter(1);
	m_X_des(2) = a(0)*(period.array().sin())(0)*(period.array().cos())(0)+m_lemniCenter(2);
}

/* Mobile Average Filter */
void KukaRMControllerRTNET::estimateF(Eigen::Vector3d& F)
{
	for(unsigned int i=0;i<CARTDOF;i++)
	{
		m_bufSum(i) = m_bufSum(i)-m_FTab(i,m_Fidx);
		m_FTab(i,m_Fidx) = F(i);
		m_bufSum(i) = m_bufSum(i)+F(i);
	}
	
	m_Fidx = m_Fidx+1;
	if(m_Fidx>m_delta-1){
		m_Fidx = 0;
		m_FBoucle = true;
	}
	
	if(m_FBoucle){
		F = m_bufSum/m_delta;
	}
	else{
		F = m_bufSum/m_Fidx;
	}
	
	for(unsigned int i=0;i<CARTDOF;i++)
	{
		if((F.array().abs())(i)>m_FSat)
		{
			if(F(i)>0.0){
				F(i) = m_FSat;
			}
			else{
				F(i) = -m_FSat;
			}
		}
	}
}

/* Adaptive Mobile Average Filter */
void KukaRMControllerRTNET::estimateAMAF(Eigen::Vector3d& F)
{
	Eigen::Vector3d er;
	Eigen::Vector3d sc;
	Eigen::Vector3d precL;
	double fast = 2/(m_fastLen+1);
	double slow = 2/(m_slowLen+1);
	
	for(unsigned int i=0;i<CARTDOF;i++)
	{
		m_bufSum(i) = m_bufSum(i)-m_FTab(i,m_Fidx);
		
		/* Memorization of the nbAxe-m_lg value */
    	if(m_FBoucle){
        	precL(i) = m_buf(i,m_Fidx);
        }
    	else{
        	precL(i) = m_buf(i,0);
    	}
    
    	/* Update the i-ieme value */
    	m_buf(i,m_Fidx) = ((F-m_prec).array().abs())(i);
    	
    	/* Update the sum */
    	m_bufSum(i) = m_bufSum(i)+m_buf(i,m_Fidx); 
    	
    	/* Calculate ER */
    	if(m_bufSum(i)==0){
        	er(i) = 0;
        }
    	else{
        	er(i) = ((F-precL).array().abs())(i)/m_bufSum(i);
    	}
    	
    	/* Calculate sc */
    	sc(i) = (er(i)*(fast-slow)+slow)*(er(i)*(fast-slow)+slow);
    
    	/* Calculate the moving average */
    	m_AMA(i) = m_AMA(i)*(1.0-sc(i))+sc(i)*F(i);
    	
    	/* Output */
    	F(i) = m_AMA(i);
    	
    	/* Update */
    	m_prec(i) = m_buf(i,m_Fidx);
	}
	
	/* Update parameters for the next use */
    m_Fidx=m_Fidx+1;
    if(m_Fidx>m_L)
    {
        m_Fidx = 1;
        m_FBoucle = 1;
    }
    if(m_fastLen<m_FL){
        m_fastLen = m_fastLen+1;
    }
    if(m_slowLen<m_SL){
        m_slowLen = m_slowLen+1;
    }
}

/* Finite Difference Method */
void KukaRMControllerRTNET::fdm(std::vector<double> const& qp_new)
{
	for(unsigned int i=0;i<LWRDOF;i++){
		m_qpp_mes[i] = (qp_new[i]-m_qp_mes[i])/m_dt;
	}
}

/* PID */
void KukaRMControllerRTNET::pidCalc()
{
	for(unsigned int i=0;i<CARTDOF;i++)
	{
		m_pid(i) = m_Kp(i)*(m_X_des(i)-m_X_mes(i))-m_Kd(i)*m_Xp_mes(i);
	}
}

/* Update the qp variables
 *
 * The main is : 	min(x) 0.5*||Cx-d||^2
 *					x = tau, qpp
 *					lb <= x <= ub 
 *					lbA <= A*x <= ubA
 *
 * It can be written as 	min(x) 0.5*x'Hx+x'g
 *				  where		H = C'*C	
 *							g = -C'd
 *
 */
void KukaRMControllerRTNET::updateQpVar(Eigen::Vector3d F)
{
	Eigen::MatrixXd C(3,14),Aeq(7,14),H(14,14);
	Eigen::Vector3d d;
	Eigen::VectorXd beq(7),LB(14),UB(14),v(7),w(7),g(14);
	unsigned const int h2 = 70;
	
	C << m_J*m_H.inverse(), Eigen::MatrixXd::Zero(3,7); 
	d << m_Xpp_mes+m_J*m_H.inverse()*m_J.transpose()*m_pid-m_Jp*m_qp_mes+m_J*m_H.inverse()*(m_c+m_g)-F;
	Aeq << Eigen::MatrixXd::Identity(7,7), -m_H;
	beq << m_c+m_g;
	v << (-m_qp_max-m_qp_mes)/h2;
	w << 2*(-m_q_max-(m_q_mes+m_qp_mes*h2))/(h2*h2);
	LB << -m_tau_max, v.array().max(w.array());
	v = (m_qp_max-m_qp_mes)/h2;
	w = 2*(m_q_max-(m_q_mes+m_qp_mes*h2))/(h2*h2);
	UB << m_tau_max, v.array().min(w.array());
	
	H << C.transpose()*C;
	g << -C.transpose()*d;
	for(unsigned int i=0;i<14;i++)
	{
    	for(unsigned int j=0;j<14;j++){
    		m_realH[i*j] = H(i,j);
    	}
    	m_realg[i] = g(i);
    	m_realub[i] = UB(i);
    	m_reallb[i] = LB(i);
    }
    
    for(unsigned int i=0;i<7;i++)
	{
    	for(unsigned int j=0;j<14;j++){
    		m_realA[i*j] = Aeq(i,j);
    	}
    	m_realubA[i] = beq(i);
    	m_reallbA[i] = beq(i);
    }
}

/* Set lemniscate angle frequency */
void KukaRMControllerRTNET::setTheta(double const Theta)
{
	m_theta = Theta;
}

/* Set lemniscate width */
void KukaRMControllerRTNET::setA(double const A)
{
	m_a = A;
}

void KukaRMControllerRTNET::setXinit(std::vector<double> const& Xinit)
{
	for(unsigned int i=0;i<CARTDOF;i++){
		m_X_des(i) = Xinit[i];
		m_lemniCenter(i) = Xinit[i];
	}
}

/* Set files names */
void KukaRMControllerRTNET::setName(std::string const& fluxErrName, std::string const& fluxTrqName)
{
	m_fluxErrName = fluxErrName;
	m_fluxTrqName = fluxTrqName;
}

/* Set and open flux for writing into files */
void KukaRMControllerRTNET::setAllowWritingData(bool writeInFile)
{
	std::string ss;	
	
	m_writeInFile = writeInFile;
	if(m_writeInFile)
	{
		if(!(m_fluxErrName==""))
		{
			ss = "/home/kuka/src/groovy_workspace/model_free_control_art/"+m_fluxErrName;
			m_fluxErr.open(ss.c_str());
    		m_fluxErr << "err_q = ["; 
    		ss = "/home/kuka/src/groovy_workspace/model_free_control_art/"+m_fluxTrqName;
    		m_fluxTrq.open(ss.c_str());
    		m_fluxTrq << "tau_mes = [..."; 
    		m_flux.open("/home/kuka/src/groovy_workspace/model_free_control_art/data.dat");
    	}
    	else{
    		std::cout << "Error ! Name of files has not been provided" << std::endl;
    	}
    }
}

void KukaRMControllerRTNET::setAMAFLength(unsigned int const& FL, unsigned int const& L, unsigned int const& SL)
{
	m_SL = SL;
	m_FL = FL;
	m_L  =  L;
	m_buf.resize(CARTDOF,m_L);
	m_buf.setZero();
}

void KukaRMControllerRTNET::setMAFDelta(unsigned int const& delta)
{
	m_delta = delta;
	m_FTab.resize(CARTDOF,m_delta);
	m_FTab.setZero();
}

void KukaRMControllerRTNET::setJointImpedance(std::vector<double> &stiffness, std::vector<double> &damping)
{
	if(stiffness.size() != LWRDOF || damping.size() != LWRDOF)
	{
		std::cout << "Wrong vector size, should be " << LWRDOF << ", " << LWRDOF << std::endl;
		return;
	}
	else{
		lwr_fri::FriJointImpedance joint_impedance_command;
		for(unsigned int i = 0; i < LWRDOF; i++){
			joint_impedance_command.stiffness[i] = stiffness[i];
			joint_impedance_command.damping[i] = damping[i];
		}
		oport_joint_impedance.write(joint_impedance_command);
	}
}

/*
 * Using this macro, only one component may live
 * in one library *and* you may *not* link this library
 * with another component library. Use
 * ORO_CREATE_COMPONENT_TYPE()
 * ORO_LIST_COMPONENT_TYPE(kuka_ming_xing_controller)
 * In case you want to link with another library that
 * already contains components.
 *
 * If you have put your component class
 * in a namespace, don't forget to add it here too:
 */
ORO_CREATE_COMPONENT(KukaRMControllerRTNET)
