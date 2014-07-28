
#include "kukafixed.h"

#include <map>
#include <vector>
#include <iostream>

#include <Eigen/Eigenvalues>

typedef  Eigen::Displacementd::AdjointMatrix  AdjointMatrix;

AdjointMatrix getAdjacencyMatrix(const Eigen::Twistd& tw)
{

    AdjointMatrix dAdj;

    dAdj << 0    ,-tw(2), tw(1), 0    , 0    , 0    ,
            tw(2), 0    ,-tw(0), 0    , 0    , 0    ,
           -tw(1), tw(0), 0    , 0    , 0    , 0    ,
            0    ,-tw(5), tw(4), 0    ,-tw(2), tw(1),
            tw(5), 0    ,-tw(3), tw(2), 0    ,-tw(0),
           -tw(4), tw(3), 0    , -tw(1), tw(0), 0    ;

    return dAdj;

}

Eigen::Matrix3d getSkewSymmetricMatrix(const Eigen::Vector3d& vec)
{

    Eigen::Matrix3d mat;

    mat << 0     ,-vec(2), vec(1),
           vec(2), 0     ,-vec(0),
          -vec(1), vec(0), 0     ;

    return mat;

}

class RevoluteJoint
{
public:
    //=================== Variables ====================//
    int                          parent;
    int                          child;
    int                          internaldof;
    int                          dof;
    Eigen::Displacementd         H_seg_joint;
    Eigen::Matrix<double,6,1>    jacobian;
    Eigen::Matrix<double,6,1>    djacobian;
    Eigen::Vector3d              axis;

    RevoluteJoint(int _parent, int _child, int _internaldof, int _dof, const Eigen::Displacementd& _H_seg_joint, const Eigen::Vector3d& _axis)
    : parent(_parent)
    , child(_child)
    , internaldof(_internaldof)
    , dof(_dof)
    , H_seg_joint(_H_seg_joint)
    , jacobian(Eigen::Matrix<double,6,1>::Zero())
    , djacobian(Eigen::Matrix<double,6,1>::Zero())
    , axis(_axis)
    {
        jacobian << axis[0], axis[1], axis[2], 0, 0, 0 ;
    };

    Eigen::Displacementd get_H_joint_child(double q_dof)
    {
        double val = q_dof/2.;
        double c = std::cos(val);
        double s = std::sin(val);
        return Eigen::Displacementd(0,0,0,c, axis[0]*s,axis[1]*s,axis[2]*s);
    }

    Eigen::Twistd get_twist(double dq_dof)
    {
        return Eigen::Twistd(axis[0]*dq_dof,axis[1]*dq_dof,axis[2]*dq_dof,0,0,0);
    }

    AdjointMatrix getdAdjointInverse(double q_dof, double dq_dof)
    {

        AdjointMatrix  iAdj = get_H_joint_child( q_dof ).inverse().adjoint();

        Eigen::Twistd  tw  = - iAdj * get_twist( dq_dof );

        return iAdj * getAdjacencyMatrix(tw);

    }


};

struct kukafixed::Pimpl
{
public:
    //=============== General Variables ================//
    int                                      nbSeg;
    Eigen::VectorXd                          actuatedDofs;
    Eigen::VectorXd                          lowerLimits;
    Eigen::VectorXd                          upperLimits;
    Eigen::VectorXd                          q;
    Eigen::VectorXd                          dq;
    Eigen::Displacementd                     Hroot;
    Eigen::Twistd                            Troot;
    Eigen::VectorXd                          alldq;
    bool                                     isFreeFlyer;
    int                                      nbDofs;
    Eigen::Matrix<double,6,1>                gravity_dtwist;

    //=============== Dynamic Variables ================//
    Eigen::MatrixXd                          M;
    Eigen::MatrixXd                          Minv;
    Eigen::MatrixXd                          B;
    Eigen::VectorXd                          n;
    Eigen::VectorXd                          l;
    Eigen::VectorXd                          g;

    //================= CoM Variables ==================//
    double                                   total_mass;
    Eigen::Vector3d                          comPosition;
    Eigen::Vector3d                          comVelocity;
    Eigen::Vector3d                          comJdotQdot;
    Eigen::Matrix<double,3,Eigen::Dynamic>   comJacobian;
    Eigen::Matrix<double,3,Eigen::Dynamic>   comJacobianDot;

    //=============== Segments Variables ===============//
    std::vector< double >                                   segMass;
    std::vector< Eigen::Vector3d >                          segCoM;
    std::vector< Eigen::Matrix<double,6,6> >                segMassMatrix;
    std::vector< Eigen::Vector3d >                          segMomentsOfInertia;
    std::vector< Eigen::Rotation3d >                        segInertiaAxes;
    std::vector< Eigen::Displacementd >                     segPosition;
    std::vector< Eigen::Twistd >                            segVelocity;
    std::vector< Eigen::Matrix<double,6,Eigen::Dynamic> >   segJacobian;
    std::vector< Eigen::Matrix<double,6,Eigen::Dynamic> >   segJdot;
    std::vector< Eigen::Matrix<double,6,Eigen::Dynamic> >   segJointJacobian;
    std::vector< Eigen::Twistd >                            segJdotQdot;

    //================ Other Variables =================//
    std::map< std::string, int >             segIndexFromName;
    std::vector< std::string >               segNameFromIndex;
    bool                                     isUpToDate;

    //================= Update Status ==================//
    bool                           M_isUpToDate;
    bool                           Minv_isUpToDate;
    bool                           n_isUpToDate;
    bool                           g_isUpToDate;
    bool                           comPosition_isUpToDate;
    bool                           comVelocity_isUpToDate;
    bool                           comJdotQdot_isUpToDate;
    bool                           comJacobian_isUpToDate;
    bool                           comJacobianDot_isUpToDate;
    std::vector< bool >            segPosition_isUpToDate;
    std::vector< bool >            segVelocity_isUpToDate;
    std::vector< bool >            segJacobian_isUpToDate;
    std::vector< bool >            segJdot_isUpToDate;
    std::vector< bool >            segJdotQdot_isUpToDate;
    std::vector< bool >            segNonLinearEffects_isUpToDate;

    //============== Recursive Variables ===============//
    std::vector< std::vector< RevoluteJoint* > >  segJoints;
    std::vector< Eigen::Matrix<double,6,6> >      segNonLinearEffects;
    Eigen::Matrix<double,6,Eigen::Dynamic>        Jroot;
    Eigen::Matrix<double,6,Eigen::Dynamic>        dJroot;

    Pimpl()
    : nbSeg(8)
    , actuatedDofs(7)
    , lowerLimits(7)
    , upperLimits(7)
    , q(7)
    , dq(7)
    , Hroot(0,0,0)
    , Troot(0,0,0,0,0,0)
    , M(7,7)
    , Minv(7,7)
    , B(7,7)
    , n(7)
    , l(7)
    , g(7)
    , total_mass(0)
    , comPosition(0,0,0)
    , comVelocity(0,0,0)
    , comJdotQdot(0,0,0)
    , comJacobian(3,7)
    , comJacobianDot(3,7)
    , segMass(8, double(0))
    , segCoM(8, Eigen::Vector3d(0,0,0))
    , segMassMatrix(8, Eigen::Matrix<double,6,6>())
    , segMomentsOfInertia(8, Eigen::Vector3d(0,0,0))
    , segInertiaAxes(8, Eigen::Rotation3d())
    , segPosition(8, Eigen::Displacementd(0,0,0))
    , segVelocity(8, Eigen::Twistd(0,0,0,0,0,0))
    , segJacobian(8, Eigen::Matrix<double,6,Eigen::Dynamic>(6,7))
    , segJdot(8, Eigen::Matrix<double,6,Eigen::Dynamic>(6,7))
    , segJointJacobian(8, Eigen::Matrix<double,6,Eigen::Dynamic>(6,7))
    , segJdotQdot(8, Eigen::Twistd(0,0,0,0,0,0))
    , segIndexFromName()
    , segNameFromIndex(8, std::string(""))
    , isUpToDate(false)
    , alldq(7)
    , isFreeFlyer(false)
    , nbDofs(7)
    , M_isUpToDate(false)
    , Minv_isUpToDate(false)
    , n_isUpToDate(false)
    , g_isUpToDate(false)
    , comPosition_isUpToDate(false)
    , comVelocity_isUpToDate(false)
    , comJdotQdot_isUpToDate(false)
    , comJacobian_isUpToDate(false)
    , comJacobianDot_isUpToDate(false)
    , segPosition_isUpToDate(8, bool(false))
    , segVelocity_isUpToDate(8, bool(false))
    , segJacobian_isUpToDate(8, bool(false))
    , segJdot_isUpToDate(8, bool(false))
    , segJdotQdot_isUpToDate(8, bool(false))
    , gravity_dtwist(Eigen::Matrix<double,6,1>::Zero())
    , segJoints(8)
    , Jroot(6,7)
    , dJroot(6,7)
    , segNonLinearEffects(8, Eigen::Matrix<double,6,6>::Zero())
    , segNonLinearEffects_isUpToDate(8)
    {
        //================== Register all Joint for the kinematic structure ==================//
        segJoints[0].push_back( new RevoluteJoint(0, 1, 0, 0, Eigen::Displacementd(0.0,0.0,0.0,1.0,0.0,0.0,0.0), Eigen::Vector3d(0,0,1)) );

        segJoints[1].push_back( new RevoluteJoint(1, 2, 1, 1, Eigen::Displacementd(0.0,0.0,0.3105,0.707106781187,0.707106781187,0.0,0.0), Eigen::Vector3d(0,0,1)) );

        segJoints[2].push_back( new RevoluteJoint(2, 3, 2, 2, Eigen::Displacementd(0.0,0.2,0.0,-0.707106781187,0.707106781187,0.0,0.0), Eigen::Vector3d(0,0,1)) );

        segJoints[3].push_back( new RevoluteJoint(3, 4, 3, 3, Eigen::Displacementd(0.0,0.0,0.2,-0.707106781187,0.707106781187,0.0,0.0), Eigen::Vector3d(0,0,1)) );

        segJoints[4].push_back( new RevoluteJoint(4, 5, 4, 4, Eigen::Displacementd(0.0,-0.195, 0.0, 0.707106781187,0.707106781187,0.0,0.0), Eigen::Vector3d(0,0,1)) );

        segJoints[5].push_back( new RevoluteJoint(5, 6, 5, 5, Eigen::Displacementd(0.0,0.06,0.195,0.707106781187,0.707106781187,0.0,0.0), Eigen::Vector3d(0,0,1)) );

        segJoints[6].push_back( new RevoluteJoint(6, 7, 6, 6, Eigen::Displacementd(0.0,0.078,0.06,-0.707106781187,0.707106781187,0.0,0.0), Eigen::Vector3d(0,0,1)) );

//        segJoints[0].push_back( new RevoluteJoint(0, 1, 0, 0, Eigen::Displacementd(0.0,0.0,0.0,1.0,0.0,0.0,0.0), Eigen::Vector3d(0,0,1)) );

//        segJoints[1].push_back( new RevoluteJoint(1, 2, 1, 1, Eigen::Displacementd(0.0,0.0,0.3105,0.707106781187,0.707106781187,0.0,0.0), Eigen::Vector3d(0,0,1)) );

//        segJoints[2].push_back( new RevoluteJoint(2, 3, 2, 2, Eigen::Displacementd(0.0,0.2,0.0,-0.707106781187,0.707106781187,0.0,0.0), Eigen::Vector3d(0,0,1)) );

//        segJoints[3].push_back( new RevoluteJoint(3, 4, 3, 3, Eigen::Displacementd(0.0,0.0,0.2,-0.707106781187,0.707106781187,0.0,0.0), Eigen::Vector3d(0,0,1)) );

//        segJoints[4].push_back( new RevoluteJoint(4, 5, 4, 4, Eigen::Displacementd(0.0,-0.195,0.707106781187,0.707106781187,0.0,0.0,0.0), Eigen::Vector3d(0,0,1)) );

//        segJoints[5].push_back( new RevoluteJoint(5, 6, 5, 5, Eigen::Displacementd(0.0,0.06,0.195,0.707106781187,0.707106781187,0.0,0.0), Eigen::Vector3d(0,0,1)) );

 //       segJoints[6].push_back( new RevoluteJoint(6, 7, 6, 6, Eigen::Displacementd(0.0,0.078,0.06,-0.707106781187,0.707106781187,0.0,0.0), Eigen::Vector3d(0,0,1)) );

        Jroot.setZero();
        dJroot.setZero();
        //================== Register all constant data other than kinematic structure ==================//
        segIndexFromName["kuka.00"] = 0;
        segNameFromIndex[0]     = "kuka.00";
        segMass[0]             = 1.5;
        segCoM[0]              = Eigen::Vector3d(0.00997877068205, -0.000481698546007, -0.254855093128);
        segMomentsOfInertia[0] = Eigen::Vector3d(0.00805112828804, 0.00845152236845, 0.00863585824764);
        segInertiaAxes[0]      = Eigen::Rotation3d(0.951651679011, 0.258168415469, -0.0562893581257, 0.156651394028);
        segMassMatrix[0]      << 0.105521428227, -9.20984846959e-05, 0.00373943706408, 5.17283020888e-18, 0.382282639692, -0.000722547819011, -9.20984846959e-05, 0.106046330919, -0.000285129190607, -0.382282639692, 1.52465930506e-20, -0.0149681560231, 0.00373943706408, -0.000285129190607, 0.00872352893124, 0.000722547819011, 0.0149681560231, 2.76979773752e-19, 1.09419716126e-17, -0.382282639692, 0.000722547819011, 1.5, 1.19211416997e-17, -4.36493018379e-17, 0.382282639692, 1.42095070858e-17, 0.0149681560231, -9.21317736729e-17, 1.5, 1.14545959523e-16, -0.000722547819011, -0.0149681560231, 3.78623727423e-19, -7.88892605755e-17, -5.30174862345e-17, 1.5;

        segIndexFromName["kuka.01"] = 1;
        segNameFromIndex[1]     = "kuka.01";
        segMass[1]             = 1.5;
        segCoM[1]              = Eigen::Vector3d(2.18269988807e-05, -0.0216238404123, -0.0791517192429);
        segMomentsOfInertia[1] = Eigen::Vector3d(0.0129545076394, 0.0128178432409, 0.00442072765916);
        segInertiaAxes[1]      = Eigen::Rotation3d(0.992539983081, 0.119380627464, 0.0056097965287, 0.0241076327054);
        segMassMatrix[1]      << 0.0230506546446, 4.09985248127e-05, -0.00013593386262, 9.0844283593e-20, 0.118727578864, -0.0324357606185, 4.09985248127e-05, 0.0217451642023, -0.000636485219866, -0.118727578864, 7.38824488243e-19, -3.27404983211e-05, -0.00013593386262, -0.000636485219866, 0.00559501652169, 0.0324357606185, 3.27404983211e-05, 7.88270036539e-20, -2.22769665128e-19, -0.118727578864, 0.0324357606185, 1.5, 1.32857117777e-18, 5.58702932009e-18, 0.118727578864, 1.29278403575e-18, 3.27404983211e-05, -1.4420735927e-17, 1.5, -6.32089866559e-17, -0.0324357606185, -3.27404983211e-05, 7.82446685026e-20, 7.28956554407e-18, -4.93854089567e-17, 1.5;

        segIndexFromName["kuka.02"] = 2;
        segNameFromIndex[2]     = "kuka.02";
        segMass[2]             = 1.5;
        segCoM[2]              = Eigen::Vector3d(-0.00100970535425, 0.0284784646483, 0.082033627936);
        segMomentsOfInertia[2] = Eigen::Vector3d(0.0132950018592, 0.0131003199685, 0.0044657387969);
        segInertiaAxes[2]      = Eigen::Rotation3d(0.991567954729, 0.126070202997, -0.00469624429149, 0.0296182437199);
        segMassMatrix[2]      << 0.0246050813717, 5.03285562013e-05, 0.000142932317751, 8.19080859995e-19, -0.123050441904, 0.0427176969725, 5.03285562013e-05, 0.0226558432223, -0.00141172174308, 0.123050441904, 9.66067546123e-19, 0.00151455803138, 0.000142932317751, -0.00141172174308, 0.00622481172849, -0.0427176969725, -0.00151455803138, -1.14534736337e-19, -8.55595920955e-19, 0.123050441904, -0.0427176969725, 1.5, -8.00853769761e-18, -8.9785492409e-19, -0.123050441904, 5.40062913213e-19, -0.00151455803138, -2.12125107959e-17, 1.5, 1.26716128909e-16, 0.0427176969725, 0.00151455803138, -3.83811804225e-21, -3.27462937409e-18, 7.2072339416e-17, 1.5;

        segIndexFromName["kuka.03"] = 3;
        segNameFromIndex[3]     = "kuka.03";
        segMass[3]             = 1.5;
        segCoM[3]              = Eigen::Vector3d(-0.000685124725844, 0.0227876276433, -0.0788064754543);
        segMomentsOfInertia[3] = Eigen::Vector3d(0.0127763513867, 0.0128011372133, 0.0043627376553);
        segInertiaAxes[3]      = Eigen::Rotation3d(0.994372066457, -0.105471539353, -0.00889191872746, -0.00456964088122);
        segMassMatrix[3]      << 0.0228686071417, 5.32028853369e-05, 5.64766769911e-05, 6.36862897217e-20, 0.118209713181, 0.0341814414649, 5.32028853369e-05, 0.0217459743213, 0.000962707553649, -0.118209713181, 1.74700545371e-21, 0.00102768708877, 5.64766769911e-05, 0.000962707553649, 0.00551626262133, -0.0341814414649, -0.00102768708877, 2.41669087763e-20, -1.94817577868e-20, -0.118209713181, -0.0341814414649, 1.5, -1.43275623028e-18, -4.57736604696e-18, 0.118209713181, 5.75982404133e-20, -0.00102768708877, -1.82239138602e-18, 1.5, 1.17364885172e-17, 0.0341814414649, 0.00102768708877, 6.05099161695e-20, -5.00935285006e-18, -1.81603863891e-18, 1.5;

        segIndexFromName["kuka.04"] = 4;
        segNameFromIndex[4]     = "kuka.04";
        segMass[4]             = 1.5;
        segCoM[4]              = Eigen::Vector3d(-0.000251206204044, 0.0270405432089, 0.0825998702666);
        segMomentsOfInertia[4] = Eigen::Vector3d(0.0132276730313, 0.0130734228464, 0.00439977729281);
        segInertiaAxes[4]      = Eigen::Rotation3d(0.991343563443, 0.129173854529, 0.00191487774184, 0.0234176807098);
        segMassMatrix[4]      << 0.0245573859768, 3.92848301246e-05, -5.1082782682e-05, -1.17525821432e-20, -0.1238998054, 0.0405608148133, 3.92848301246e-05, 0.0227394327974, -0.00120378365379, 0.1238998054, 4.18990141315e-19, 0.000376809306065, -5.1082782682e-05, -0.00120378365379, 0.00606603234525, -0.0405608148133, -0.000376809306065, 1.61995051162e-20, 6.10181359378e-19, 0.1238998054, -0.0405608148133, 1.5, -2.30975296804e-17, 1.79316874934e-18, -0.1238998054, -7.6757066889e-19, -0.000376809306065, -2.83868269196e-17, 1.5, -1.07498645402e-16, 0.0405608148133, 0.000376809306065, 5.61688723148e-20, -2.53262851229e-19, -8.06646416329e-17, 1.5;

        segIndexFromName["kuka.05"] = 5;
        segNameFromIndex[5]     = "kuka.05";
        segMass[5]             = 1.5;
        segCoM[5]              = Eigen::Vector3d(-7.2951660344e-05, 0.03430055852, -0.0828364494194);
        segMomentsOfInertia[5] = Eigen::Vector3d(0.0124407882151, 0.0124710412724, 0.00438005667842);
        segInertiaAxes[5]      = Eigen::Rotation3d(0.980059747892, -0.198267427159, -0.0127957403804, 0.0030309930321);
        segMassMatrix[5]      << 0.0244928283852, 8.60478297e-05, 0.000186080559021, 4.87255702908e-19, 0.124254674129, 0.05145083778, 8.60478297e-05, 0.0215423549773, 0.00136643367477, -0.124254674129, 1.07996700775e-20, 0.000109427490516, 0.000186080559021, 0.00136643367477, 0.00737193577082, -0.05145083778, -0.000109427490516, 1.57548128189e-19, -1.58395161137e-19, -0.124254674129, -0.05145083778, 1.5, -2.96800344718e-18, -1.50195882207e-17, 0.124254674129, -2.02229116157e-20, -0.000109427490516, -2.28614192464e-18, 1.5, 2.54245409448e-17, 0.05145083778, 0.000109427490516, 2.74544554029e-19, -1.54193877718e-17, 6.85215773011e-17, 1.5;

        segIndexFromName["kuka.06"] = 6;
        segNameFromIndex[6]     = "kuka.06";
        segMass[6]             = 0.75;
        segCoM[6]              = Eigen::Vector3d(0.00332657851002, -0.000250606397711, 0.000625058111516);
        segMomentsOfInertia[6] = Eigen::Vector3d(0.00284169077177, 0.00298256896714, 0.00296892557443);
        segInertiaAxes[6]      = Eigen::Rotation3d(0.97793407857, 0.00668809838671, -0.206811726513, -0.0287943933321);
        segMassMatrix[6]      << 0.00286329239677, 8.21664995199e-06, -4.84746567004e-05, 8.11960489282e-21, -0.000468793583637, -0.000187954798283, 8.21664995199e-06, 0.00299066976421, 3.4957933272e-06, 0.000468793583637, 2.65442258568e-21, -0.00249493388252, -4.84746567004e-05, 3.4957933272e-06, 0.00295650259105, 0.000187954798283, 0.00249493388252, -6.90861247604e-21, 7.54388718648e-22, 0.000468793583637, 0.000187954798283, 0.75, -6.10176065422e-18, 8.34835672814e-18, -0.000468793583637, -7.91198255907e-22, 0.00249493388252, -5.70982262766e-18, 0.75, -5.22995199326e-18, -0.000187954798283, -0.00249493388252, -8.24533634593e-21, -5.61074624261e-18, -6.88256621291e-18, 0.75;

        segIndexFromName["kuka.07"] = 7;
        segNameFromIndex[7]     = "kuka.07";
        segMass[7]             = 0.075;
        segCoM[7]              = Eigen::Vector3d(-1.05088630008e-05, 0.000493914659019, 0.0625078296156);
        segMomentsOfInertia[7] = Eigen::Vector3d(5.81502716972e-05, 5.84472283068e-05, 9.68983240391e-05);
        segInertiaAxes[7]      = Eigen::Rotation3d(0.952020561258, -0.0024010296426, -0.00517544935602, -0.305980882937);
        segMassMatrix[7]      << 0.000351314252074, 1.38490605053e-07, -2.75855087e-07, 1.55811291645e-19, -0.00468808722117, 3.70435994264e-05, 1.38490605053e-07, 0.000351390911357, -2.01600865809e-06, 0.00468808722117, -2.21622753359e-19, 7.88164725062e-07, -2.75855087e-07, -2.01600865809e-06, 9.69115844182e-05, -3.70435994264e-05, -7.88164725062e-07, -3.34491160201e-23, 2.4051458637e-19, 0.00468808722117, -3.70435994264e-05, 0.075, -3.64247824685e-18, 4.95514274144e-20, -0.00468808722117, -1.92177243203e-19, -7.88164725062e-07, -4.78070854822e-18, 0.075, -6.68097237147e-20, 3.70435994264e-05, 7.88164725062e-07, -1.19165706996e-23, 8.77737891592e-20, -4.51574440005e-20, 0.075;

        total_mass        = 9.825;
        actuatedDofs      = Eigen::VectorXd::Constant(nbDofs, 1);
        lowerLimits      << -10.0, -10.0, -10.0, -10.0, -10.0, -10.0, -10.0;
        upperLimits      << 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0;
        gravity_dtwist   << 0, 0 ,0, 0, 0, -9.80665;
        B.setZero();
        B.diagonal()      = Eigen::VectorXd::Constant(nbDofs, 0.001);
        outdateModel();
        l.setZero();

    };

    ~Pimpl()
    {
        for (int i=0; i<nbSegments(); ++i)
        {
            while (!segJoints[i].empty())
            {
                delete segJoints[i].back();
                segJoints[i].pop_back();
            }
        }
    }

    int nbSegments()
    {
        return nbSeg;
    }

    const Eigen::VectorXd& getActuatedDofs()
    {
        return actuatedDofs;
    }

    const Eigen::VectorXd& getJointLowerLimits()
    {
        return lowerLimits;
    }

    const Eigen::VectorXd& getJointUpperLimits()
    {
        return upperLimits;
    }

    const Eigen::VectorXd& getJointPositions()
    {
        return q;
    }

    const Eigen::VectorXd& getJointVelocities()
    {
        return dq;
    }

    const Eigen::Displacementd& getFreeFlyerPosition()
    {
        return Hroot;
    }

    const Eigen::Twistd& getFreeFlyerVelocity()
    {
        return Troot;
    }

    const Eigen::MatrixXd& getInertiaMatrix()
    {
        if (!M_isUpToDate)
        {
            updateModel();
            M.setZero();
            for (int i=0; i<8; i++)
            {
                M.noalias() += getSegmentJacobian(i).transpose() * getSegmentMassMatrix(i) * getSegmentJacobian(i);
            }

            M_isUpToDate = true;
        }
        return M;
    }

    const Eigen::MatrixXd& getInertiaMatrixInverse()
    {
        if (!Minv_isUpToDate)
        {
            updateModel();
//            Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es( getInertiaMatrix() );
//            const Eigen::VectorXd invdiag = es.eigenvalues().array().inverse();
//            const Eigen::MatrixXd& R      = es.eigenvectors();
//            Minv = R * invdiag.asDiagonal() * R.transpose();

            Minv = getInertiaMatrix().inverse();

            Minv_isUpToDate = true;
        }
        return Minv;
    }

    const Eigen::MatrixXd& getDampingMatrix()
    {
        return B;
    }

    const Eigen::VectorXd& getNonLinearTerms()
    {
        if (!n_isUpToDate)
        {
            updateModel();
            n.setZero();
            for (int i=0; i<8; i++)
            {
                n.noalias() += getSegmentJacobian(i).transpose() * ( getSegmentMassMatrix(i) * getSegmentJdotQdot(i) + getSegmentNonLinearEffectsMatrix(i) * getSegmentVelocity(i) ) ;
            }

            n_isUpToDate = true;
        }
        return n;
    }

    const Eigen::VectorXd& getLinearTerms()
    {
        return l;
    }

    const Eigen::VectorXd& getGravityTerms()
    {
        if (!g_isUpToDate)
        {
            updateModel();
            g.setZero();
            for (int i=0; i<8; i++)
            {
                g -= ( getSegmentJacobian(i).transpose() * getSegmentMassMatrix(i) * getSegmentPosition(i).inverse().adjoint() ) * gravity_dtwist;
            }

            g_isUpToDate = true;
        }
        return g;
    }

    double getMass()
    {
        return total_mass;
    }

    const Eigen::Vector3d& getCoMPosition()
    {
        if (!comPosition_isUpToDate)
        {
            updateModel();
            comPosition.setZero();
            for (int i=0; i<8; i++)
            {
                comPosition.noalias() += (getSegmentPosition(i).getTranslation() + getSegmentPosition(i).getRotation() * getSegmentCoM(i)) * (getSegmentMass(i)/total_mass);
            }

            comPosition_isUpToDate = true;
        }
        return comPosition;
    }

    const Eigen::Vector3d& getCoMVelocity()
    {
        if (!comVelocity_isUpToDate)
        {
            updateModel();
            comVelocity = getCoMJacobian() * alldq;

            comVelocity_isUpToDate = true;
        }
        return comVelocity;
    }

    const Eigen::Vector3d& getCoMJdotQdot()
    {
        if (!comJdotQdot_isUpToDate)
        {
            updateModel();
            comJdotQdot = getCoMJacobianDot() * alldq;

            comJdotQdot_isUpToDate = true;
        }
        return comJdotQdot;
    }

    const Eigen::Matrix<double,3,Eigen::Dynamic>& getCoMJacobian()
    {
        if (!comJacobian_isUpToDate)
        {
            updateModel();
            comJacobian.setZero();
            for (int i=0; i<8; i++)
            {
                Eigen::Matrix<double,6,7> tmpJ = Eigen::Displacementd(getSegmentCoM(i), getSegmentPosition(i).getRotation().inverse()).inverse().adjoint() * getSegmentJacobian(i);
                comJacobian.noalias()              += ( tmpJ.bottomRows<3>() ) * (getSegmentMass(i)/ total_mass);
            }

            comJacobian_isUpToDate = true;
        }
        return comJacobian;
    }

    const Eigen::Matrix<double,3,Eigen::Dynamic>& getCoMJacobianDot()
    {
        if (!comJacobianDot_isUpToDate)
        {
            updateModel();
            comJacobianDot.setZero();
            for (int i=0; i<8; i++)
            {
                Eigen::Twistd T_comR_seg_comR       = getSegmentVelocity(i);
                T_comR_seg_comR.getLinearVelocity().setZero();
                AdjointMatrix Ad_comR_seg           = Eigen::Displacementd(getSegmentCoM(i), getSegmentPosition(i).getRotation().inverse()).inverse().adjoint();
                AdjointMatrix dAd_comR_seg          = Ad_comR_seg * getAdjacencyMatrix( T_comR_seg_comR );
                Eigen::Matrix<double,6,7> tmpJ = Ad_comR_seg * getSegmentJdot(i)  +  dAd_comR_seg * getSegmentJacobian(i);
                comJacobianDot.noalias()           += ( tmpJ.bottomRows<3>() ) * (getSegmentMass(i)/ total_mass);
            }

            comJacobianDot_isUpToDate = true;
        }
        return comJacobianDot;
    }

    double getSegmentMass(int index)
    {
        return segMass[index];
    }

    const Eigen::Vector3d& getSegmentCoM(int index)
    {
        return segCoM[index];
    }

    const Eigen::Matrix<double,6,6>& getSegmentMassMatrix(int index)
    {
        return segMassMatrix[index];
    }

    const Eigen::Vector3d& getSegmentMomentsOfInertia(int index)
    {
        return segMomentsOfInertia[index];
    }

    const Eigen::Rotation3d& getSegmentInertiaAxes(int index)
    {
        return segInertiaAxes[index];
    }

    const Eigen::Displacementd& getSegmentPosition(int index)
    {
        if (!segPosition_isUpToDate[index])
        {
            updateModel();


            segPosition_isUpToDate[index] = true;
        }
        return segPosition[index];
    }

    const Eigen::Twistd& getSegmentVelocity(int index)
    {
        if (!segVelocity_isUpToDate[index])
        {
            updateModel();


            segVelocity_isUpToDate[index] = true;
        }
        return segVelocity[index];
    }

    const Eigen::Matrix<double,6,Eigen::Dynamic>& getSegmentJacobian(int index)
    {
        if (!segJacobian_isUpToDate[index])
        {
            updateModel();


            segJacobian_isUpToDate[index] = true;
        }
        return segJacobian[index];
    }

    const Eigen::Matrix<double,6,Eigen::Dynamic>& getSegmentJdot(int index)
    {
        if (!segJdot_isUpToDate[index])
        {
            updateModel();


            segJdot_isUpToDate[index] = true;
        }
        return segJdot[index];
    }

    const Eigen::Matrix<double,6,Eigen::Dynamic>& getJointJacobian(int index)
    {
        segJointJacobian[index] = getSegmentJacobian(index);
        return segJointJacobian[index];
    }

    const Eigen::Twistd& getSegmentJdotQdot(int index)
    {
        if (!segJdotQdot_isUpToDate[index])
        {
            updateModel();


            segJdotQdot_isUpToDate[index] = true;
        }
        return segJdotQdot[index];
    }

    void doSetJointPositions(const Eigen::VectorXd& _q)
    {
        q = _q;
        outdateModel();
    }

    void doSetJointVelocities(const Eigen::VectorXd& _dq)
    {
        dq    = _dq;
        alldq =  dq;
        outdateModel();
    }

    void doSetFreeFlyerPosition(const Eigen::Displacementd& _Hroot)
    {
        Hroot = _Hroot;
        outdateModel();
    }

    void doSetFreeFlyerVelocity(const Eigen::Twistd& _Troot)
    {
        Troot = _Troot;
        outdateModel();
    }

    int doGetSegmentIndex(const std::string& name)
    {
        try
        {
            return segIndexFromName.at(name);
        }
        catch(std::out_of_range)
        {
            std::stringstream ss;
            ss << "[kukafixed::doGetSegmentIndex]: The input segment name '"+name+"' does not exist in this robot; possible key:\n";
            for (std::map< std::string, int >::iterator it = segIndexFromName.begin(); it != segIndexFromName.end(); ++it)
                ss << it->first <<"\n";
            throw std::out_of_range(ss.str());
        }

    }

    const std::string& doGetSegmentName(int index)
    {
        return segNameFromIndex.at(index);
    }

    void outdateModel()
    {
        isUpToDate                     = false;
        M_isUpToDate                   = false;
        Minv_isUpToDate                = false;
        n_isUpToDate                   = false;
        g_isUpToDate                   = false;
        comPosition_isUpToDate         = false;
        comVelocity_isUpToDate         = false;
        comJdotQdot_isUpToDate         = false;
        comJacobian_isUpToDate         = false;
        comJacobianDot_isUpToDate      = false;
        for (int i=0; i<8; i++)
        {
            segPosition_isUpToDate[i]  = false;
            segVelocity_isUpToDate[i]  = false;
            segJacobian_isUpToDate[i]  = false;
            segJdot_isUpToDate[i]      = false;
            segJdotQdot_isUpToDate[i]  = false;
        }
    }

    void updateModel()
    {
        if (!isUpToDate)
        {
            recursiveUpdateModel(0, Hroot, Troot, Jroot, dJroot); //0 because this index must be the root link of the robot
            isUpToDate = true;
        }
    }

    void recursiveUpdateModel(int segIdx, const Eigen::Displacementd& H_0_seg, const Eigen::Twistd& T_s_0_s, const Eigen::Matrix<double,6,Eigen::Dynamic>& J_s_0_s, const Eigen::Matrix<double,6,Eigen::Dynamic>& dJ_s_0_s)
    {

        segPosition[segIdx] = H_0_seg;
        segVelocity[segIdx] = T_s_0_s;
        segJacobian[segIdx] = J_s_0_s;
        segJdot[segIdx]     = dJ_s_0_s;

        segJdotQdot[segIdx] = dJ_s_0_s * alldq;

        Eigen::Matrix3d wx                  = getSkewSymmetricMatrix( T_s_0_s.getAngularVelocity() );
        Eigen::Matrix3d rx                  = getSkewSymmetricMatrix( segCoM[segIdx] );
        Eigen::Matrix<double, 6, 6 > _segN  = Eigen::Matrix<double, 6, 6 >::Zero();
        _segN.topLeftCorner<3,3>()          = wx;
        _segN.bottomRightCorner<3,3>()      = wx;
        _segN.topRightCorner<3,3>()         = rx * wx - wx * rx;

        segNonLinearEffects[segIdx]         = _segN * segMassMatrix[segIdx];

        // register update status
        segPosition_isUpToDate[segIdx]         = true;
        segVelocity_isUpToDate[segIdx]         = true;
        segJacobian_isUpToDate[segIdx]         = true;
        segJdot_isUpToDate[segIdx]             = true;
        segJdotQdot_isUpToDate[segIdx]         = true;
        segNonLinearEffects_isUpToDate[segIdx] = true;


        // -----> recursive computation of children
        for (int jidx=0; jidx < segJoints[segIdx].size(); jidx++)
        {
            int cdof         = segJoints[segIdx][jidx]->dof;
            int cinternaldof = segJoints[segIdx][jidx]->internaldof;

            double q_dof  =  q[cinternaldof];
            double dq_dof = dq[cinternaldof];


            // DATA //
            Eigen::Displacementd  H_seg_child   =  segJoints[segIdx][jidx]->H_seg_joint * segJoints[segIdx][jidx]->get_H_joint_child( q_dof ) ;
            AdjointMatrix         Ad_child_seg  =  H_seg_child.inverse().adjoint();
            AdjointMatrix         Ad_joint_seg  =  segJoints[segIdx][jidx]->H_seg_joint.inverse().adjoint();
            Eigen::Twistd         T_c_j_c       =  segJoints[segIdx][jidx]->get_twist( dq_dof );
            AdjointMatrix         dAd_child_seg =  segJoints[segIdx][jidx]->getdAdjointInverse( q_dof, dq_dof ) * Ad_joint_seg;


            //                    H_0_child    = H_0_parent * H_parent_child
            //                                 = H_0_parent * H_parent_joint * H_joint_child(q_i)
            Eigen::Displacementd  H_0_child    = H_0_seg * H_seg_child;

            // T_child/0(child)    =  T_child/parent(child) + T_parent/0(child)
            //                     =  T_child/parent(child) + Ad_child_parent * T_parent/0(parent)
            //                     =  T_child/joint(child) + T_joint/parent(child)[=0 because same link] + Ad_child_parent * T_parent/0(parent)
            //                     =  T_child/joint(child)(dq_i) + Ad_child_parent * T_parent/0(parent)
            Eigen::Twistd  T_c_0_c =  T_c_j_c  +  Ad_child_seg * T_s_0_s;

            // J_child/0(child)  =  (... same reasoning as for T_child/0(child) ...)
            //                   =  [ 0 0 0 ... J_child/joint(child) ... 0 0 ]  +  Ad_child_parent * J_parent/0(parent)
            Eigen::Matrix<double,6,Eigen::Dynamic> J_c_0_c  = Eigen::Matrix<double,6,Eigen::Dynamic>::Zero(6, nbDofs);
            J_c_0_c.col(cdof)    = segJoints[segIdx][jidx]->jacobian;
            J_c_0_c.noalias()   += Ad_child_seg * J_s_0_s;

            // dJ_child/0(child)  =  (... same reasoning as for T_child/0(child) ...)
            //                    =  [ 0 0 0 ... dJ_child/joint(child) ... 0 0 ]  +  dAd_child_parent * J_parent/0(parent) + Ad_child_parent * dJ_parent/0(parent)
            Eigen::Matrix<double,6,Eigen::Dynamic> dJ_c_0_c = Eigen::Matrix<double,6,Eigen::Dynamic>::Zero(6, nbDofs);
            dJ_c_0_c.col(cdof)    = segJoints[segIdx][jidx]->djacobian;
            dJ_c_0_c.noalias()   += dAd_child_seg * J_s_0_s + Ad_child_seg * dJ_s_0_s;


            recursiveUpdateModel(segJoints[segIdx][jidx]->child, H_0_child, T_c_0_c, J_c_0_c, dJ_c_0_c);
        }

    }

    const Eigen::Matrix<double,6,6>& getSegmentNonLinearEffectsMatrix(int index)
    {
        if (!segNonLinearEffects_isUpToDate[index])
        {
            updateModel();


            segNonLinearEffects_isUpToDate[index] = true;
        }
        return segNonLinearEffects[index];
    }


};

kukafixed::kukafixed(const std::string& robotName)
: orc::Model(robotName, 7, false)
, pimpl( new Pimpl() )
{

};

kukafixed::~kukafixed()
{

}

int kukafixed::nbSegments() const
{
    return pimpl->nbSegments();
}

const Eigen::VectorXd& kukafixed::getActuatedDofs() const
{
    return pimpl->getActuatedDofs();
}

const Eigen::VectorXd& kukafixed::getJointLowerLimits() const
{
    return pimpl->getJointLowerLimits();
}

const Eigen::VectorXd& kukafixed::getJointUpperLimits() const
{
    return pimpl->getJointUpperLimits();
}

const Eigen::VectorXd& kukafixed::getJointPositions() const
{
    return pimpl->getJointPositions();
}

const Eigen::VectorXd& kukafixed::getJointVelocities() const
{
    return pimpl->getJointVelocities();
}

const Eigen::Displacementd& kukafixed::getFreeFlyerPosition() const
{
    return pimpl->getFreeFlyerPosition();
}

const Eigen::Twistd& kukafixed::getFreeFlyerVelocity() const
{
    return pimpl->getFreeFlyerVelocity();
}

const Eigen::MatrixXd& kukafixed::getInertiaMatrix() const
{
    return pimpl->getInertiaMatrix();
}

const Eigen::MatrixXd& kukafixed::getInertiaMatrixInverse() const
{
    return pimpl->getInertiaMatrixInverse();
}

const Eigen::MatrixXd& kukafixed::getDampingMatrix() const
{
    return pimpl->getDampingMatrix();
}

const Eigen::VectorXd& kukafixed::getNonLinearTerms() const
{
    return pimpl->getNonLinearTerms();
}

const Eigen::VectorXd& kukafixed::getLinearTerms() const
{
    return pimpl->getLinearTerms();
}

const Eigen::VectorXd& kukafixed::getGravityTerms() const
{
    return pimpl->getGravityTerms();
}

double kukafixed::getMass() const
{
    return pimpl->getMass();
}

const Eigen::Vector3d& kukafixed::getCoMPosition() const
{
    return pimpl->getCoMPosition();
}

const Eigen::Vector3d& kukafixed::getCoMVelocity() const
{
    return pimpl->getCoMVelocity();
}

const Eigen::Vector3d& kukafixed::getCoMJdotQdot() const
{
    return pimpl->getCoMJdotQdot();
}

const Eigen::Matrix<double,3,Eigen::Dynamic>& kukafixed::getCoMJacobian() const
{
    return pimpl->getCoMJacobian();
}

const Eigen::Matrix<double,3,Eigen::Dynamic>& kukafixed::getCoMJacobianDot() const
{
    return pimpl->getCoMJacobianDot();
}

double kukafixed::getSegmentMass(int index) const
{
    return pimpl->getSegmentMass(index);
}

const Eigen::Vector3d& kukafixed::getSegmentCoM(int index) const
{
    return pimpl->getSegmentCoM(index);
}

const Eigen::Matrix<double,6,6>& kukafixed::getSegmentMassMatrix(int index) const
{
    return pimpl->getSegmentMassMatrix(index);
}

const Eigen::Vector3d& kukafixed::getSegmentMomentsOfInertia(int index) const
{
    return pimpl->getSegmentMomentsOfInertia(index);
}

const Eigen::Rotation3d& kukafixed::getSegmentInertiaAxes(int index) const
{
    return pimpl->getSegmentInertiaAxes(index);
}

const Eigen::Displacementd& kukafixed::getSegmentPosition(int index) const
{
    return pimpl->getSegmentPosition(index);
}

const Eigen::Twistd& kukafixed::getSegmentVelocity(int index) const
{
    return pimpl->getSegmentVelocity(index);
}

const Eigen::Matrix<double,6,Eigen::Dynamic>& kukafixed::getSegmentJacobian(int index) const
{
    return pimpl->getSegmentJacobian(index);
}

const Eigen::Matrix<double,6,Eigen::Dynamic>& kukafixed::getSegmentJdot(int index) const
{
    return pimpl->getSegmentJdot(index);
}

const Eigen::Matrix<double,6,Eigen::Dynamic>& kukafixed::getJointJacobian(int index) const
{
    return pimpl->getJointJacobian(index);
}

const Eigen::Twistd& kukafixed::getSegmentJdotQdot(int index) const
{
    return pimpl->getSegmentJdotQdot(index);
}

void kukafixed::doSetJointPositions(const Eigen::VectorXd& q)
{
    pimpl->doSetJointPositions(q);
}

void kukafixed::doSetJointVelocities(const Eigen::VectorXd& dq)
{
    pimpl->doSetJointVelocities(dq);
}

void kukafixed::doSetFreeFlyerPosition(const Eigen::Displacementd& Hroot)
{
    pimpl->doSetFreeFlyerPosition(Hroot);
}

void kukafixed::doSetFreeFlyerVelocity(const Eigen::Twistd& Troot)
{
    pimpl->doSetFreeFlyerVelocity(Troot);
}

int kukafixed::doGetSegmentIndex(const std::string& name) const
{
    return pimpl->doGetSegmentIndex(name);
}

const std::string& kukafixed::doGetSegmentName(int index) const
{
    return pimpl->doGetSegmentName(index);
}

void kukafixed::printAllData() const
{
    std::cout<<"nbSeg:\n";
    std::cout<<nbSegments()<<"\n";

    std::cout<<"actuatedDofs:\n";
    std::cout<<getActuatedDofs()<<"\n";

    std::cout<<"lowerLimits:\n";
    std::cout<<getJointLowerLimits()<<"\n";

    std::cout<<"upperLimits:\n";
    std::cout<<getJointUpperLimits()<<"\n";

    std::cout<<"q:\n";
    std::cout<<getJointPositions()<<"\n";

    std::cout<<"dq:\n";
    std::cout<<getJointVelocities()<<"\n";

    std::cout<<"Hroot:\n";
    std::cout<<getFreeFlyerPosition()<<"\n";

    std::cout<<"Troot:\n";
    std::cout<<getFreeFlyerVelocity()<<"\n";

    std::cout<<"total_mass:\n";
    std::cout<<getMass()<<"\n";

    std::cout<<"comPosition:\n";
    std::cout<<getCoMPosition()<<"\n";

    std::cout<<"comVelocity:\n";
    std::cout<<getCoMVelocity()<<"\n";

    std::cout<<"comJdotQdot:\n";
    std::cout<<getCoMJdotQdot()<<"\n";

    std::cout<<"comJacobian:\n";
    std::cout<<getCoMJacobian()<<"\n";

    std::cout<<"comJacobianDot:\n";
    std::cout<<getCoMJacobianDot()<<"\n";

    std::cout<<"M:\n";
    std::cout<<getInertiaMatrix()<<"\n";

    std::cout<<"Minv:\n";
    std::cout<<getInertiaMatrixInverse()<<"\n";

    std::cout<<"B:\n";
    std::cout<<getDampingMatrix()<<"\n";

    std::cout<<"n:\n";
    std::cout<<getNonLinearTerms()<<"\n";

    std::cout<<"l:\n";
    std::cout<<getLinearTerms()<<"\n";

    std::cout<<"g:\n";
    std::cout<<getGravityTerms()<<"\n";


    for (int idx=0; idx<nbSegments(); idx++)
    {
        std::cout<<"segMass "<<idx<<":\n";
        std::cout<<getSegmentMass(idx)<<"\n";

        std::cout<<"segCoM "<<idx<<":\n";
        std::cout<<getSegmentCoM(idx)<<"\n";

        std::cout<<"segMassMatrix "<<idx<<":\n";
        std::cout<<getSegmentMassMatrix(idx)<<"\n";

        std::cout<<"segMomentsOfInertia "<<idx<<":\n";
        std::cout<<getSegmentMomentsOfInertia(idx)<<"\n";

        std::cout<<"segInertiaAxes "<<idx<<":\n";
        std::cout<<getSegmentInertiaAxes(idx)<<"\n";

        std::cout<<"segPosition "<<idx<<":\n";
        std::cout<<getSegmentPosition(idx)<<"\n";

        std::cout<<"segVelocity "<<idx<<":\n";
        std::cout<<getSegmentVelocity(idx)<<"\n";

        std::cout<<"segJacobian "<<idx<<":\n";
        std::cout<<getSegmentJacobian(idx)<<"\n";

        std::cout<<"segJdot "<<idx<<":\n";
        std::cout<<getSegmentJdot(idx)<<"\n";

        std::cout<<"segJointJacobian "<<idx<<":\n";
        std::cout<<getJointJacobian(idx)<<"\n";

        std::cout<<"segJdotQdot "<<idx<<":\n";
        std::cout<<getSegmentJdotQdot(idx)<<"\n";

    }

}





