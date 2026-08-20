#include "iLQRSolver.hpp"

void buildAcrobotModel(pinocchio::Model& model)
{
        double m1 = 1.0, m2 = 1.0;
        double l1 = 1.0, l2 = 1.0;
        double J1_val = 1.0, J2_val = 1.0;

        model.name = "Acrobot";

        /// Junta 1 ///
        pinocchio::JointIndex joint1_id = model.addJoint(
                0,                          // ID do nó pai, que é o universo aqui
                pinocchio::JointModelRZ(),  // O tipo de junta
                pinocchio::SE3::Identity(), // Na origem do mundo
                "shoulder"                  // Nome da junta
        );

        Eigen::Vector3d com1(l1, 0.0, 0.0); // Centro de massa do elo 1
        Eigen::Matrix3d I1 = Eigen::Matrix3d::Zero(); //Matriz de inércia do elo 1
        I1(2, 2) = J1_val;

        pinocchio::Inertia inertia1(m1, com1, I1);

        model.appendBodyToJoint(joint1_id, inertia1, pinocchio::SE3::Identity());

        /// Junta 2///
        pinocchio::SE3 joint2_placement(
                Eigen::Matrix3d::Identity(), // Sem rotação inicial relativa
                Eigen::Vector3d(l1, 0.0, 0.0) // Deslocado l1 no eixo X
        );


        pinocchio::JointIndex joint2_id = model.addJoint(
                joint1_id,
                pinocchio::JointModelRZ(),
                joint2_placement,
                "elbow"
        );

        Eigen::Vector3d com2(l2, 0.0, 0.0);
        Eigen::Matrix3d I2 = Eigen::Matrix3d::Zero();
        I2(2, 2) = J2_val;

        pinocchio::Inertia inertia2(m2, com2, I2);

        model.appendBodyToJoint(joint2_id, inertia2, pinocchio::SE3::Identity());


        model.gravity.linear(Eigen::Vector3d(0,-9.81,0));
}

int main()
{
	pinocchio::Model model;
	buildAcrobotModel(model);

	Eigen::Matrix<double, 2, 1> S;
	S << 1.0, 
	     0.0;

	double h = 0.1;
	double Tfinal = 10.0;
	int Nt = (int)(Tfinal/h) + 1;

	iLQRSolver ilqr(model, S, h, Nt);

	Eigen::Vector4d x0(-M_PI/2.0, 0.0, 0.0, 0.0);
        Eigen::Vector4d xgoal(M_PI/2.0, 0.0, 0.0, 0.0);

	Eigen::Matrix4d Q = Eigen::Vector4d(1.0, 1.0, 0.1, 0.1).asDiagonal();
	Eigen::Matrix<double, 1, 1> R(0.01);
	Eigen::Matrix4d Qn = 100.0 * Eigen::MatrixXd::Identity(4, 4);

	Eigen::Vector2d c(1.0, 1.0);

	ilqr.setCosts(Q, R, Qn);
	ilqr.setBoundaryConditions(x0, xgoal);
	ilqr.setViscousDamping(c);

	ilqr.solve();

	auto xtraj = ilqr.getXTraj();
	auto utraj = ilqr.getUTraj();

	std::cout << "Estado final alcancado:\n" << xtraj.back().transpose() << std::endl;

        std::ofstream arquivo("trajetoria.csv");
        arquivo << "theta1,theta2,omega1,omega2,torque\n";
    
         for (int i = 0; i < Nt - 1; ++i) {
                arquivo << xtraj[i][0] << "," << xtraj[i][1] << ","
                        << xtraj[i][2] << "," << xtraj[i][3] << ","
                        << utraj[i](0, 0) << "\n";
        }
    
        arquivo << xtraj.back()[0] << "," << xtraj.back()[1] << ","
                    << xtraj.back()[2] << "," << xtraj.back()[3] << ",0.0\n";
                
        arquivo.close();
        std::cout << "Trajetoria salva com sucesso em 'trajetoria.csv'!" << std::endl;

	return 0;
}

