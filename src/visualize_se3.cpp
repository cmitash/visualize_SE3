#include <iostream>

#include <pcl/point_types.h>
#include <pcl/point_cloud.h>
#include <pcl/common/transforms.h>
#include <pcl/io/ply_io.h>
#include <pcl/features/moment_of_inertia_estimation.h>
#include <pcl/visualization/cloud_viewer.h>

#define POINT_SIZE 0.002
#define LINE_SIZE 0.001
#define GT_POINT_SIZE 0.002
enum modes {KUFFNER, SE3};

void toEulerianAngle(Eigen::Quaternionf q, Eigen::Vector3f& eulAngles){
	// roll (x-axis rotation)
	double sinr = +2.0 * (q.w() * q.x() + q.y() * q.z());
	double cosr = +1.0 - 2.0 * (q.x() * q.x() + q.y() * q.y());
	eulAngles[0] = atan2(sinr, cosr);

	// pitch (y-axis rotation)
	double sinp = +2.0 * (q.w() * q.y() - q.z() * q.x());
	if (fabs(sinp) >= 1)
		eulAngles[1] = copysign(M_PI / 2, sinp); // use 90 degrees if out of range
	else
		eulAngles[1] = asin(sinp);

	// yaw (z-axis rotation)
	double siny = +2.0 * (q.w() * q.z() + q.x() * q.y());
	double cosy = +1.0 - 2.0 * (q.y() * q.y() + q.z() * q.z());  
	eulAngles[2] = atan2(siny, cosy);
}

void get_rotation_4d(float qw, float qx, float qy, float qz, pcl::PointXYZ &point_start, pcl::PointXYZ &point_end) {
	float angle;
	Eigen::Vector3f axis;
	Eigen::Matrix3f rotm;
	Eigen::Vector3f tangent_line;
	Eigen::Vector3f end_point;
	Eigen::Vector3f eulAngles;

	// get axis angle representation
	angle = 2 * acos(qw);
	double s = sqrt(1-qw*qw);

	if (s < 0.001) { // test to avoid divide by zero, s is always positive due to sqrt
		// if s close to zero then direction of axis not important
		axis[0] = qx;
		axis[1] = qy;
		axis[2] = qz;
	} else {
		axis[0] = qx / s;
		axis[1] = qy / s;
		axis[2] = qz / s;
	}

	toEulerianAngle(Eigen::Quaternionf(qw, qx, qy, qz) ,eulAngles);
	rotm = Eigen::AngleAxisf(angle, axis);

	std::cout << "\n\n\neuler angle: " << eulAngles[0]*180/M_PI << " " << eulAngles[1]*180/M_PI << " " << eulAngles[2]*180/M_PI << "\n"; 
	std::cout << "axis: " << axis[0] << " " << axis[1] << " " << axis[2] << "\nangle: " << angle; 
	
	// ambiguity in quaternions/axis-angle
	// if(axis[2] < 0){
	// 	axis = -axis;
	// 	angle = -angle;
	// }

	tangent_line = axis.cross(Eigen::Vector3f(0,0,1));
	end_point = axis + LINE_SIZE*tangent_line;
	end_point = rotm*end_point;

	point_start.x = axis[0];
	point_start.y = axis[1];
	point_start.z = axis[2];

	point_end.x = end_point[0];
	point_end.y = end_point[1];
	point_end.z = end_point[2];
}

void get_custom_se3(float tx, float ty, float tz, float qw, float qx, float qy, float qz, 
	pcl::PointXYZ &point_start, pcl::PointXYZ &point_end_x, pcl::PointXYZ &point_end_y, Eigen::Vector3f &z_axis) {
	Eigen::Quaternionf q;
  	q.w() = qw;
	q.x() = qx;
	q.y() = qy;
	q.z() = qz;
	Eigen::Matrix3f rotMat;
	rotMat = q.toRotationMatrix();

	Eigen::Vector3f x_axis = rotMat.col(0);
	Eigen::Vector3f y_axis = rotMat.col(1);
	z_axis = rotMat.col(2);

	point_start.x = tx;
	point_start.y = ty;
	point_start.z = tz;

	point_end_x.x = point_start.x + LINE_SIZE*x_axis[0];
	point_end_x.y = point_start.y + LINE_SIZE*x_axis[1];
	point_end_x.z = point_start.z + LINE_SIZE*x_axis[2];

	point_end_y.x = point_start.x + LINE_SIZE*y_axis[0];
	point_end_y.y = point_start.y + LINE_SIZE*y_axis[1];
	point_end_y.z = point_start.z + LINE_SIZE*y_axis[2];
}

void custom_se3_vizualization(ifstream &pFile, std::vector<boost::shared_ptr<pcl::visualization::PCLVisualizer> > &viewers) {
	float tx, ty, tz, qw, qx, qy, qz;
  	pcl::PointXYZ point_start, point_end_x, point_end_y;
  	Eigen::Vector3f z_axis;
  	boost::shared_ptr<pcl::visualization::PCLVisualizer> viewer_se3 (new pcl::visualization::PCLVisualizer ("se3"));
	int num_poses = 0;

	viewer_se3->setBackgroundColor (0, 0, 0);
	viewer_se3->addCoordinateSystem (0.1);
	viewer_se3->initCameraParameters ();

	// get ground truth pose
	pFile >> tx >> ty >> tz >> qw >> qx >> qy >> qz;
  	get_custom_se3(tx, ty, tz, qw, qx, qy, qz, point_start, point_end_x, point_end_y, z_axis);
  	pcl::ModelCoefficients coeffs;
  	coeffs.values.clear ();
	coeffs.values.push_back (tx);
	coeffs.values.push_back (ty);
	coeffs.values.push_back (tz);
	coeffs.values.push_back (LINE_SIZE*z_axis[0]);
	coeffs.values.push_back (LINE_SIZE*z_axis[1]);
	coeffs.values.push_back (LINE_SIZE*z_axis[2]);
	coeffs.values.push_back (15.0);
  	viewer_se3->addCone(coeffs, "cone_gt");

	// viewer_se3->addLine (point_start, point_end_x, 0.0f, 0.0f, 1.0f, "gt_x");
 //  	viewer_se3->addLine (point_start, point_end_y, 0.0f, 0.0f, 0.0f, "gt_y");

  	while(pFile >> tx >> ty >> tz >> qw >> qx >> qy >> qz) {
  		get_custom_se3(tx, ty, tz, qw, qx, qy, qz, point_start, point_end_x, point_end_y, z_axis);
  		pcl::ModelCoefficients coeffs;
	  	coeffs.values.clear ();
		coeffs.values.push_back (tx);
		coeffs.values.push_back (ty);
		coeffs.values.push_back (tz);
		coeffs.values.push_back (LINE_SIZE*z_axis[0]);
		coeffs.values.push_back (LINE_SIZE*z_axis[1]);
		coeffs.values.push_back (LINE_SIZE*z_axis[2]);
		coeffs.values.push_back (15.0);
	  	viewer_se3->addCone(coeffs, std::to_string(num_poses));
  		// viewer_se3->addLine (point_start, point_end_x, 1.0f, 0.0f, 0.0f, std::to_string(num_poses) + "_x");
  		// viewer_se3->addLine (point_start, point_end_y, 0.0f, 1.0f, 0.0f, std::to_string(num_poses) + "_y");
  		num_poses++;
  	}

  	viewers.push_back(viewer_se3);

  	std::cout << "visualizing " << num_poses << " poses." << std::endl;	
}

// SO(3) is visualized as oriented arrows on a sphere
void kuffner_visualization(ifstream &pFile, std::vector<boost::shared_ptr<pcl::visualization::PCLVisualizer> > &viewers) {
	float tx, ty, tz, qw, qx, qy, qz;
  	pcl::PointXYZ point_start, point_end;
  	boost::shared_ptr<pcl::visualization::PCLVisualizer> viewer_translation (new pcl::visualization::PCLVisualizer ("translation"));
	boost::shared_ptr<pcl::visualization::PCLVisualizer> viewer_rotation (new pcl::visualization::PCLVisualizer ("rotation"));
	int num_poses = 0;

	viewer_translation->setBackgroundColor (0, 0, 0);
	viewer_translation->addCoordinateSystem (0.1);
	viewer_translation->initCameraParameters ();

	viewer_rotation->setBackgroundColor (0, 0, 0);
	viewer_rotation->addCoordinateSystem (0.1);
	viewer_rotation->initCameraParameters ();

	// get ground truth pose
	pFile >> tx >> ty >> tz >> qw >> qx >> qy >> qz;
  	get_rotation_4d(qw, qx, qy, qz, point_start, point_end);
  	viewer_translation->addSphere(pcl::PointXYZ(tx, ty,tz), GT_POINT_SIZE, 0.0f, 1.0f, 0.0f, "gt");
  	viewer_rotation->addArrow (point_end, point_start, 0.0f, 1.0f, 0.0f, false, "gt");

  	while(pFile >> tx >> ty >> tz >> qw >> qx >> qy >> qz) {
  		get_rotation_4d(qw, qx, qy, qz, point_start, point_end);

  		viewer_translation->addSphere(pcl::PointXYZ(tx, ty,tz), POINT_SIZE, std::to_string(num_poses));
  		viewer_rotation->addArrow (point_end, point_start, 1.0f, 0.0f, 0.0f, false, std::to_string(num_poses));
  		num_poses++;
  	}

  	viewers.push_back(viewer_translation);
	viewers.push_back(viewer_rotation);

  	std::cout << "visualizing " << num_poses << " poses." << std::endl;
}

int main(int argc, char* argv[]) {

	if(argc < 3){
		std::cout << "<mode of operation><filepath for a file with se3 configurations>\n";
		exit(-1);
	}

	ifstream pFile;
  	std::vector<boost::shared_ptr<pcl::visualization::PCLVisualizer> > viewers;
  	modes mode;
  	
  	
  	mode = modes(atoi(argv[1]));
  	pFile.open (argv[2], std::ofstream::in);
	if(!pFile) {
		cerr << "Can't read the input file !!!\n";
		exit(-1);
	}

	switch(mode) {
		case KUFFNER: 
			kuffner_visualization(pFile, viewers);
			break;
		case SE3:
			custom_se3_vizualization(pFile, viewers);
			break;
		default:
			break;
	}
	
	if(viewers.size() == 0) {
		std::cout << "no viewers found!!!\n";
		exit(-1);
	}

	while(!viewers[0]->wasStopped()) {
		for(auto viewer_it:viewers)
			viewer_it->spinOnce (100);
	    boost::this_thread::sleep (boost::posix_time::microseconds (100000));
	}

	pFile.close();
}           