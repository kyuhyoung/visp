
//! \example tutorial-hand-eye-calibration.cpp
#include <visp3/vision/vpHandEyeCalibration.h>
#include <vector>
//#include "hand-eye-calibration-without-visp.h"



#define PI  3.14159265358979

double deg2rad(double degree)
{
    return degree * PI / 180;
    }

//vector<string> split_string_by_delimiter(const string& strin, const string& delimiter)
std::vector<std::string> split_string_by_delimiter(const std::string& strin, std::string const& delims)
{
    //std::string const delims{ " .,:;!?" };
    std::vector<std::string> li_str;
    size_t beg, pos = 0;
    while ((beg = strin.find_first_not_of(delims, pos)) != std::string::npos)
    {
        pos = strin.find_first_of(delims, beg + 1);
        li_str.push_back(strin.substr(beg, pos - beg));
        //std::cout << li_str.back() << std::endl;
    }
    return li_str;
}    

bool is_only_number(const std::string& s)
{
    bool is_number = true;
    int n_pure_num = 0, n_dot = 0;
    if(s.empty()) is_number = false;
    else
    {
        for(auto it = s.begin(); it != s.end(); it++)
        {
            //std::cout << "*it : " << *it << std::endl;
            if(s.begin() == it)
            {
                if(!(std::isdigit(*it) || '-' == *it || '.' == *it)) {    is_number = false;  break; }
                if(std::isdigit(*it)) n_pure_num++;
                else
                {
                    if('.' == *it) n_dot++;
                }     
            }
            else
            {
                if('.' == *it) {
                    if(++n_dot > 1) { is_number = false;  break; }   
                }     
                else
                {
                    if(std::isdigit(*it)) n_pure_num++;
                    else { is_number = false;  break; }
                }     
            }
        }
    }
    if(0 == n_pure_num) is_number = false;
    //std::cout << "s : " << s << ", is_number : " << is_number << std::endl; 
    //exit(0);
    return is_number;
}

bool make_yaml_from_txt(std::string& prefix_base2gripper_camera2board, std::string prefix_yaml_base2gripper, std::string& prefix_yaml_camera2board, int id_start, int ndata)
{
    bool succeed = true;
    for(int id = id_start; id < id_start + ndata; id++)
    {
        std::vector<double> li_num;
        std::string path_txt = prefix_base2gripper_camera2board + std::to_string(id) + ".txt";
        std::cout << "Converting " << path_txt << " into YAML." << std::endl;
        std::ifstream input(path_txt);
        if(input.is_open()) {
            std::string line, delimiter = " =[],;";
            std::string const delims{ delimiter };
            while(std::getline(input, line)) {
                // using printf() in all tests for consistency
                //printf("%s\n", line.c_str());
                //if(line.empty() || '-' == line.at(0)) continue;
                if(line.empty()) continue;
                std::vector<std::string> li_str = split_string_by_delimiter(line, delims);

                for(auto str : li_str) 
                {
                    //str = "";
                    //str = "+45.7";
                    //str = "-45.7";
                    //str = "-0.457";
                    //str = "0.45.7";
                    //str = "457.";
                    //str = ".457";
                    if(is_only_number(str)) li_num.push_back(stof(str));
                }
            }
            input.close();
            //std::cout << "li_num.size() : " << li_num.size() << std::endl;
            //exit(0);
        }
        if(12 != li_num.size()) { succeed = false;    break;}
        double 
        cPo_tx = li_num[0], 
        cPo_ty = li_num[1], 
        cPo_tz = li_num[2], 
        cPo_rx = li_num[3], 
        cPo_ry = li_num[4], 
        cPo_rz = li_num[5],
        fPe_tx = li_num[6], 
        fPe_ty = li_num[7], 
        fPe_tz = li_num[8], 
        fPe_rx = li_num[9], 
        fPe_ry = li_num[10], 
        fPe_rz = li_num[11];
        
        if(cPo_tx > 10 || cPo_ty > 10 || cPo_tz > 10) { cPo_tx /= 1000.0;   cPo_ty /= 1000.0;   cPo_tz /= 1000.0; }
        
        if(fPe_tx > 10 || fPe_ty > 10 || fPe_tz > 10) { fPe_tx /= 1000.0;   fPe_ty /= 1000.0;   fPe_tz /= 1000.0; }
              
        if(cPo_rx > PI || cPo_ry > PI || cPo_rz > PI) { cPo_rx = deg2rad(cPo_rx);   cPo_ry = deg2rad(cPo_ry);   cPo_rz = deg2rad(cPo_rz); }
        if(fPe_rx > PI || fPe_ry > PI || fPe_rz > PI) { fPe_rx = deg2rad(fPe_rx);   fPe_ry = deg2rad(fPe_ry);   fPe_rz = deg2rad(fPe_rz); }

        vpPoseVector cPo(cPo_tx, cPo_ty, cPo_tz, cPo_rx, cPo_ry, cPo_rz), fPe(fPe_tx, fPe_ty, fPe_tz, fPe_rx, fPe_ry, fPe_rz);
        std::string path_cPo = prefix_yaml_camera2board + std::to_string(id) + ".yaml", path_fPe = prefix_yaml_base2gripper + std::to_string(id) + ".yaml";
        bool is_cPo_saved = cPo.saveYAML(path_cPo, cPo), is_fPe_saved = fPe.saveYAML(path_fPe, fPe);
        if(!(is_cPo_saved && is_fPe_saved))
        {
            if(!is_cPo_saved) std::cout << "Could NOT save cPo at : " << path_cPo << std::endl;
            if(!is_fPe_saved) std::cout << "Could NOT save fPe at : " << path_fPe << std::endl;
            exit(0);
        }     
        else
        {
            std::cout << "Saved cPo at : " << path_cPo << std::endl;
            std::cout << "Saved fPe at : " << path_fPe << std::endl;
        }
    }

    return succeed;
}

#ifdef HE_WO_VISP

Mat rot_mat_2_quaternion_2(const Mat& m33) {

    double qx, qy, qz, qw,  
        m00 = m33.at<double>(0, 0), m01 = m33.at<double>(0, 1), m02 = m33.at<double>(0, 2),
        m10 = m33.at<double>(1, 0), m11 = m33.at<double>(1, 1), m12 = m33.at<double>(1, 2),
        m20 = m33.at<double>(2, 0), m21 = m33.at<double>(2, 1), m22 = m33.at<double>(2, 2);

    float tr = m00 + m11 + m22;

    if (tr > 0) { 
        float S = sqrt(tr+1.0) * 2; // S=4*qw 
    qw = 0.25 * S;
      qx = (m21 - m12) / S;
        qy = (m02 - m20) / S; 
          qz = (m10 - m01) / S; 
          } else if ((m00 > m11)&(m00 > m22)) { 
            float S = sqrt(1.0 + m00 - m11 - m22) * 2; // S=4*qx 
              qw = (m21 - m12) / S;
                qx = 0.25 * S;
                  qy = (m01 + m10) / S; 
                    qz = (m02 + m20) / S; 
                    } else if (m11 > m22) { 
                      float S = sqrt(1.0 + m11 - m00 - m22) * 2; // S=4*qy
                        qw = (m02 - m20) / S;
                          qx = (m01 + m10) / S; 
                            qy = 0.25 * S;
                              qz = (m12 + m21) / S; 
                              } else { 
                                float S = sqrt(1.0 + m22 - m00 - m11) * 2; // S=4*qz
                                  qw = (m10 - m01) / S;
                                    qx = (m02 + m20) / S;
                                      qy = (m12 + m21) / S;
                                        qz = 0.25 * S;
                                        }
    
    Mat Q = (Mat_<double>(4, 1) << qx, qy, qz, qw);
    //Mat Q = (Mat_<double>(4, 1) << qw, qx, qy, qz);
    return Q;
}    

//void getQuaternion(Mat R, double Q[])
Mat rot_mat_2_quaternion(const Mat& R)
{
    Mat Q(4, 1, CV_64F);
    double trace = R.at<double>(0,0) + R.at<double>(1,1) + R.at<double>(2,2);
     
         if (trace > 0.0) 
             {
                     double s = sqrt(trace + 1.0);
                             //Q[3] = (s * 0.5);
                             Q.at<double>(3) = (s * 0.5);
                                     s = 0.5 / s;
                                             //Q[0] = ((R.at<double>(2,1) - R.at<double>(1,2)) * s);
                                             Q.at<double>(0) = ((R.at<double>(2,1) - R.at<double>(1,2)) * s);
                                                     //Q[1] = ((R.at<double>(0,2) - R.at<double>(2,0)) * s);
                                                     Q.at<double>(1) = ((R.at<double>(0,2) - R.at<double>(2,0)) * s);
                                                             //Q[2] = ((R.at<double>(1,0) - R.at<double>(0,1)) * s);
                                                             Q.at<double>(2) = ((R.at<double>(1,0) - R.at<double>(0,1)) * s);
                                                                 } 
                                                                     
                                                                         else 
                                                                             {
                                                                                     int i = R.at<double>(0,0) < R.at<double>(1,1) ? (R.at<double>(1,1) < R.at<double>(2,2) ? 2 : 1) : (R.at<double>(0,0) < R.at<double>(2,2) ? 2 : 0); 
                                                                                             int j = (i + 1) % 3;  
                                                                                                     int k = (i + 2) % 3;

                                                                                                             double s = sqrt(R.at<double>(i, i) - R.at<double>(j,j) - R.at<double>(k,k) + 1.0);
                                                                                                                     //Q[i] = s * 0.5;
                                                                                                                     Q.at<double>(i) = s * 0.5;
                                                                                                                             s = 0.5 / s;

                                                                                                                                     //Q[3] = (R.at<double>(k,j) - R.at<double>(j,k)) * s;
                                                                                                                                     Q.at<double>(3) = (R.at<double>(k,j) - R.at<double>(j,k)) * s;
                                                                                                                                             //Q[j] = (R.at<double>(j,i) + R.at<double>(i,j)) * s;
                                                                                                                                             Q.at<double>(j) = (R.at<double>(j,i) + R.at<double>(i,j)) * s;
                                                                                                                                                     //Q[k] = (R.at<double>(k,i) + R.at<double>(i,k)) * s;
                                                                                                                                                     Q.at<double>(k) = (R.at<double>(k,i) + R.at<double>(i,k)) * s;
                                                                                                                                                         }
                                                                                                                                                         return Q;
                                                                                                                                                         }

Mat vpHomo2Mat(const vpHomogeneousMatrix& vpHomo)
{
    int n_row = vpHomo.getRows(), n_col = vpHomo.getCols();//, CV_32F);
    Mat matHomo(n_row, n_col, CV_64F);
    double *data = (double *)matHomo.data;
    for(int iR = 0; iR < n_row; iR++)
    {
        int offset = iR * n_col;
        for(int iC = 0; iC < n_col; iC++)
        {
            data[offset + iC] = vpHomo[iR][iC];
        }
    }
    //std::cout << "vpHomo : " << std::endl << vpHomo << std::endl << "matHomo : " << std::endl << matHomo << std::endl; exit(0);
    return matHomo;
}
#endif  //  HE_WO_VISP

int main(int argc, char *argv[])
{
    unsigned int id_start = 0, ndata = 0;
    bool is_floyd = false;
    std::string prefix_base2gripper_camera2board, fn_gripper2camera, prefix_yaml_base2gripper, prefix_yaml_camera2board;
    for (int i = 1; i < argc; i++) {
        if (std::string(argv[i]) == "--ndata" && i + 1 < argc) {
            ndata = atoi(argv[i + 1]);
        }
        
        else if (std::string(argv[i]) == "--id_start" && i + 1 < argc) {
            id_start = atoi(argv[i + 1]);
        }

        else if (std::string(argv[i]) == "--fn_gripper2camera" && i + 1 < argc) 
        {
            fn_gripper2camera = argv[i + 1];
            //std::cout << "fn_gripper2camera : " << fn_gripper2camera << std::endl;  exit(0);
        }
        else if (std::string(argv[i]) == "--prefix_yaml_base2gripper" && i + 1 < argc) 
        {
            prefix_yaml_base2gripper = argv[i + 1];
        }
        else if (std::string(argv[i]) == "--prefix_yaml_camera2board" && i + 1 < argc) 
        {
            prefix_yaml_camera2board = argv[i + 1];
        }
        
        else if (std::string(argv[i]) == "--floyd_txt" && i + 1 < argc) 
        {
            is_floyd = true;    prefix_base2gripper_camera2board = argv[i + 1];
        }
       
        else if (std::string(argv[i]) == "--help" || std::string(argv[i]) == "-h") {
            std::cout << argv[0] << " [--ndata <number of data to process>] "
                "[--help] [-h]" << std::endl;
            return EXIT_SUCCESS;
        }
    }
    if (ndata == 0) {
        std::cout << "Number of data to process not specified" << std::endl;
        std::cout << argv[0] << " --help" << std::endl;
        return EXIT_SUCCESS;
    }
    std::vector<vpHomogeneousMatrix> cMo(ndata);
    std::vector<vpHomogeneousMatrix> wMe(ndata);
    vpHomogeneousMatrix eMc;

    
    if(is_floyd) 
    {   
        prefix_yaml_camera2board = prefix_base2gripper_camera2board + "cPo_";
        prefix_yaml_base2gripper = prefix_base2gripper_camera2board + "fPe_";

        make_yaml_from_txt(prefix_base2gripper_camera2board, prefix_yaml_base2gripper, prefix_yaml_camera2board, id_start, ndata);
    }
    //for (unsigned int i = id_start; i <= ndata; i++) {
    for (unsigned int idx = 0; idx < ndata; idx++) {
        int i = idx + id_start;
        std::ostringstream ss_fPe, ss_cPo;
        ss_fPe << prefix_yaml_base2gripper << i << ".yaml"; ss_cPo << prefix_yaml_camera2board << i << ".yaml";  
        std::cout << "Use fPe=" << ss_fPe.str() << ", cPo=" << ss_cPo.str() << std::endl;

        vpPoseVector wPe;
        if (wPe.loadYAML(ss_fPe.str(), wPe) == false) {
            std::cout << "Unable to read data from: " << ss_fPe.str() << std::endl;
            return EXIT_FAILURE;
        }
        //wMe[i - 1] = vpHomogeneousMatrix(wPe);
        wMe[idx] = vpHomogeneousMatrix(wPe);

        vpPoseVector cPo;
        if (cPo.loadYAML(ss_cPo.str(), cPo)  == false) {
            std::cout << "Unable to read data from: " << ss_cPo.str() << std::endl;
            return EXIT_FAILURE;
        }
        //cMo[i-1] = vpHomogeneousMatrix(cPo);
        cMo[idx] = vpHomogeneousMatrix(cPo);
    }

    int ret = vpHandEyeCalibration::calibrate(cMo, wMe, eMc);
#ifdef HE_WO_VISP
    int n_pose = cMo.size();
    std::vector<Mat> li_cMo(n_pose), li_wMe(n_pose);
    Mat eMc_wo_visp;
    for(int iP = 0; iP < n_pose; iP++) 
    {
        li_cMo[iP] = vpHomo2Mat(cMo[iP]);
        //std::cout << "li_cMo[iP] " << std::endl << li_cMo[iP] << std::endl;    exit(0);
        li_wMe[iP] = vpHomo2Mat(wMe[iP]);
    }     
    int ret2 = vpHandEyeCalibration::calibrate_wo_visp(li_cMo, li_wMe, eMc_wo_visp, 0);

    if (ret2 == 0) {
        std::cout << std::endl << "** Hand-eye calibration succeed : w/o visp" << std::endl;
        std::cout << std::endl << "** Hand-eye (eMc) transformation estimated:" << std::endl;
        std::cout << eMc_wo_visp << std::endl;
        //Mat
        std::pair<Mat, Mat> rot_tra = vpHandEyeCalibration::split_homogeneous_transform_matrix_into_rotation_and_translation(eMc_wo_visp, 0);
        Mat eTRc, ePc, eRc = rot_tra.first, eTc = rot_tra.second;
        Rodrigues(eRc, ePc);
        eTRc.push_back(eTc);    eTRc.push_back(ePc);
        std::cout << "** Corresponding pose vector: " << eTRc.t() << std::endl;
        //vpThetaUVector erc(eMc.getRotationMatrix());
        std::cout << std::endl << "** Translation [m]: " << eMc_wo_visp.at<double>(0, 3) << " " << eMc_wo_visp.at<double>(1, 3) << " " << eMc_wo_visp.at<double>(2, 3) << std::endl;
        std::cout << "** Rotation (theta-u representation) [rad]: " << ePc.t() << std::endl;
        std::cout << "** Rotation (theta-u representation) [deg]: " << vpHandEyeCalibration::rad2deg(ePc.at<double>(0)) << " " << vpHandEyeCalibration::rad2deg(ePc.at<double>(1)) << " " << vpHandEyeCalibration::rad2deg(ePc.at<double>(2)) << std::endl;
        //vpQuaternionVector quaternion(eMc.getRotationMatrix());
        Mat eQc = rot_mat_2_quaternion(eRc);
        std::cout << "** Rotation (quaternion representation) [rad]: " << eQc.t() << std::endl;
        Mat eQc2 = rot_mat_2_quaternion_2(eRc);
        std::cout << "** Rotation (quaternion representation) [rad] 2 : " << eQc2.t() << std::endl;


#if 0
        // save eMc
        //std::ofstream file_eMc("eMc.txt");
        std::ofstream file_eMc(fn_gripper2camera + ".txt");
        eMc.save(file_eMc);
        vpPoseVector pose_vec(eMc);
        //std::string output_filename("eMc.yaml");
        std::string output_filename(fn_gripper2camera + ".yaml");
        std::cout << std::endl << "Save transformation matrix eMc as a vpPoseVector in " << output_filename << std::endl;
        pose_vec.saveYAML(output_filename, pose_vec);
#endif  //  0        
    }
    else {
        std::cout << std::endl << "** Hand-eye calibration failed : w/o visp" << std::endl;
        std::cout << std::endl << "Check your input data and ensure they are covering the half sphere over the chessboard." << std::endl;
        std::cout << std::endl << "See https://visp-doc.inria.fr/doxygen/visp-daily/tutorial-calibration-extrinsic.html" << std::endl;
    }


#endif  //  HE_WO_VISP    

    if (ret == 0) {
        std::cout << std::endl << "** Hand-eye calibration succeed : ori" << std::endl;
        std::cout << std::endl << "** Hand-eye (eMc) transformation estimated:" << std::endl;
        std::cout << eMc << std::endl;
        std::cout << "** Corresponding pose vector: " << vpPoseVector(eMc).t() << std::endl;

        vpThetaUVector erc(eMc.getRotationMatrix());
        std::cout << std::endl << "** Translation [m]: " << eMc[0][3] << " " << eMc[1][3] << " " << eMc[2][3] << std::endl;
        std::cout << "** Rotation (theta-u representation) [rad]: " << erc.t() << std::endl;
        std::cout << "** Rotation (theta-u representation) [deg]: " << vpMath::deg(erc[0]) << " " << vpMath::deg(erc[1]) << " " << vpMath::deg(erc[2]) << std::endl;
        vpQuaternionVector quaternion(eMc.getRotationMatrix());
        std::cout << "** Rotation (quaternion representation) [rad]: " << quaternion.t() << std::endl;

        // save eMc
        //std::ofstream file_eMc("eMc.txt");
        std::ofstream file_eMc(fn_gripper2camera + ".txt");
        eMc.save(file_eMc);

        vpPoseVector pose_vec(eMc);
        //std::string output_filename("eMc.yaml");
        std::string output_filename(fn_gripper2camera + ".yaml");
        std::cout << std::endl << "Save transformation matrix eMc as a vpPoseVector in " << output_filename << std::endl;
        pose_vec.saveYAML(output_filename, pose_vec);
    }
    else {
        std::cout << std::endl << "** Hand-eye calibration failed : ori" << std::endl;
        std::cout << std::endl << "Check your input data and ensure they are covering the half sphere over the chessboard." << std::endl;
        std::cout << std::endl << "See https://visp-doc.inria.fr/doxygen/visp-daily/tutorial-calibration-extrinsic.html" << std::endl;
    }

    return EXIT_SUCCESS;
}

