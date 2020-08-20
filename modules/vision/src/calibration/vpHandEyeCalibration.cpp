/****************************************************************************
 *
 * ViSP, open source Visual Servoing Platform software.
 * Copyright (C) 2005 - 2019 by Inria. All rights reserved.
 *
 * This software is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 * See the file LICENSE.txt at the root directory of this source
 * distribution for additional information about the GNU GPL.
 *
 * For using ViSP with software that can not be combined with the GNU
 * GPL, please contact Inria about acquiring a ViSP Professional
 * Edition License.
 *
 * See http://visp.inria.fr for more information.
 *
 * This software was developed at:
 * Inria Rennes - Bretagne Atlantique
 * Campus Universitaire de Beaulieu
 * 35042 Rennes Cedex
 * France
 *
 * If you have questions regarding the use of this file, please contact
 * Inria at visp@inria.fr
 *
 * This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
 * WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
 *
 * Description:
 * Hand-eye calibration.
 *
 * Authors:
 * Francois Chaumette
 * Fabien Spindler
 *
 *****************************************************************************/

#include <cmath>  // std::fabs
#include <limits> // numeric_limits

#include <visp3/vision/vpHandEyeCalibration.h>

#define DEBUG_LEVEL1 0
#define DEBUG_LEVEL2 0

/*!
  \brief Compute the distances of the data to the mean obtained.

  \param[in] cMo : Vector of homogeneous matrices representing the transformation
  between the camera and the scene.
  \param[in] rMe : Vector of homogeneous matrices representing the transformation
  between the effector (where the camera is fixed) and the reference
  coordinates (base of the manipulator). Must be the same size as cMo.
  \param[in] eMc : Homogeneous matrix between the effector and the camera.
*/
void vpHandEyeCalibration::calibrationVerifrMo(const std::vector<vpHomogeneousMatrix> &cMo, const std::vector<vpHomogeneousMatrix> &rMe, const vpHomogeneousMatrix &eMc)
{
  unsigned int nbPose = (unsigned int) cMo.size();
  std::vector<vpTranslationVector> rTo(nbPose);
  std::vector<vpRotationMatrix> rRo(nbPose);

  for (unsigned int i = 0; i < nbPose; i++) {
    vpHomogeneousMatrix rMo = rMe[i] * eMc * cMo[i];
    rRo[i] = rMo.getRotationMatrix();
    rTo[i] = rMo.getTranslationVector();
  }
  vpRotationMatrix meanRot = vpRotationMatrix::mean(rRo);
  vpTranslationVector meanTrans = vpTranslationVector::mean(rTo);

#if DEBUG_LEVEL2
  {
    std::cout << "Mean  " << std::endl;
    std::cout << "Translation: " << meanTrans.t() << std::endl;
    vpThetaUVector P(meanRot);
    std::cout << "Rotation : theta (deg) = " << vpMath::deg(sqrt(P.sumSquare())) << " Matrice : " << std::endl << meanRot  << std::endl;
    std::cout << "theta U (deg): " << vpMath::deg(P[0]) << " " << vpMath::deg(P[1]) << " " << vpMath::deg(P[2]) << std::endl;
  }
#endif

  // standard deviation, rotational part
  double resRot = 0.0;
  for (unsigned int i = 0; i < nbPose; i++) {
    vpRotationMatrix R = meanRot.t() * rRo[i]; // Rm^T  Ri
    vpThetaUVector P(R);
    // theta = Riemannian distance d(Rm,Ri)
    double theta = sqrt(P.sumSquare());
    std::cout << "Distance theta between rMo(" << i << ") and mean (deg) = " << vpMath::deg(theta) << std::endl;
    // Euclidean distance d(Rm,Ri) not used
    // theta = 2.0*sqrt(2.0)*sin(theta/2.0);
    resRot += theta*theta;
  }
  resRot = sqrt(resRot/nbPose);
  std::cout << "Mean residual rMo(" << nbPose << ") - rotation (deg) = " << vpMath::deg(resRot) << std::endl;
  // standard deviation, translational part
  double resTrans = 0.0;
  for (unsigned int i = 0; i < nbPose; i++) {
    vpColVector errTrans = ((vpColVector) rTo[i]) - meanTrans;
    resTrans += errTrans.sumSquare();
    std::cout << "Distance d between rMo(" << i << ") and mean (m) = " << sqrt(errTrans.sumSquare()) << std::endl;
  }
  resTrans = sqrt(resTrans/nbPose);
  std::cout << "Mean residual rMo(" << nbPose << ") - translation (m) = " << resTrans << std::endl;
  double resPos = (resRot*resRot + resTrans*resTrans)*nbPose;
  resPos = sqrt(resPos/(2*nbPose));
  std::cout << "Mean residual rMo(" << nbPose << ") - global = " << resPos << std::endl;
}

/*!
  \brief Compute the rotation part (eRc) of hand-eye pose by solving a
  Procrustes problem [... (theta u)_e ...] = eRc [ ... (theta u)_c ...]

  \param[in] cMo : Vector of homogeneous matrices representing the transformation
  between the camera and the scene (input)
  \param[in] rMe : Vector of homogeneous matrices representing the transformation
  between the effector (where the camera is fixed) and the reference
  coordinates (base of the manipulator) (input). Must be the same size as cMo.
  \param[out] eRc : Rotation matrix  between the effector and the camera (output)
*/
int vpHandEyeCalibration::calibrationRotationProcrustes(const std::vector<vpHomogeneousMatrix> &cMo, const std::vector<vpHomogeneousMatrix> &rMe,vpRotationMatrix &eRc)
{
  // Method by solving the orthogonal Procrustes problem
  // [... (theta u)_e ...] = eRc [ ... (theta u)_c ...]
  // similar to E^T = eRc C^T below

  vpMatrix Et,Ct;
  vpMatrix A;
  unsigned int k = 0;
  unsigned int nbPose = (unsigned int) cMo.size();

  // for all couples ij
  for (unsigned int i = 0; i < nbPose; i++) {
    vpRotationMatrix rRei, ciRo;
    rMe[i].extract(rRei);
    cMo[i].extract(ciRo);
    // std::cout << "rMei: " << std::endl << rMe[i] << std::endl;

    for (unsigned int j = 0; j < nbPose; j++) {
      if (j > i) // we don't use two times same couples...
      {
        vpRotationMatrix rRej, cjRo;
        rMe[j].extract(rRej);
        cMo[j].extract(cjRo);
        // std::cout << "rMej: " << std::endl << rMe[j] << std::endl;

        vpRotationMatrix ejRei = rRej.t() * rRei;
        vpThetaUVector ejPei(ejRei);
        vpColVector xe = ejPei;

        vpRotationMatrix cjRci = cjRo * ciRo.t();
        vpThetaUVector cjPci(cjRci);
        vpColVector xc = cjPci;

        if (k == 0) {
          Et = xe.t();
          Ct = xc.t();
        } else {
          Et.stack(xe.t());
          Ct.stack(xc.t());
        }
        k++;
      }
    }
  }
  // std::cout << "Et "  << std::endl << Et << std::endl;
  // std::cout << "Ct "  << std::endl << Ct << std::endl;

  // R obtained from the SVD of (E C^T) with all singular values equal to 1
  A = Et.t() * Ct;
  vpMatrix M, U, V;
  vpColVector sv;
  int rank = A.pseudoInverse(M, sv, 1e-6, U, V);
  if (rank != 3) return -1;
  A = U * V.t();
  eRc = vpRotationMatrix(A);
#if 0
    vpThetaUVector ePc(eRc);
    std::cout << "Rotation from Procrustes method : ori" << std::endl;
    std::cout << "theta U (deg): " << vpMath::deg(ePc[0]) << " " << vpMath::deg(ePc[1]) << " " << vpMath::deg(ePc[2]) << std::endl;
    // Residual
    vpMatrix residual;
    residual = A * Ct.t() - Et.t();
    //  std::cout << "Residual: " << std::endl << residual << std::endl;
    double res = sqrt(residual.sumSquare()/(residual.getRows()*residual.getCols()));
    printf("Mean residual (rotation) = %lf\n",res);
#endif  //  1

#if DEBUG_LEVEL2
  {
    vpThetaUVector ePc(eRc);
    std::cout << "Rotation from Procrustes method " << std::endl;
    std::cout << "theta U (deg): " << vpMath::deg(ePc[0]) << " " << vpMath::deg(ePc[1]) << " " << vpMath::deg(ePc[2]) << std::endl;
    // Residual
    vpMatrix residual;
    residual = A * Ct.t() - Et.t();
    //  std::cout << "Residual: " << std::endl << residual << std::endl;
    double res = sqrt(residual.sumSquare()/(residual.getRows()*residual.getCols()));
    printf("Mean residual (rotation) = %lf\n",res);
  }
#endif
  return 0;
}

/*!
  \brief Compute the rotation part (eRc) of hand-eye pose by solving a
  linear system using R. Tsai and R. Lorenz method
  \cite Tsai89a.

  \param[in] cMo : Vector of homogeneous matrices representing the transformation
  between the camera and the scene (input)
  \param[in] rMe : Vector of homogeneous matrices representing the transformation
  between the effector (where the camera is fixed) and the reference
  coordinates (base of the manipulator) (input). Must be the same size as cMo.
  \param[out] eRc : Rotation matrix  between the effector and the camera (output)
*/
int vpHandEyeCalibration::calibrationRotationTsai(const std::vector<vpHomogeneousMatrix> &cMo, const std::vector<vpHomogeneousMatrix> &rMe,vpRotationMatrix &eRc)
{
  vpMatrix A;
  vpColVector B;
  unsigned int nbPose = (unsigned int) cMo.size();
  unsigned int k = 0;
  // for all couples ij
  for (unsigned int i = 0; i < nbPose; i++) {
    vpRotationMatrix rRei, ciRo;
    rMe[i].extract(rRei);
    cMo[i].extract(ciRo);
    // std::cout << "rMei: " << std::endl << rMe[i] << std::endl;

    for (unsigned int j = 0; j < nbPose; j++) {
      if (j > i) // we don't use two times same couples...
      {
        vpRotationMatrix rRej, cjRo;
        rMe[j].extract(rRej);
        cMo[j].extract(cjRo);
        // std::cout << "rMej: " << std::endl << rMe[j] << std::endl;

        vpRotationMatrix ejRei = rRej.t() * rRei;
        vpThetaUVector ejPei(ejRei);

        vpRotationMatrix cjRci = cjRo * ciRo.t();
        vpThetaUVector cjPci(cjRci);
        // std::cout << "theta U (camera) " << cjPci.t() << std::endl;

        vpMatrix As;
        vpColVector b(3);

        As = vpColVector::skew(vpColVector(ejPei) + vpColVector(cjPci));

        b =  (vpColVector)cjPci - (vpColVector) ejPei; // A.40

        if (k == 0) {
          A = As;
          B = b;
        } else {
          A = vpMatrix::stack(A, As);
          B = vpColVector::stack(B, b);
        }
        k++;
      }
    }
  }
#if DEBUG_LEVEL2
  {
    std::cout << "Tsai method: system A X = B "  << std::endl;
    std::cout << "A "  << std::endl << A << std::endl;
    std::cout << "B "  << std::endl << B << std::endl;
  }
#endif
  vpMatrix Ap;
  // the linear system A x = B is solved
  // using x = A^+ B

  int rank = A.pseudoInverse(Ap);
  if (rank != 3) return -1;

  vpColVector x = Ap * B;

  // extraction of theta U

  // x = tan(theta/2) U
  double norm =  x.sumSquare();
  double c = 1 / sqrt(1 + norm);  // cos(theta/2)
  double alpha = acos(c);         // theta/2
  norm = 2.0*c/vpMath::sinc(alpha);  // theta / tan(theta/2)
  for (unsigned int i = 0; i < 3; i++) x[i] *= norm;

  // Building of the rotation matrix eRc
  vpThetaUVector xP(x[0], x[1], x[2]);
  eRc = vpRotationMatrix(xP);
#if 0
    std::cout << "Rotation from Tsai method : ori" << std::endl;
    std::cout << "theta U (deg): " << vpMath::deg(x[0]) << " " << vpMath::deg(x[1]) << " " << vpMath::deg(x[2]) << std::endl;
    // Residual
    for (unsigned int i = 0; i < 3; i++) x[i] /= norm; /* original x */
    vpColVector residual;
    residual = A*x-B;
    // std::cout << "Residual: " << std::endl << residual << std::endl;
    double res = sqrt(residual.sumSquare()/residual.getRows());
    printf("Mean residual (rotation) = %lf\n",res);
#endif  //  1
#if DEBUG_LEVEL2
  {
    std::cout << "Rotation from Tsai method" << std::endl;
    std::cout << "theta U (deg): " << vpMath::deg(x[0]) << " " << vpMath::deg(x[1]) << " " << vpMath::deg(x[2]) << std::endl;
    // Residual
    for (unsigned int i = 0; i < 3; i++) x[i] /= norm; /* original x */
    vpColVector residual;
    residual = A*x-B;
    // std::cout << "Residual: " << std::endl << residual << std::endl;
    double res = sqrt(residual.sumSquare()/residual.getRows());
    printf("Mean residual (rotation) = %lf\n",res);
  }
#endif
  return 0;
}

/*!
  \brief Old ViSP implementation for computing the rotation part (eRc)
  of hand-eye pose by solving a linear system using R. Tsai and R. Lorenz method
  \cite Tsai89a.

  \param[in] cMo : Vector of homogeneous matrices representing the transformation
  between the camera and the scene (input)
  \param[in] rMe : Vector of homogeneous matrices representing the transformation
  between the effector (where the camera is fixed) and the reference
  coordinates (base of the manipulator) (input). Must be the same size as cMo.
  \param[ou] eRc : Rotation matrix  between the effector and the camera (output)
*/
int vpHandEyeCalibration::calibrationRotationTsaiOld(const std::vector<vpHomogeneousMatrix> &cMo, const std::vector<vpHomogeneousMatrix> &rMe,vpRotationMatrix &eRc)
{
  unsigned int nbPose = (unsigned int) cMo.size();
  vpMatrix A;
  vpColVector B;
  vpColVector x;
  unsigned int k = 0;
  // for all couples ij
  for (unsigned int i = 0; i < nbPose; i++) {
    vpRotationMatrix rRei, ciRo;
    rMe[i].extract(rRei);
    cMo[i].extract(ciRo);
#if 0
    std::cout << std::endl << "i : " << i << std::endl;
    std::cout << "rRei ori : " << std::endl << rRei << std::endl;
    std::cout << "ciRo ori : " << std::endl << ciRo << std::endl;
#endif  //  1
    for (unsigned int j = 0; j < nbPose; j++) {
      if (j > i) { // we don't use two times same couples...
        vpRotationMatrix rRej, cjRo;
        rMe[j].extract(rRej);
        cMo[j].extract(cjRo);
        // std::cout << "rMej: " << std::endl << rMe[j] << std::endl;

        vpRotationMatrix rReij = rRej.t() * rRei;

        vpRotationMatrix cijRo = cjRo * ciRo.t();

        vpThetaUVector rPeij(rReij);

        double theta = sqrt(rPeij[0] * rPeij[0] + rPeij[1] * rPeij[1] + rPeij[2] * rPeij[2]);
#if 0
        std::cout << std::endl << "calibrationRotationTsaiOld ori" << std::endl;
        std::cout << i << " " << j << " " << "ejRei: " << std::endl << rReij << std::endl;
        std::cout << "theta (robot) " << theta << std::endl;
        std::cout << "theta U (robot) " << rPeij << std::endl;
        std::cout << "cjRci: " << std::endl << cijRo.t() << std::endl;
#endif  //  1
        for (unsigned int m = 0; m < 3; m++) {
          rPeij[m] = rPeij[m] * vpMath::sinc(theta / 2);
        }

        vpThetaUVector cijPo(cijRo);
        theta = sqrt(cijPo[0] * cijPo[0] + cijPo[1] * cijPo[1] + cijPo[2] * cijPo[2]);
        for (unsigned int m = 0; m < 3; m++) {
          cijPo[m] = cijPo[m] * vpMath::sinc(theta / 2);
        }

        // std::cout << "theta (camera) " << theta << std::endl;
        // std::cout << "theta U (camera) " << cijPo.t() << std::endl;

        vpMatrix As;
        vpColVector b(3);

        As = vpColVector::skew(vpColVector(rPeij) + vpColVector(cijPo));

        b = (vpColVector)cijPo - (vpColVector)rPeij; // A.40

        if (k == 0) {
          A = As;
          B = b;
        } else {
          A = vpMatrix::stack(A, As);
          B = vpColVector::stack(B, b);
        }
        k++;
      }
    }
  }

    //std::cout << "A ori "  << std::endl << A << std::endl << "B ori "  << std::endl << B << std::endl;

  // the linear system is defined
  // x = AtA^-1AtB is solved
  vpMatrix AtA = A.AtA();

  vpMatrix Ap;
  int rank = AtA.pseudoInverse(Ap, 1e-6);
  if (rank != 3) return -1;

  x = Ap * A.t() * B;
  vpColVector x2 = x; /* pour calcul residu */
  //std::cout << "x2 ori : " << std::endl << x2 << std::endl << std::endl;

  //     {
  //       // Residual
  //       vpColVector residual;
  //       residual = A*x-B;
  //       std::cout << "Residual: " << std::endl << residual << std::endl;

  //       double res = 0;
  //       for (int i=0; i < residual.getRows(); i++)
  // 	res += residual[i]*residual[i];
  //       res = sqrt(res/residual.getRows());
  //       printf("Mean residual = %lf\n",res);
  //     }

  // extraction of theta and U
  double theta;
  double d = x.sumSquare();
  for (unsigned int i = 0; i < 3; i++)
    x[i] = 2 * x[i] / sqrt(1 + d);
  //std::cout << "x *= 2.0 / sqrt(1.0 + d) ori : " << std::endl << x << std::endl << std::endl;
  theta = sqrt(x.sumSquare()) / 2;
  theta = 2 * asin(theta);
  // if (theta !=0)
  if (std::fabs(theta) > std::numeric_limits<double>::epsilon()) {
    for (unsigned int i = 0; i < 3; i++)
      x[i] *= theta / (2 * sin(theta / 2));
  } else
    x = 0;

  // Building of the rotation matrix eRc
  vpThetaUVector xP(x[0], x[1], x[2]);
  eRc = vpRotationMatrix(xP);

#if 0
    std::cout << "Rotation from Old Tsai method ori" << std::endl;
    std::cout << "theta U (deg): " << vpMath::deg(x[0]) << " " << vpMath::deg(x[1]) << " " << vpMath::deg(x[2]) << std::endl;
    // Residual
    vpColVector residual;
    residual = A*x2-B;
    // std::cout << "Residual: " << std::endl << residual << std::endl;
    double res = sqrt(residual.sumSquare()/residual.getRows());
    printf("Mean residual (rotation) = %lf\n",res);
#endif

#if DEBUG_LEVEL2
  {
    std::cout << "Rotation from Old Tsai method" << std::endl;
    std::cout << "theta U (deg): " << vpMath::deg(x[0]) << " " << vpMath::deg(x[1]) << " " << vpMath::deg(x[2]) << std::endl;
    // Residual
    vpColVector residual;
    residual = A*x2-B;
    // std::cout << "Residual: " << std::endl << residual << std::endl;
    double res = sqrt(residual.sumSquare()/residual.getRows());
    printf("Mean residual (rotation) = %lf\n",res);
  }
#endif
  return 0;
}

/*!
  \brief Compute the translation part (eTc) of hand-eye pose by solving a
  linear system (see for instance R. Tsai and R. Lorenz method)
  \cite Tsai89a.

  \param[in] cMo : Vector of homogeneous matrices representing the transformation
  between the camera and the scene (input)
  \param[in] rMe : Vector of homogeneous matrices representing the transformation
  between the effector (where the camera is fixed) and the reference
  coordinates (base of the manipulator) (input). Must be the same size as cMo.
  \param[out] eRc : Rotation matrix  between the effector and the camera (input)
  \param[out] eTc : Translation  between the effector and the camera (output)
*/
int vpHandEyeCalibration::calibrationTranslation(const std::vector<vpHomogeneousMatrix> &cMo,
                                                 const std::vector<vpHomogeneousMatrix> &rMe,
                                                 vpRotationMatrix &eRc,
                                                 vpTranslationVector &eTc)
{
  vpMatrix I3(3,3);
  I3.eye();
  unsigned int k = 0;
  unsigned int nbPose = (unsigned int)cMo.size();
  vpMatrix A(3*nbPose,3);
  vpColVector B(3*nbPose);
  // Building of the system for the translation estimation
  // for all couples ij
  for (unsigned int i = 0; i < nbPose; i++) {
    for (unsigned int j = 0; j < nbPose; j++) {
      if (j > i) { // we don't use two times same couples...
        vpHomogeneousMatrix ejMei = rMe[j].inverse() * rMe[i];
        vpHomogeneousMatrix cjMci = cMo[j] * cMo[i].inverse();

        vpRotationMatrix ejRei, cjRci;
        vpTranslationVector ejTei, cjTci;

        ejMei.extract(ejRei);
        ejMei.extract(ejTei);

        cjMci.extract(cjRci);
        cjMci.extract(cjTci);

        vpMatrix a = vpMatrix(ejRei) - I3;
        vpTranslationVector b = eRc * cjTci - ejTei;

        if (k == 0) {
          A = a;
          B = b;
        } else {
          A = vpMatrix::stack(A, a);
          B = vpColVector::stack(B, b);
        }
        k++;
      }
    }
  }

  // the linear system A x = B is solved
  // using x = A^+ B
  vpMatrix Ap;
  int rank = A.pseudoInverse(Ap);
  if (rank != 3) return -1;

  vpColVector x = Ap * B;
  eTc = (vpTranslationVector) x;
#if 0
    printf("New Hand-eye calibration ori : ");
    std::cout << "Translation: " << eTc[0] << " " << eTc[1] << " " << eTc[2] << std::endl;
    // residual
    vpColVector residual;
    residual = A*x-B;
    // std::cout << "Residual: " << std::endl << residual << std::endl;
    double res = sqrt(residual.sumSquare()/residual.getRows());
    printf("Mean residual (translation) = %lf\n",res);
#endif  //  1

#if DEBUG_LEVEL2
  {
    printf("New Hand-eye calibration : ");
    std::cout << "Translation: " << eTc[0] << " " << eTc[1] << " " << eTc[2] << std::endl;
    // residual
    vpColVector residual;
    residual = A*x-B;
    // std::cout << "Residual: " << std::endl << residual << std::endl;
    double res = sqrt(residual.sumSquare()/residual.getRows());
    printf("Mean residual (translation) = %lf\n",res);
  }
#endif
  return 0;
}

/*!
  \brief Old method to compute the translation part (eTc) of hand-eye pose
  by solving a linear system (see for instance R. Tsai and R. Lorenz method)
  \cite Tsai89a.

  \param[in] cMo : vector of homogeneous matrices representing the transformation
  between the camera and the scene (input)
  \param[in] rMe : vector of homogeneous matrices representing the transformation
  between the effector (where the camera is fixed) and the reference
  coordinates (base of the manipulator) (input). Must be the same size as cMo.
  \param[out] eRc : rotation matrix  between the effector and the camera (input)
  \param[out] eTc : translation  between the effector and the camera (output)
*/
int vpHandEyeCalibration::calibrationTranslationOld(const std::vector<vpHomogeneousMatrix> &cMo,
                                                    const std::vector<vpHomogeneousMatrix> &rMe,
                                                    vpRotationMatrix &eRc,
                                                    vpTranslationVector &eTc)
{
  vpMatrix A;
  vpColVector B;
  // Building of the system for the translation estimation
  // for all couples ij
  vpRotationMatrix I3;
  I3.eye();
  int k = 0;
  unsigned int nbPose = (unsigned int)cMo.size();

  for (unsigned int i = 0; i < nbPose; i++) {
    vpRotationMatrix rRei, ciRo;
    vpTranslationVector rTei, ciTo;
    rMe[i].extract(rRei);
    cMo[i].extract(ciRo);
    rMe[i].extract(rTei);
    cMo[i].extract(ciTo);

    for (unsigned int j = 0; j < nbPose; j++) {
      if (j > i) // we don't use two times same couples...
      {

        vpRotationMatrix rRej, cjRo;
        rMe[j].extract(rRej);
        cMo[j].extract(cjRo);

        vpTranslationVector rTej, cjTo;
        rMe[j].extract(rTej);
        cMo[j].extract(cjTo);

        vpRotationMatrix rReij = rRej.t() * rRei;

        vpTranslationVector rTeij = rTej + (-rTei);

        rTeij = rRej.t() * rTeij;

        vpMatrix a = vpMatrix(rReij) - vpMatrix(I3);

        vpTranslationVector b;
        b = eRc * cjTo - rReij * eRc * ciTo + rTeij;

        if (k == 0) {
          A = a;
          B = b;
        } else {
          A = vpMatrix::stack(A, a);
          B = vpColVector::stack(B, b);
        }
        k++;
      }
    }
  }

  // the linear system is solved
  // x = AtA^-1AtB is solved
  vpMatrix AtA = A.AtA();
  vpMatrix Ap;
  vpColVector AeTc;
  int rank = AtA.pseudoInverse(Ap, 1e-6);
  if (rank != 3) return -1;

  AeTc = Ap * A.t() * B;
  eTc = (vpTranslationVector) AeTc;
#if 0
    printf("Old Hand-eye calibration ori : ");
    std::cout << "Translation: " << eTc[0] << " " << eTc[1] << " " << eTc[2] << std::endl;

    // residual
    vpColVector residual;
    residual = A*AeTc-B;
    // std::cout << "Residual: " << std::endl << residual << std::endl;
    double res = 0;
    for (unsigned int i=0; i < residual.getRows(); i++)
      res += residual[i]*residual[i];
    res = sqrt(res/residual.getRows());
    printf("Mean residual (translation) = %lf\n",res);
#endif  //  1

#if DEBUG_LEVEL2
  {
    printf("Old Hand-eye calibration : ");
    std::cout << "Translation: " << eTc[0] << " " << eTc[1] << " " << eTc[2] << std::endl;

    // residual
    vpColVector residual;
    residual = A*AeTc-B;
    // std::cout << "Residual: " << std::endl << residual << std::endl;
    double res = 0;
    for (unsigned int i=0; i < residual.getRows(); i++)
      res += residual[i]*residual[i];
    res = sqrt(res/residual.getRows());
    printf("Mean residual (translation) = %lf\n",res);
  }
#endif
  return 0;
}

/*!
  \brief Compute the set of errors minimised by VVS.

  \param[in] cMo : vector of homogeneous matrices representing the transformation
  between the camera and the scene.
  \param[in] rMe : vector of homogeneous matrices representing the transformation
  between the effector (where the camera is fixed) and the reference
  coordinates (base of the manipulator). Must be the same size as cMo.
  \param[in] eMc : homogeneous matrix between the effector and the camera (input)
  \param[out] err: set of errors minimised by VVS (3 for rotation, 3 for translation, etc.) (output)
*/

double vpHandEyeCalibration::calibrationErrVVS(const std::vector<vpHomogeneousMatrix> &cMo, const std::vector<vpHomogeneousMatrix> &rMe, const vpHomogeneousMatrix &eMc, vpColVector &errVVS)
{
  unsigned int nbPose = (unsigned int) cMo.size();
  vpMatrix I3(3,3);
  I3.eye();
  vpRotationMatrix eRc;
  vpTranslationVector eTc;
  eMc.extract(eRc);
  eMc.extract(eTc);

  unsigned int k = 0;
  for (unsigned int i = 0; i < nbPose; i++) {
    for (unsigned int j = 0; j < nbPose; j++) {
      if (j > i) // we don't use two times same couples...
      {
        vpColVector s(3);

        vpHomogeneousMatrix ejMei = rMe[j].inverse() * rMe[i];
        vpHomogeneousMatrix cjMci = cMo[j] * cMo[i].inverse();

        vpRotationMatrix ejRei, cjRci;
        vpTranslationVector ejTei, cjTci;

        ejMei.extract(ejRei);
        vpThetaUVector ejPei(ejRei);
        ejMei.extract(ejTei);

        cjMci.extract(cjRci);
        vpThetaUVector cjPci(cjRci);
        cjMci.extract(cjTci);
        // terms due to rotation
        s = vpMatrix(eRc) * vpColVector(cjPci) - vpColVector(ejPei);
        if (k == 0) {
          errVVS = s;
        } else {
          errVVS = vpColVector::stack(errVVS, s);
        }
        k++;
        // terms due to translation
        s = (vpMatrix(ejRei) - I3) * eTc - eRc * cjTci + ejTei;
        errVVS = vpColVector::stack(errVVS, s);
      } // enf if i > j
    } // end for j
  } // end for i

  double resRot, resTrans, resPos;
  resRot = resTrans = resPos = 0.0;
  //std::cout << "errVVS.size() : " << errVVS.size() << std::endl;
  for (unsigned int i=0; i < (unsigned int) errVVS.size() ; i += 6)
  {
    resRot += errVVS[i]*errVVS[i];
    resRot += errVVS[i+1]*errVVS[i+1];
    resRot += errVVS[i+2]*errVVS[i+2];
    resTrans += errVVS[i+3]*errVVS[i+3];
    resTrans += errVVS[i+4]*errVVS[i+4];
    resTrans += errVVS[i+5]*errVVS[i+5];
  }
  resPos = resRot + resTrans;
  resRot = sqrt(resRot*2/errVVS.size());
  resTrans = sqrt(resTrans*2/errVVS.size());
  resPos = sqrt(resPos/errVVS.size());
#if 0
    printf("\ncalibrationErrVVS ori\n");
    printf("Mean VVS residual - rotation (deg) = %lf\n",vpMath::deg(resRot));
    printf("Mean VVS residual - translation = %lf\n",resTrans);
    printf("Mean VVS residual - global = %lf\n",resPos);
#endif  //  HE_WO_VISP

#if DEBUG_LEVEL1
  {
    printf("Mean VVS residual - rotation (deg) = %lf\n",vpMath::deg(resRot));
    printf("Mean VVS residual - translation = %lf\n",resTrans);
    printf("Mean VVS residual - global = %lf\n",resPos);
  }
#endif
  return resPos;
}

/*!
  \brief Hand-Eye Calibration by VVS.

  \param[in] cMo : vector of homogeneous matrices representing the transformation
  between the camera and the scene.
  \param[in] rMe : vector of homogeneous matrices representing the transformation
  between the effector (where the camera is fixed) and the reference
  coordinates (base of the manipulator). Must be the same size as cMo.
  \param[in,out] eMc : homogeneous matrix between the effector and the camera.
*/
#define NB_ITER_MAX 30

int vpHandEyeCalibration::calibrationVVS(const std::vector<vpHomogeneousMatrix> &cMo,
                                         const std::vector<vpHomogeneousMatrix> &rMe,
                                         vpHomogeneousMatrix &eMc)
{
  unsigned int it = 0;
  double res = 1.0;
  unsigned int nbPose = (unsigned int) cMo.size();
  vpColVector err;
  vpMatrix L;
  vpMatrix I3(3,3);
  I3.eye();
  vpRotationMatrix eRc;
  vpTranslationVector eTc;
  eMc.extract(eRc);
  eMc.extract(eTc);

  /* FC : on recalcule 2 fois tous les ejMei et cjMci a chaque iteration
    alors qu'ils sont constants. Ce serait sans doute mieux de les
    calculer une seule fois et de les stocker. Pourraient alors servir
    dans les autres fonctions HandEye. A voir si vraiment interessant vu la
    combinatoire. Idem pour les theta u */
  while ((res > 1e-7) && (it < NB_ITER_MAX))
  {
    /* compute s - s^* */
    vpHandEyeCalibration::calibrationErrVVS(cMo, rMe, eMc, err);
    /* compute L_s */
    unsigned int k = 0;
    for (unsigned int i = 0; i < nbPose; i++) {
      for (unsigned int j = 0; j < nbPose; j++) {
        if (j > i) // we don't use two times same couples...
        {
          vpMatrix Ls(3,6),Lv(3,3),Lw(3,3);

          vpHomogeneousMatrix ejMei = rMe[j].inverse() * rMe[i];
          vpHomogeneousMatrix cjMci = cMo[j] * cMo[i].inverse();

          vpRotationMatrix ejRei;
          ejMei.extract(ejRei);
          vpThetaUVector cjPci(cjMci);

          vpTranslationVector cjTci;

          cjMci.extract(cjTci);
          // terms due to rotation
          //Lv.diag(0.0); //
          Lv = 0.0;
          Lw = -vpMatrix(eRc) * vpColVector::skew(vpColVector(cjPci));
          for (unsigned int m=0;m<3;m++)
            for (unsigned int n=0;n<3;n++)
            {
              Ls[m][n] = Lv[m][n];
              Ls[m][n+3] = Lw[m][n];
            }
          if (k == 0) {
            L = Ls;
          } else {
            L = vpMatrix::stack(L,Ls);
          }
          k++;
          // terms due to translation
          Lv = (vpMatrix(ejRei) - I3) * vpMatrix(eRc);
          Lw =  vpMatrix(eRc) * vpColVector::skew(vpColVector(cjTci));
          for (unsigned int m=0;m<3;m++)
            for (unsigned int n=0;n<3;n++)
            {
              Ls[m][n] = Lv[m][n];
              Ls[m][n+3] = Lw[m][n];
            }
          L = vpMatrix::stack(L,Ls);

        } // enf if i > j
      } // end for j
    } // end for i
    double lambda = 0.9;
    vpMatrix Lp;
    int rank = L.pseudoInverse(Lp);
    if (rank != 6) return -1;

    vpColVector e = Lp * err;
    vpColVector v = - e * lambda;
    //  std::cout << "e: "  << e.t() << std::endl;
    eMc = eMc * vpExponentialMap::direct(v);
    eMc.extract(eRc);
    eMc.extract(eTc);
    res = sqrt(v.sumSquare()/v.getRows());
#if DEBUG_LEVEL2
    std::cout << "it : " << it << ",\t res : " << res << std::endl; 
#endif    //  DEBUG_LEVEL2
    it++;
  } // end while
#if DEBUG_LEVEL2
  {
    printf("Iteration number for NL hand-eye minimisation ori : %d / %d\n", it, NB_ITER_MAX);
    vpThetaUVector ePc(eRc);
    std::cout << "theta U (deg): " << vpMath::deg(ePc[0]) << " " << vpMath::deg(ePc[1]) << " " << vpMath::deg(ePc[2]) << std::endl;
    std::cout << "Translation: " << eTc[0] << " " << eTc[1] << " " << eTc[2] << std::endl;
    // Residual
    double res = err.sumSquare();
    res = sqrt(res/err.getRows());
    printf("Mean residual (rotation+translation) = %lf\n",res);
  }
#endif    //  DEBUG_LEVEL2
  if (it == NB_ITER_MAX) return 1;  // VVS has not converged before NB_ITER_MAX
  else return 0;
}

#undef NB_ITER_MAX

#define HE_I 0
#define HE_TSAI_OROT 1
#define HE_TSAI_ORNT 2
#define HE_TSAI_NROT 3
#define HE_TSAI_NRNT 4
#define HE_PROCRUSTES_OT 5
#define HE_PROCRUSTES_NT 6

/*!
  Compute extrinsic camera parameters : the constant transformation from
  the effector to the camera frames (eMc).

  \param[in] cMo : vector of homogeneous matrices representing the transformation
  between the camera and the scene.
  \param[in] rMe : vector of homogeneous matrices representing the transformation
  between the effector (where the camera is fixed) and the reference
  coordinates (base of the manipulator). Must be the same size as cMo.
  \param[out] eMc : homogeneous matrix representing the transformation
  between the effector and the camera (output)

  \return 0 if calibration succeed, -1 if the system is not full rank, 1 if the algorithm doesn't converge.
*/
int vpHandEyeCalibration::calibrate(const std::vector<vpHomogeneousMatrix> &cMo,
                                    const std::vector<vpHomogeneousMatrix> &rMe, vpHomogeneousMatrix &eMc)
{
  if (cMo.size() != rMe.size())
    throw vpException(vpException::dimensionError, "cMo and rMe have different sizes");

  vpRotationMatrix eRc;
  vpTranslationVector eTc;
  vpColVector errVVS;
  double resPos;

  /* initialisation of eMc to I in case all other methods fail */
  eMc.eye();
  resPos = vpHandEyeCalibration::calibrationErrVVS(cMo, rMe, eMc, errVVS);
  double vmin = resPos;  // will serve to determine the best method
#if DEBUG_LEVEL1
  int He_method = HE_I;  // will serve to know which is the best method
#endif
  vpHomogeneousMatrix eMcMin = eMc;  // best initial estimation for VSS
  // Method using Old Tsai implementation
  int err = vpHandEyeCalibration::calibrationRotationTsaiOld(cMo, rMe, eRc);
  if (err != 0) printf("\n Problem in solving Hand-Eye Rotation by Old Tsai method \n");
  else
  {
    eMc.insert(eRc);
    err = vpHandEyeCalibration::calibrationTranslationOld(cMo, rMe, eRc, eTc);
    if (err != 0) printf("\n Problem in solving Hand-Eye Translation by Old Tsai method after Old Tsai method for Rotation\n");
    else
    {
      eMc.insert(eTc);
#if DEBUG_LEVEL1
      {
        printf("\nRotation by (old) Tsai, old implementation for translation ori\n");
        vpHandEyeCalibration::calibrationVerifrMo(cMo, rMe, eMc);
      }
#endif
      resPos = vpHandEyeCalibration::calibrationErrVVS(cMo, rMe, eMc, errVVS);
      if (resPos < vmin)
      {
        vmin = resPos;
        eMcMin = eMc;
#if DEBUG_LEVEL1
        He_method = HE_TSAI_OROT;
#endif
      }
    }
    err = vpHandEyeCalibration::calibrationTranslation(cMo, rMe, eRc, eTc);
    if (err != 0) printf("\n Problem in solving Hand-Eye Translation after Old Tsai method for Rotation\n");
    else
    {
      eMc.insert(eTc);
#if 0
        printf("\nRotation by (old) Tsai, new implementation for translation : ori\n");
        vpHandEyeCalibration::calibrationVerifrMo(cMo, rMe, eMc);
#endif  //  1

#if DEBUG_LEVEL1
      {
        printf("\nRotation by (old) Tsai, new implementation for translation\n");
        vpHandEyeCalibration::calibrationVerifrMo(cMo, rMe, eMc);
      }
#endif
      resPos = vpHandEyeCalibration::calibrationErrVVS(cMo, rMe, eMc, errVVS);
      if (resPos < vmin)
      {
        vmin = resPos;
        eMcMin = eMc;
#if DEBUG_LEVEL1
        He_method = HE_TSAI_ORNT;
#endif
      }
    }
  }
  // First method using Tsai formulation
  err = vpHandEyeCalibration::calibrationRotationTsai(cMo, rMe, eRc);
  if (err != 0) printf("\n Problem in solving Hand-Eye Rotation by Tsai method \n");
  else
  {
    eMc.insert(eRc);
    err = vpHandEyeCalibration::calibrationTranslationOld(cMo, rMe, eRc, eTc);
    if (err != 0) printf("\n Problem in solving Hand-Eye Translation by Old Tsai method after Tsai method for Rotation\n");
    else
    {
      eMc.insert(eTc);
#if DEBUG_LEVEL1
      {
        printf("\nRotation by Tsai, old implementation for translation : ori\n");
        vpHandEyeCalibration::calibrationVerifrMo(cMo, rMe, eMc);
      }
#endif
      resPos = vpHandEyeCalibration::calibrationErrVVS(cMo, rMe, eMc, errVVS);
      if (resPos < vmin)
      {
        vmin = resPos;
        eMcMin = eMc;
#if DEBUG_LEVEL1
        He_method = HE_TSAI_NROT;
#endif
      }
    }
    err = vpHandEyeCalibration::calibrationTranslation(cMo, rMe, eRc, eTc);
    if (err != 0) printf("\n Problem in solving Hand-Eye Translation after Tsai method for Rotation \n");
    else
    {
      eMc.insert(eTc);

#if 0
    printf("\nRotation by Tsai, new implementation for translation : ori\n");
    vpHandEyeCalibration::calibrationVerifrMo(cMo, rMe, eMc);
#endif //   1

#if DEBUG_LEVEL1
      {
        printf("\nRotation by Tsai, new implementation for translation\n");
        vpHandEyeCalibration::calibrationVerifrMo(cMo, rMe, eMc);
      }
#endif
      resPos = vpHandEyeCalibration::calibrationErrVVS(cMo, rMe, eMc, errVVS);
      if (resPos < vmin)
      {
        vmin = resPos;
        eMcMin = eMc;
#if DEBUG_LEVEL1
        He_method = HE_TSAI_NRNT;
#endif
      }
    }
  }




  // Second method by solving the orthogonal Procrustes problem
  err = vpHandEyeCalibration::calibrationRotationProcrustes(cMo, rMe, eRc);
  if (err != 0) printf("\n Problem in solving Hand-Eye Rotation by Procrustes method \n");
  else
  {
    eMc.insert(eRc);
    err = vpHandEyeCalibration::calibrationTranslationOld(cMo, rMe, eRc, eTc);
    if (err != 0) printf("\n Problem in solving Hand-Eye Translation by Old Tsai method after Procrustes method for Rotation\n");
    else
    {
      eMc.insert(eTc);
#if DEBUG_LEVEL1
      {
        printf("\nRotation by Procrustes, old implementation for translation : ori\n");
        vpHandEyeCalibration::calibrationVerifrMo(cMo, rMe, eMc);
      }
#endif
      resPos = vpHandEyeCalibration::calibrationErrVVS(cMo, rMe, eMc, errVVS);
      if (resPos < vmin)
      {
        vmin = resPos;
        eMcMin = eMc;
#if DEBUG_LEVEL1
        He_method = HE_PROCRUSTES_OT;
#endif
      }
    }
    err = vpHandEyeCalibration::calibrationTranslation(cMo, rMe, eRc, eTc);
    if(err != 0) printf("\n Problem in solving Hand-Eye Translation after Procrustes method for Rotation\n");
    else
    {
      eMc.insert(eTc);
#if DEBUG_LEVEL1
      {
        printf("\nRotation by Procrustes, new implementation for translation : ori\n");
        vpHandEyeCalibration::calibrationVerifrMo(cMo, rMe, eMc);
      }
#endif
      resPos = vpHandEyeCalibration::calibrationErrVVS(cMo, rMe, eMc, errVVS);
      if (resPos < vmin)
      {
        eMcMin = eMc;
#if DEBUG_LEVEL1
        vmin = resPos;
        He_method = HE_PROCRUSTES_NT;
#endif
      }
    }
  }

  /* determination of the best method in case at least one succeeds */
  eMc = eMcMin;
#if 0
    std::cout << "AAA : ori" << std::endl;
    vpThetaUVector ePc(eMc);
    std::cout << "theta U (deg): " << vpMath::deg(ePc[0]) << " " << vpMath::deg(ePc[1]) << " " << vpMath::deg(ePc[2]) << std::endl;
    std::cout << "Translation: " << eMc[0][3] << " " << eMc[1][3] << " " << eMc[2][3] << std::endl;
#endif  //  1

#if DEBUG_LEVEL1
  {
    if (He_method == HE_I) printf("Best method : I !!!, vmin = %lf\n",vmin);
    if (He_method == HE_TSAI_OROT) printf("Best method : TSAI_OROT, vmin = %lf\n",vmin);
    if (He_method == HE_TSAI_ORNT) printf("Best method : TSAI_ORNT, vmin = %lf\n",vmin);
    if (He_method == HE_TSAI_NROT) printf("Best method : TSAI_NROT, vmin = %lf\n",vmin);
    if (He_method == HE_TSAI_NRNT) printf("Best method : TSAI_NRNT, vmin = %lf\n",vmin);
    if (He_method == HE_PROCRUSTES_OT) printf("Best method : PROCRUSTES_OT, vmin = %lf\n",vmin);
    if (He_method == HE_PROCRUSTES_NT) printf("Best method : PROCRUSTES_NT, vmin = %lf\n",vmin);
    vpThetaUVector ePc(eMc);
    std::cout << "theta U (deg): " << vpMath::deg(ePc[0]) << " " << vpMath::deg(ePc[1]) << " " << vpMath::deg(ePc[2]) << std::endl;
    std::cout << "Translation: " << eMc[0][3] << " " << eMc[1][3] << " " << eMc[2][3] << std::endl;
  }
#endif

  // Non linear iterative minimization to estimate simultaneouslty eRc and eTc
  err = vpHandEyeCalibration::calibrationVVS(cMo, rMe, eMc);
  // FC : err : 0 si tout OK, -1 si pb de rang, 1 si pas convergence
  if (err != 0) printf("\n Problem in solving Hand-Eye Calibration by VVS \n");
  else
  {
    printf("\nRotation and translation after VVS : ori\n");
    vpHandEyeCalibration::calibrationVerifrMo(cMo, rMe, eMc);
  }
  return err;
}

#ifdef HE_WO_VISP

#define PI              3.14159265
#define RAD_MIN_SINC    1.0e-8
#define RAD_MIN_MC      2.5e-4
#define NB_ITER_MAX     30

void cout_indented(int n_space, const std::string& str)
{
    if(n_space >= 0) std::cout << std::string(n_space * 2, ' ') << str << std::endl;
    }


std::string mat_type_2_str(int type, int n_sp) 
{
    cout_indented(n_sp, "mat_type_2_str");
        std::string r;
            uchar depth = type & CV_MAT_DEPTH_MASK;
                uchar chans = 1 + (type >> CV_CN_SHIFT);
                    switch ( depth ) {
                            case CV_8U:  r = "8U"; break;
                                    case CV_8S:  r = "8S"; break;
                                            case CV_16U: r = "16U"; break;
                                                    case CV_16S: r = "16S"; break;
                                                            case CV_32S: r = "32S"; break;
                                                                    case CV_32F: r = "32F"; break;
                                                                            case CV_64F: r = "64F"; break;
                                                                                    default:     r = "User"; break;
                                                                                        }
                                                                                            r += "C";
                                                                                                r += (chans + '0');
                                                                                                    return r;
                                                                                                    }



void print_mat_type(const Mat& mat, int n_sp)
{
    cout_indented(n_sp, "print_mat_type");
        cout_indented(n_sp + 1, "mat type : " + mat_type_2_str(mat.type(), n_sp + 1));      
        }

double vpHandEyeCalibration::rad2deg(double radian)
{
    return radian * 180.0 / PI;
    }



double sum_of_squared(const cv::Mat& m)
{
    double sqrt_sos = cv::norm(m, cv::NORM_L2);
        return sqrt_sos * sqrt_sos;
        }

double sqrt_of_mean_square(const cv::Mat& m)
{
    double sos = sum_of_squared(m);
        double ms = sos / m.total();
            return sqrt(ms);        
            }

// Checks if a matrix is a valid rotation matrix.
bool is_rotation_matrix(const cv::Mat &R)
{
    cv::Mat Rt;
        cv::transpose(R, Rt);
            cv::Mat shouldBeIdentity = Rt * R;
                cv::Mat I = Mat::eye(3, 3, shouldBeIdentity.type());    
                    return  cv::norm(I, shouldBeIdentity) < 1e-6;    
                    }


cv::Mat combine_rotation_translation_into_homogeneous_matrix(const cv::Mat& rot, const cv::Mat& tra)
{
    cv::Mat rot_mat;
        bool is_rot_valid = true;
            int n_elem_rot = rot.rows * rot.cols;
                if(3 == n_elem_rot)
                    {
                            cv::Rodrigues(rot, rot_mat);
                                }
                                    else if(9 == n_elem_rot)
                                        {
                                                if(is_rotation_matrix(rot)) rot_mat = rot;
                                                        else is_rot_valid = false;
                                                            }
                                                                else is_rot_valid = false;
                                                                    if(!is_rot_valid) { std::cout << "The given rotation matirx is NOT a real rotation matrix." << std::endl;        exit(0); }
                                                                        cv::Mat mat_homo = cv::Mat::eye(4, 4, CV_64F);
                                                                            rot_mat.copyTo(mat_homo(cv::Rect(0, 0, 3, 3)));
                                                                                tra.copyTo(mat_homo(cv::Rect(3, 0, 1, 3)));
                                                                                    return mat_homo;
                                                                                    }


std::pair<cv::Mat, cv::Mat> vpHandEyeCalibration::split_homogeneous_transform_matrix_into_rotation_and_translation(const Mat& mat_homo, int n_sp)
{
    //cout_indented(n_sp, "split_homogeneous_transform_matrix_into_rotation_and_translation START");
    if(!(mat_homo.rows >= 3 && mat_homo.cols >= 3 && mat_homo.cols == mat_homo.rows))
        {
                std::cout << "The given matrix is NOT homogeneous transformation matrix." << std::endl; exit(0);
                    }
                        int dim_rot = mat_homo.rows - 1;
                            cv::Mat mat_rot = mat_homo(cv::Rect(0, 0, dim_rot, dim_rot)).clone();
                                cv::Mat vec_tra = mat_homo(cv::Rect(dim_rot, 0, 1, dim_rot)).clone();
    //cout_indented(n_sp, "split_homogeneous_transform_matrix_into_rotation_and_translation END");
                                    return std::pair<cv::Mat, cv::Mat>(mat_rot, vec_tra);
                                    }


double mcosc_wo_visp(double cosx, double x)
{
  //if (fabs(x) < ang_min_mc)
  if (fabs(x) < RAD_MIN_MC)
      return 0.5;
        else
            return ((1.0 - cosx) / x / x);
            }



double msinc_wo_visp(double sinx, double x)
{
  //if (fabs(x) < ang_min_mc)
  if (fabs(x) < RAD_MIN_MC)
      return (1. / 6.0);
        else
            return ((1.0 - sinx / x) / x / x);
            }

double sinc_wo_visp(double sinx, double rad)
{
  if (fabs(rad) < RAD_MIN_SINC)
      return 1.0;
        else
            return sinx / rad;
            }



double sinc_wo_visp(double rad)
{
  if (fabs(rad) < RAD_MIN_SINC)
      return 1.0;
        else
            return sin(rad) / rad;
            }


void svdOpenCV_wo_visp(const Mat& A, Mat& U, Mat &w, Mat &V)
{
    std::cout << "VISP_HAVE_OPENCV_VERSION : " << VISP_HAVE_OPENCV_VERSION << std::endl;
    cv::Mat m;
    A.convertTo(m, CV_64F);
    //int rows = U.rows, cols = U.cols;
    //cv::Mat m(rows, cols, CV_64F, this->data);
    cv::SVD opencvSVD(m);
    V = opencvSVD.vt.t();
    w = opencvSVD.w;
    U = opencvSVD.u;
    return;
}





//void compute_pseudo_inverse(const vpMatrix &U, const vpColVector &sv, const vpMatrix &V, unsigned int nrows_orig, unsigned int ncols_orig, double svThreshold, vpMatrix &Ap, unsigned int &rank, vpMatrix &imA, vpMatrix &imAt, vpMatrix &kerAt)
void compute_pseudo_inverse_wo_visp(const Mat &U, const Mat &sv, const Mat &V, unsigned int nrows_orig, unsigned int ncols_orig, double svThreshold, Mat &Ap, unsigned int &rank, Mat &imA, Mat &imAt, Mat &kerAt)
{
    //Ap.resize(ncols_orig, nrows_orig, true);
    Ap = Mat::zeros(ncols_orig, nrows_orig, CV_64F);
    // compute the highest singular value and the rank of h
    //double maxsv = fabs(sv[0]);
    double maxsv = sv.at<double>(0);
    rank = 0;
    for(unsigned int i = 0; i < ncols_orig; i++) {
        //if (fabs(sv[i]) > maxsv * svThreshold) {
        if(fabs(sv.at<double>(i)) > maxsv * svThreshold) {
            rank++;
        }
        for(unsigned int j = 0; j < nrows_orig; j++) {
            //      Ap[i][j] = 0.0;
            for(unsigned int k = 0; k < ncols_orig; k++) {
                //if (fabs(sv[k]) > maxsv * svThreshold) {
                double sv_k = sv.at<double>(k);
                if(fabs(sv_k) > maxsv * svThreshold) {
                    //Ap[i][j] += V[i][k] * U[j][k] / sv[k];
                    Ap.at<double>(i, j) += V.at<double>(i, k) * U.at<double>(j, k) / sv_k;
                }
            }
        }
    }
    // Compute im(A) and im(At)
    //imA.resize(nrows_orig, rank);
    //imAt.resize(ncols_orig, rank);
    //for (unsigned int i = 0; i < nrows_orig; i++) {
    //    for (unsigned int j = 0; j < rank; j++) {
    //        imA[i][j] = U[i][j];
    //    }
    //}
    U(Rect(0, 0, rank, nrows_orig)).copyTo(imA);  
    //for (unsigned int i = 0; i < ncols_orig; i++) {
    //    for (unsigned int j = 0; j < rank; j++) {
    //        imAt[i][j] = V[i][j];
    //    }
    //}
    V(Rect(0, 0, rank, ncols_orig)).copyTo(imAt);  

    //kerAt.resize(ncols_orig - rank, ncols_orig);
    if(rank != ncols_orig) {
        kerAt = Mat::zeros(ncols_orig - rank, ncols_orig, CV_64F);
        for (unsigned int j = 0, k = 0; j < ncols_orig; j++) {
            // // if( v.col(j) in kernel and non zero )
            //if((fabs(sv[j]) <= maxsv * svThreshold) && (std::fabs(V.getCol(j).sumSquare()) > std::numeric_limits<double>::epsilon())) {
            if(fabs(sv.at<double>(j)) <= maxsv * svThreshold) {
                double sos = sum_of_squared(V.col(j)); 
                if(sos > std::numeric_limits<double>::epsilon())
                {
                    //for (unsigned int i = 0; i < V.getRows(); i++) {
                    for(unsigned int i = 0; i < V.rows; i++) {
                        //kerAt[k][i] = V[i][j];
                        kerAt.at<double>(k, i) = V.at<double>(i, j);
                    }
                    k++;
                }
            }
        }
    }     
}



void compute_pseudo_inverse_wo_visp(const Mat &a, const Mat &sv, const Mat &v, unsigned int nrows, unsigned int ncols, unsigned int nrows_orig, unsigned int ncols_orig, double svThreshold, Mat &Ap, unsigned int &rank)
{
    Mat a1 = Mat::zeros(ncols, nrows, a.type());
    // compute the highest singular value and the rank of h
    double maxsv = 0;
    for (unsigned int i = 0; i < ncols; i++) {
        double sv_abs = fabs(sv.at<double>(i));
        if(sv_abs > maxsv) maxsv = sv_abs;
    }
    rank = 0;
    for (unsigned int i = 0; i < ncols; i++) {
        if(fabs(sv.at<double>(i)) > maxsv * svThreshold) rank++;
        for(unsigned int j = 0; j < nrows; j++) {
            for (unsigned int k = 0; k < ncols; k++) {
                if(fabs(sv.at<double>(k)) > maxsv * svThreshold) {
                    a1.at<double>(i, j) += v.at<double>(i, k) * a.at<double>(j, k) / sv.at<double>(k);
                }
            }
        }
    }
    Ap = nrows_orig >= ncols_orig ? a1 : a1.t();
    return;
}


//unsigned int vpMatrix::pseudoInverseOpenCV(vpMatrix &Ap, vpColVector &sv, double svThreshold, vpMatrix &imA,vpMatrix &imAt, vpMatrix &kerA) const
unsigned int pseudoInverse_wo_visp(const Mat& A, Mat &Ap, Mat &sv, double svThreshold, Mat &imA, Mat &imAt, Mat &kerA) 
{
    //std::cout << "VISP_HAVE_OPENCV : " << VISP_HAVE_OPENCV_VERSION << std::endl << std::endl;
    unsigned int rank, nrows = A.rows, ncols = A.cols;
    //vpMatrix U, V;  vpColVector sv_;
    Mat A_dummy, U, V, sv_;
    if(nrows < ncols) {
        //U.resize(ncols, ncols, true);   sv.resize(nrows, false);
        A_dummy = Mat::zeros(ncols, ncols, A.type());  sv = Mat::zeros(nrows, 1, A.type());
    } 
    else {
        //U.resize(nrows, ncols, false);  sv.resize(ncols, false);
        A_dummy = Mat::zeros(nrows, ncols, A.type());  sv = Mat::zeros(ncols, 1, A.type());
    }
    //U.insert(*this, 0, 0);
    A.copyTo(A_dummy(Rect(0, 0, ncols, nrows)));
    //U.svdOpenCV(sv_, V);
    svdOpenCV_wo_visp(A_dummy, U, sv_, V);
    //compute_pseudo_inverse(U, sv_, V, nrows, ncols, svThreshold, Ap, rank, imA, imAt, kerA);
    compute_pseudo_inverse_wo_visp(U, sv_, V, nrows, ncols, svThreshold, Ap, rank, imA, imAt, kerA);
    // Remove singular values equal to to that correspond to the lines of 0
    // introduced when m < n
    //for(unsigned int i = 0; i < sv.size(); i++) sv[i] = sv_[i];
    sv_(Rect(0, 0, 1, sv.rows)).copyTo(sv);
    return rank;
}


//unsigned int vpMatrix::pseudoInverseOpenCV(vpMatrix &Ap, double svThreshold) const
unsigned int pseudoInverse_wo_visp(const Mat& A, Mat &Ap, double svThreshold)
{
    //std::cout << "VISP_HAVE_OPENCV : " << VISP_HAVE_OPENCV_VERSION << std::endl << std::endl;
    unsigned int rank, nrows, ncols, nrows_orig = A.rows, ncols_orig = A.cols;
    Ap = Mat(ncols_orig, nrows_orig, A.type());
    if (nrows_orig >= ncols_orig) {
        nrows = nrows_orig;
        ncols = ncols_orig;
    } 
    else {
        nrows = ncols_orig;
        ncols = nrows_orig;
    }
    Mat U, sv, V, A_dummy = nrows_orig >= ncols_orig ? A.clone() : A.t();
    svdOpenCV_wo_visp(A_dummy, U, sv, V);
    compute_pseudo_inverse_wo_visp(U, sv, V, nrows, ncols, nrows_orig, ncols_orig, svThreshold, Ap, rank);
    return rank;
}


template<typename T>
Mat
skewMat( const Mat_<T> &x )
{
  Mat_<T> skew(3,3);
    skew <<   0 , -x(2),  x(1),
              x(2),    0 , -x(0),
                       -x(1),  x(0),    0;

                         return std::move(skew);
                         }

                         Mat
                         skew_symmetric( InputArray _x )
                         {
                           const Mat x = _x.getMat();
                             const int depth = x.depth();
                               CV_Assert( x.size() == Size(3,1) || x.size() == Size(1,3) );
                                 CV_Assert( depth == CV_32F || depth == CV_64F );

                                   Mat skewMatrix;
                                     if( depth == CV_32F )
                                       {
                                           skewMatrix = skewMat<float>(x);
                                             }
                                               else if( depth == CV_64F )
                                                 {
                                                     skewMatrix = skewMat<double>(x);
                                                       }
                                                         else
                                                           {
                                                               //CV_Error(CV_StsBadArg, "The DataType must be CV_32F or CV_64F");
                                                                 }

                                                                   return skewMatrix;
                                                                   }



int vpHandEyeCalibration::calibrationRotationTsaiOld_wo_visp(const std::vector<Mat> &cMo, const std::vector<Mat> &rMe, Mat& eRc, int n_sp)
{
    cout_indented(n_sp, "calibrationRotationTsaiOld_wo_visp START");   
    unsigned int nbPose = (unsigned int) cMo.size();
    int n_combi = nbPose * (nbPose - 1) / 2; 
    Mat A(n_combi * 3, 3, CV_64F), B(n_combi * 3, 1, CV_64F), x;
    //unsigned int k = 0;
    // for all couples ij
    std::vector<Mat> li_rRe(nbPose), li_cRo(nbPose);
    for(unsigned int i = 0; i < nbPose; i++) {
        //std::cout << "rMe[i]: " << std::endl << rMe[i] << std::endl;
        //std::cout << "cMo[i]: " << std::endl << cMo[i] << std::endl;
        std::pair<Mat, Mat> rot_tra_re = split_homogeneous_transform_matrix_into_rotation_and_translation(rMe[i], n_sp + 1);
        li_rRe[i] = rot_tra_re.first;
        std::pair<Mat, Mat> rot_tra_co = split_homogeneous_transform_matrix_into_rotation_and_translation(cMo[i], n_sp + 1);
        li_cRo[i] = rot_tra_co.first;
    }
    for (unsigned int k = 0, i = 0; i < nbPose; i++) {
        Mat rRei = li_rRe[i], ciRo = li_cRo[i];
#if 0
        std::cout << std::endl << "i : " << i << std::endl;
        std::cout << "rRei w/o visp : " << std::endl << rRei << std::endl;
        std::cout << "ciRo w/o visp : " << std::endl << ciRo << std::endl;
#endif  //  1
        for (unsigned int j = i + 1; j < nbPose; j++, k += 3) {
            Mat cijPo, rPeij, rRej = li_rRe[j], cjRo = li_cRo[j];

            Mat rReij = rRej.t() * rRei, cijRo = cjRo * ciRo.t();
            Rodrigues(rReij, rPeij);
            double theta_re = norm(rPeij);
#if 0            
            std::cout << std::endl << "calibrationRotationTsaiOld w/o visp" << std::endl;
            std::cout << i << " " << j << " " << "ejRei: " << std::endl << rReij << std::endl;
            std::cout << "theta (robot) " << theta_re << std::endl;
            std::cout << "theta U (robot) " << rPeij << std::endl;
            std::cout << "cjRci: " << std::endl << cijRo.t() << std::endl;
#endif  //  1            
            rPeij *= sinc_wo_visp(theta_re / 2.0);

            Rodrigues(cijRo, cijPo);
            double theta_co = norm(cijPo);
            cijPo *= sinc_wo_visp(theta_co / 2.0);
            // std::cout << "theta (camera) " << theta_co << std::endl;
            // std::cout << "theta U (camera) " << cijPo.t() << std::endl;
            //Mat As = cv::sfm::skew(rPeij + cijPo);
            Mat As = skew_symmetric(rPeij + cijPo);
            Mat b = cijPo - rPeij;
            As.copyTo(A(Rect(0, k, 3, 3)));
            b.copyTo(B(Rect(0, k, 1, 3))); 
        }
    }     
  
    //std::cout << "A w/o visp"  << std::endl << A << std::endl << "B w/o visp"  << std::endl << B << std::endl;

    // the linear system is defined
    // x = AtA^-1AtB is solved
    Mat Ap, AtA = A.t() * A;
    int rank = pseudoInverse_wo_visp(AtA, Ap, 1e-6);
    if(rank != 3) return -1;
    //print_mat_type(Ap, n_sp + 1);   print_mat_type(A.t(), n_sp + 1);    print_mat_type(B, n_sp + 1);
    x = Ap * A.t() * B;
    Mat x2 = x.clone(); /* pour calcul residu */
    //std::cout << "x2 wo visp : " << std::endl << x2 << std::endl << std::endl;
    // extraction of theta and U
    double d = sum_of_squared(x);
    //std::cout << "DDD tsai rot old" << std::endl << std::endl;
    x *= 2.0 / sqrt(1.0 + d); 
    //std::cout << "x *= 2.0 / sqrt(1.0 + d) wo visp : " << std::endl << x << std::endl << std::endl;
    //double d = x.sumSquare();
    //for (unsigned int i = 0; i < 3; i++) x[i] = 2 * x[i] / sqrt(1 + d);
    //double theta = sqrt(x.sumSquare()) / 2;
    double theta = cv::norm(x, cv::NORM_L2) / 2.0;
    theta = 2.0 * asin(theta);
  // if (theta !=0)
    if(std::fabs(theta) > std::numeric_limits<double>::epsilon()) {
        //for(unsigned int i = 0; i < 3; i++) x[i] *= theta / (2 * sin(theta / 2));
        x *= theta / (2.0 * sin(theta / 2.0));
    }
    else x = Mat::zeros(x.size(), x.type());

    // Building of the rotation matrix eRc
    //vpThetaUVector xP(x[0], x[1], x[2]);
    //eRc = vpRotationMatrix(xP);
    Rodrigues(x, eRc);

#if 0
    std::cout << "Rotation from Old Tsai method w/o visp" << std::endl;
    std::cout << "theta U (deg): " << rad2deg(x.at<double>(0)) << " " << rad2deg(x.at<double>(1)) << " " << rad2deg(x.at<double>(2)) << std::endl;
    // Residual
    Mat residual = A * x2 - B;
    // std::cout << "Residual: " << std::endl << residual << std::endl;
    // double res = sqrt(residual.sumSquare()/residual.getRows());
    double res = sqrt_of_mean_square(residual);
    printf("Mean residual (rotation) = %lf\n", res);
    //exit(0);
#endif

#if DEBUG_LEVEL2

    std::cout << "Rotation from Old Tsai method" << std::endl;
    std::cout << "theta U (deg): " << rad2deg(x.at<double>(0)) << " " << rad2deg(x.at<double>(1)) << " " << rad2deg(x.at<double>(2)) << std::endl;
    // Residual
    Mat residual = A * x2 - B;
    // std::cout << "Residual: " << std::endl << residual << std::endl;
    // double res = sqrt(residual.sumSquare()/residual.getRows());
    double res = sqrt(sum_of_squared(residual) / float(residual.rows));
    printf("Mean residual (rotation) = %lf\n", res);

#endif
    cout_indented(n_sp, "calibrationRotationTsaiOld_wo_visp END");   
    return 0;
}        

double vpHandEyeCalibration::calibrationErrVVS_wo_visp(const std::vector<Mat> &cMo, const std::vector<Mat> &rMe, const Mat &eMc, Mat& errVVS, int n_sp)
{
    cout_indented(n_sp, "calibrationErrVVS_wo_visp START");
    //errVVS.release();
    unsigned int nbPose = (unsigned int) cMo.size();
    int n_combi = nbPose * (nbPose - 1) / 2; 
    //std::cout << "n_combi : " << n_combi << std::endl;  exit(0);
    Mat errVVSRot(3, n_combi, CV_64F), errVVSTra(3, n_combi, CV_64F), I3 = Mat::eye(3, 3, CV_64F);
    errVVS = Mat::zeros(n_combi * 6, 1, CV_64F);
    std::pair<Mat, Mat> rot_tra = split_homogeneous_transform_matrix_into_rotation_and_translation(eMc, n_sp + 1);
    Mat eRc = rot_tra.first, eTc = rot_tra.second;
   
    //unsigned int k = 0;
    for (unsigned int k = 0, i = 0; i < nbPose; i++) {
        for (unsigned int j = i + 1; j < nbPose; j++, k++) {
            //print_mat_type(rMe[i], n_sp + 2);
            //print_mat_type(cMo[j], n_sp + 2);
            //Mat s(3, 1, CV_32F);
            Mat ejMei = rMe[j].inv() * rMe[i];
            Mat cjMci = cMo[j] * cMo[i].inv();

            std::pair<Mat, Mat> rot_tra_e = split_homogeneous_transform_matrix_into_rotation_and_translation(ejMei, n_sp + 2);
            Mat ejPei, ejRei = rot_tra_e.first, ejTei = rot_tra_e.second;
            Rodrigues(ejRei, ejPei);  

            std::pair<Mat, Mat> rot_tra_c = split_homogeneous_transform_matrix_into_rotation_and_translation(cjMci, n_sp + 2);
            Mat cjPci, cjRci = rot_tra_c.first, cjTci = rot_tra_c.second;
            Rodrigues(cjRci, cjPci);  

            // terms due to rotation
            //cout_indented(n_sp + 2, "AAA");
            //print_mat_type(eRc, n_sp + 2);
            //print_mat_type(cjMci, n_sp + 2);
            //print_mat_type(ejMei, n_sp + 2);
            //print_mat_type(cjRci, n_sp + 2);
            //print_mat_type(ejRei, n_sp + 2);
            //print_mat_type(cjPci, n_sp + 2);
            //print_mat_type(ejPei, n_sp + 2);
            Mat s_rot = eRc * cjPci - ejPei;
            //cout_indented(n_sp + 2, "BBB");
            s_rot.copyTo(errVVSRot(Rect(k, 0, 1, 3))); 
            s_rot.copyTo(errVVS(Rect(0, k * 6, 1, 3))); 
            //cout_indented(n_sp + 2, "CCC");
            // terms due to translation
            Mat s_tra = (ejRei - I3) * eTc - eRc * cjTci + ejTei;
            //cout_indented(n_sp + 2, "DDD");
            s_tra.copyTo(errVVSTra(Rect(k, 0, 1, 3))); 
            s_tra.copyTo(errVVS(Rect(0, k * 6 + 3, 1, 3))); 
            //cout_indented(n_sp + 2, "EEE");
#if 0         
            std::cout << "k : " << k << std::endl;
            std::cout << "ejRei : " << std::endl << ejRei << std::endl;
            std::cout << "ejTei : " << std::endl << ejTei << std::endl;
            std::cout << "s_rot : " << std::endl << s_rot << std::endl;
            std::cout << "s_tra : " << std::endl << s_tra << std::endl << std::endl;
            exit(0);
#endif            
        }
    } 
    //std::cout << "errVVSRot : " << std::endl << errVVSRot << std::endl;  exit(0); 
    Mat mat_dp_rot = errVVSRot.t() * errVVSRot, mat_dp_tra = errVVSTra.t() * errVVSTra;
    double sum_of_sq_err_rot = sum((errVVSRot.t() * errVVSRot).diag())[0], 
            sum_of_sq_err_tra = sum((errVVSTra.t() * errVVSTra).diag())[0];
    double resRot = sqrt(sum_of_sq_err_rot / (n_combi * 3)), 
            resTrans = sqrt(sum_of_sq_err_tra / (n_combi * 3)), 
            resPos = sqrt((sum_of_sq_err_rot + sum_of_sq_err_tra) / (n_combi * 6));
#if 0   
    printf("\ncalibrationErrVVS w/o visp\n");
    printf("Mean VVS residual - rotation (deg) = %lf\n", rad2deg(resRot));
    printf("Mean VVS residual - translation = %lf\n", resTrans);
    printf("Mean VVS residual - global = %lf\n", resPos);
#endif    
#if DEBUG_LEVEL1
    printf("Mean VVS residual - rotation (deg) = %lf\n", rad2deg(resRot));
    printf("Mean VVS residual - translation = %lf\n", resTrans);
    printf("Mean VVS residual - global = %lf\n", resPos);   
#endif
    cout_indented(n_sp, "calibrationErrVVS_wo_visp END");
    return resPos;   
}

int vpHandEyeCalibration::calibrationTranslationOld_wo_visp(const std::vector<Mat> &cMo, const std::vector<Mat> &rMe, Mat &eRc, Mat &eTc, int n_sp)
{
    cout_indented(n_sp, "calibrationTranslationOld_wo_visp START");
    //vpMatrix A; vpColVector B;
    // Building of the system for the translation estimation
    // for all couples ij
    //vpRotationMatrix I3;    I3.eye();
    //int k = 0;
    unsigned int nbPose = (unsigned int)cMo.size();
    int n_combi = nbPose * (nbPose - 1) / 2; 

    Mat I3 = Mat::eye(3, 3, CV_64F), A(n_combi * 3, 3, CV_64F), B(n_combi * 3, 1, CV_64F);
    std::vector<Mat> li_rRe(nbPose), li_cRo(nbPose), li_rTe(nbPose), li_cTo(nbPose);
    for(unsigned int i = 0; i < nbPose; i++) {
        std::pair<Mat, Mat> rot_tra_re = split_homogeneous_transform_matrix_into_rotation_and_translation(rMe[i], n_sp + 2);
        li_rRe[i] = rot_tra_re.first;   li_rTe[i] = rot_tra_re.second;   
        std::pair<Mat, Mat> rot_tra_co = split_homogeneous_transform_matrix_into_rotation_and_translation(cMo[i], n_sp + 2);
        li_cRo[i] = rot_tra_co.first;   li_cTo[i] = rot_tra_co.second;
    }

    for (unsigned int k = 0, i = 0; i < nbPose; i++) {
        Mat rRei = li_rRe[i], ciRo = li_cRo[i], rTei = li_rTe[i], ciTo = li_cTo[i];
        for (unsigned int j = i + 1; j < nbPose; j++, k += 3) {
            //cout_indented(n_sp + 1, "k : " + std::to_string(k));
            Mat cijPo, rPeij, rRej = li_rRe[j], cjRo = li_cRo[j], rTej = li_rTe[j], cjTo = li_cTo[j];

            //cout_indented(n_sp + 2, "AAA");
            Mat rReij = rRej.t() * rRei;
            //cout_indented(n_sp + 2, "BBB");
            Mat rTeij = rRej.t() * (rTej - rTei);
            //cout_indented(n_sp + 2, "CCC");
            Mat a = rReij - I3;
            //cout_indented(n_sp + 2, "DDD");
            Mat b = eRc * cjTo - rReij * eRc * ciTo + rTeij;
            //cout_indented(n_sp + 2, "EEE");
            a.copyTo(A(Rect(0, k, 3, 3)));
            b.copyTo(B(Rect(0, k, 1, 3))); 
            //cout_indented(n_sp + 2, "FFF");
        }
    }
    // the linear system is solved
    // x = AtA^-1AtB is solved
    Mat Ap, AtA = A.t() * A;
    //vpColVector AeTc;
    //int rank = AtA.pseudoInverse(Ap, 1e-6);
    int rank = pseudoInverse_wo_visp(AtA, Ap, 1e-6);
    if (rank != 3) return -1;
    eTc = Ap * A.t() * B;
#if 0
    printf("Old Hand-eye calibration w/o visp : ");
    std::cout << "Translation: " << eTc.at<double>(0) << " " << eTc.at<double>(1) << " " << eTc.at<double>(2) << std::endl;
    // residual
    Mat residual = A * eTc - B;
    // std::cout << "Residual: " << std::endl << residual << std::endl;
    double res = sqrt_of_mean_square(residual);
    printf("Mean residual (translation) = %lf\n",res);
    //exit(0);
#endif  //  1

#if DEBUG_LEVEL2
  {
    printf("Old Hand-eye calibration w/o visp : ");
    std::cout << "Translation: " << eTc.at<double>(0) << " " << eTc.at<double>(1) << " " << eTc.at<double>(2) << std::endl;
    // residual
    Mat residual = A * AeTc - B;
    // std::cout << "Residual: " << std::endl << residual << std::endl;
    double res = sqrt_of_mean_square(residual);
    printf("Mean residual (translation) = %lf\n",res);

  }
#endif
    cout_indented(n_sp, "calibrationTranslationOld_wo_visp END");
    return 0;
}

Mat mean_traslation_vector(const std::vector<Mat> li_tra)
{
    int iT, n_tra = li_tra.size();
    Mat meanT = li_tra[0].clone();
    for(iT = 1; iT < n_tra; iT++)
    {
        meanT += li_tra[iT];
    }
    meanT /= n_tra;
    return meanT;
}

Mat mean_rotation_matrix_moakher(const std::vector<Mat> &vec_R)
{
    int n_rot = vec_R.size();
    Mat meanR = Mat::zeros(3, 3, CV_64F);
    for(size_t i = 0; i < vec_R.size(); i++) meanR += vec_R[i];
    meanR /= static_cast<double>(n_rot);
    // Euclidean mean of the rotation matrix following Moakher's method (SIAM 2002)
    //vpMatrix M, U, V;   vpColVector sv;
    //meanR.pseudoInverse(M, sv, 1e-6, U, V);
    Mat M, U, V, sv, kerA;
    pseudoInverse_wo_visp(meanR, M, sv, 1e-6, U, V, kerA);
    //double det = sv[0]*sv[1]*sv[2];
    double det = sv.at<double>(0) * sv.at<double>(1) * sv.at<double>(2);
    if(det > 0) {
        meanR = U * V.t();
    }
    else {
        //vpMatrix D(3, 3);   D = 0.0;    D[0][0] = D[1][1] = 1.0; D[2][2] = -1;
        Mat D = Mat::eye(3, 3, CV_64F); D.at<double>(2, 2) = -1;
        meanR = U * D * V.t();
    }
    //R = meanR;
    //return R;
    return meanR;
}

void vpHandEyeCalibration::calibrationVerifrMo_wo_visp(const std::vector<Mat> &cMo, const std::vector<Mat> &rMe, const Mat &eMc, int n_sp)
{
    cout_indented(n_sp, "calibrationVerifrMo_wo_visp START");
    unsigned int nbPose = (unsigned int) cMo.size();
    std::vector<Mat> rTo(nbPose), rRo(nbPose);
    for (unsigned int i = 0; i < nbPose; i++) {
        Mat rMo = rMe[i] * eMc * cMo[i];
        std::pair<Mat, Mat> rot_tra_ro = split_homogeneous_transform_matrix_into_rotation_and_translation(rMo, n_sp + 1);
        rRo[i] = rot_tra_ro.first;  rTo[i] = rot_tra_ro.second;
    }
    Mat meanRot = mean_rotation_matrix_moakher(rRo), meanTrans = mean_traslation_vector(rTo);
#if 0
    std::cout << "Mean w/o visp" << std::endl;
    std::cout << "Translation: " << meanTrans.t() << std::endl;
    //vpThetaUVector P(meanRot);
    Mat P;  Rodrigues(meanRot, P);
    double rad_mean = norm(P);
    //std::cout << "Rotation : theta (deg) = " <<  vpMath::deg(sqrt(P.sumSquare())) << " Matrice : " << std::endl << meanRot  << std::endl;
    std::cout << "Rotation : theta (deg) = " << rad2deg(rad_mean) << " Matrice : " << std::endl << meanRot  << std::endl;
    //std::cout << "theta U (deg): " << vpMath::deg(P[0]) << " " << vpMath::deg(P[1]) << " " << vpMath::deg(P[2]) << std::endl;
    std::cout << "theta U (deg): " << rad2deg(P.at<double>(0)) << " " << rad2deg(P.at<double>(1)) << " " << rad2deg(P.at<double>(2)) << std::endl;
#endif  //  1

#if DEBUG_LEVEL2
  {
    std::cout << "Mean  " << std::endl;
    std::cout << "Translation: " << meanTrans.t() << std::endl;
    vpThetaUVector P(meanRot);
    std::cout << "Rotation : theta (deg) = " << vpMath::deg(sqrt(P.sumSquare())) << " Matrice : " << std::endl << meanRot  << std::endl;
    std::cout << "theta U (deg): " << vpMath::deg(P[0]) << " " << vpMath::deg(P[1]) << " " << vpMath::deg(P[2]) << std::endl;
  }
#endif

  // standard deviation, rotational part
    double resRot = 0.0;
    for(unsigned int i = 0; i < nbPose; i++) {
        //vpRotationMatrix R = meanRot.t() * rRo[i]; // Rm^T  Ri
        Mat R = meanRot.t() * rRo[i]; // Rm^T  Ri
        //vpThetaUVectorP(R);
        Mat P;  Rodrigues(R, P);
        // theta = Riemannian distance d(Rm,Ri)
        //double theta = sqrt(P.sumSquare());
        double theta = norm(P);
        //std::cout << "Distance theta between rMo(" << i << ") and mean (deg) = " << vpMath::deg(theta) << std::endl;
        std::cout << "Distance theta between rMo(" << i << ") and mean (deg) = " << rad2deg(theta) << std::endl;
    // Euclidean distance d(Rm,Ri) not used
    // theta = 2.0*sqrt(2.0)*sin(theta/2.0);
        resRot += theta * theta;
    }
    resRot = sqrt(resRot / (double)nbPose);
    //std::cout << "Mean residual rMo(" << nbPose << ") - rotation (deg) = " << vpMath::deg(resRot) << std::endl;
    std::cout << "Mean residual rMo(" << nbPose << ") - rotation (deg) = " << rad2deg(resRot) << std::endl;
    // standard deviation, translational part
    double resTrans = 0.0;
    for(unsigned int i = 0; i < nbPose; i++) {
        //vpColVector errTrans = ((vpColVector) rTo[i]) - meanTrans;
        //resTrans += errTrans.sumSquare();
        double norrm = norm(rTo[i], meanTrans, NORM_L2);
        resTrans += norrm * norrm;
        std::cout << "Distance d between rMo(" << i << ") and mean (m) = " << norrm << std::endl;
    }
    resTrans = sqrt(resTrans / (double)nbPose);
    std::cout << "Mean residual rMo(" << nbPose << ") - translation (m) = " << resTrans << std::endl;
    double resPos = (resRot * resRot + resTrans * resTrans) * nbPose;
    resPos = sqrt(resPos / (2 * nbPose));
    std::cout << "Mean residual rMo(" << nbPose << ") - global = " << resPos << std::endl;
    cout_indented(n_sp, "calibrationVerifrMo_wo_visp END");
}

int vpHandEyeCalibration::calibrationTranslation_wo_visp(const std::vector<Mat> &cMo, const std::vector<Mat> &rMe, Mat &eRc, Mat &eTc, int n_sp)
{
    cout_indented(n_sp, "calibrationTranslation_wo_visp START");
    //vpMatrix I3(3, 3);  I3.eye();
    //unsigned int k = 0;
    unsigned int nbPose = (unsigned int)cMo.size();
    int n_combi = nbPose * (nbPose - 1) / 2; 
    Mat I3 = Mat::eye(3, 3, CV_64F), A(n_combi * 3, 3, CV_64F), B(n_combi * 3, 1, CV_64F);
    //vpMatrix A(3*nbPose,3);    vpColVector B(3*nbPose);
    // Building of the system for the translation estimation
    // for all couples ij
    for (unsigned int k = 0, i = 0; i < nbPose; i++) {
//        for (unsigned int j = 0; j < nbPose; j++) {
//            if (j > i) { // we don't use two times same couples...
        for (unsigned int j = i + 1; j < nbPose; j++, k += 3) {
            //vpHomogeneousMatrix ejMei = rMe[j].inverse() * rMe[i];
            //vpHomogeneousMatrix cjMci = cMo[j] * cMo[i].inverse();
            Mat ejMei = rMe[j].inv() * rMe[i], cjMci = cMo[j] * cMo[i].inv();
            //vpRotationMatrix ejRei, cjRci;  vpTranslationVector ejTei, cjTci;
            //ejMei.extract(ejRei);   ejMei.extract(ejTei);
            //cjMci.extract(cjRci);   cjMci.extract(cjTci);
            std::pair<Mat, Mat> rot_tra_ee = split_homogeneous_transform_matrix_into_rotation_and_translation(ejMei, n_sp + 2);
            std::pair<Mat, Mat> rot_tra_cc = split_homogeneous_transform_matrix_into_rotation_and_translation(cjMci, n_sp + 2);
            Mat ejRei = rot_tra_ee.first, ejTei = rot_tra_ee.second, 
                cjRci = rot_tra_cc.first, cjTci = rot_tra_cc.second;
            //vpMatrix a = vpMatrix(ejRei) - I3;
            //vpTranslationVector b = eRc * cjTci - ejTei;
            Mat a = ejRei - I3, b = eRc * cjTci - ejTei;
            a.copyTo(A(Rect(0, k, 3, 3)));
            b.copyTo(B(Rect(0, k, 1, 3))); 
            //if (k == 0) {
            //  A = a;
            //  B = b;
            //} 
            //else {
            //    A = vpMatrix::stack(A, a);
            //    B = vpColVector::stack(B, b);
            //}
            //k++;
        }
    }
    //// the linear system A x = B is solved
    //// using x = A^+ B
    //vpMatrix Ap;
    //int rank = A.pseudoInverse(Ap);
    Mat Ap;
    int rank = pseudoInverse_wo_visp(A, Ap, 1e-6); 
    if(rank != 3) return -1;
    //vpColVector x = Ap * B;
    //eTc = (vpTranslationVector) x;
    eTc = Ap * B;
#if 0
    printf("New Hand-eye calibration w/o visp : ");
    std::cout << "Translation: " << eTc.at<double>(0) << " " << eTc.at<double>(1) << " " << eTc.at<double>(2) << std::endl;
    // residual
    //vpColVector residual;
    //residual = A*x-B;
    Mat residual = A * eTc - B;
    // std::cout << "Residual: " << std::endl << residual << std::endl;
    //double res = sqrt(residual.sumSquare()/residual.getRows());
    double res = sqrt_of_mean_square(residual);
    printf("Mean residual (translation) = %lf\n",res);
#endif  //  1
#if DEBUG_LEVEL2
  {
    printf("New Hand-eye calibration : ");
    std::cout << "Translation: " << eTc[0] << " " << eTc[1] << " " << eTc[2] << std::endl;
    // residual
    vpColVector residual;
    residual = A*x-B;
    // std::cout << "Residual: " << std::endl << residual << std::endl;
    double res = sqrt(residual.sumSquare()/residual.getRows());
    printf("Mean residual (translation) = %lf\n",res);
  }
#endif
    cout_indented(n_sp, "calibrationTranslation_wo_visp END");
  return 0;
}


int vpHandEyeCalibration::calibrationRotationTsai_wo_visp(const std::vector<Mat> &cMo, const std::vector<Mat> &rMe, Mat &eRc, int n_sp)
{
    cout_indented(n_sp, "calibrationRotationTsai_wo_visp START");   
    //vpMatrix A;
    //vpColVector B;
    unsigned int nbPose = (unsigned int) cMo.size();
    int n_combi = nbPose * (nbPose - 1) / 2; 
    Mat A(n_combi * 3, 3, CV_64F), B(n_combi * 3, 1, CV_64F);
    //unsigned int k = 0;
    // for all couples ij
    std::vector<Mat> li_rRe(nbPose), li_cRo(nbPose);
    for(unsigned int i = 0; i < nbPose; i++) {
        //vpRotationMatrix rRei, ciRo;
        //rMe[i].extract(rRei);        cMo[i].extract(ciRo);
        std::pair<Mat, Mat> rot_tra_re = split_homogeneous_transform_matrix_into_rotation_and_translation(rMe[i], n_sp + 1);
        li_rRe[i] = rot_tra_re.first;
        std::pair<Mat, Mat> rot_tra_co = split_homogeneous_transform_matrix_into_rotation_and_translation(cMo[i], n_sp + 1);
        li_cRo[i] = rot_tra_co.first;
        // std::cout << "rMei: " << std::endl << rMe[i] << std::endl;
    }

    for (unsigned int k = 0, i = 0; i < nbPose; i++) 
    {
        Mat rRei = li_rRe[i], ciRo = li_cRo[i];
        //for (unsigned int j = 0; j < nbPose; j++) {
        //    if (j > i) // we don't use two times same couples...
        for (unsigned int j = i + 1; j < nbPose; j++, k += 3) 
        {
            Mat cjPci, ejPei, rRej = li_rRe[j], cjRo = li_cRo[j];
            //vpRotationMatrix rRej, cjRo;
            //rMe[j].extract(rRej);            cMo[j].extract(cjRo);
            //// std::cout << "rMej: " << std::endl << rMe[j] << std::endl;
            //vpRotationMatrix ejRei = rRej.t() * rRei;
            //vpThetaUVector ejPei(ejRei);
            Mat ejRei = rRej.t() * rRei;
            Rodrigues(ejRei, ejPei); 

            //vpRotationMatrix cjRci = cjRo * ciRo.t();
            //vpThetaUVector cjPci(cjRci);
            Mat cjRci = cjRo * ciRo.t();
            Rodrigues(cjRci, cjPci);
            // std::cout << "theta U (camera) " << cjPci.t() << std::endl;

            //vpMatrix As;    vpColVector b(3);
            //As = vpColVector::skew(vpColVector(ejPei) + vpColVector(cjPci));
            //b =  (vpColVector)cjPci - (vpColVector) ejPei; // A.40
            Mat As = skew_symmetric(ejPei + cjPci), 
                b = cjPci - ejPei;
                
            As.copyTo(A(Rect(0, k, 3, 3)));
            b.copyTo(B(Rect(0, k, 1, 3))); 
            //if (k == 0) {
            //  A = As;
            //  B = b;
            //} else {
            //  A = vpMatrix::stack(A, As);
            //  B = vpColVector::stack(B, b);
            //}
            //k++;
        }
    }
#if DEBUG_LEVEL2
    {
        std::cout << "Tsai method: system A X = B "  << std::endl;
        std::cout << "A "  << std::endl << A << std::endl;
        std::cout << "B "  << std::endl << B << std::endl;
    }
#endif
    //// the linear system A x = B is solved
    //// using x = A^+ B
    //vpMatrix Ap;
    //int rank = A.pseudoInverse(Ap);
    Mat Ap;
    int rank = pseudoInverse_wo_visp(A, Ap, 1e-6);
    if(rank != 3) return -1;
    //vpColVector x = Ap * B;
    Mat x = Ap * B;
    // extraction of theta U
    // x = tan(theta/2) U
    //double norm =  x.sumSquare();
    double norm =  sum_of_squared(x);
    double c = 1.0 / sqrt(1.0 + norm);  // cos(theta/2)
    double alpha = acos(c);         // theta/2
    //norm = 2.0*c/vpMath::sinc(alpha);  // theta / tan(theta/2)
    norm = 2.0 * c / sinc_wo_visp(alpha);  // theta / tan(theta/2)
    //for (unsigned int i = 0; i < 3; i++) x[i] *= norm;
    Mat x_ori = x.clone();
    x *= norm;
    //// Building of the rotation matrix eRc
    //vpThetaUVector xP(x[0], x[1], x[2]);
    //eRc = vpRotationMatrix(xP);
    Rodrigues(x, eRc);
#if 0
    std::cout << "Rotation from Tsai method : w/o visp" << std::endl;
    std::cout << "theta U (deg): " << rad2deg(x.at<double>(0)) << " " << rad2deg(x.at<double>(1)) << " " << rad2deg(x.at<double>(2)) << std::endl;
    // Residual
    //for (unsigned int i = 0; i < 3; i++) x[i] /= norm; /* original x */
    Mat residual = A * x_ori - B;
    // std::cout << "Residual: " << std::endl << residual << std::endl;
    // double res = sqrt(residual.sumSquare()/residual.getRows());
    double res = sqrt_of_mean_square(residual);
    printf("Mean residual (rotation) = %lf\n", res);
    //exit(0);
#endif

#if DEBUG_LEVEL2
    {
        std::cout << "Rotation from Tsai method" << std::endl;
        std::cout << "theta U (deg): " << vpMath::deg(x[0]) << " " << vpMath::deg(x[1]) << " " << vpMath::deg(x[2]) << std::endl;
    // Residual
        for (unsigned int i = 0; i < 3; i++) x[i] /= norm; /* original x */
        vpColVector residual;
        residual = A*x-B;
        // std::cout << "Residual: " << std::endl << residual << std::endl;
        double res = sqrt(residual.sumSquare()/residual.getRows());
        printf("Mean residual (rotation) = %lf\n",res);
    }
#endif
    cout_indented(n_sp, "calibrationRotationTsai_wo_visp END");   
    return 0;
}


int vpHandEyeCalibration::calibrationRotationProcrustes_wo_visp(const std::vector<Mat> &cMo, const std::vector<Mat> &rMe, Mat &eRc, int n_sp)
{
    cout_indented(n_sp, "calibrationRotationProcrustes_wo_visp START");   
    // Method by solving the orthogonal Procrustes problem
    // [... (theta u)_e ...] = eRc [ ... (theta u)_c ...]
    // similar to E^T = eRc C^T below

    //vpMatrix Et, Ct;
    //vpMatrix A;
    //unsigned int k = 0;
    unsigned int nbPose = (unsigned int) cMo.size();
    int n_combi = nbPose * (nbPose - 1) / 2; 
    Mat Et(n_combi, 3, CV_64F), Ct(n_combi, 3, CV_64F);
    std::vector<Mat> li_rRe(nbPose), li_cRo(nbPose);
    for(unsigned int i = 0; i < nbPose; i++) {
        std::pair<Mat, Mat> rot_tra_re = split_homogeneous_transform_matrix_into_rotation_and_translation(rMe[i], n_sp + 1);
        li_rRe[i] = rot_tra_re.first;
        std::pair<Mat, Mat> rot_tra_co = split_homogeneous_transform_matrix_into_rotation_and_translation(cMo[i], n_sp + 1);
        li_cRo[i] = rot_tra_co.first;
    }

    // for all couples ij
    for (unsigned int k = 0, i = 0; i < nbPose; i++) {
        //vpRotationMatrix rRei, ciRo;
        //rMe[i].extract(rRei);   cMo[i].extract(ciRo);
        //// std::cout << "rMei: " << std::endl << rMe[i] << std::endl;
        Mat rRei = li_rRe[i], ciRo = li_cRo[i];
        for(unsigned int j = i + 1; j < nbPose; j++, k++) {
            //if (j > i) // we don't use two times same couples...
            //{
            //vpRotationMatrix rRej, cjRo;
            //rMe[j].extract(rRej);   cMo[j].extract(cjRo);
            //// std::cout << "rMej: " << std::endl << rMe[j] << std::endl;
            Mat cjPci, ejPei, rRej = li_rRe[j], cjRo = li_cRo[j];
            //vpRotationMatrix ejRei = rRej.t() * rRei;
            Mat ejRei = rRej.t() * rRei, cjRci = cjRo * ciRo.t();
            //vpThetaUVector ejPei(ejRei);
            Rodrigues(ejRei, ejPei);
            //vpColVector xe = ejPei;
            //vpRotationMatrix cjRci = cjRo * ciRo.t();
            //vpThetaUVector cjPci(cjRci);
            Rodrigues(cjRci, cjPci);
            //vpColVector xc = cjPci;

            //(ejPei.t()).copyTo(Et(Rect(0, k, 3, 1)));   (cjPci.t()).copyTo(Ct(Rect(0, k, 3, 1)));
            ejPei = ejPei.t();  ejPei.copyTo(Et(Rect(0, k, 3, 1)));  
            cjPci = cjPci.t();  cjPci.copyTo(Ct(Rect(0, k, 3, 1)));
            //if (k == 0) {
            //    Et = xe.t();
            //    Ct = xc.t();
            //} 
            //else {
            //    Et.stack(xe.t());
            //    Ct.stack(xc.t());
            //}
            //k++;
        }
    }
    //}
    // std::cout << "Et "  << std::endl << Et << std::endl;
    // std::cout << "Ct "  << std::endl << Ct << std::endl;

    // R obtained from the SVD of (E C^T) with all singular values equal to 1
    //A = Et.t() * Ct;
    //vpMatrix M, U, V;
    //vpColVector sv;
    //int rank = A.pseudoInverse(M, sv, 1e-6, U, V);
    Mat M, U, V, sv, kerA, A = Et.t() * Ct;
    int rank = pseudoInverse_wo_visp(A, M, sv, 1e-6, U, V, kerA);
    if (rank != 3) return -1;
    //A = U * V.t();
    //eRc = vpRotationMatrix(A);
    eRc = U * V.t();
#if 0
    //vpThetaUVector ePc(eRc);
    Mat ePc;    Rodrigues(eRc, ePc);
    std::cout << "Rotation from Procrustes method : w/o visp" << std::endl;
    std::cout << "theta U (deg): " << rad2deg(ePc.at<double>(0)) << " " << rad2deg(ePc.at<double>(1)) << " " << rad2deg(ePc.at<double>(2)) << std::endl;
    // Residual
    //vpMatrix residual;      residual = A * Ct.t() - Et.t();
    Mat residual = eRc * Ct.t() - Et.t();
    //  std::cout << "Residual: " << std::endl << residual << std::endl;
    //double res = sqrt(residual.sumSquare()/(residual.getRows()*residual.getCols()));
    double res = sqrt_of_mean_square(residual);
    printf("Mean residual (rotation) = %lf\n",res);
#endif  //  1

#if DEBUG_LEVEL2
  {
    vpThetaUVector ePc(eRc);
    std::cout << "Rotation from Procrustes method " << std::endl;
    std::cout << "theta U (deg): " << vpMath::deg(ePc[0]) << " " << vpMath::deg(ePc[1]) << " " << vpMath::deg(ePc[2]) << std::endl;
    // Residual
    vpMatrix residual;
    residual = A * Ct.t() - Et.t();
    //  std::cout << "Residual: " << std::endl << residual << std::endl;
    double res = sqrt(residual.sumSquare()/(residual.getRows()*residual.getCols()));
    printf("Mean residual (rotation) = %lf\n",res);
  }
#endif
    cout_indented(n_sp, "calibrationRotationProcrustes_wo_visp END");   
    return 0;
}




        
//  vpHomogeneousMatrix vpExponentialMap::direct(const vpColVector &v) { return vpExponentialMap::direct(v, 1.0); }
//vpHomogeneousMatrix vpExponentialMap::direct(const vpColVector &v, const double &delta_t)
Mat exponential_map_direct(const Mat &v, const double &delta_t, int n_sp)
{
    cout_indented(n_sp, "exponential_map_direct START");   
    //if(v.size() != 6) {
    //    throw(vpException(vpException::dimensionError, "Cannot compute direct exponential map from a %d-dim velocity vector. Should be 6-dim.", v.size()));
    //}
    if(!(6 == v.rows && 1 == v.cols)) throw std::length_error("Cannot compute direct exponential map from a " + std::to_string(v.rows) + " by " + std::to_string(v.cols) + ". Should be 6 by 1.");

    double theta, si, co, sinc, mcosc, msinc;
    //vpThetaUVector u;
    //vpRotationMatrix rd;
    //vpTranslationVector dt;
    //vpColVector v_dt = v * delta_t;
    Mat rd, u, dt(3, 1, CV_64F), v_dt = v * delta_t;
    //u[0] = v_dt[3]; u[1] = v_dt[4]; u[2] = v_dt[5];
    v_dt(Rect(0, 3, 1, 3)).copyTo(u);
    
    //rd.buildFrom(u);
    Rodrigues(u, rd);

    
    //theta = sqrt(u[0] * u[0] + u[1] * u[1] + u[2] * u[2]);
    theta = norm(u, NORM_L2);
    
    si = sin(theta);
    co = cos(theta);
    //sinc = vpMath::sinc(si, theta);
    sinc = sinc_wo_visp(si, theta);
    //mcosc = vpMath::mcosc(co, theta);
    mcosc = mcosc_wo_visp(co, theta);
    //msinc = vpMath::msinc(si, theta);
    msinc = msinc_wo_visp(si, theta);

    double  u_0 =       u.at<double>(0),    u_1 =       u.at<double>(1),    u_2 =       u.at<double>(2), 
            v_dt_0 =    v_dt.at<double>(0), v_dt_1 =    v_dt.at<double>(1), v_dt_2 =    v_dt.at<double>(2);


    //dt[0] = v_dt[0] * (sinc + u[0] * u[0] * msinc) + 
    //        v_dt[1] * (u[0] * u[1] * msinc - u[2] * mcosc) + 
    //        v_dt[2] * (u[0] * u[2] * msinc + u[1] * mcosc);
    dt.at<double>(0) =  v_dt_0 * (sinc + u_0 * u_0 * msinc) +
                        v_dt_1 * (u_0 * u_1 * msinc - u_2 * mcosc) + 
                        v_dt_2 * (u_0 * u_2 * msinc + u_1 * mcosc);
    //dt[1] =   v_dt[0] * (u[0] * u[1] * msinc + u[2] * mcosc) + 
    //          v_dt[1] * (sinc + u[1] * u[1] * msinc) +
    //          v_dt[2] * (u[1] * u[2] * msinc - u[0] * mcosc);
    dt.at<double>(1) =  v_dt_0 * (u_0 * u_1 * msinc + u_2 * mcosc) + 
                        v_dt_1 * (sinc + u_1 * u_1 * msinc) +
                        v_dt_2 * (u_1 * u_2 * msinc - u_0 * mcosc);

    //dt[2] = v_dt[0] * (u[0] * u[2] * msinc - u[1] * mcosc) + 
    //        v_dt[1] * (u[1] * u[2] * msinc + u[0] * mcosc) +
    //        v_dt[2] * (sinc + u[2] * u[2] * msinc);

    dt.at<double>(2) =  v_dt_0 * (u_0 * u_2 * msinc - u_1 * mcosc) +
                        v_dt_1 * (u_1 * u_2 * msinc + u_0 * mcosc) +
                        v_dt_2 * (sinc + u_2 * u_2 * msinc);

    //vpHomogeneousMatrix Delta;
    //Delta.insert(rd);
    //Delta.insert(dt);
    Mat Delta = combine_rotation_translation_into_homogeneous_matrix(rd, dt);
#if 0
    //if(0) // test new version wrt old version
    if(1) // test new version wrt old version
    {
        // old version
        unsigned int i, j;
        double s;
        // double u[3];
        //  vpRotationMatrix rd ;
        //  vpTranslationVector dt ;
        s = sqrt(v_dt[3] * v_dt[3] + v_dt[4] * v_dt[4] + v_dt[5] * v_dt[5]);
        if (s > 1.e-15) {
            for (i = 0; i < 3; i++)
                u[i] = v_dt[3 + i] / s;
            double sinu = sin(s);
            double cosi = cos(s);
            double mcosi = 1 - cosi;
            rd[0][0] = cosi + mcosi * u[0] * u[0];
            rd[0][1] = -sinu * u[2] + mcosi * u[0] * u[1];
            rd[0][2] = sinu * u[1] + mcosi * u[0] * u[2];
            rd[1][0] = sinu * u[2] + mcosi * u[1] * u[0];
            rd[1][1] = cosi + mcosi * u[1] * u[1];
            rd[1][2] = -sinu * u[0] + mcosi * u[1] * u[2];
            rd[2][0] = -sinu * u[1] + mcosi * u[2] * u[0];
            rd[2][1] = sinu * u[0] + mcosi * u[2] * u[1];
            rd[2][2] = cosi + mcosi * u[2] * u[2];
            dt[0] = v_dt[0] * (sinu / s + u[0] * u[0] * (1 - sinu / s)) + 
                    v_dt[1] * (u[0] * u[1] * (1 - sinu / s) - u[2] * mcosi / s) +
                    v_dt[2] * (u[0] * u[2] * (1 - sinu / s) + u[1] * mcosi / s);
            dt[1] = v_dt[0] * (u[0] * u[1] * (1 - sinu / s) + u[2] * mcosi / s) +
                    v_dt[1] * (sinu / s + u[1] * u[1] * (1 - sinu / s)) +
                    v_dt[2] * (u[1] * u[2] * (1 - sinu / s) - u[0] * mcosi / s);
            dt[2] = v_dt[0] * (u[0] * u[2] * (1 - sinu / s) - u[1] * mcosi / s) +
                    v_dt[1] * (u[1] * u[2] * (1 - sinu / s) + u[0] * mcosi / s) +
                    v_dt[2] * (sinu / s + u[2] * u[2] * (1 - sinu / s));
        } 
        else {
            for (i = 0; i < 3; i++) {
                for (j = 0; j < 3; j++)
                    rd[i][j] = 0.0;
                rd[i][i] = 1.0;
                dt[i] = v_dt[i];
            }
        }
        // end old version

        // Test of the new version
        vpHomogeneousMatrix Delta_old;
        Delta_old.insert(rd);
        Delta_old.insert(dt);
        int pb = 0;
        for (i = 0; i < 4; i++) {
            for (j = 0; j < 4; j++)
                if (fabs(Delta[i][j] - Delta_old[i][j]) > 1.e-5)
                    pb = 1;
        }
        if (pb == 1) {
            printf("pb vpHomogeneousMatrix::expMap\n");
            std::cout << " Delta : " << std::endl << Delta << std::endl;
            std::cout << " Delta_old : " << std::endl << Delta_old << std::endl;
        }
        // end of the test
    }
#endif  //  0    
    cout_indented(n_sp, "exponential_map_direct END");   
    return Delta;
}





int vpHandEyeCalibration::calibrationVVS_wo_visp(const std::vector<Mat> &cMo, const std::vector<Mat> &rMe, Mat &eMc, int n_sp)
{
    cout_indented(n_sp, "calibrationVVS_wo_visp START");   
    unsigned int it = 0;
    double res = 1.0;
    unsigned int nbPose = (unsigned int) cMo.size();
    int n_combi = nbPose * (nbPose - 1) / 2; 
    //vpColVector err;
    //vpMatrix L;
    //vpMatrix I3(3,3);   I3.eye();
    Mat err, L = Mat::zeros(n_combi * 6, 6, CV_64F), I3 = Mat::eye(3, 3, CV_64F);
    //vpRotationMatrix eRc;   vpTranslationVector eTc;
    //eMc.extract(eRc);   eMc.extract(eTc);
    std::pair<Mat, Mat> rot_tra = split_homogeneous_transform_matrix_into_rotation_and_translation(eMc, n_sp + 1);
    Mat eRc = rot_tra.first, eTc = rot_tra.second;
  /* FC : on recalcule 2 fois tous les ejMei et cjMci a chaque iteration
    alors qu'ils sont constants. Ce serait sans doute mieux de les
    calculer une seule fois et de les stocker. Pourraient alors servir
    dans les autres fonctions HandEye. A voir si vraiment interessant vu la
    combinatoire. Idem pour les theta u */
    while ((res > 1e-7) && (it < NB_ITER_MAX))
    {
        /* compute s - s^* */
        //vpHandEyeCalibration::calibrationErrVVS(cMo, rMe, eMc, err);
        vpHandEyeCalibration::calibrationErrVVS_wo_visp(cMo, rMe, eMc, err, n_sp + 1);
        /* compute L_s */
        //unsigned int k = 0;
        //for (unsigned int i = 0; i < nbPose; i++) {
        //   for (unsigned int j = 0; j < nbPose; j++) {
        //        if (j > i) // we don't use two times same couples...
        //        {
        
        for(unsigned int k = 0, i = 0; i < nbPose; i++) {
            for(unsigned int j = i + 1; j < nbPose; j++, k += 6) {
                //vpMatrix Ls(3,6), Lv(3,3), Lw(3,3);
                //vpHomogeneousMatrix ejMei = rMe[j].inverse() * rMe[i];
                //vpHomogeneousMatrix cjMci = cMo[j] * cMo[i].inverse();
                Mat ejMei = rMe[j].inv() * rMe[i], cjMci = cMo[j] * cMo[i].inv();
                //vpRotationMatrix ejRei;
                //ejMei.extract(ejRei);
                //vpThetaUVector cjPci(cjMci);
                //vpTranslationVector cjTci;
                //cjMci.extract(cjTci);
                std::pair<Mat, Mat> rot_tra_ee = split_homogeneous_transform_matrix_into_rotation_and_translation(ejMei, n_sp + 2);
                std::pair<Mat, Mat> rot_tra_cc = split_homogeneous_transform_matrix_into_rotation_and_translation(cjMci, n_sp + 2);
                Mat cjPci, ejRei = rot_tra_ee.first, cjRci = rot_tra_cc.first, cjTci = rot_tra_cc.second;
                Rodrigues(cjRci, cjPci);
                // terms due to rotation
                //Lv.diag(0.0); //
                //Lv = 0.0;
                //Lw = -vpMatrix(eRc) * vpColVector::skew(vpColVector(cjPci));
                Mat Lw_rot = -eRc * skew_symmetric(cjPci); 
                //for (unsigned int m=0;m<3;m++)
                //    for (unsigned int n=0;n<3;n++)
                //    {
                //        Ls[m][n] = Lv[m][n];
                //        Ls[m][n+3] = Lw[m][n];
                //    }
                Lw_rot.copyTo(L(Rect(3, k, 3, 3)));     
                //if (k == 0) {
                //    L = Ls;
                //} 
                //else {
                //    L = vpMatrix::stack(L,Ls);
                //}
                //k++;
                // terms due to translation
                //Lv = (vpMatrix(ejRei) - I3) * vpMatrix(eRc);
                //Lw =  vpMatrix(eRc) * vpColVector::skew(vpColVector(cjTci));
                Mat Lv_tra = (ejRei - I3) * eRc, Lw_tra = eRc * skew_symmetric(cjTci);
                //for (unsigned int m=0;m<3;m++)
                //    for (unsigned int n=0;n<3;n++)
                //    {
                //        Ls[m][n] = Lv[m][n];
                //        Ls[m][n+3] = Lw[m][n];
                //    }
                //L = vpMatrix::stack(L,Ls);
                Lv_tra.copyTo(L(Rect(0, k + 3, 3, 3)));     
                Lw_tra.copyTo(L(Rect(3, k + 3, 3, 3)));     

                //} // enf if i > j
            } // end for j
        } // end for i
        double lambda = 0.9;
        //vpMatrix Lp;
        //int rank = L.pseudoInverse(Lp);
        Mat Lp;
        int rank = pseudoInverse_wo_visp(L, Lp, 1e-6);
        if (rank != 6) return -1;

        //vpColVector e = Lp * err;
        //vpColVector v = - e * lambda;

        Mat e = Lp * err;
        Mat v = -lambda * e;

        //  std::cout << "e: "  << e.t() << std::endl;
        //eMc = eMc * vpExponentialMap::direct(v);
        //  vpHomogeneousMatrix vpExponentialMap::direct(const vpColVector &v) { return vpExponentialMap::direct(v, 1.0); }
        eMc = eMc * exponential_map_direct(v, 1.0, n_sp + 1);
        //eMc.extract(eRc);        eMc.extract(eTc);
        std::pair<Mat, Mat> rot_tra = split_homogeneous_transform_matrix_into_rotation_and_translation(eMc, n_sp + 1);
        eRc = rot_tra.first, eTc = rot_tra.second;
        //res = sqrt(v.sumSquare()/v.getRows());
        res = sqrt_of_mean_square(v);
#if 0
        std::cout << "it : " << it << ",\t res : " << res << std::endl; 
#endif  //  1
        it++;
    } // end while

#if 0
    printf("Iteration number for NL hand-eye minimisation w/o visp : %d / %d\n", it, NB_ITER_MAX);
    //vpThetaUVector ePc(eRc);
    Mat ePc;    Rodrigues(eRc, ePc);
    std::cout << "theta U (deg): " << rad2deg(ePc.at<double>(0)) << " " << rad2deg(ePc.at<double>(1)) << " " << rad2deg(ePc.at<double>(2)) << std::endl;
    std::cout << "Translation: " << eTc.at<double>(0) << " " << eTc.at<double>(1) << " " << eTc.at<double>(2) << std::endl;
    // Residual
    //double res = err.sumSquare();   res = sqrt(res/err.getRows());
    res = sqrt_of_mean_square(err);
    printf("Mean residual (rotation+translation) = %lf\n", res);
#endif  //  1

#if DEBUG_LEVEL2
  {
    printf(" Iteration number for NL hand-eye minimisation : %d\n",it);
    vpThetaUVector ePc(eRc);
    std::cout << "theta U (deg): " << vpMath::deg(ePc[0]) << " " << vpMath::deg(ePc[1]) << " " << vpMath::deg(ePc[2]) << std::endl;
    std::cout << "Translation: " << eTc[0] << " " << eTc[1] << " " << eTc[2] << std::endl;
    // Residual
    double res = err.sumSquare();
    res = sqrt(res/err.getRows());
    printf("Mean residual (rotation+translation) = %lf\n",res);
  }
#endif
    cout_indented(n_sp, "calibrationVVS_wo_visp END");   
    if (it == NB_ITER_MAX) return 1;  // VVS has not converged before NB_ITER_MAX
    else return 0;
}




int vpHandEyeCalibration::calibrate_wo_visp(const std::vector<Mat> &cMo, const std::vector<Mat> &rMe, Mat &eMc, int n_sp)
{
    cout_indented(n_sp, "calibrate_wo_visp START");
    if (cMo.size() != rMe.size()) throw vpException(vpException::dimensionError, "cMo and rMe have different sizes");

    //Mat eRc(3, 3, CV_64F), eTc(3, 1, CV_64F), errVVS(3, 1, CV_64F);
    Mat eRc, eTc, errVVS;
    eMc = Mat::eye(4, 4, CV_64F);
    double resPos;
    /* initialisation of eMc to I in case all other methods fail */
    resPos = vpHandEyeCalibration::calibrationErrVVS_wo_visp(cMo, rMe, eMc, errVVS, n_sp + 1);
    double vmin = resPos;  // will serve to determine the best method
#if DEBUG_LEVEL1
    int He_method = HE_I;  // will serve to know which is the best method
#endif
    Mat eMcMin = eMc.clone();  // best initial estimation for VSS
    // Method using Old Tsai implementation
    int err = vpHandEyeCalibration::calibrationRotationTsaiOld_wo_visp(cMo, rMe, eRc, n_sp + 1);
    if (err != 0) printf("\n Problem in solving Hand-Eye Rotation by Old Tsai method \n");
    else
    {
        eRc.copyTo(eMc(Rect(0, 0, 3, 3)));
        err = vpHandEyeCalibration::calibrationTranslationOld_wo_visp(cMo, rMe, eRc, eTc, n_sp + 1);
        if (err != 0) printf("\n Problem in solving Hand-Eye Translation by Old Tsai method after Old Tsai method for Rotation\n");

        else
        {
            //eMc.insert(eTc);
            eTc.copyTo(eMc(Rect(3, 0, 1, 3)));
#if DEBUG_LEVEL1
            printf("\nRotation by (old) Tsai, old implementation for translation w/o visp\n");
            vpHandEyeCalibration::calibrationVerifrMo_wo_visp(cMo, rMe, eMc, n_sp + 1);
            //exit(0);
#endif
#if 1
            //resPos = vpHandEyeCalibration::calibrationErrVVS_wo_visp(cMo, rMe, eMc, errVVS);
            resPos = vpHandEyeCalibration::calibrationErrVVS_wo_visp(cMo, rMe, eMc, errVVS, n_sp + 1);
            if (resPos < vmin)
            {
                vmin = resPos;
                eMc.copyTo(eMcMin);
#if DEBUG_LEVEL1
                He_method = HE_TSAI_OROT;
#endif  //  DEBUG_LEVEL1
            }
#endif  //  1            
        }  
        err = vpHandEyeCalibration::calibrationTranslation_wo_visp(cMo, rMe, eRc, eTc, n_sp + 1);

        if (err != 0) printf("\n Problem in solving Hand-Eye Translation after Old Tsai method for Rotation\n");
        else
        {
            //eMc.insert(eTc);
            eTc.copyTo(eMc(Rect(3, 0, 1, 3)));
#if 0            
            printf("\nRotation by (old) Tsai, new implementation for translation : w/o visp\n");
            vpHandEyeCalibration::calibrationVerifrMo_wo_visp(cMo, rMe, eMc, n_sp + 1);
            //exit(0);
#endif  //  1            

#if DEBUG_LEVEL1
            {
                printf("\nRotation by (old) Tsai, new implementation for translation\n");
                vpHandEyeCalibration::calibrationVerifrMo(cMo, rMe, eMc);
            }
#endif  //  DEBUG_LEVEL1
            //resPos = vpHandEyeCalibration::calibrationErrVVS_wo_visp(cMo, rMe, eMc, errVVS);
            resPos = vpHandEyeCalibration::calibrationErrVVS_wo_visp(cMo, rMe, eMc, errVVS, n_sp + 1);
            if (resPos < vmin)
            {
                vmin = resPos;
                //eMcMin = eMc;
                eMc.copyTo(eMcMin);
#if DEBUG_LEVEL1
                He_method = HE_TSAI_ORNT;
#endif  //  DEBUG_LEVEL1
            }
        }
    }


    // First method using Tsai formulation
    err = vpHandEyeCalibration::calibrationRotationTsai_wo_visp(cMo, rMe, eRc, n_sp + 1);
  
  
    if (err != 0) printf("\n Problem in solving Hand-Eye Rotation by Tsai method \n");
    else
    {
        //eMc.insert(eRc);
        eRc.copyTo(eMc(Rect(0, 0, 3, 3)));
        //err = vpHandEyeCalibration::calibrationTranslationOld(cMo, rMe, eRc, eTc);
        err = vpHandEyeCalibration::calibrationTranslationOld_wo_visp(cMo, rMe, eRc, eTc, n_sp + 1);
        if (err != 0) printf("\n Problem in solving Hand-Eye Translation by Old Tsai method after Tsai method for Rotation\n");
        else
        {
            //eMc.insert(eTc);
            eTc.copyTo(eMc(Rect(3, 0, 1, 3)));
#if DEBUG_LEVEL1
            {
                printf("\nRotation by Tsai, old implementation for translation : w/o visp\n");
                //vpHandEyeCalibration::calibrationVerifrMo(cMo, rMe, eMc);
                vpHandEyeCalibration::calibrationVerifrMo_wo_visp(cMo, rMe, eMc, n_sp + 1);
                //exit(0);
            }
#endif  //  DEBUG_LEVEL1
            //resPos = vpHandEyeCalibration::calibrationErrVVS(cMo, rMe, eMc, errVVS);
            resPos = vpHandEyeCalibration::calibrationErrVVS_wo_visp(cMo, rMe, eMc, errVVS, n_sp + 1);
            if (resPos < vmin)
            {
                vmin = resPos;
                //eMcMin = eMc;
                eMc.copyTo(eMcMin);
#if DEBUG_LEVEL1
                He_method = HE_TSAI_NROT;
#endif  //  DEBUG_LEVEL1
            }
        }
        //err = vpHandEyeCalibration::calibrationTranslation(cMo, rMe, eRc, eTc);
        err = vpHandEyeCalibration::calibrationTranslation_wo_visp(cMo, rMe, eRc, eTc, n_sp + 1);
        if (err != 0) printf("\n Problem in solving Hand-Eye Translation after Tsai method for Rotation \n");
        else
        {
            //eMc.insert(eTc);
            eTc.copyTo(eMc(Rect(3, 0, 1, 3)));
#if 0
            printf("\nRotation by Tsai, new implementation for translation : w/o visp\n");
            vpHandEyeCalibration::calibrationVerifrMo_wo_visp(cMo, rMe, eMc, n_sp + 1);
            //exit(0);
#endif  //  1

#if DEBUG_LEVEL1
            {
                printf("\nRotation by Tsai, new implementation for translation\n");
                vpHandEyeCalibration::calibrationVerifrMo(cMo, rMe, eMc);
            }
#endif  //  DEBUG_LEVEL1
            resPos = vpHandEyeCalibration::calibrationErrVVS_wo_visp(cMo, rMe, eMc, errVVS, n_sp + 1);
            if (resPos < vmin)
            {
                vmin = resPos;
                //eMcMin = eMc;
                eMc.copyTo(eMcMin);
#if DEBUG_LEVEL1
                He_method = HE_TSAI_NRNT;
#endif  //  DEBUG_LEVEL1
            }
        }
    }
    // Second method by solving the orthogonal Procrustes problem
    
    err = vpHandEyeCalibration::calibrationRotationProcrustes_wo_visp(cMo, rMe, eRc, n_sp + 1);
    if(err != 0) printf("\n Problem in solving Hand-Eye Rotation by Procrustes method \n");
    else
    {
        //eMc.insert(eRc);
        eRc.copyTo(eMc(Rect(0, 0, 3, 3)));
        //err = vpHandEyeCalibration::calibrationTranslationOld(cMo, rMe, eRc, eTc);
        err = vpHandEyeCalibration::calibrationTranslationOld_wo_visp(cMo, rMe, eRc, eTc, n_sp + 1);
        if(err != 0) printf("\n Problem in solving Hand-Eye Translation by Old Tsai method after Procrustes method for Rotation\n");
        else
        {
            //eMc.insert(eTc);
            eTc.copyTo(eMc(Rect(3, 0, 1, 3)));
#if 0
            printf("\nRotation by Procrustes, old implementation for translation : w/o visp\n");
            vpHandEyeCalibration::calibrationVerifrMo_wo_visp(cMo, rMe, eMc, n_sp + 1);
#endif  //  1
#if DEBUG_LEVEL1
            {
                printf("\nRotation by Procrustes, old implementation for translation\n");
                vpHandEyeCalibration::calibrationVerifrMo(cMo, rMe, eMc);
            }
#endif
            resPos = vpHandEyeCalibration::calibrationErrVVS_wo_visp(cMo, rMe, eMc, errVVS, n_sp + 1);
            if(resPos < vmin)
            {
                vmin = resPos;
                //eMcMin = eMc;
                eMc.copyTo(eMcMin);
#if DEBUG_LEVEL1
                He_method = HE_PROCRUSTES_OT;
#endif  //  DEBUG_LEVEL1
            }
        }
        err = vpHandEyeCalibration::calibrationTranslation_wo_visp(cMo, rMe, eRc, eTc, n_sp + 1);
        if(err != 0) printf("\n Problem in solving Hand-Eye Translation after Procrustes method for Rotation\n");
        else
        {
            //eMc.insert(eTc);
            eTc.copyTo(eMc(Rect(3, 0, 1, 3)));
#if 0        
            printf("\nRotation by Procrustes, new implementation for translation : w/o visp\n");
            vpHandEyeCalibration::calibrationVerifrMo_wo_visp(cMo, rMe, eMc, n_sp + 1);
#endif  //  1

#if DEBUG_LEVEL1
            {
                printf("\nRotation by Procrustes, new implementation for translation\n");
                vpHandEyeCalibration::calibrationVerifrMo(cMo, rMe, eMc);
            }
#endif
            resPos = vpHandEyeCalibration::calibrationErrVVS_wo_visp(cMo, rMe, eMc, errVVS, n_sp + 1);
            if (resPos < vmin)
            {
                vmin = resPos;
                //eMcMin = eMc;
                eMc.copyTo(eMcMin);
#if DEBUG_LEVEL1
                He_method = HE_PROCRUSTES_NT;
#endif
            }
        }
    }

  /* determination of the best method in case at least one succeeds */
    //eMc = eMcMin;
    eMcMin.copyTo(eMc);

#if 0
    std::cout << "AAA : w/o visp" << std::endl;
    //vpThetaUVector ePc(eMc);
    Mat ePc;    Rodrigues(eMc(Rect(0, 0, 3, 3)), ePc);
    std::cout << "theta U (deg): " << rad2deg(ePc.at<double>(0)) << " " << rad2deg(ePc.at<double>(1)) << " " << rad2deg(ePc.at<double>(2)) << std::endl;
    std::cout << "Translation: " << eMc.at<double>(0, 3) << " " << eMc.at<double>(1, 3) << " " << eMc.at<double>(2, 3) << std::endl;
#endif  //  1

#if DEBUG_LEVEL1
  {
    if (He_method == HE_I) printf("Best method : I !!!, vmin = %lf\n",vmin);
    if (He_method == HE_TSAI_OROT) printf("Best method : TSAI_OROT, vmin = %lf\n",vmin);
    if (He_method == HE_TSAI_ORNT) printf("Best method : TSAI_ORNT, vmin = %lf\n",vmin);
    if (He_method == HE_TSAI_NROT) printf("Best method : TSAI_NROT, vmin = %lf\n",vmin);
    if (He_method == HE_TSAI_NRNT) printf("Best method : TSAI_NRNT, vmin = %lf\n",vmin);
    if (He_method == HE_PROCRUSTES_OT) printf("Best method : PROCRUSTES_OT, vmin = %lf\n",vmin);
    if (He_method == HE_PROCRUSTES_NT) printf("Best method : PROCRUSTES_NT, vmin = %lf\n",vmin);
    vpThetaUVector ePc(eMc);
    std::cout << "theta U (deg): " << vpMath::deg(ePc[0]) << " " << vpMath::deg(ePc[1]) << " " << vpMath::deg(ePc[2]) << std::endl;
    std::cout << "Translation: " << eMc[0][3] << " " << eMc[1][3] << " " << eMc[2][3] << std::endl;
  }
#endif

    // Non linear iterative minimization to estimate simultaneouslty eRc and eTc
    err = vpHandEyeCalibration::calibrationVVS_wo_visp(cMo, rMe, eMc, n_sp + 1);
    // FC : err : 0 si tout OK, -1 si pb de rang, 1 si pas convergence
    if (err != 0) printf("\n Problem in solving Hand-Eye Calibration by VVS \n");
    else
    {
        printf("\nRotation and translation after VVS : w/o visp\n");
        vpHandEyeCalibration::calibrationVerifrMo_wo_visp(cMo, rMe, eMc, n_sp + 1);
    }
    cout_indented(n_sp, "calibrate_wo_visp END");
    return err;
}
#endif  //  HE_WO_VISP

#undef HE_I
#undef HE_TSAI_OROT
#undef HE_TSAI_ORNT
#undef HE_TSAI_NROT
#undef HE_TSAI_NRNT
#undef HE_PROCRUSTES_OT
#undef HE_PROCRUSTES_NT

#undef DEBUG_LEVEL1
#undef DEBUG_LEVEL2
