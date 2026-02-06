
// X->tautau
map<string, py::array_t<double>> kinfit_3pr_3pr(
						py::array_t<double> pt_1_vec,
						py::array_t<double> eta_1_vec,
						py::array_t<double> phi_1_vec,
						py::array_t<double> mass_1_vec,
						// 2-nd tau
						py::array_t<double> pt_2_vec,
						py::array_t<double> eta_2_vec,
						py::array_t<double> phi_2_vec,
						py::array_t<double> mass_2_vec,
						// MET and covariance
						py::array_t<double> met_pt_vec,
						py::array_t<double> met_phi_vec,
						py::array_t<double> metcov_xx_vec,
						py::array_t<double> metcov_xy_vec,
						py::array_t<double> metcov_yy_vec,
						// SV vector
						py::array_t<double> sv_x_vec,
						py::array_t<double> sv_y_vec,
						py::array_t<double> sv_z_vec,
						// SV covariance
						py::array_t<double> sv_xx_vec,
						py::array_t<double> sv_xy_vec,
						py::array_t<double> sv_yy_vec,
						py::array_t<double> sv_xz_vec,
						py::array_t<double> sv_yz_vec,
						py::array_t<double> sv_zz_vec,
						// SV vector
						py::array_t<double> sv2_x_vec,
						py::array_t<double> sv2_y_vec,
						py::array_t<double> sv2_z_vec,
						// SV covariance
						py::array_t<double> sv2_xx_vec,
						py::array_t<double> sv2_xy_vec,
						py::array_t<double> sv2_yy_vec,
						py::array_t<double> sv2_xz_vec,
						py::array_t<double> sv2_yz_vec,
						py::array_t<double> sv2_zz_vec,
						// steering parameters
						double mX) {

  unsigned int N = pt_1_vec.size();

  double delta=1.0/1.15;
  double reg_order=6.0;
  double massConstraint = mX > 0.;
  
  // Casting of variables
  auto pt_1 = pt_1_vec.unchecked<1>();
  auto eta_1 = eta_1_vec.unchecked<1>();
  auto phi_1 = phi_1_vec.unchecked<1>();
  auto mass_1 = mass_1_vec.unchecked<1>();
  
  auto pt_2 = pt_2_vec.unchecked<1>();
  auto eta_2 = eta_2_vec.unchecked<1>();
  auto phi_2 = phi_2_vec.unchecked<1>();
  auto mass_2 = mass_2_vec.unchecked<1>();

  auto met_pt  = met_pt_vec.unchecked<1>();
  auto met_phi = met_phi_vec.unchecked<1>();
  auto metcov_xx = metcov_xx_vec.unchecked<1>();
  auto metcov_xy = metcov_xy_vec.unchecked<1>();
  auto metcov_yy = metcov_yy_vec.unchecked<1>();

  auto sv_x = sv_x_vec.unchecked<1>();
  auto sv_y = sv_y_vec.unchecked<1>();
  auto sv_z = sv_z_vec.unchecked<1>();
  
  auto sv2_x = sv2_x_vec.unchecked<1>();
  auto sv2_y = sv2_y_vec.unchecked<1>();
  auto sv2_z = sv2_z_vec.unchecked<1>();
  
  auto sv_xx = sv_xx_vec.unchecked<1>();
  auto sv_xy = sv_xy_vec.unchecked<1>();
  auto sv_xz = sv_xz_vec.unchecked<1>();
  auto sv_yy = sv_yy_vec.unchecked<1>();
  auto sv_yz = sv_yz_vec.unchecked<1>();
  auto sv_zz = sv_zz_vec.unchecked<1>();

  auto sv2_xx = sv2_xx_vec.unchecked<1>();
  auto sv2_xy = sv2_xy_vec.unchecked<1>();
  auto sv2_xz = sv2_xz_vec.unchecked<1>();
  auto sv2_yy = sv2_yy_vec.unchecked<1>();
  auto sv2_yz = sv2_yz_vec.unchecked<1>();
  auto sv2_zz = sv2_zz_vec.unchecked<1>();

  // 
  
  py::buffer_info buffer = pt_1_vec.request();
  auto px_1_vec = py::array_t<double>(buffer.size);
  auto py_1_vec = py::array_t<double>(buffer.size);
  auto pz_1_vec = py::array_t<double>(buffer.size);
  auto px_2_vec = py::array_t<double>(buffer.size);
  auto py_2_vec = py::array_t<double>(buffer.size);
  auto pz_2_vec = py::array_t<double>(buffer.size);
  auto chi2_vec = py::array_t<double>(buffer.size);
  
  py::buffer_info px_1_buffer = px_1_vec.request();
  py::buffer_info py_1_buffer = py_1_vec.request();
  py::buffer_info pz_1_buffer = pz_1_vec.request();
  py::buffer_info px_2_buffer = px_2_vec.request();
  py::buffer_info py_2_buffer = py_2_vec.request();
  py::buffer_info pz_2_buffer = pz_2_vec.request();
  py::buffer_info chi2_buffer = chi2_vec.request();

  
  double * px_1_opt = static_cast<double *>(px_1_buffer.ptr);
  double * py_1_opt = static_cast<double *>(py_1_buffer.ptr);
  double * pz_1_opt = static_cast<double *>(pz_1_buffer.ptr);
  double * px_2_opt = static_cast<double *>(px_2_buffer.ptr);
  double * py_2_opt = static_cast<double *>(py_2_buffer.ptr);
  double * pz_2_opt = static_cast<double *>(pz_2_buffer.ptr);
  double * chi2_opt = static_cast<double *>(chi2_buffer.ptr);
  
  for (unsigned int i=0; i<N; ++i) {  

    // 4-vectors of taus
    // First tau decays to a1
    double m_1 = mass_1(i);
    double px_1 = pt_1(i)*cos(phi_1(i));
    double py_1 = pt_1(i)*sin(phi_1(i));
    double pz_1 = pt_1(i)*sinh(eta_1(i));
    double E_1 = sqrt(px_1*px_1+py_1*py_1+pz_1*pz_1+m_1*m_1);
    double ptot_1 = sqrt(px_1*px_1+py_1*py_1+pz_1*pz_1);
    
    double m_2 = mass_2(i);
    double px_2 = pt_2(i)*cos(phi_2(i));
    double py_2 = pt_2(i)*sin(phi_2(i));
    double pz_2 = pt_2(i)*sinh(eta_2(i));
    double E_2 = sqrt(px_2*px_2+py_2*py_2+pz_2*pz_2+m_2*m_2);
    double ptot_2 = sqrt(px_2*px_2+py_2*py_2+pz_2*pz_2);
    
    // avoiding lorentzvectors from ROOT and awkward
    double px_vis = px_1 + px_2;
    double py_vis = py_1 + py_2;
    double pz_vis = pz_1 + pz_2;
    double E_vis  = E_1 + E_2;

    double met_x = met_pt(i)*cos(met_phi(i));
    double met_y = met_pt(i)*sin(met_phi(i));
      
    double p_vis = sqrt(px_vis*px_vis+
			py_vis*py_vis+
			pz_vis*pz_vis);
    
    double m_vis = sqrt(E_vis*E_vis-p_vis*p_vis);
    double x1_min = min(double(1.0), (m_1/m_tau)*(m_1/m_tau));
    double x2_min = min(double(1.0), (m_2/m_tau)*(m_2/m_tau));
    
    // invert met covariance matrix, calculate determinant
    double metcovinv_xx = metcov_yy(i);
    double metcovinv_yy = metcov_xx(i);
    double metcovinv_xy = -metcov_xy(i);
    double metcovinv_yx = -metcov_xy(i);
    double metcovinv_det = (metcovinv_xx*metcovinv_yy -
			    metcovinv_yx*metcovinv_xy);
    if (metcovinv_det<1e-12) { 
      printf("Warning! Ill-conditioned MET covariance at event index %1i\n",i);
      printf("Diagonilizing MET covariance\n");
      metcovinv_xx = metcov_yy(i);
      metcovinv_yy = metcov_xx(i);
      metcovinv_det = metcovinv_xx*metcovinv_yy;
    }

    // SV
    double svcovinv[3][3];
    double sv_unit[3];

    double sv_mag = sqrt(sv_x(i)*sv_x(i)+sv_y(i)*sv_y(i)+sv_z(i)*sv_z(i));
    
    sv_unit[0] = sv_x(i)/sv_mag;
    sv_unit[1] = sv_y(i)/sv_mag;
    sv_unit[2] = sv_z(i)/sv_mag;
    
    double sv_matrix[3][3]; 
    /*      
    sv_matrix[0][0] = 1.0e+4*sv_xx(i);
    sv_matrix[0][1] = 1.0e+4*sv_xy(i);
    sv_matrix[0][2] = 1.0e+4*sv_xz(i);
    
    sv_matrix[1][0] = 1.0e+4*sv_xy(i);
    sv_matrix[1][1] = 1.0e+4*sv_yy(i);
    sv_matrix[1][2] = 1.0e+4*sv_yz(i);

    sv_matrix[2][0] = 1.0e+4*sv_xz(i);
    sv_matrix[2][1] = 1.0e+4*sv_yz(i);
    sv_matrix[2][2] = 1.0e+4*sv_zz(i);
    
    double svcovinv_det = invertedMatrix(sv_matrix,svcovinv);
    */
    
    //    if (svcovinv_det<1E-12) {
    //      printf("Warning! Ill-conditioned SV covariance at event index %1i\n",i);
    //      printf("Diagonilizing SV covariance\n");
    sv_matrix[0][1] = 0.;
    sv_matrix[1][0] = 0.;
    sv_matrix[0][2] = 0.;
    sv_matrix[2][0] = 0.;
    sv_matrix[1][2] = 0.;
    sv_matrix[2][1] = 0.;
    sv_matrix[0][0] = sv_xx(i);
    sv_matrix[1][1] = sv_yy(i);
    sv_matrix[2][2] = sv_zz(i);
    //    double svcovinv_det = invertedMatrix(sv_matrix,svcovinv);
    //    }
      
    double nz[3] = {px_1/ptot_1, py_1/ptot_1, pz_1/ptot_1};    
    double ny[3] = {1.,1.,1.};
    Cross(nz,sv_unit,ny);
    double norm = sqrt(ny[0]*ny[0]+ny[1]*ny[1]+ny[2]*ny[2]);
    ny[0] /= norm;
    ny[1] /= norm;
    ny[2] /= norm;
    double nx[3] = {1.,1.,1.};
    Cross(ny,nz,nx);
    norm = sqrt(nx[0]*nx[0]+nx[1]*nx[1]+nx[2]*nx[2]);
    nx[0] /= norm;
    nx[1] /= norm;
    nx[2] /= norm;

    // ***********
    // *** SV2 ***
    // ***********
    
    double svcovinv2[3][3];
    double sv2_unit[3];

    double sv2_mag = sqrt(sv2_x(i)*sv2_x(i)+sv2_y(i)*sv2_y(i)+sv2_z(i)*sv2_z(i));
    
    sv2_unit[0] = sv2_x(i)/sv2_mag;
    sv2_unit[1] = sv2_y(i)/sv2_mag;
    sv2_unit[2] = sv2_z(i)/sv2_mag;
    
    double sv2_matrix[3][3]; 

    /*
    sv2_matrix[0][0] = 1.0e+4*sv2_xx(i);
    sv2_matrix[0][1] = 1.0e+4*sv2_xy(i);
    sv2_matrix[0][2] = 1.0e+4*sv2_xz(i);
    
    sv2_matrix[1][0] = 1.0e+4*sv2_xy(i);
    sv2_matrix[1][1] = 1.0e+4*sv2_yy(i);
    sv2_matrix[1][2] = 1.0e+4*sv2_yz(i);

    sv2_matrix[2][0] = 1.0e+4*sv2_xz(i);
    sv2_matrix[2][1] = 1.0e+4*sv2_yz(i);
    sv2_matrix[2][2] = 1.0e+4*sv2_zz(i);
    
    double svcovinv2_det = invertedMatrix(sv2_matrix,svcovinv2);
    */
    //    if (svcovinv2_det<1E-12) {
    //      printf("Warning! Ill-conditioned SV covariance at event index %1i\n",i);
    //      printf("Diagonilizing SV covariance\n");
    sv2_matrix[0][1] = 0.;
    sv2_matrix[1][0] = 0.;
    sv2_matrix[0][2] = 0.;
    sv2_matrix[2][0] = 0.;
    sv2_matrix[1][2] = 0.;
    sv2_matrix[2][1] = 0.;
    sv2_matrix[0][0] = sv2_xx(i);
    sv2_matrix[1][1] = sv2_yy(i);
    sv2_matrix[2][2] = sv2_zz(i);
    //    double svcovinv2_det = invertedMatrix(sv2_matrix,svcovinv2);
    //    }
      
    double nz2[3] = {px_2/ptot_2, py_2/ptot_2, pz_2/ptot_2};    
    double ny2[3];
    Cross(nz2,sv2_unit,ny2);
    double norm2 = sqrt(ny2[0]*ny2[0]+ny2[1]*ny2[1]+ny2[2]*ny2[2]);
    ny2[0] /= norm2;
    ny2[1] /= norm2;
    ny2[2] /= norm2;
    double nx2[3];
    Cross(ny2,nz2,nx2);
    norm2 = sqrt(nx2[0]*nx2[0]+nx2[1]*nx2[1]+nx2[2]*nx2[2]);
    nx2[0] /= norm2;
    nx2[1] /= norm2;
    nx2[2] /= norm2;
    
    // Perform chi2 scan
    
    double min_chi2 = 1.0e+12;
    
    double thetaGJmax = ThetaGJmax(ptot_1,m_1);    
    double theta_range = thetaGJmax;  
    
    // scan over theta 
    unsigned int n_theta = 50;
    double dtheta = theta_range/float(n_theta); 

    for (unsigned int i1 = 0; i1 <= n_theta; ++i1) {
      
      double theta_scan = dtheta*double(i1);
      
      double cosTheta = cos(theta_scan);
      double sinTheta = sin(theta_scan);
      double difference[3] = {0.,0.,0.};
      double direction[3] = {0.,0.,0.};
      LinearSum(cosTheta,sinTheta,nz,nx,direction);
      double a = 1.0;
      double b = -1.0;
      LinearSum(a,b,direction,sv_unit,difference);
      double chi2_sv = sv_mag*sv_mag*(
				     difference[0]*difference[0]/sv_matrix[0][0]+
				     difference[1]*difference[1]/sv_matrix[1][1]+
				     difference[2]*difference[2]/sv_matrix[2][2]
				     );

      // solutions for momentum  
      vector<double> solutions = TauMomGJ(E_1,cosTheta,m_1);
      
      for (auto solution : solutions) {

	if (solution<0.0) continue;
	
	double pscan_1 = solution;
	double x1 = ptot_1/pscan_1;
	if ((x1 > 1.0) || (x1 < x1_min)) continue;
	
	if (massConstraint) {
	  double x2 = m_vis*m_vis/(x1*mX*mX);
	  double pscan_2 = ptot_2/x2;

	  if ((x2 < x2_min) || (x2 > 1.0)) continue; 

	  // test weighted four-vectors
	  double nu_test_px = px_1/x1 + px_2/x2 - px_vis;
	  double nu_test_py = py_1/x1 + py_2/x2 - py_vis;
	  
	  // calculate MET transfer function 
	  double residual_x = met_x - nu_test_px;
	  double residual_y = met_y - nu_test_py;
	  double chi2_met = (residual_x*(metcovinv_xx*residual_x + 
					 metcovinv_xy*residual_y) +
			     residual_y*(metcovinv_yx*residual_x +
					 metcovinv_yy*residual_y));
	  chi2_met /= metcovinv_det;

	  double theta_scan2 = ThetaGJ(ptot_2,pscan_2,m_2);

	  double cosTheta2 = cos(theta_scan2);
	  double sinTheta2 = sin(theta_scan2);
	  double difference2[3] = {0.,0.,0.};
	  double direction2[3] = {0.,0.,0.};
	  LinearSum(cosTheta2,sinTheta2,nz2,nx2,direction2);
	  double a2 = 1.0;
	  double b2 = -1.0;
	  LinearSum(a2,b2,direction2,sv2_unit,difference2);
	  double chi2_sv2 = sv2_mag*sv2_mag*(
					     difference2[0]*difference2[0]/sv2_matrix[0][0]+
					     difference2[1]*difference2[1]/sv2_matrix[1][1]+
					     difference2[2]*difference2[2]/sv2_matrix[2][2]
					     );
	
	  double chi2 = chi2_met + chi2_sv + chi2_sv2;

	//	  std::cout << "p1  " << pscan_1 <<  "  chi2_met = " << chi2_met << "  chi2_sv = " << chi2_sv << std::endl;
	
	  if (chi2 < min_chi2) {
	    min_chi2 = chi2;
	    px_1_opt[i] = pscan_1*direction[0];
	    py_1_opt[i] = pscan_1*direction[1];
	    pz_1_opt[i] = pscan_1*direction[2];
	    px_2_opt[i] = pscan_2*direction2[0];
	    py_2_opt[i] = pscan_2*direction2[1];
	    pz_2_opt[i] = pscan_2*direction2[2];
	  }
	}
	else {
	  for (int j=0; j<50; ++j) {
	    double x2 = 0.02 + double(j)*0.02;	    
	    if ((x2 < x2_min) || (x2 > 1.0)) continue;
	    double pscan_2 = ptot_2/x2;

	    double test_mass = m_vis/sqrt(x1*x2);
	
	    // calculate mass likelihood integral 
	    double m_shift = test_mass * delta;

	    if (m_shift < m_vis)
	      continue;

	    double x2_mini = max((m_2/m_tau)*(m_2/m_tau), 
				(m_vis/m_shift)*(m_vis/m_shift));
	    double x2_maxi = min(double(1.0), (m_vis/m_shift)*(m_vis/m_shift)/x1_min);

	    if (x2_maxi < x2_mini)
	      continue;

	    double J = 2*m_vis*m_vis * pow(m_shift, -reg_order);
	    double I_x2 = log(x2_maxi) - log(x2_mini);
	    double I_tot = I_x2;
	    double logL_mass = 10.0 + log(J) + log(I_tot);
	    
	    // test weighted four-vectors
	    double nu_test_px = px_1/x1 + px_2/x2 - px_vis;
	    double nu_test_py = py_1/x1 + py_2/x2 - py_vis;
	  
	    // calculate MET transfer function 
	    double residual_x = met_x - nu_test_px;
	    double residual_y = met_y - nu_test_py;
	    double chi2_met = (residual_x*(metcovinv_xx*residual_x + 
					   metcovinv_xy*residual_y) +
			       residual_y*(metcovinv_yx*residual_x +
					   metcovinv_yy*residual_y));
	    chi2_met /= metcovinv_det;

	    double theta_scan2 = ThetaGJ(ptot_2,pscan_2,m_2);

	    double cosTheta2 = cos(theta_scan2);
	    double sinTheta2 = sin(theta_scan2);
	    double difference2[3] = {0.,0.,0.};
	    double direction2[3] = {0.,0.,0.};
	    LinearSum(cosTheta2,sinTheta2,nz2,nx2,direction2);
	    double a2 = 1.0;
	    double b2 = -1.0;
	    LinearSum(a2,b2,direction2,sv2_unit,difference2);
	    double chi2_sv2 = sv2_mag*sv2_mag*(
					       difference2[0]*difference2[0]/sv2_matrix[0][0]+
					       difference2[1]*difference2[1]/sv2_matrix[1][1]+
					       difference2[2]*difference2[2]/sv2_matrix[2][2]
					       );
	    double chi2 = 0.5*chi2_met + 0.5*chi2_sv + 0.5*chi2_sv2 - logL_mass;
	
	    if (chi2 < min_chi2) {
	      min_chi2 = chi2;
	      px_1_opt[i] = pscan_1*direction[0];
	      py_1_opt[i] = pscan_1*direction[1];
	      pz_1_opt[i] = pscan_1*direction[2];
	      px_2_opt[i] = pscan_2*direction2[0];
	      py_2_opt[i] = pscan_2*direction2[1];
	      pz_2_opt[i] = pscan_2*direction2[2];
	    }
	  }
	}
      }
    }
    chi2_opt[i] = min_chi2;
  }

  map<string, py::array_t<double> > results = {
    {"chi2",chi2_vec},
    {"px_1",px_1_vec},
    {"py_1",py_1_vec},
    {"pz_1",pz_1_vec},
    {"px_2",px_2_vec},
    {"py_2",py_2_vec},
    {"pz_2",pz_2_vec},
  };
  
  return results;
    
}


