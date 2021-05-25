/**
 * 2D Serial-link robot class, modeled after Corke's toolbox
 */
class SerialLink2 {

  int n_; ///< number of links
  Link2[] links_; /// The links
  Pose2 base_, tool_; ///< arbitrary base and tool transformations

  /// Constructors with base and tool constants
  SerialLink2(int n, Link2[] links, Pose2 base, Pose2 tool) {
    n_ = n;
    links_ = links;
    base_ = base;
    tool_ = tool;
  }

  /// Constructor with links only, only references are copied
  SerialLink2(int n, Link2[] links) {
    this(n, links, new Pose2(0, 0, 0), new Pose2(0, 0, 0));
  }

  /// Forward kinematics, code similar to Robot toolbox
  Pose2 fkine(double[] q) {
    Pose2 t = base_;
    for (int i=0;i<n_;i++) {
      Pose2 Ti = links_[i].T(q[i]);
      t = t.compose(Ti);
    }
    t = t.compose(tool_);
    return t;
  }

  /// Analytical Jacobian for three-link manipulator
  static double[][] jacobian3(double q[]) {
    double[][] J = new double[3][3];

    double alpha = q[0]+q[1], beta = alpha + q[2];
    double sin1 = Math.sin(q[0]), cos1 = Math.cos(q[0]);
    double cosa = Math.cos(alpha), sina = Math.sin(alpha);
    double cosb = Math.cos(beta), sinb = Math.sin(beta);

    // x
    J[0][0] = -3.5*sin1 - 3.5*sina - 2*sinb;
    J[0][1] =            -3.5*sina - 2*sinb;
    J[0][2] =                       -2*sinb;

    // y
    J[1][0] = 3.5*cos1 + 3.5*cosa + 2*cosb;
    J[1][1] =            3.5*cosa + 2*cosb;
    J[1][2] =                       2*cosb;

    // theta
    J[2][0] = 1.0;
    J[2][1] = 1.0;
    J[2][2] = 1.0;

    return J;
  }

  // DO in-place QR factorization of A, also computes Rdiag
  static void inPlaceQR(double[][] A, int m, int n, double[] Rdiag) {

    for (int k = 0; k < n; k++) {
      // Compute 2-norm of k-th column without under/overflow.
      double nrm = 0;
      for (int i = k; i < m; i++) {
        nrm = Math.hypot(nrm, A[i][k]);
      }

      if (nrm != 0.0) {
        // Form k-th Householder vector.
        if (A[k][k] < 0) {
          nrm = -nrm;
        }
        for (int i = k; i < m; i++) {
          A[i][k] /= nrm;
        }
        A[k][k] += 1.0;

        // Apply transformation to remaining columns.
        for (int j = k+1; j < n; j++) {
          double s = 0.0;
          for (int i = k; i < m; i++) {
            s += A[i][k]*A[i][j];
          }
          s = -s/A[k][k];
          for (int i = k; i < m; i++) {
            A[i][j] += s*A[i][k];
          }
        }
      }
      Rdiag[k] = -nrm;
    }
  }


  // Solve for RHS b, given QR factorization A
  static void solveRHS(double[][] A, double[] Rdiag, int m, int n, double[] b) {

    // Compute Y = transpose(Q)*B
    for (int k = 0; k < n; k++) {
      double s = 0.0;
      for (int i = k; i < m; i++) {
        s += A[i][k]*b[i];
      }
      s = -s/A[k][k];
      for (int i = k; i < m; i++) {
        b[i] += s*A[i][k];
      }
    }

    // Solve R*X = Y;
    for (int k = n-1; k >= 0; k--) {
      b[k] /= Rdiag[k];
      for (int i = 0; i < k; i++) {
        b[i] -= b[k]*A[i][k];
      }
    }
  }

  /// Inverse Jacobian for three-link manipulator
  /// Fairly elaborate, using QR, because that's the code I had around.
  static double[][] inverseJacobian3(double q[]) {
    // Perform QR on Jacobian
    double[][] A = jacobian3(q);
    double[] Rdiag = new double[3];
    inPlaceQR(A, 3, 3, Rdiag);

    // Recover the inverse by using as RHS the columns of identity matrix
    double[] b = new double[3];
    double[][] invJ = new double[3][3];

    b[0] = 1;
    b[1] = 0;
    b[2] = 0;
    solveRHS(A, Rdiag, 3, 3, b);
    invJ[0][0] = b[0];
    invJ[1][0] = b[1];
    invJ[2][0] = b[2];

    b[0] = 0;
    b[1] = 1;
    b[2] = 0;
    solveRHS(A, Rdiag, 3, 3, b);
    invJ[0][1] = b[0];
    invJ[1][1] = b[1];
    invJ[2][1] = b[2];

    b[0] = 0;
    b[1] = 0;
    b[2] = 1;
    solveRHS(A, Rdiag, 3, 3, b);
    invJ[0][2] = b[0];
    invJ[1][2] = b[1];
    invJ[2][2] = b[2];

    return invJ;
  }

  static double[] proportionalJointAngleController(double[] q, double[] qd, double Kp) {
    ///////////////////////////////////////////////////////////////////////
    ///WHAT DOES THIS FUNCTION DO: 
    //This function takes in current joint angles (q), and desired destination joint angles (qd) of
    //a 3 revolute joint robot, as well as a  gain factor Kp, and calculates and returns the next 
    //joint angle position/waypoint to get to, in order to follow joint-space motion control.
    
    //INPUTS :  
    //q[] is an array containing the current joint angles indexed in ascending order of the joints starting from the base joint angle
    //qd[] is an array containing the desired end joint angles of the robot, in the same form as q
    //Kp is the gain factor. an example value used for Kp is 0.05.
    
    //OUTPUTS:
    //res[] is an array that is the same size as q[], and is the next joint angle position/waypoint
    //to get to, in order to follow joint-space motion control
//
    ///////////////////////////////////////////////////////////////////////
    
    //copying q current joint angles to another array, res. We will return this array, res.
    double[] res = new double[q.length];
    for (int i = 0; i < res.length; i++) {
      res[i] = q[i];
    }
    
    ///////////////////////////////////////////////////////////////////////
    ///STUDENT CODE START: THIS IS WHERE THE STUDENT JOINT ANGLE TRAJECTORY CONTROL IS IMPLEMENTED
    //TO DO FOR STUDENT. 
    ///////////////////////////////////////////////////////////////////////
    
    //calculate error in joint angles and remember to make angles within PI and -PI
    for (int i = 0; i < q.length; i++) {
      // qt+1 = qt + Kp(qd - qt
      double error = qd[i] - res[i];
      error = (error + Math.PI) % (2 * Math.PI) - Math.PI;
      res[i] = res[i] + Kp*error;
    }
    //do pi control for angles. Use equation from slide 14 of https://dellaert.github.io/20S-3630/Slides/L11_Jacobian_Control.pdf
    //also use as reference https://dellaert.github.io/20S-3630/notes/3-manipulators.pdf
    
    ///////////////////////////////////////////////////////////////////////
    ///STUDENT CODE END: THIS IS WHERE THE STUDENT JOINT ANGLE TRAJECTORY CONTROL IS IMPLEMENTED
    //TO DO FOR STUDENT. 
    ///////////////////////////////////////////////////////////////////////
    return res;
  }

  static double[] proportionalCartesianController(double xd, double yd, double theta_d, double[] q, Pose2 tool, double Kp) {
    ///////////////////////////////////////////////////////////////////////
    ///WHAT DOES THIS FUNCTION DO: 
    //This function takes in current and desired x,y, theta, as well as a  gain factor Kp, of
    //a 3 revolute joint robot end effector, and calculates and returns the next 
    //joint angle position/waypoint to get to, in order to follow cartesian motion control.
    
    //INPUTS :  
    //xd is the desired x position of the end effector
    //yd is the desired y position of the end effector
    //theta_d is the desired theta of the end effector, which is the theta of the last joint
    //q[] is an array containing the current joint angles indexed in ascending order of the joints starting from the base joint angle
    //tool is a variable of type Pose2 that contains the current pose of the end effector.
    //Use tool.x(), tool.y(), and tool.theta() to get current pose of end effector
    //Kp is the gain factor. an example value used for Kp is 0.05.
    
    //OUTPUTS:
    //res[] is an array that is the same size as q[], and is the next joint angle position/waypoint
    //to get to, in order to follow cartesian motion control
//
    ///////////////////////////////////////////////////////////////////////    
    
    //copying q current joint angles to another array, res. We will return this array, res.
    double[] res = new double[q.length];
    for (int i = 0; i < res.length; i++) {
      res[i] = q[i];
    }

    ///////////////////////////////////////////////////////////////////////
    ///STUDENT CODE START: THIS IS WHERE THE STUDENT CARTESIAN TRAJECTORY CONTROL IS IMPLEMENTED
    //TO DO FOR STUDENT. 
    ///////////////////////////////////////////////////////////////////////
    
    // 1. Calculate Inverse Jacobian USING SerialLink2.inverseJacobian FUNCTION
    double[][] invJac = inverseJacobian3(q);
    // 2.Calculate position error
    double xError = xd - tool.x();
    double yError = yd - tool.y();
    double thetaError = theta_d - tool.theta();
    thetaError = (thetaError + Math.PI) % (2* Math.PI) - Math.PI;
    double[] Et = {xError, yError, thetaError};
    // 3. Do PI Control on position. Use equation from slide 21 of https://dellaert.github.io/20S-3630/Slides/L11_Jacobian_Control.pdf
    //also use as reference https://dellaert.github.io/20S-3630/notes/3-manipulators.pdf
    for (int i = 0; i < res.length; i++) {
      double deltaQ = Kp * (invJac[i][0]*Et[0] + invJac[i][1] * Et[1] + invJac[i][2] * Et[2]);
      res[i] = res[i] + deltaQ;
    }
    ///////////////////////////////////////////////////////////////////////
    ///STUDENT CODE END: THIS IS WHERE THE STUDENT CARTESIAN TRAJECTORY CONTROL IS IMPLEMENTED
    //TO DO FOR STUDENT. 
    ///////////////////////////////////////////////////////////////////////
    
    return res;
  }

}
