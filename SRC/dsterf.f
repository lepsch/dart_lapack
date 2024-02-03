      SUBROUTINE DSTERF( N, D, E, INFO )

*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                INFO, N;
      // ..
      // .. Array Arguments ..
      double             D( * ), E( * );
      // ..

*  =====================================================================

      // .. Parameters ..
      double             ZERO, ONE, TWO, THREE;
      const              ZERO = 0.0D0, ONE = 1.0D0, TWO = 2.0D0, THREE = 3.0D0 ;
      int                MAXIT;
      const              MAXIT = 30 ;
      // ..
      // .. Local Scalars ..
      int                I, ISCALE, JTOT, L, L1, LEND, LENDSV, LSV, M, NMAXIT;
      double             ALPHA, ANORM, BB, C, EPS, EPS2, GAMMA, OLDC, OLDGAM, P, R, RT1, RT2, RTE, S, SAFMAX, SAFMIN, SIGMA, SSFMAX, SSFMIN, RMAX;
      // ..
      // .. External Functions ..
      double             DLAMCH, DLANST, DLAPY2;
      // EXTERNAL DLAMCH, DLANST, DLAPY2
      // ..
      // .. External Subroutines ..
      // EXTERNAL DLAE2, DLASCL, DLASRT, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, SIGN, SQRT
      // ..
      // .. Executable Statements ..

      // Test the input parameters.

      INFO = 0

      // Quick return if possible

      if ( N.LT.0 ) {
         INFO = -1
         xerbla('DSTERF', -INFO );
         RETURN
      }
      if (N.LE.1) RETURN;

      // Determine the unit roundoff for this environment.

      EPS = DLAMCH( 'E' )
      EPS2 = EPS**2
      SAFMIN = DLAMCH( 'S' )
      SAFMAX = ONE / SAFMIN
      SSFMAX = SQRT( SAFMAX ) / THREE
      SSFMIN = SQRT( SAFMIN ) / EPS2
      RMAX = DLAMCH( 'O' )

      // Compute the eigenvalues of the tridiagonal matrix.

      NMAXIT = N*MAXIT
      SIGMA = ZERO
      JTOT = 0

      // Determine where the matrix splits and choose QL or QR iteration
      // for each block, according to whether top or bottom diagonal
      // element is smaller.

      L1 = 1

      } // 10
      if (L1.GT.N) GO TO 170       IF( L1.GT.1 ) E( L1-1 ) = ZERO;
      for (M = L1; M <= N - 1; M++) { // 20
         if ( ABS( E( M ) ).LE.( SQRT( ABS( D( M ) ) )*SQRT( ABS( D( M+ 1 ) ) ) )*EPS ) {
            E( M ) = ZERO
            GO TO 30
         }
      } // 20
      M = N

      } // 30
      L = L1
      LSV = L
      LEND = M
      LENDSV = LEND
      L1 = M + 1
      if (LEND.EQ.L) GO TO 10;

      // Scale submatrix in rows and columns L to LEND

      ANORM = DLANST( 'M', LEND-L+1, D( L ), E( L ) )
      ISCALE = 0
      if (ANORM.EQ.ZERO) GO TO 10;
      if ( (ANORM.GT.SSFMAX) ) {
         ISCALE = 1
         dlascl('G', 0, 0, ANORM, SSFMAX, LEND-L+1, 1, D( L ), N, INFO );
         dlascl('G', 0, 0, ANORM, SSFMAX, LEND-L, 1, E( L ), N, INFO );
      } else if ( ANORM.LT.SSFMIN ) {
         ISCALE = 2
         dlascl('G', 0, 0, ANORM, SSFMIN, LEND-L+1, 1, D( L ), N, INFO );
         dlascl('G', 0, 0, ANORM, SSFMIN, LEND-L, 1, E( L ), N, INFO );
      }

      for (I = L; I <= LEND - 1; I++) { // 40
         E( I ) = E( I )**2
      } // 40

      // Choose between QL and QR iteration

      if ( ABS( D( LEND ) ).LT.ABS( D( L ) ) ) {
         LEND = LSV
         L = LENDSV
      }

      if ( LEND.GE.L ) {

         // QL Iteration

         // Look for small subdiagonal element.

         } // 50
         if ( L.NE.LEND ) {
            for (M = L; M <= LEND - 1; M++) { // 60
               IF( ABS( E( M ) ).LE.EPS2*ABS( D( M )*D( M+1 ) ) ) GO TO 70
            } // 60
         }
         M = LEND

         } // 70
         if (M.LT.LEND) E( M ) = ZERO;
         P = D( L )
         if (M.EQ.L) GO TO 90;

         // If remaining matrix is 2 by 2, use DLAE2 to compute its
         // eigenvalues.

         if ( M.EQ.L+1 ) {
            RTE = SQRT( E( L ) )
            dlae2(D( L ), RTE, D( L+1 ), RT1, RT2 );
            D( L ) = RT1
            D( L+1 ) = RT2
            E( L ) = ZERO
            L = L + 2
            if (L.LE.LEND) GO TO 50;
            GO TO 150
         }

         if (JTOT.EQ.NMAXIT) GO TO 150;
         JTOT = JTOT + 1

         // Form shift.

         RTE = SQRT( E( L ) )
         SIGMA = ( D( L+1 )-P ) / ( TWO*RTE )
         R = DLAPY2( SIGMA, ONE )
         SIGMA = P - ( RTE / ( SIGMA+SIGN( R, SIGMA ) ) )

         C = ONE
         S = ZERO
         GAMMA = D( M ) - SIGMA
         P = GAMMA*GAMMA

         // Inner loop

         DO 80 I = M - 1, L, -1
            BB = E( I )
            R = P + BB
            if (I.NE.M-1) E( I+1 ) = S*R;
            OLDC = C
            C = P / R
            S = BB / R
            OLDGAM = GAMMA
            ALPHA = D( I )
            GAMMA = C*( ALPHA-SIGMA ) - S*OLDGAM
            D( I+1 ) = OLDGAM + ( ALPHA-GAMMA )
            if ( C.NE.ZERO ) {
               P = ( GAMMA*GAMMA ) / C
            } else {
               P = OLDC*BB
            }
         } // 80

         E( L ) = S*P
         D( L ) = SIGMA + GAMMA
         GO TO 50

         // Eigenvalue found.

         } // 90
         D( L ) = P

         L = L + 1
         if (L.LE.LEND) GO TO 50;
         GO TO 150

      } else {

         // QR Iteration

         // Look for small superdiagonal element.

         } // 100
         DO 110 M = L, LEND + 1, -1
            IF( ABS( E( M-1 ) ).LE.EPS2*ABS( D( M )*D( M-1 ) ) ) GO TO 120
         } // 110
         M = LEND

         } // 120
         if (M.GT.LEND) E( M-1 ) = ZERO;
         P = D( L )
         if (M.EQ.L) GO TO 140;

         // If remaining matrix is 2 by 2, use DLAE2 to compute its
         // eigenvalues.

         if ( M.EQ.L-1 ) {
            RTE = SQRT( E( L-1 ) )
            dlae2(D( L ), RTE, D( L-1 ), RT1, RT2 );
            D( L ) = RT1
            D( L-1 ) = RT2
            E( L-1 ) = ZERO
            L = L - 2
            if (L.GE.LEND) GO TO 100;
            GO TO 150
         }

         if (JTOT.EQ.NMAXIT) GO TO 150;
         JTOT = JTOT + 1

         // Form shift.

         RTE = SQRT( E( L-1 ) )
         SIGMA = ( D( L-1 )-P ) / ( TWO*RTE )
         R = DLAPY2( SIGMA, ONE )
         SIGMA = P - ( RTE / ( SIGMA+SIGN( R, SIGMA ) ) )

         C = ONE
         S = ZERO
         GAMMA = D( M ) - SIGMA
         P = GAMMA*GAMMA

         // Inner loop

         for (I = M; I <= L - 1; I++) { // 130
            BB = E( I )
            R = P + BB
            if (I.NE.M) E( I-1 ) = S*R;
            OLDC = C
            C = P / R
            S = BB / R
            OLDGAM = GAMMA
            ALPHA = D( I+1 )
            GAMMA = C*( ALPHA-SIGMA ) - S*OLDGAM
            D( I ) = OLDGAM + ( ALPHA-GAMMA )
            if ( C.NE.ZERO ) {
               P = ( GAMMA*GAMMA ) / C
            } else {
               P = OLDC*BB
            }
         } // 130

         E( L-1 ) = S*P
         D( L ) = SIGMA + GAMMA
         GO TO 100

         // Eigenvalue found.

         } // 140
         D( L ) = P

         L = L - 1
         if (L.GE.LEND) GO TO 100;
         GO TO 150

      }

      // Undo scaling if necessary

      } // 150
      if (ISCALE.EQ.1) CALL DLASCL( 'G', 0, 0, SSFMAX, ANORM, LENDSV-LSV+1, 1, D( LSV ), N, INFO )       IF( ISCALE.EQ.2 ) CALL DLASCL( 'G', 0, 0, SSFMIN, ANORM, LENDSV-LSV+1, 1, D( LSV ), N, INFO );

      // Check for no convergence to an eigenvalue after a total
      // of N*MAXIT iterations.

      if (JTOT.LT.NMAXIT) GO TO 10;
      for (I = 1; I <= N - 1; I++) { // 160
         IF( E( I ).NE.ZERO ) INFO = INFO + 1
      } // 160
      GO TO 180

      // Sort eigenvalues in increasing order.

      } // 170
      dlasrt('I', N, D, INFO );

      } // 180
      RETURN

      // End of DSTERF

      }
