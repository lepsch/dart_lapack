      SUBROUTINE CLAEIN( RIGHTV, NOINIT, N, H, LDH, W, V, B, LDB, RWORK, EPS3, SMLNUM, INFO );

*  -- LAPACK auxiliary routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      bool               NOINIT, RIGHTV;
      int                INFO, LDB, LDH, N;
      REAL               EPS3, SMLNUM;
      COMPLEX            W;
      // ..
      // .. Array Arguments ..
      REAL               RWORK( * );
      COMPLEX            B( LDB, * ), H( LDH, * ), V( * );
      // ..

*  =====================================================================

      // .. Parameters ..
      REAL               ONE, TENTH;
      const              ONE = 1.0, TENTH = 1.0e-1 ;
      COMPLEX            ZERO;
      const              ZERO = ( 0.0, 0.0 ) ;
      // ..
      // .. Local Scalars ..
      String             NORMIN, TRANS;
      int                I, IERR, ITS, J;
      REAL               GROWTO, NRMSML, ROOTN, RTEMP, SCALE, VNORM;
      COMPLEX            CDUM, EI, EJ, TEMP, X;
      // ..
      // .. External Functions ..
      int                ICAMAX;
      REAL               SCASUM, SCNRM2;
      COMPLEX            CLADIV;
      // EXTERNAL ICAMAX, SCASUM, SCNRM2, CLADIV
      // ..
      // .. External Subroutines ..
      // EXTERNAL CLATRS, CSSCAL
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, AIMAG, MAX, REAL, SQRT
      // ..
      // .. Statement Functions ..
      REAL               CABS1;
      // ..
      // .. Statement Function definitions ..
      CABS1( CDUM ) = ABS( REAL( CDUM ) ) + ABS( AIMAG( CDUM ) );
      // ..
      // .. Executable Statements ..

      INFO = 0;

      // GROWTO is the threshold used in the acceptance test for an
      // eigenvector.

      ROOTN = SQRT( REAL( N ) );
      GROWTO = TENTH / ROOTN;
      NRMSML = MAX( ONE, EPS3*ROOTN )*SMLNUM;

      // Form B = H - W*I (except that the subdiagonal elements are not
      // stored).

      for (J = 1; J <= N; J++) { // 20
         for (I = 1; I <= J - 1; I++) { // 10
            B( I, J ) = H( I, J );
         } // 10
         B( J, J ) = H( J, J ) - W;
      } // 20

      if ( NOINIT ) {

         // Initialize V.

         for (I = 1; I <= N; I++) { // 30
            V( I ) = EPS3;
         } // 30
      } else {

         // Scale supplied initial vector.

         VNORM = SCNRM2( N, V, 1 );
         csscal(N, ( EPS3*ROOTN ) / MAX( VNORM, NRMSML ), V, 1 );
      }

      if ( RIGHTV ) {

         // LU decomposition with partial pivoting of B, replacing zero
         // pivots by EPS3.

         for (I = 1; I <= N - 1; I++) { // 60
            EI = H( I+1, I );
            if ( CABS1( B( I, I ) ) < CABS1( EI ) ) {

               // Interchange rows and eliminate.

               X = CLADIV( B( I, I ), EI );
               B( I, I ) = EI;
               for (J = I + 1; J <= N; J++) { // 40
                  TEMP = B( I+1, J );
                  B( I+1, J ) = B( I, J ) - X*TEMP;
                  B( I, J ) = TEMP;
               } // 40
            } else {

               // Eliminate without interchange.

               IF( B( I, I ) == ZERO ) B( I, I ) = EPS3;
               X = CLADIV( EI, B( I, I ) );
               if ( X != ZERO ) {
                  for (J = I + 1; J <= N; J++) { // 50
                     B( I+1, J ) = B( I+1, J ) - X*B( I, J );
                  } // 50
               }
            }
         } // 60
         IF( B( N, N ) == ZERO ) B( N, N ) = EPS3;

         TRANS = 'N';

      } else {

         // UL decomposition with partial pivoting of B, replacing zero
         // pivots by EPS3.

         DO 90 J = N, 2, -1;
            EJ = H( J, J-1 );
            if ( CABS1( B( J, J ) ) < CABS1( EJ ) ) {

               // Interchange columns and eliminate.

               X = CLADIV( B( J, J ), EJ );
               B( J, J ) = EJ;
               for (I = 1; I <= J - 1; I++) { // 70
                  TEMP = B( I, J-1 );
                  B( I, J-1 ) = B( I, J ) - X*TEMP;
                  B( I, J ) = TEMP;
               } // 70
            } else {

               // Eliminate without interchange.

               IF( B( J, J ) == ZERO ) B( J, J ) = EPS3;
               X = CLADIV( EJ, B( J, J ) );
               if ( X != ZERO ) {
                  for (I = 1; I <= J - 1; I++) { // 80
                     B( I, J-1 ) = B( I, J-1 ) - X*B( I, J );
                  } // 80
               }
            }
         } // 90
         IF( B( 1, 1 ) == ZERO ) B( 1, 1 ) = EPS3;

         TRANS = 'C';

      }

      NORMIN = 'N';
      for (ITS = 1; ITS <= N; ITS++) { // 110

         // Solve U*x = scale*v for a right eigenvector
           // or U**H *x = scale*v for a left eigenvector,
         // overwriting x on v.

         clatrs('Upper', TRANS, 'Nonunit', NORMIN, N, B, LDB, V, SCALE, RWORK, IERR );
         NORMIN = 'Y';

         // Test for sufficient growth in the norm of v.

         VNORM = SCASUM( N, V, 1 );
         if (VNORM >= GROWTO*SCALE) GO TO 120;

         // Choose new orthogonal starting vector and try again.

         RTEMP = EPS3 / ( ROOTN+ONE );
         V( 1 ) = EPS3;
         for (I = 2; I <= N; I++) { // 100
            V( I ) = RTEMP;
         } // 100
         V( N-ITS+1 ) = V( N-ITS+1 ) - EPS3*ROOTN;
      } // 110

      // Failure to find eigenvector in N iterations.

      INFO = 1;

      } // 120

      // Normalize eigenvector.

      I = ICAMAX( N, V, 1 );
      csscal(N, ONE / CABS1( V( I ) ), V, 1 );

      RETURN;

      // End of CLAEIN

      }
