      SUBROUTINE DDRGVX( NSIZE, THRESH, NIN, NOUT, A, LDA, B, AI, BI, ALPHAR, ALPHAI, BETA, VL, VR, ILO, IHI, LSCALE, RSCALE, S, DTRU, DIF, DIFTRU, WORK, LWORK, IWORK, LIWORK, RESULT, BWORK, INFO );

// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                IHI, ILO, INFO, LDA, LIWORK, LWORK, NIN, NOUT, NSIZE;
      double             THRESH;
      // ..
      // .. Array Arguments ..
      bool               BWORK( * );
      int                IWORK( * );
      double             A( LDA, * ), AI( LDA, * ), ALPHAI( * ), ALPHAR( * ), B( LDA, * ), BETA( * ), BI( LDA, * ), DIF( * ), DIFTRU( * ), DTRU( * ), LSCALE( * ), RESULT( 4 ), RSCALE( * ), S( * ), VL( LDA, * ), VR( LDA, * ), WORK( * );
      // ..

// =====================================================================

      // .. Parameters ..
      double             ZERO, ONE, TEN, TNTH, HALF;
      const              ZERO = 0.0, ONE = 1.0, TEN = 1.0e+1, TNTH = 1.0e-1, HALF = 0.5 ;
      // ..
      // .. Local Scalars ..
      int                I, IPTYPE, IWA, IWB, IWX, IWY, J, LINFO, MAXWRK, MINWRK, N, NERRS, NMAX, NPTKNT, NTESTT;
      double             ABNORM, ANORM, BNORM, RATIO1, RATIO2, THRSH2, ULP, ULPINV;
      // ..
      // .. Local Arrays ..
      double             WEIGHT( 5 );
      // ..
      // .. External Functions ..
      int                ILAENV;
      double             DLAMCH, DLANGE;
      // EXTERNAL ILAENV, DLAMCH, DLANGE
      // ..
      // .. External Subroutines ..
      // EXTERNAL ALASVM, DGET52, DGGEVX, DLACPY, DLATM6, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, MAX, SQRT
      // ..
      // .. Executable Statements ..

      // Check for errors

      INFO = 0;

      NMAX = 5;

      if ( NSIZE < 0 ) {
         INFO = -1;
      } else if ( THRESH < ZERO ) {
         INFO = -2;
      } else if ( NIN <= 0 ) {
         INFO = -3;
      } else if ( NOUT <= 0 ) {
         INFO = -4;
      } else if ( LDA < 1 || LDA < NMAX ) {
         INFO = -6;
      } else if ( LIWORK < NMAX+6 ) {
         INFO = -26;
      }

      // Compute workspace
       // (Note: Comments in the code beginning "Workspace:" describe the
        // minimal amount of workspace needed at that point in the code,
        // as well as the preferred amount for good performance.
        // NB refers to the optimal block size for the immediately
        // following subroutine, as returned by ILAENV.)

      MINWRK = 1;
      if ( INFO == 0 && LWORK >= 1 ) {
         MINWRK = 2*NMAX*NMAX + 12*NMAX + 16;
         MAXWRK = 6*NMAX + NMAX*ILAENV( 1, 'DGEQRF', ' ', NMAX, 1, NMAX, 0 );
         MAXWRK = MAX( MAXWRK, 2*NMAX*NMAX+12*NMAX+16 );
         WORK( 1 ) = MAXWRK;
      }

      if (LWORK < MINWRK) INFO = -24;

      if ( INFO != 0 ) {
         xerbla('DDRGVX', -INFO );
         return;
      }

      N = 5;
      ULP = DLAMCH( 'P' );
      ULPINV = ONE / ULP;
      THRSH2 = TEN*THRESH;
      NERRS = 0;
      NPTKNT = 0;
      NTESTT = 0;

      if (NSIZE == 0) GO TO 90;

      // Parameters used for generating test matrices.

      WEIGHT( 1 ) = TNTH;
      WEIGHT( 2 ) = HALF;
      WEIGHT( 3 ) = ONE;
      WEIGHT( 4 ) = ONE / WEIGHT( 2 );
      WEIGHT( 5 ) = ONE / WEIGHT( 1 );

      for (IPTYPE = 1; IPTYPE <= 2; IPTYPE++) { // 80
         for (IWA = 1; IWA <= 5; IWA++) { // 70
            for (IWB = 1; IWB <= 5; IWB++) { // 60
               for (IWX = 1; IWX <= 5; IWX++) { // 50
                  for (IWY = 1; IWY <= 5; IWY++) { // 40

                     // generated a test matrix pair

                     dlatm6(IPTYPE, 5, A, LDA, B, VR, LDA, VL, LDA, WEIGHT( IWA ), WEIGHT( IWB ), WEIGHT( IWX ), WEIGHT( IWY ), DTRU, DIFTRU );

                     // Compute eigenvalues/eigenvectors of (A, B).
                     // Compute eigenvalue/eigenvector condition numbers
                     // using computed eigenvectors.

                     dlacpy('F', N, N, A, LDA, AI, LDA );
                     dlacpy('F', N, N, B, LDA, BI, LDA );

                     dggevx('N', 'V', 'V', 'B', N, AI, LDA, BI, LDA, ALPHAR, ALPHAI, BETA, VL, LDA, VR, LDA, ILO, IHI, LSCALE, RSCALE, ANORM, BNORM, S, DIF, WORK, LWORK, IWORK, BWORK, LINFO );
                     if ( LINFO != 0 ) {
                        RESULT( 1 ) = ULPINV;
                        WRITE( NOUT, FMT = 9999 )'DGGEVX', LINFO, N, IPTYPE;
                        GO TO 30;
                     }

                     // Compute the norm(A, B)

                     dlacpy('Full', N, N, AI, LDA, WORK, N );
                     dlacpy('Full', N, N, BI, LDA, WORK( N*N+1 ), N );
                     ABNORM = DLANGE( 'Fro', N, 2*N, WORK, N, WORK );

                     // Tests (1) and (2)

                     RESULT( 1 ) = ZERO;
                     dget52( true , N, A, LDA, B, LDA, VL, LDA, ALPHAR, ALPHAI, BETA, WORK, RESULT( 1 ) );
                     if ( RESULT( 2 ) > THRESH ) {
                        WRITE( NOUT, FMT = 9998 )'Left', 'DGGEVX', RESULT( 2 ), N, IPTYPE, IWA, IWB, IWX, IWY;
                     }

                     RESULT( 2 ) = ZERO;
                     dget52( false , N, A, LDA, B, LDA, VR, LDA, ALPHAR, ALPHAI, BETA, WORK, RESULT( 2 ) );
                     if ( RESULT( 3 ) > THRESH ) {
                        WRITE( NOUT, FMT = 9998 )'Right', 'DGGEVX', RESULT( 3 ), N, IPTYPE, IWA, IWB, IWX, IWY;
                     }

                     // Test (3)

                     RESULT( 3 ) = ZERO;
                     for (I = 1; I <= N; I++) { // 10
                        if ( S( I ) == ZERO ) {
                           IF( DTRU( I ) > ABNORM*ULP ) RESULT( 3 ) = ULPINV;
                        } else if ( DTRU( I ) == ZERO ) {
                           IF( S( I ) > ABNORM*ULP ) RESULT( 3 ) = ULPINV;
                        } else {
                           WORK( I ) = MAX( ABS( DTRU( I ) / S( I ) ), ABS( S( I ) / DTRU( I ) ) );
                           RESULT( 3 ) = MAX( RESULT( 3 ), WORK( I ) );
                        }
                     } // 10

                     // Test (4)

                     RESULT( 4 ) = ZERO;
                     if ( DIF( 1 ) == ZERO ) {
                        IF( DIFTRU( 1 ) > ABNORM*ULP ) RESULT( 4 ) = ULPINV;
                     } else if ( DIFTRU( 1 ) == ZERO ) {
                        IF( DIF( 1 ) > ABNORM*ULP ) RESULT( 4 ) = ULPINV;
                     } else if ( DIF( 5 ) == ZERO ) {
                        IF( DIFTRU( 5 ) > ABNORM*ULP ) RESULT( 4 ) = ULPINV;
                     } else if ( DIFTRU( 5 ) == ZERO ) {
                        IF( DIF( 5 ) > ABNORM*ULP ) RESULT( 4 ) = ULPINV;
                     } else {
                        RATIO1 = MAX( ABS( DIFTRU( 1 ) / DIF( 1 ) ), ABS( DIF( 1 ) / DIFTRU( 1 ) ) )                         RATIO2 = MAX( ABS( DIFTRU( 5 ) / DIF( 5 ) ), ABS( DIF( 5 ) / DIFTRU( 5 ) ) );
                        RESULT( 4 ) = MAX( RATIO1, RATIO2 );
                     }

                     NTESTT = NTESTT + 4;

                     // Print out tests which fail.

                     for (J = 1; J <= 4; J++) { // 20
                        if ( ( RESULT( J ) >= THRSH2 && J >= 4 ) || ( RESULT( J ) >= THRESH && J <= 3 ) ) {

                        // If this is the first test to fail,
                        // print a header to the data file.

                           if ( NERRS == 0 ) {
                              WRITE( NOUT, FMT = 9997 )'DXV';

                           // Print out messages for built-in examples

                           // Matrix types

                              WRITE( NOUT, FMT = 9995 );
                              WRITE( NOUT, FMT = 9994 );
                              WRITE( NOUT, FMT = 9993 );

                           // Tests performed

                              WRITE( NOUT, FMT = 9992 )'''', 'transpose', '''';

                           }
                           NERRS = NERRS + 1;
                           if ( RESULT( J ) < 10000.0 ) {
                              WRITE( NOUT, FMT = 9991 )IPTYPE, IWA, IWB, IWX, IWY, J, RESULT( J );
                           } else {
                              WRITE( NOUT, FMT = 9990 )IPTYPE, IWA, IWB, IWX, IWY, J, RESULT( J );
                           }
                        }
                     } // 20

                     } // 30

                  } // 40
               } // 50
            } // 60
         } // 70
      } // 80

      GO TO 150;

      } // 90

      // Read in data from file to check accuracy of condition estimation
      // Read input data until N=0

      READ( NIN, FMT = *, END = 150 )N;
      if (N == 0) GO TO 150;
      for (I = 1; I <= N; I++) { // 100
         READ( NIN, FMT = * )( A( I, J ), J = 1, N );
      } // 100
      for (I = 1; I <= N; I++) { // 110
         READ( NIN, FMT = * )( B( I, J ), J = 1, N );
      } // 110
      READ( NIN, FMT = * )( DTRU( I ), I = 1, N );
      READ( NIN, FMT = * )( DIFTRU( I ), I = 1, N );

      NPTKNT = NPTKNT + 1;

      // Compute eigenvalues/eigenvectors of (A, B).
      // Compute eigenvalue/eigenvector condition numbers
      // using computed eigenvectors.

      dlacpy('F', N, N, A, LDA, AI, LDA );
      dlacpy('F', N, N, B, LDA, BI, LDA );

      dggevx('N', 'V', 'V', 'B', N, AI, LDA, BI, LDA, ALPHAR, ALPHAI, BETA, VL, LDA, VR, LDA, ILO, IHI, LSCALE, RSCALE, ANORM, BNORM, S, DIF, WORK, LWORK, IWORK, BWORK, LINFO );

      if ( LINFO != 0 ) {
         RESULT( 1 ) = ULPINV;
         WRITE( NOUT, FMT = 9987 )'DGGEVX', LINFO, N, NPTKNT;
         GO TO 140;
      }

      // Compute the norm(A, B)

      dlacpy('Full', N, N, AI, LDA, WORK, N );
      dlacpy('Full', N, N, BI, LDA, WORK( N*N+1 ), N );
      ABNORM = DLANGE( 'Fro', N, 2*N, WORK, N, WORK );

      // Tests (1) and (2)

      RESULT( 1 ) = ZERO;
      dget52( true , N, A, LDA, B, LDA, VL, LDA, ALPHAR, ALPHAI, BETA, WORK, RESULT( 1 ) );
      if ( RESULT( 2 ) > THRESH ) {
         WRITE( NOUT, FMT = 9986 )'Left', 'DGGEVX', RESULT( 2 ), N, NPTKNT;
      }

      RESULT( 2 ) = ZERO;
      dget52( false , N, A, LDA, B, LDA, VR, LDA, ALPHAR, ALPHAI, BETA, WORK, RESULT( 2 ) );
      if ( RESULT( 3 ) > THRESH ) {
         WRITE( NOUT, FMT = 9986 )'Right', 'DGGEVX', RESULT( 3 ), N, NPTKNT;
      }

      // Test (3)

      RESULT( 3 ) = ZERO;
      for (I = 1; I <= N; I++) { // 120
         if ( S( I ) == ZERO ) {
            IF( DTRU( I ) > ABNORM*ULP ) RESULT( 3 ) = ULPINV;
         } else if ( DTRU( I ) == ZERO ) {
            IF( S( I ) > ABNORM*ULP ) RESULT( 3 ) = ULPINV;
         } else {
            WORK( I ) = MAX( ABS( DTRU( I ) / S( I ) ), ABS( S( I ) / DTRU( I ) ) );
            RESULT( 3 ) = MAX( RESULT( 3 ), WORK( I ) );
         }
      } // 120

      // Test (4)

      RESULT( 4 ) = ZERO;
      if ( DIF( 1 ) == ZERO ) {
         IF( DIFTRU( 1 ) > ABNORM*ULP ) RESULT( 4 ) = ULPINV;
      } else if ( DIFTRU( 1 ) == ZERO ) {
         IF( DIF( 1 ) > ABNORM*ULP ) RESULT( 4 ) = ULPINV;
      } else if ( DIF( 5 ) == ZERO ) {
         IF( DIFTRU( 5 ) > ABNORM*ULP ) RESULT( 4 ) = ULPINV;
      } else if ( DIFTRU( 5 ) == ZERO ) {
         IF( DIF( 5 ) > ABNORM*ULP ) RESULT( 4 ) = ULPINV;
      } else {
         RATIO1 = MAX( ABS( DIFTRU( 1 ) / DIF( 1 ) ), ABS( DIF( 1 ) / DIFTRU( 1 ) ) )          RATIO2 = MAX( ABS( DIFTRU( 5 ) / DIF( 5 ) ), ABS( DIF( 5 ) / DIFTRU( 5 ) ) );
         RESULT( 4 ) = MAX( RATIO1, RATIO2 );
      }

      NTESTT = NTESTT + 4;

      // Print out tests which fail.

      for (J = 1; J <= 4; J++) { // 130
         if ( RESULT( J ) >= THRSH2 ) {

            // If this is the first test to fail,
            // print a header to the data file.

            if ( NERRS == 0 ) {
               WRITE( NOUT, FMT = 9997 )'DXV';

               // Print out messages for built-in examples

               // Matrix types

               WRITE( NOUT, FMT = 9996 );

               // Tests performed

               WRITE( NOUT, FMT = 9992 )'''', 'transpose', '''';

            }
            NERRS = NERRS + 1;
            if ( RESULT( J ) < 10000.0 ) {
               WRITE( NOUT, FMT = 9989 )NPTKNT, N, J, RESULT( J );
            } else {
               WRITE( NOUT, FMT = 9988 )NPTKNT, N, J, RESULT( J );
            }
         }
      } // 130

      } // 140

      GO TO 90;
      } // 150

      // Summary

      alasvm('DXV', NOUT, NERRS, NTESTT, 0 );

      WORK( 1 ) = MAXWRK;

      return;

 9999 FORMAT( ' DDRGVX: ', A, ' returned INFO=', I6, '.', / 9X, 'N=', I6, ', JTYPE=', I6, ')' );

 9998 FORMAT( ' DDRGVX: ', A, ' Eigenvectors from ', A, ' incorrectly ', 'normalized.', / ' Bits of error=', 0P, G10.3, ',', 9X, 'N=', I6, ', JTYPE=', I6, ', IWA=', I5, ', IWB=', I5, ', IWX=', I5, ', IWY=', I5 );

 9997 FORMAT( / 1X, A3, ' -- Real Expert Eigenvalue/vector', ' problem driver' );

 9996 FORMAT( ' Input Example' );

 9995 FORMAT( ' Matrix types: ', / );

 9994 FORMAT( ' TYPE 1: Da is diagonal, Db is identity, ', / '     A = Y^(-H) Da X^(-1), B = Y^(-H) Db X^(-1) ', / '     YH and X are left and right eigenvectors. ', / );

 9993 FORMAT( ' TYPE 2: Da is quasi-diagonal, Db is identity, ', / '     A = Y^(-H) Da X^(-1), B = Y^(-H) Db X^(-1) ', / '     YH and X are left and right eigenvectors. ', / );

 9992 FORMAT( / ' Tests performed:  ', / 4X, ' a is alpha, b is beta, l is a left eigenvector, ', / 4X, ' r is a right eigenvector and ', A, ' means ', A, '.', / ' 1 = max | ( b A - a B )', A, ' l | / const.', / ' 2 = max | ( b A - a B ) r | / const.', / ' 3 = max ( Sest/Stru, Stru/Sest ) ', ' over all eigenvalues', / ' 4 = max( DIFest/DIFtru, DIFtru/DIFest ) ', ' over the 1st and 5th eigenvectors', / );

 9991 FORMAT( ' Type=', I2, ',', ' IWA=', I2, ', IWB=', I2, ', IWX=', I2, ', IWY=', I2, ', result ', I2, ' is', 0P, F8.2 );
 9990 FORMAT( ' Type=', I2, ',', ' IWA=', I2, ', IWB=', I2, ', IWX=', I2, ', IWY=', I2, ', result ', I2, ' is', 1P, D10.3 );
 9989 FORMAT( ' Input example #', I2, ', matrix order=', I4, ',', ' result ', I2, ' is', 0P, F8.2 );
 9988 FORMAT( ' Input example #', I2, ', matrix order=', I4, ',', ' result ', I2, ' is', 1P, D10.3 );
 9987 FORMAT( ' DDRGVX: ', A, ' returned INFO=', I6, '.', / 9X, 'N=', I6, ', Input example #', I2, ')' );

 9986 FORMAT( ' DDRGVX: ', A, ' Eigenvectors from ', A, ' incorrectly ', 'normalized.', / ' Bits of error=', 0P, G10.3, ',', 9X, 'N=', I6, ', Input Example #', I2, ')' );


      // End of DDRGVX

      }
