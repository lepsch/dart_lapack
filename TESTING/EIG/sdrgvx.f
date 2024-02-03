      SUBROUTINE SDRGVX( NSIZE, THRESH, NIN, NOUT, A, LDA, B, AI, BI, ALPHAR, ALPHAI, BETA, VL, VR, ILO, IHI, LSCALE, RSCALE, S, STRU, DIF, DIFTRU, WORK, LWORK, IWORK, LIWORK, RESULT, BWORK, INFO )

*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                IHI, ILO, INFO, LDA, LIWORK, LWORK, NIN, NOUT, NSIZE;
      REAL               THRESH
      // ..
      // .. Array Arguments ..
      bool               BWORK( * );
      int                IWORK( * );
      REAL               A( LDA, * ), AI( LDA, * ), ALPHAI( * ), ALPHAR( * ), B( LDA, * ), BETA( * ), BI( LDA, * ), DIF( * ), DIFTRU( * ), LSCALE( * ), RESULT( 4 ), RSCALE( * ), S( * ), STRU( * ), VL( LDA, * ), VR( LDA, * ), WORK( * )
      // ..

*  =====================================================================

      // .. Parameters ..
      REAL               ZERO, ONE, TEN, TNTH, HALF
      const              ZERO = 0.0E+0, ONE = 1.0E+0, TEN = 1.0E+1, TNTH = 1.0E-1, HALF = 0.5D+0 ;
      // ..
      // .. Local Scalars ..
      int                I, IPTYPE, IWA, IWB, IWX, IWY, J, LINFO, MAXWRK, MINWRK, N, NERRS, NMAX, NPTKNT, NTESTT;
      REAL               ABNORM, ANORM, BNORM, RATIO1, RATIO2, THRSH2, ULP, ULPINV;
      // ..
      // .. Local Arrays ..
      REAL               WEIGHT( 5 )
      // ..
      // .. External Functions ..
      int                ILAENV;
      REAL               SLAMCH, SLANGE
      // EXTERNAL ILAENV, SLAMCH, SLANGE
      // ..
      // .. External Subroutines ..
      // EXTERNAL ALASVM, SGET52, SGGEVX, SLACPY, SLATM6, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, MAX, SQRT
      // ..
      // .. Executable Statements ..

      // Check for errors

      INFO = 0

      NMAX = 5

      if ( NSIZE.LT.0 ) {
         INFO = -1
      } else if ( THRESH.LT.ZERO ) {
         INFO = -2
      } else if ( NIN.LE.0 ) {
         INFO = -3
      } else if ( NOUT.LE.0 ) {
         INFO = -4
      } else if ( LDA.LT.1 .OR. LDA.LT.NMAX ) {
         INFO = -6
      } else if ( LIWORK.LT.NMAX+6 ) {
         INFO = -26
      }

      // Compute workspace
       // (Note: Comments in the code beginning "Workspace:" describe the
        // minimal amount of workspace needed at that point in the code,
        // as well as the preferred amount for good performance.
        // NB refers to the optimal block size for the immediately
        // following subroutine, as returned by ILAENV.)

      MINWRK = 1
      if ( INFO.EQ.0 .AND. LWORK.GE.1 ) {
         MINWRK = 2*NMAX*NMAX + 12*NMAX + 16
         MAXWRK = 6*NMAX + NMAX*ILAENV( 1, 'SGEQRF', ' ', NMAX, 1, NMAX, 0 )
         MAXWRK = MAX( MAXWRK, 2*NMAX*NMAX+12*NMAX+16 )
         WORK( 1 ) = MAXWRK
      }

      IF( LWORK.LT.MINWRK ) INFO = -24

      if ( INFO.NE.0 ) {
         xerbla('SDRGVX', -INFO );
         RETURN
      }

      N = 5
      ULP = SLAMCH( 'P' )
      ULPINV = ONE / ULP
      THRSH2 = TEN*THRESH
      NERRS = 0
      NPTKNT = 0
      NTESTT = 0

      IF( NSIZE.EQ.0 ) GO TO 90

      // Parameters used for generating test matrices.

      WEIGHT( 1 ) = TNTH
      WEIGHT( 2 ) = HALF
      WEIGHT( 3 ) = ONE
      WEIGHT( 4 ) = ONE / WEIGHT( 2 )
      WEIGHT( 5 ) = ONE / WEIGHT( 1 )

      DO 80 IPTYPE = 1, 2
         DO 70 IWA = 1, 5
            DO 60 IWB = 1, 5
               DO 50 IWX = 1, 5
                  DO 40 IWY = 1, 5

                     // generated a test matrix pair

                     slatm6(IPTYPE, 5, A, LDA, B, VR, LDA, VL, LDA, WEIGHT( IWA ), WEIGHT( IWB ), WEIGHT( IWX ), WEIGHT( IWY ), STRU, DIFTRU );

                     // Compute eigenvalues/eigenvectors of (A, B).
                     // Compute eigenvalue/eigenvector condition numbers
                     // using computed eigenvectors.

                     slacpy('F', N, N, A, LDA, AI, LDA );
                     slacpy('F', N, N, B, LDA, BI, LDA );

                     sggevx('N', 'V', 'V', 'B', N, AI, LDA, BI, LDA, ALPHAR, ALPHAI, BETA, VL, LDA, VR, LDA, ILO, IHI, LSCALE, RSCALE, ANORM, BNORM, S, DIF, WORK, LWORK, IWORK, BWORK, LINFO );
                     if ( LINFO.NE.0 ) {
                        RESULT( 1 ) = ULPINV
                        WRITE( NOUT, FMT = 9999 )'SGGEVX', LINFO, N, IPTYPE
                        GO TO 30
                     }

                     // Compute the norm(A, B)

                     slacpy('Full', N, N, AI, LDA, WORK, N );
                     slacpy('Full', N, N, BI, LDA, WORK( N*N+1 ), N );
                     ABNORM = SLANGE( 'Fro', N, 2*N, WORK, N, WORK )

                     // Tests (1) and (2)

                     RESULT( 1 ) = ZERO
                     sget52(.TRUE., N, A, LDA, B, LDA, VL, LDA, ALPHAR, ALPHAI, BETA, WORK, RESULT( 1 ) );
                     if ( RESULT( 2 ).GT.THRESH ) {
                        WRITE( NOUT, FMT = 9998 )'Left', 'SGGEVX', RESULT( 2 ), N, IPTYPE, IWA, IWB, IWX, IWY
                     }

                     RESULT( 2 ) = ZERO
                     sget52(.FALSE., N, A, LDA, B, LDA, VR, LDA, ALPHAR, ALPHAI, BETA, WORK, RESULT( 2 ) );
                     if ( RESULT( 3 ).GT.THRESH ) {
                        WRITE( NOUT, FMT = 9998 )'Right', 'SGGEVX', RESULT( 3 ), N, IPTYPE, IWA, IWB, IWX, IWY
                     }

                     // Test (3)

                     RESULT( 3 ) = ZERO
                     DO 10 I = 1, N
                        if ( S( I ).EQ.ZERO ) {
                           IF( STRU( I ).GT.ABNORM*ULP ) RESULT( 3 ) = ULPINV
                        } else if ( STRU( I ).EQ.ZERO ) {
                           IF( S( I ).GT.ABNORM*ULP ) RESULT( 3 ) = ULPINV
                        } else {
                           WORK( I ) = MAX( ABS( STRU( I ) / S( I ) ), ABS( S( I ) / STRU( I ) ) )
                           RESULT( 3 ) = MAX( RESULT( 3 ), WORK( I ) )
                        }
   10                CONTINUE

                     // Test (4)

                     RESULT( 4 ) = ZERO
                     if ( DIF( 1 ).EQ.ZERO ) {
                        IF( DIFTRU( 1 ).GT.ABNORM*ULP ) RESULT( 4 ) = ULPINV
                     } else if ( DIFTRU( 1 ).EQ.ZERO ) {
                        IF( DIF( 1 ).GT.ABNORM*ULP ) RESULT( 4 ) = ULPINV
                     } else if ( DIF( 5 ).EQ.ZERO ) {
                        IF( DIFTRU( 5 ).GT.ABNORM*ULP ) RESULT( 4 ) = ULPINV
                     } else if ( DIFTRU( 5 ).EQ.ZERO ) {
                        IF( DIF( 5 ).GT.ABNORM*ULP ) RESULT( 4 ) = ULPINV
                     } else {
                        RATIO1 = MAX( ABS( DIFTRU( 1 ) / DIF( 1 ) ), ABS( DIF( 1 ) / DIFTRU( 1 ) ) )                         RATIO2 = MAX( ABS( DIFTRU( 5 ) / DIF( 5 ) ), ABS( DIF( 5 ) / DIFTRU( 5 ) ) )
                        RESULT( 4 ) = MAX( RATIO1, RATIO2 )
                     }

                     NTESTT = NTESTT + 4

                     // Print out tests which fail.

                     DO 20 J = 1, 4
                        if ( ( RESULT( J ).GE.THRSH2 .AND. J.GE.4 ) .OR. ( RESULT( J ).GE.THRESH .AND. J.LE.3 ) ) {

                        // If this is the first test to fail,
                        // print a header to the data file.

                           if ( NERRS.EQ.0 ) {
                              WRITE( NOUT, FMT = 9997 )'SXV'

                           // Print out messages for built-in examples

                           // Matrix types

                              WRITE( NOUT, FMT = 9995 )
                              WRITE( NOUT, FMT = 9994 )
                              WRITE( NOUT, FMT = 9993 )

                           // Tests performed

                              WRITE( NOUT, FMT = 9992 )'''', 'transpose', ''''

                           }
                           NERRS = NERRS + 1
                           if ( RESULT( J ).LT.10000.0 ) {
                              WRITE( NOUT, FMT = 9991 )IPTYPE, IWA, IWB, IWX, IWY, J, RESULT( J )
                           } else {
                              WRITE( NOUT, FMT = 9990 )IPTYPE, IWA, IWB, IWX, IWY, J, RESULT( J )
                           }
                        }
   20                CONTINUE

   30                CONTINUE

   40             CONTINUE
   50          CONTINUE
   60       CONTINUE
   70    CONTINUE
   80 CONTINUE

      GO TO 150

   90 CONTINUE

      // Read in data from file to check accuracy of condition estimation
      // Read input data until N=0

      READ( NIN, FMT = *, END = 150 )N
      IF( N.EQ.0 ) GO TO 150
      DO 100 I = 1, N
         READ( NIN, FMT = * )( A( I, J ), J = 1, N )
  100 CONTINUE
      DO 110 I = 1, N
         READ( NIN, FMT = * )( B( I, J ), J = 1, N )
  110 CONTINUE
      READ( NIN, FMT = * )( STRU( I ), I = 1, N )
      READ( NIN, FMT = * )( DIFTRU( I ), I = 1, N )

      NPTKNT = NPTKNT + 1

      // Compute eigenvalues/eigenvectors of (A, B).
      // Compute eigenvalue/eigenvector condition numbers
      // using computed eigenvectors.

      slacpy('F', N, N, A, LDA, AI, LDA );
      slacpy('F', N, N, B, LDA, BI, LDA );

      sggevx('N', 'V', 'V', 'B', N, AI, LDA, BI, LDA, ALPHAR, ALPHAI, BETA, VL, LDA, VR, LDA, ILO, IHI, LSCALE, RSCALE, ANORM, BNORM, S, DIF, WORK, LWORK, IWORK, BWORK, LINFO );

      if ( LINFO.NE.0 ) {
         RESULT( 1 ) = ULPINV
         WRITE( NOUT, FMT = 9987 )'SGGEVX', LINFO, N, NPTKNT
         GO TO 140
      }

      // Compute the norm(A, B)

      slacpy('Full', N, N, AI, LDA, WORK, N );
      slacpy('Full', N, N, BI, LDA, WORK( N*N+1 ), N );
      ABNORM = SLANGE( 'Fro', N, 2*N, WORK, N, WORK )

      // Tests (1) and (2)

      RESULT( 1 ) = ZERO
      sget52(.TRUE., N, A, LDA, B, LDA, VL, LDA, ALPHAR, ALPHAI, BETA, WORK, RESULT( 1 ) );
      if ( RESULT( 2 ).GT.THRESH ) {
         WRITE( NOUT, FMT = 9986 )'Left', 'SGGEVX', RESULT( 2 ), N, NPTKNT
      }

      RESULT( 2 ) = ZERO
      sget52(.FALSE., N, A, LDA, B, LDA, VR, LDA, ALPHAR, ALPHAI, BETA, WORK, RESULT( 2 ) );
      if ( RESULT( 3 ).GT.THRESH ) {
         WRITE( NOUT, FMT = 9986 )'Right', 'SGGEVX', RESULT( 3 ), N, NPTKNT
      }

      // Test (3)

      RESULT( 3 ) = ZERO
      DO 120 I = 1, N
         if ( S( I ).EQ.ZERO ) {
            IF( STRU( I ).GT.ABNORM*ULP ) RESULT( 3 ) = ULPINV
         } else if ( STRU( I ).EQ.ZERO ) {
            IF( S( I ).GT.ABNORM*ULP ) RESULT( 3 ) = ULPINV
         } else {
            WORK( I ) = MAX( ABS( STRU( I ) / S( I ) ), ABS( S( I ) / STRU( I ) ) )
            RESULT( 3 ) = MAX( RESULT( 3 ), WORK( I ) )
         }
  120 CONTINUE

      // Test (4)

      RESULT( 4 ) = ZERO
      if ( DIF( 1 ).EQ.ZERO ) {
         IF( DIFTRU( 1 ).GT.ABNORM*ULP ) RESULT( 4 ) = ULPINV
      } else if ( DIFTRU( 1 ).EQ.ZERO ) {
         IF( DIF( 1 ).GT.ABNORM*ULP ) RESULT( 4 ) = ULPINV
      } else if ( DIF( 5 ).EQ.ZERO ) {
         IF( DIFTRU( 5 ).GT.ABNORM*ULP ) RESULT( 4 ) = ULPINV
      } else if ( DIFTRU( 5 ).EQ.ZERO ) {
         IF( DIF( 5 ).GT.ABNORM*ULP ) RESULT( 4 ) = ULPINV
      } else {
         RATIO1 = MAX( ABS( DIFTRU( 1 ) / DIF( 1 ) ), ABS( DIF( 1 ) / DIFTRU( 1 ) ) )          RATIO2 = MAX( ABS( DIFTRU( 5 ) / DIF( 5 ) ), ABS( DIF( 5 ) / DIFTRU( 5 ) ) )
         RESULT( 4 ) = MAX( RATIO1, RATIO2 )
      }

      NTESTT = NTESTT + 4

      // Print out tests which fail.

      DO 130 J = 1, 4
         if ( RESULT( J ).GE.THRSH2 ) {

            // If this is the first test to fail,
            // print a header to the data file.

            if ( NERRS.EQ.0 ) {
               WRITE( NOUT, FMT = 9997 )'SXV'

               // Print out messages for built-in examples

               // Matrix types

               WRITE( NOUT, FMT = 9996 )

               // Tests performed

               WRITE( NOUT, FMT = 9992 )'''', 'transpose', ''''

            }
            NERRS = NERRS + 1
            if ( RESULT( J ).LT.10000.0 ) {
               WRITE( NOUT, FMT = 9989 )NPTKNT, N, J, RESULT( J )
            } else {
               WRITE( NOUT, FMT = 9988 )NPTKNT, N, J, RESULT( J )
            }
         }
  130 CONTINUE

  140 CONTINUE

      GO TO 90
  150 CONTINUE

      // Summary

      alasvm('SXV', NOUT, NERRS, NTESTT, 0 );

      WORK( 1 ) = MAXWRK

      RETURN

 9999 FORMAT( ' SDRGVX: ', A, ' returned INFO=', I6, '.', / 9X, 'N=', I6, ', JTYPE=', I6, ')' )

 9998 FORMAT( ' SDRGVX: ', A, ' Eigenvectors from ', A, ' incorrectly ', 'normalized.', / ' Bits of error=', 0P, G10.3, ',', 9X, 'N=', I6, ', JTYPE=', I6, ', IWA=', I5, ', IWB=', I5, ', IWX=', I5, ', IWY=', I5 )

 9997 FORMAT( / 1X, A3, ' -- Real Expert Eigenvalue/vector', ' problem driver' )

 9996 FORMAT( ' Input Example' )

 9995 FORMAT( ' Matrix types: ', / )

 9994 FORMAT( ' TYPE 1: Da is diagonal, Db is identity, ', / '     A = Y^(-H) Da X^(-1), B = Y^(-H) Db X^(-1) ', / '     YH and X are left and right eigenvectors. ', / )

 9993 FORMAT( ' TYPE 2: Da is quasi-diagonal, Db is identity, ', / '     A = Y^(-H) Da X^(-1), B = Y^(-H) Db X^(-1) ', / '     YH and X are left and right eigenvectors. ', / )

 9992 FORMAT( / ' Tests performed:  ', / 4X, ' a is alpha, b is beta, l is a left eigenvector, ', / 4X, ' r is a right eigenvector and ', A, ' means ', A, '.', / ' 1 = max | ( b A - a B )', A, ' l | / const.', / ' 2 = max | ( b A - a B ) r | / const.', / ' 3 = max ( Sest/Stru, Stru/Sest ) ', ' over all eigenvalues', / ' 4 = max( DIFest/DIFtru, DIFtru/DIFest ) ', ' over the 1st and 5th eigenvectors', / )

 9991 FORMAT( ' Type=', I2, ',', ' IWA=', I2, ', IWB=', I2, ', IWX=', I2, ', IWY=', I2, ', result ', I2, ' is', 0P, F8.2 )
 9990 FORMAT( ' Type=', I2, ',', ' IWA=', I2, ', IWB=', I2, ', IWX=', I2, ', IWY=', I2, ', result ', I2, ' is', 1P, E10.3 )
 9989 FORMAT( ' Input example #', I2, ', matrix order=', I4, ',', ' result ', I2, ' is', 0P, F8.2 )
 9988 FORMAT( ' Input example #', I2, ', matrix order=', I4, ',', ' result ', I2, ' is', 1P, E10.3 )
 9987 FORMAT( ' SDRGVX: ', A, ' returned INFO=', I6, '.', / 9X, 'N=', I6, ', Input example #', I2, ')' )

 9986 FORMAT( ' SDRGVX: ', A, ' Eigenvectors from ', A, ' incorrectly ', 'normalized.', / ' Bits of error=', 0P, G10.3, ',', 9X, 'N=', I6, ', Input Example #', I2, ')' )


      // End of SDRGVX

      }
