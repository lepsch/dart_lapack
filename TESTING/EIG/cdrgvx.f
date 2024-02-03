      SUBROUTINE CDRGVX( NSIZE, THRESH, NIN, NOUT, A, LDA, B, AI, BI, ALPHA, BETA, VL, VR, ILO, IHI, LSCALE, RSCALE, S, STRU, DIF, DIFTRU, WORK, LWORK, RWORK, IWORK, LIWORK, RESULT, BWORK, INFO )
*
*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      int                IHI, ILO, INFO, LDA, LIWORK, LWORK, NIN, NOUT, NSIZE;
      REAL               THRESH
*     ..
*     .. Array Arguments ..
      bool               BWORK( * );
      int                IWORK( * );
      REAL               DIF( * ), DIFTRU( * ), LSCALE( * ), RESULT( 4 ), RSCALE( * ), RWORK( * ), S( * ), STRU( * )       COMPLEX            A( LDA, * ), AI( LDA, * ), ALPHA( * ), B( LDA, * ), BETA( * ), BI( LDA, * ), VL( LDA, * ), VR( LDA, * ), WORK( * )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      REAL               ZERO, ONE, TEN, TNTH, HALF
      PARAMETER          ( ZERO = 0.0E+0, ONE = 1.0E+0, TEN = 1.0E+1, TNTH = 1.0E-1, HALF = 0.5E+0 )
*     ..
*     .. Local Scalars ..
      int                I, IPTYPE, IWA, IWB, IWX, IWY, J, LINFO, MAXWRK, MINWRK, N, NERRS, NMAX, NPTKNT, NTESTT       REAL               ABNORM, ANORM, BNORM, RATIO1, RATIO2, THRSH2, ULP, ULPINV;
*     ..
*     .. Local Arrays ..
      COMPLEX            WEIGHT( 5 )
*     ..
*     .. External Functions ..
      int                ILAENV;
      REAL               CLANGE, SLAMCH
      // EXTERNAL ILAENV, CLANGE, SLAMCH
*     ..
*     .. External Subroutines ..
      // EXTERNAL ALASVM, CGET52, CGGEVX, CLACPY, CLATM6, XERBLA
*     ..
*     .. Intrinsic Functions ..
      // INTRINSIC ABS, CMPLX, MAX, SQRT
*     ..
*     .. Executable Statements ..
*
*     Check for errors
*
      INFO = 0
*
      NMAX = 5
*
      IF( NSIZE.LT.0 ) THEN
         INFO = -1
      ELSE IF( THRESH.LT.ZERO ) THEN
         INFO = -2
      ELSE IF( NIN.LE.0 ) THEN
         INFO = -3
      ELSE IF( NOUT.LE.0 ) THEN
         INFO = -4
      ELSE IF( LDA.LT.1 .OR. LDA.LT.NMAX ) THEN
         INFO = -6
      ELSE IF( LIWORK.LT.NMAX+2 ) THEN
         INFO = -26
      END IF
*
*     Compute workspace
*      (Note: Comments in the code beginning "Workspace:" describe the
*       minimal amount of workspace needed at that point in the code,
*       as well as the preferred amount for good performance.
*       NB refers to the optimal block size for the immediately
*       following subroutine, as returned by ILAENV.)
*
      MINWRK = 1
      IF( INFO.EQ.0 .AND. LWORK.GE.1 ) THEN
         MINWRK = 2*NMAX*( NMAX+1 )
         MAXWRK = NMAX*( 1+ILAENV( 1, 'CGEQRF', ' ', NMAX, 1, NMAX, 0 ) )
         MAXWRK = MAX( MAXWRK, 2*NMAX*( NMAX+1 ) )
         WORK( 1 ) = MAXWRK
      END IF
*
      IF( LWORK.LT.MINWRK ) INFO = -23
*
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'CDRGVX', -INFO )
         RETURN
      END IF
*
      N = 5
      ULP = SLAMCH( 'P' )
      ULPINV = ONE / ULP
      THRSH2 = TEN*THRESH
      NERRS = 0
      NPTKNT = 0
      NTESTT = 0
*
      IF( NSIZE.EQ.0 ) GO TO 90
*
*     Parameters used for generating test matrices.
*
      WEIGHT( 1 ) = CMPLX( TNTH, ZERO )
      WEIGHT( 2 ) = CMPLX( HALF, ZERO )
      WEIGHT( 3 ) = ONE
      WEIGHT( 4 ) = ONE / WEIGHT( 2 )
      WEIGHT( 5 ) = ONE / WEIGHT( 1 )
*
      DO 80 IPTYPE = 1, 2
         DO 70 IWA = 1, 5
            DO 60 IWB = 1, 5
               DO 50 IWX = 1, 5
                  DO 40 IWY = 1, 5
*
*                    generated a pair of test matrix
*
                     CALL CLATM6( IPTYPE, 5, A, LDA, B, VR, LDA, VL, LDA, WEIGHT( IWA ), WEIGHT( IWB ), WEIGHT( IWX ), WEIGHT( IWY ), STRU, DIFTRU )
*
*                    Compute eigenvalues/eigenvectors of (A, B).
*                    Compute eigenvalue/eigenvector condition numbers
*                    using computed eigenvectors.
*
                     CALL CLACPY( 'F', N, N, A, LDA, AI, LDA )
                     CALL CLACPY( 'F', N, N, B, LDA, BI, LDA )
*
                     CALL CGGEVX( 'N', 'V', 'V', 'B', N, AI, LDA, BI, LDA, ALPHA, BETA, VL, LDA, VR, LDA, ILO, IHI, LSCALE, RSCALE, ANORM, BNORM, S, DIF, WORK, LWORK, RWORK, IWORK, BWORK, LINFO )
                     IF( LINFO.NE.0 ) THEN
                        WRITE( NOUT, FMT = 9999 )'CGGEVX', LINFO, N, IPTYPE, IWA, IWB, IWX, IWY
                        GO TO 30
                     END IF
*
*                    Compute the norm(A, B)
*
                     CALL CLACPY( 'Full', N, N, AI, LDA, WORK, N )
                     CALL CLACPY( 'Full', N, N, BI, LDA, WORK( N*N+1 ), N )
                     ABNORM = CLANGE( 'Fro', N, 2*N, WORK, N, RWORK )
*
*                    Tests (1) and (2)
*
                     RESULT( 1 ) = ZERO
                     CALL CGET52( .TRUE., N, A, LDA, B, LDA, VL, LDA, ALPHA, BETA, WORK, RWORK, RESULT( 1 ) )
                     IF( RESULT( 2 ).GT.THRESH ) THEN
                        WRITE( NOUT, FMT = 9998 )'Left', 'CGGEVX', RESULT( 2 ), N, IPTYPE, IWA, IWB, IWX, IWY
                     END IF
*
                     RESULT( 2 ) = ZERO
                     CALL CGET52( .FALSE., N, A, LDA, B, LDA, VR, LDA, ALPHA, BETA, WORK, RWORK, RESULT( 2 ) )
                     IF( RESULT( 3 ).GT.THRESH ) THEN
                        WRITE( NOUT, FMT = 9998 )'Right', 'CGGEVX', RESULT( 3 ), N, IPTYPE, IWA, IWB, IWX, IWY
                     END IF
*
*                    Test (3)
*
                     RESULT( 3 ) = ZERO
                     DO 10 I = 1, N
                        IF( S( I ).EQ.ZERO ) THEN
                           IF( STRU( I ).GT.ABNORM*ULP ) RESULT( 3 ) = ULPINV
                        ELSE IF( STRU( I ).EQ.ZERO ) THEN
                           IF( S( I ).GT.ABNORM*ULP ) RESULT( 3 ) = ULPINV
                        ELSE
                           RWORK( I ) = MAX( ABS( STRU( I ) / S( I ) ), ABS( S( I ) / STRU( I ) ) )
                           RESULT( 3 ) = MAX( RESULT( 3 ), RWORK( I ) )
                        END IF
   10                CONTINUE
*
*                    Test (4)
*
                     RESULT( 4 ) = ZERO
                     IF( DIF( 1 ).EQ.ZERO ) THEN
                        IF( DIFTRU( 1 ).GT.ABNORM*ULP ) RESULT( 4 ) = ULPINV
                     ELSE IF( DIFTRU( 1 ).EQ.ZERO ) THEN
                        IF( DIF( 1 ).GT.ABNORM*ULP ) RESULT( 4 ) = ULPINV
                     ELSE IF( DIF( 5 ).EQ.ZERO ) THEN
                        IF( DIFTRU( 5 ).GT.ABNORM*ULP ) RESULT( 4 ) = ULPINV
                     ELSE IF( DIFTRU( 5 ).EQ.ZERO ) THEN
                        IF( DIF( 5 ).GT.ABNORM*ULP ) RESULT( 4 ) = ULPINV
                     ELSE
                        RATIO1 = MAX( ABS( DIFTRU( 1 ) / DIF( 1 ) ), ABS( DIF( 1 ) / DIFTRU( 1 ) ) )                         RATIO2 = MAX( ABS( DIFTRU( 5 ) / DIF( 5 ) ), ABS( DIF( 5 ) / DIFTRU( 5 ) ) )
                        RESULT( 4 ) = MAX( RATIO1, RATIO2 )
                     END IF
*
                     NTESTT = NTESTT + 4
*
*                    Print out tests which fail.
*
                     DO 20 J = 1, 4
                        IF( ( RESULT( J ).GE.THRSH2 .AND. J.GE.4 ) .OR. ( RESULT( J ).GE.THRESH .AND. J.LE.3 ) ) THEN
*
*                       If this is the first test to fail,
*                       print a header to the data file.
*
                           IF( NERRS.EQ.0 ) THEN
                              WRITE( NOUT, FMT = 9997 )'CXV'
*
*                          Print out messages for built-in examples
*
*                          Matrix types
*
                              WRITE( NOUT, FMT = 9995 )
                              WRITE( NOUT, FMT = 9994 )
                              WRITE( NOUT, FMT = 9993 )
*
*                          Tests performed
*
                              WRITE( NOUT, FMT = 9992 )'''', 'transpose', ''''
*
                           END IF
                           NERRS = NERRS + 1
                           IF( RESULT( J ).LT.10000.0 ) THEN
                              WRITE( NOUT, FMT = 9991 )IPTYPE, IWA, IWB, IWX, IWY, J, RESULT( J )
                           ELSE
                              WRITE( NOUT, FMT = 9990 )IPTYPE, IWA, IWB, IWX, IWY, J, RESULT( J )
                           END IF
                        END IF
   20                CONTINUE
*
   30                CONTINUE
*
   40             CONTINUE
   50          CONTINUE
   60       CONTINUE
   70    CONTINUE
   80 CONTINUE
*
      GO TO 150
*
   90 CONTINUE
*
*     Read in data from file to check accuracy of condition estimation
*     Read input data until N=0
*
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
*
      NPTKNT = NPTKNT + 1
*
*     Compute eigenvalues/eigenvectors of (A, B).
*     Compute eigenvalue/eigenvector condition numbers
*     using computed eigenvectors.
*
      CALL CLACPY( 'F', N, N, A, LDA, AI, LDA )
      CALL CLACPY( 'F', N, N, B, LDA, BI, LDA )
*
      CALL CGGEVX( 'N', 'V', 'V', 'B', N, AI, LDA, BI, LDA, ALPHA, BETA, VL, LDA, VR, LDA, ILO, IHI, LSCALE, RSCALE, ANORM, BNORM, S, DIF, WORK, LWORK, RWORK, IWORK, BWORK, LINFO )
*
      IF( LINFO.NE.0 ) THEN
         WRITE( NOUT, FMT = 9987 )'CGGEVX', LINFO, N, NPTKNT
         GO TO 140
      END IF
*
*     Compute the norm(A, B)
*
      CALL CLACPY( 'Full', N, N, AI, LDA, WORK, N )
      CALL CLACPY( 'Full', N, N, BI, LDA, WORK( N*N+1 ), N )
      ABNORM = CLANGE( 'Fro', N, 2*N, WORK, N, RWORK )
*
*     Tests (1) and (2)
*
      RESULT( 1 ) = ZERO
      CALL CGET52( .TRUE., N, A, LDA, B, LDA, VL, LDA, ALPHA, BETA, WORK, RWORK, RESULT( 1 ) )
      IF( RESULT( 2 ).GT.THRESH ) THEN
         WRITE( NOUT, FMT = 9986 )'Left', 'CGGEVX', RESULT( 2 ), N, NPTKNT
      END IF
*
      RESULT( 2 ) = ZERO
      CALL CGET52( .FALSE., N, A, LDA, B, LDA, VR, LDA, ALPHA, BETA, WORK, RWORK, RESULT( 2 ) )
      IF( RESULT( 3 ).GT.THRESH ) THEN
         WRITE( NOUT, FMT = 9986 )'Right', 'CGGEVX', RESULT( 3 ), N, NPTKNT
      END IF
*
*     Test (3)
*
      RESULT( 3 ) = ZERO
      DO 120 I = 1, N
         IF( S( I ).EQ.ZERO ) THEN
            IF( STRU( I ).GT.ABNORM*ULP ) RESULT( 3 ) = ULPINV
         ELSE IF( STRU( I ).EQ.ZERO ) THEN
            IF( S( I ).GT.ABNORM*ULP ) RESULT( 3 ) = ULPINV
         ELSE
            RWORK( I ) = MAX( ABS( STRU( I ) / S( I ) ), ABS( S( I ) / STRU( I ) ) )
            RESULT( 3 ) = MAX( RESULT( 3 ), RWORK( I ) )
         END IF
  120 CONTINUE
*
*     Test (4)
*
      RESULT( 4 ) = ZERO
      IF( DIF( 1 ).EQ.ZERO ) THEN
         IF( DIFTRU( 1 ).GT.ABNORM*ULP ) RESULT( 4 ) = ULPINV
      ELSE IF( DIFTRU( 1 ).EQ.ZERO ) THEN
         IF( DIF( 1 ).GT.ABNORM*ULP ) RESULT( 4 ) = ULPINV
      ELSE IF( DIF( 5 ).EQ.ZERO ) THEN
         IF( DIFTRU( 5 ).GT.ABNORM*ULP ) RESULT( 4 ) = ULPINV
      ELSE IF( DIFTRU( 5 ).EQ.ZERO ) THEN
         IF( DIF( 5 ).GT.ABNORM*ULP ) RESULT( 4 ) = ULPINV
      ELSE
         RATIO1 = MAX( ABS( DIFTRU( 1 ) / DIF( 1 ) ), ABS( DIF( 1 ) / DIFTRU( 1 ) ) )          RATIO2 = MAX( ABS( DIFTRU( 5 ) / DIF( 5 ) ), ABS( DIF( 5 ) / DIFTRU( 5 ) ) )
         RESULT( 4 ) = MAX( RATIO1, RATIO2 )
      END IF
*
      NTESTT = NTESTT + 4
*
*     Print out tests which fail.
*
      DO 130 J = 1, 4
         IF( RESULT( J ).GE.THRSH2 ) THEN
*
*           If this is the first test to fail,
*           print a header to the data file.
*
            IF( NERRS.EQ.0 ) THEN
               WRITE( NOUT, FMT = 9997 )'CXV'
*
*              Print out messages for built-in examples
*
*              Matrix types
*
               WRITE( NOUT, FMT = 9996 )
*
*              Tests performed
*
               WRITE( NOUT, FMT = 9992 )'''', 'transpose', ''''
*
            END IF
            NERRS = NERRS + 1
            IF( RESULT( J ).LT.10000.0 ) THEN
               WRITE( NOUT, FMT = 9989 )NPTKNT, N, J, RESULT( J )
            ELSE
               WRITE( NOUT, FMT = 9988 )NPTKNT, N, J, RESULT( J )
            END IF
         END IF
  130 CONTINUE
*
  140 CONTINUE
*
      GO TO 90
  150 CONTINUE
*
*     Summary
*
      CALL ALASVM( 'CXV', NOUT, NERRS, NTESTT, 0 )
*
      WORK( 1 ) = MAXWRK
*
      RETURN
*
 9999 FORMAT( ' CDRGVX: ', A, ' returned INFO=', I6, '.', / 9X, 'N=',
     $      I6, ', JTYPE=', I6, ')' )
*
 9998 FORMAT( ' CDRGVX: ', A, ' Eigenvectors from ', A, ' incorrectly ',
     $      'normalized.', / ' Bits of error=', 0P, G10.3, ',', 9X,
     $      'N=', I6, ', JTYPE=', I6, ', IWA=', I5, ', IWB=', I5,
     $      ', IWX=', I5, ', IWY=', I5 )
*
 9997 FORMAT( / 1X, A3, ' -- Complex Expert Eigenvalue/vector',
     $      ' problem driver' )
*
 9996 FORMAT( 'Input Example' )
*
 9995 FORMAT( ' Matrix types: ', / )
*
 9994 FORMAT( ' TYPE 1: Da is diagonal, Db is identity, ',
     $      / '     A = Y^(-H) Da X^(-1), B = Y^(-H) Db X^(-1) ',
     $      / '     YH and X are left and right eigenvectors. ', / )
*
 9993 FORMAT( ' TYPE 2: Da is quasi-diagonal, Db is identity, ',
     $      / '     A = Y^(-H) Da X^(-1), B = Y^(-H) Db X^(-1) ',
     $      / '     YH and X are left and right eigenvectors. ', / )
*
 9992 FORMAT( / ' Tests performed:  ', / 4X,
     $      ' a is alpha, b is beta, l is a left eigenvector, ', / 4X,
     $      ' r is a right eigenvector and ', A, ' means ', A, '.',
     $      / ' 1 = max | ( b A - a B )', A, ' l | / const.',
     $      / ' 2 = max | ( b A - a B ) r | / const.',
     $      / ' 3 = max ( Sest/Stru, Stru/Sest ) ',
     $      ' over all eigenvalues', /
     $      ' 4 = max( DIFest/DIFtru, DIFtru/DIFest ) ',
     $      ' over the 1st and 5th eigenvectors', / )
*
 9991 FORMAT( ' Type=', I2, ',', ' IWA=', I2, ', IWB=', I2, ', IWX=',
     $      I2, ', IWY=', I2, ', result ', I2, ' is', 0P, F8.2 )
*
 9990 FORMAT( ' Type=', I2, ',', ' IWA=', I2, ', IWB=', I2, ', IWX=',
     $      I2, ', IWY=', I2, ', result ', I2, ' is', 1P, E10.3 )
*
 9989 FORMAT( ' Input example #', I2, ', matrix order=', I4, ',',
     $      ' result ', I2, ' is', 0P, F8.2 )
*
 9988 FORMAT( ' Input example #', I2, ', matrix order=', I4, ',',
     $      ' result ', I2, ' is', 1P, E10.3 )
*
 9987 FORMAT( ' CDRGVX: ', A, ' returned INFO=', I6, '.', / 9X, 'N=',
     $      I6, ', Input example #', I2, ')' )
*
 9986 FORMAT( ' CDRGVX: ', A, ' Eigenvectors from ', A, ' incorrectly ',
     $      'normalized.', / ' Bits of error=', 0P, G10.3, ',', 9X,
     $      'N=', I6, ', Input Example #', I2, ')' )
*
*     End of CDRGVX
*
      END
