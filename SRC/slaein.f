      SUBROUTINE SLAEIN( RIGHTV, NOINIT, N, H, LDH, WR, WI, VR, VI, B, LDB, WORK, EPS3, SMLNUM, BIGNUM, INFO )

*  -- LAPACK auxiliary routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      bool               NOINIT, RIGHTV;
      int                INFO, LDB, LDH, N;
      REAL               BIGNUM, EPS3, SMLNUM, WI, WR
      // ..
      // .. Array Arguments ..
      REAL               B( LDB, * ), H( LDH, * ), VI( * ), VR( * ), WORK( * )
      // ..

*  =====================================================================

      // .. Parameters ..
      REAL               ZERO, ONE, TENTH
      const              ZERO = 0.0E+0, ONE = 1.0E+0, TENTH = 1.0E-1 ;
      // ..
      // .. Local Scalars ..
      String             NORMIN, TRANS;
      int                I, I1, I2, I3, IERR, ITS, J;
      REAL               ABSBII, ABSBJJ, EI, EJ, GROWTO, NORM, NRMSML, REC, ROOTN, SCALE, TEMP, VCRIT, VMAX, VNORM, W, W1, X, XI, XR, Y
      // ..
      // .. External Functions ..
      int                ISAMAX;
      REAL               SASUM, SLAPY2, SNRM2
      // EXTERNAL ISAMAX, SASUM, SLAPY2, SNRM2
      // ..
      // .. External Subroutines ..
      // EXTERNAL SLADIV, SLATRS, SSCAL
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, MAX, REAL, SQRT
      // ..
      // .. Executable Statements ..

      INFO = 0

      // GROWTO is the threshold used in the acceptance test for an
      // eigenvector.

      ROOTN = SQRT( REAL( N ) )
      GROWTO = TENTH / ROOTN
      NRMSML = MAX( ONE, EPS3*ROOTN )*SMLNUM

      // Form B = H - (WR,WI)*I (except that the subdiagonal elements and
      // the imaginary parts of the diagonal elements are not stored).

      DO 20 J = 1, N
         DO 10 I = 1, J - 1
            B( I, J ) = H( I, J )
   10    CONTINUE
         B( J, J ) = H( J, J ) - WR
   20 CONTINUE

      if ( WI.EQ.ZERO ) {

         // Real eigenvalue.

         if ( NOINIT ) {

            // Set initial vector.

            DO 30 I = 1, N
               VR( I ) = EPS3
   30       CONTINUE
         } else {

            // Scale supplied initial vector.

            VNORM = SNRM2( N, VR, 1 )
            sscal(N, ( EPS3*ROOTN ) / MAX( VNORM, NRMSML ), VR, 1 );
         }

         if ( RIGHTV ) {

            // LU decomposition with partial pivoting of B, replacing zero
            // pivots by EPS3.

            DO 60 I = 1, N - 1
               EI = H( I+1, I )
               if ( ABS( B( I, I ) ).LT.ABS( EI ) ) {

                  // Interchange rows and eliminate.

                  X = B( I, I ) / EI
                  B( I, I ) = EI
                  DO 40 J = I + 1, N
                     TEMP = B( I+1, J )
                     B( I+1, J ) = B( I, J ) - X*TEMP
                     B( I, J ) = TEMP
   40             CONTINUE
               } else {

                  // Eliminate without interchange.

                  IF( B( I, I ).EQ.ZERO ) B( I, I ) = EPS3
                  X = EI / B( I, I )
                  if ( X.NE.ZERO ) {
                     DO 50 J = I + 1, N
                        B( I+1, J ) = B( I+1, J ) - X*B( I, J )
   50                CONTINUE
                  }
               }
   60       CONTINUE
            IF( B( N, N ).EQ.ZERO ) B( N, N ) = EPS3

            TRANS = 'N'

         } else {

            // UL decomposition with partial pivoting of B, replacing zero
            // pivots by EPS3.

            DO 90 J = N, 2, -1
               EJ = H( J, J-1 )
               if ( ABS( B( J, J ) ).LT.ABS( EJ ) ) {

                  // Interchange columns and eliminate.

                  X = B( J, J ) / EJ
                  B( J, J ) = EJ
                  DO 70 I = 1, J - 1
                     TEMP = B( I, J-1 )
                     B( I, J-1 ) = B( I, J ) - X*TEMP
                     B( I, J ) = TEMP
   70             CONTINUE
               } else {

                  // Eliminate without interchange.

                  IF( B( J, J ).EQ.ZERO ) B( J, J ) = EPS3
                  X = EJ / B( J, J )
                  if ( X.NE.ZERO ) {
                     DO 80 I = 1, J - 1
                        B( I, J-1 ) = B( I, J-1 ) - X*B( I, J )
   80                CONTINUE
                  }
               }
   90       CONTINUE
            IF( B( 1, 1 ).EQ.ZERO ) B( 1, 1 ) = EPS3

            TRANS = 'T'

         }

         NORMIN = 'N'
         DO 110 ITS = 1, N

            // Solve U*x = scale*v for a right eigenvector
              // or U**T*x = scale*v for a left eigenvector,
            // overwriting x on v.

            slatrs('Upper', TRANS, 'Nonunit', NORMIN, N, B, LDB, VR, SCALE, WORK, IERR );
            NORMIN = 'Y'

            // Test for sufficient growth in the norm of v.

            VNORM = SASUM( N, VR, 1 )
            IF( VNORM.GE.GROWTO*SCALE ) GO TO 120

            // Choose new orthogonal starting vector and try again.

            TEMP = EPS3 / ( ROOTN+ONE )
            VR( 1 ) = EPS3
            DO 100 I = 2, N
               VR( I ) = TEMP
  100       CONTINUE
            VR( N-ITS+1 ) = VR( N-ITS+1 ) - EPS3*ROOTN
  110    CONTINUE

         // Failure to find eigenvector in N iterations.

         INFO = 1

  120    CONTINUE

         // Normalize eigenvector.

         I = ISAMAX( N, VR, 1 )
         sscal(N, ONE / ABS( VR( I ) ), VR, 1 );
      } else {

         // Complex eigenvalue.

         if ( NOINIT ) {

            // Set initial vector.

            DO 130 I = 1, N
               VR( I ) = EPS3
               VI( I ) = ZERO
  130       CONTINUE
         } else {

            // Scale supplied initial vector.

            NORM = SLAPY2( SNRM2( N, VR, 1 ), SNRM2( N, VI, 1 ) )
            REC = ( EPS3*ROOTN ) / MAX( NORM, NRMSML )
            sscal(N, REC, VR, 1 );
            sscal(N, REC, VI, 1 );
         }

         if ( RIGHTV ) {

            // LU decomposition with partial pivoting of B, replacing zero
            // pivots by EPS3.

            // The imaginary part of the (i,j)-th element of U is stored in
            // B(j+1,i).

            B( 2, 1 ) = -WI
            DO 140 I = 2, N
               B( I+1, 1 ) = ZERO
  140       CONTINUE

            DO 170 I = 1, N - 1
               ABSBII = SLAPY2( B( I, I ), B( I+1, I ) )
               EI = H( I+1, I )
               if ( ABSBII.LT.ABS( EI ) ) {

                  // Interchange rows and eliminate.

                  XR = B( I, I ) / EI
                  XI = B( I+1, I ) / EI
                  B( I, I ) = EI
                  B( I+1, I ) = ZERO
                  DO 150 J = I + 1, N
                     TEMP = B( I+1, J )
                     B( I+1, J ) = B( I, J ) - XR*TEMP
                     B( J+1, I+1 ) = B( J+1, I ) - XI*TEMP
                     B( I, J ) = TEMP
                     B( J+1, I ) = ZERO
  150             CONTINUE
                  B( I+2, I ) = -WI
                  B( I+1, I+1 ) = B( I+1, I+1 ) - XI*WI
                  B( I+2, I+1 ) = B( I+2, I+1 ) + XR*WI
               } else {

                  // Eliminate without interchanging rows.

                  if ( ABSBII.EQ.ZERO ) {
                     B( I, I ) = EPS3
                     B( I+1, I ) = ZERO
                     ABSBII = EPS3
                  }
                  EI = ( EI / ABSBII ) / ABSBII
                  XR = B( I, I )*EI
                  XI = -B( I+1, I )*EI
                  DO 160 J = I + 1, N
                     B( I+1, J ) = B( I+1, J ) - XR*B( I, J ) + XI*B( J+1, I )
                     B( J+1, I+1 ) = -XR*B( J+1, I ) - XI*B( I, J )
  160             CONTINUE
                  B( I+2, I+1 ) = B( I+2, I+1 ) - WI
               }

               // Compute 1-norm of offdiagonal elements of i-th row.

               WORK( I ) = SASUM( N-I, B( I, I+1 ), LDB ) + SASUM( N-I, B( I+2, I ), 1 )
  170       CONTINUE
            IF( B( N, N ).EQ.ZERO .AND. B( N+1, N ).EQ.ZERO ) B( N, N ) = EPS3
            WORK( N ) = ZERO

            I1 = N
            I2 = 1
            I3 = -1
         } else {

            // UL decomposition with partial pivoting of conjg(B),
            // replacing zero pivots by EPS3.

            // The imaginary part of the (i,j)-th element of U is stored in
            // B(j+1,i).

            B( N+1, N ) = WI
            DO 180 J = 1, N - 1
               B( N+1, J ) = ZERO
  180       CONTINUE

            DO 210 J = N, 2, -1
               EJ = H( J, J-1 )
               ABSBJJ = SLAPY2( B( J, J ), B( J+1, J ) )
               if ( ABSBJJ.LT.ABS( EJ ) ) {

                  // Interchange columns and eliminate

                  XR = B( J, J ) / EJ
                  XI = B( J+1, J ) / EJ
                  B( J, J ) = EJ
                  B( J+1, J ) = ZERO
                  DO 190 I = 1, J - 1
                     TEMP = B( I, J-1 )
                     B( I, J-1 ) = B( I, J ) - XR*TEMP
                     B( J, I ) = B( J+1, I ) - XI*TEMP
                     B( I, J ) = TEMP
                     B( J+1, I ) = ZERO
  190             CONTINUE
                  B( J+1, J-1 ) = WI
                  B( J-1, J-1 ) = B( J-1, J-1 ) + XI*WI
                  B( J, J-1 ) = B( J, J-1 ) - XR*WI
               } else {

                  // Eliminate without interchange.

                  if ( ABSBJJ.EQ.ZERO ) {
                     B( J, J ) = EPS3
                     B( J+1, J ) = ZERO
                     ABSBJJ = EPS3
                  }
                  EJ = ( EJ / ABSBJJ ) / ABSBJJ
                  XR = B( J, J )*EJ
                  XI = -B( J+1, J )*EJ
                  DO 200 I = 1, J - 1
                     B( I, J-1 ) = B( I, J-1 ) - XR*B( I, J ) + XI*B( J+1, I )
                     B( J, I ) = -XR*B( J+1, I ) - XI*B( I, J )
  200             CONTINUE
                  B( J, J-1 ) = B( J, J-1 ) + WI
               }

               // Compute 1-norm of offdiagonal elements of j-th column.

               WORK( J ) = SASUM( J-1, B( 1, J ), 1 ) + SASUM( J-1, B( J+1, 1 ), LDB )
  210       CONTINUE
            IF( B( 1, 1 ).EQ.ZERO .AND. B( 2, 1 ).EQ.ZERO ) B( 1, 1 ) = EPS3
            WORK( 1 ) = ZERO

            I1 = 1
            I2 = N
            I3 = 1
         }

         DO 270 ITS = 1, N
            SCALE = ONE
            VMAX = ONE
            VCRIT = BIGNUM

            // Solve U*(xr,xi) = scale*(vr,vi) for a right eigenvector,
              // or U**T*(xr,xi) = scale*(vr,vi) for a left eigenvector,
            // overwriting (xr,xi) on (vr,vi).

            DO 250 I = I1, I2, I3

               if ( WORK( I ).GT.VCRIT ) {
                  REC = ONE / VMAX
                  sscal(N, REC, VR, 1 );
                  sscal(N, REC, VI, 1 );
                  SCALE = SCALE*REC
                  VMAX = ONE
                  VCRIT = BIGNUM
               }

               XR = VR( I )
               XI = VI( I )
               if ( RIGHTV ) {
                  DO 220 J = I + 1, N
                     XR = XR - B( I, J )*VR( J ) + B( J+1, I )*VI( J )
                     XI = XI - B( I, J )*VI( J ) - B( J+1, I )*VR( J )
  220             CONTINUE
               } else {
                  DO 230 J = 1, I - 1
                     XR = XR - B( J, I )*VR( J ) + B( I+1, J )*VI( J )
                     XI = XI - B( J, I )*VI( J ) - B( I+1, J )*VR( J )
  230             CONTINUE
               }

               W = ABS( B( I, I ) ) + ABS( B( I+1, I ) )
               if ( W.GT.SMLNUM ) {
                  if ( W.LT.ONE ) {
                     W1 = ABS( XR ) + ABS( XI )
                     if ( W1.GT.W*BIGNUM ) {
                        REC = ONE / W1
                        sscal(N, REC, VR, 1 );
                        sscal(N, REC, VI, 1 );
                        XR = VR( I )
                        XI = VI( I )
                        SCALE = SCALE*REC
                        VMAX = VMAX*REC
                     }
                  }

                  // Divide by diagonal element of B.

                  sladiv(XR, XI, B( I, I ), B( I+1, I ), VR( I ), VI( I ) );
                  VMAX = MAX( ABS( VR( I ) )+ABS( VI( I ) ), VMAX )
                  VCRIT = BIGNUM / VMAX
               } else {
                  DO 240 J = 1, N
                     VR( J ) = ZERO
                     VI( J ) = ZERO
  240             CONTINUE
                  VR( I ) = ONE
                  VI( I ) = ONE
                  SCALE = ZERO
                  VMAX = ONE
                  VCRIT = BIGNUM
               }
  250       CONTINUE

            // Test for sufficient growth in the norm of (VR,VI).

            VNORM = SASUM( N, VR, 1 ) + SASUM( N, VI, 1 )
            IF( VNORM.GE.GROWTO*SCALE ) GO TO 280

            // Choose a new orthogonal starting vector and try again.

            Y = EPS3 / ( ROOTN+ONE )
            VR( 1 ) = EPS3
            VI( 1 ) = ZERO

            DO 260 I = 2, N
               VR( I ) = Y
               VI( I ) = ZERO
  260       CONTINUE
            VR( N-ITS+1 ) = VR( N-ITS+1 ) - EPS3*ROOTN
  270    CONTINUE

         // Failure to find eigenvector in N iterations

         INFO = 1

  280    CONTINUE

         // Normalize eigenvector.

         VNORM = ZERO
         DO 290 I = 1, N
            VNORM = MAX( VNORM, ABS( VR( I ) )+ABS( VI( I ) ) )
  290    CONTINUE
         sscal(N, ONE / VNORM, VR, 1 );
         sscal(N, ONE / VNORM, VI, 1 );

      }

      RETURN

      // End of SLAEIN

      }
