      SUBROUTINE CLAED8( K, N, QSIZ, Q, LDQ, D, RHO, CUTPNT, Z, DLAMBDA, Q2, LDQ2, W, INDXP, INDX, INDXQ, PERM, GIVPTR, GIVCOL, GIVNUM, INFO )

*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                CUTPNT, GIVPTR, INFO, K, LDQ, LDQ2, N, QSIZ;
      REAL               RHO
      // ..
      // .. Array Arguments ..
      int                GIVCOL( 2, * ), INDX( * ), INDXP( * ), INDXQ( * ), PERM( * );
      REAL               D( * ), DLAMBDA( * ), GIVNUM( 2, * ), W( * ), Z( * );
      COMPLEX            Q( LDQ, * ), Q2( LDQ2, * )
      // ..

*  =====================================================================

      // .. Parameters ..
      REAL               MONE, ZERO, ONE, TWO, EIGHT
      const              MONE = -1.0E0, ZERO = 0.0E0, ONE = 1.0E0, TWO = 2.0E0, EIGHT = 8.0E0 ;
      // ..
      // .. Local Scalars ..
      int                I, IMAX, J, JLAM, JMAX, JP, K2, N1, N1P1, N2;
      REAL               C, EPS, S, T, TAU, TOL
      // ..
      // .. External Functions ..
      int                ISAMAX;
      REAL               SLAMCH, SLAPY2
      // EXTERNAL ISAMAX, SLAMCH, SLAPY2
      // ..
      // .. External Subroutines ..
      // EXTERNAL CCOPY, CLACPY, CSROT, SCOPY, SLAMRG, SSCAL, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, MAX, MIN, SQRT
      // ..
      // .. Executable Statements ..

      // Test the input parameters.

      INFO = 0

      if ( N.LT.0 ) {
         INFO = -2
      } else if ( QSIZ.LT.N ) {
         INFO = -3
      } else if ( LDQ.LT.MAX( 1, N ) ) {
         INFO = -5
      } else if ( CUTPNT.LT.MIN( 1, N ) .OR. CUTPNT.GT.N ) {
         INFO = -8
      } else if ( LDQ2.LT.MAX( 1, N ) ) {
         INFO = -12
      }
      if ( INFO.NE.0 ) {
         xerbla('CLAED8', -INFO );
         RETURN
      }

      // Need to initialize GIVPTR to O here in case of quick exit
      // to prevent an unspecified code behavior (usually sigfault)
      // when IWORK array on entry to *stedc is not zeroed
      // (or at least some IWORK entries which used in *laed7 for GIVPTR).

      GIVPTR = 0

      // Quick return if possible

      IF( N.EQ.0 ) RETURN

      N1 = CUTPNT
      N2 = N - N1
      N1P1 = N1 + 1

      if ( RHO.LT.ZERO ) {
         sscal(N2, MONE, Z( N1P1 ), 1 );
      }

      // Normalize z so that norm(z) = 1

      T = ONE / SQRT( TWO )
      for (J = 1; J <= N; J++) { // 10
         INDX( J ) = J
      } // 10
      sscal(N, T, Z, 1 );
      RHO = ABS( TWO*RHO )

      // Sort the eigenvalues into increasing order

      DO 20 I = CUTPNT + 1, N
         INDXQ( I ) = INDXQ( I ) + CUTPNT
      } // 20
      for (I = 1; I <= N; I++) { // 30
         DLAMBDA( I ) = D( INDXQ( I ) )
         W( I ) = Z( INDXQ( I ) )
      } // 30
      I = 1
      J = CUTPNT + 1
      slamrg(N1, N2, DLAMBDA, 1, 1, INDX );
      for (I = 1; I <= N; I++) { // 40
         D( I ) = DLAMBDA( INDX( I ) )
         Z( I ) = W( INDX( I ) )
      } // 40

      // Calculate the allowable deflation tolerance

      IMAX = ISAMAX( N, Z, 1 )
      JMAX = ISAMAX( N, D, 1 )
      EPS = SLAMCH( 'Epsilon' )
      TOL = EIGHT*EPS*ABS( D( JMAX ) )

      // If the rank-1 modifier is small enough, no more needs to be done
      // -- except to reorganize Q so that its columns correspond with the
      // elements in D.

      if ( RHO*ABS( Z( IMAX ) ).LE.TOL ) {
         K = 0
         for (J = 1; J <= N; J++) { // 50
            PERM( J ) = INDXQ( INDX( J ) )
            ccopy(QSIZ, Q( 1, PERM( J ) ), 1, Q2( 1, J ), 1 );
         } // 50
         clacpy('A', QSIZ, N, Q2( 1, 1 ), LDQ2, Q( 1, 1 ), LDQ );
         RETURN
      }

      // If there are multiple eigenvalues then the problem deflates.  Here
      // the number of equal eigenvalues are found.  As each equal
      // eigenvalue is found, an elementary reflector is computed to rotate
      // the corresponding eigensubspace so that the corresponding
      // components of Z are zero in this new basis.

      K = 0
      K2 = N + 1
      for (J = 1; J <= N; J++) { // 60
         if ( RHO*ABS( Z( J ) ).LE.TOL ) {

            // Deflate due to small z component.

            K2 = K2 - 1
            INDXP( K2 ) = J
            IF( J.EQ.N ) GO TO 100
         } else {
            JLAM = J
            GO TO 70
         }
      } // 60
      } // 70
      J = J + 1
      IF( J.GT.N ) GO TO 90
      if ( RHO*ABS( Z( J ) ).LE.TOL ) {

         // Deflate due to small z component.

         K2 = K2 - 1
         INDXP( K2 ) = J
      } else {

         // Check if eigenvalues are close enough to allow deflation.

         S = Z( JLAM )
         C = Z( J )

         // Find sqrt(a**2+b**2) without overflow or
         // destructive underflow.

         TAU = SLAPY2( C, S )
         T = D( J ) - D( JLAM )
         C = C / TAU
         S = -S / TAU
         if ( ABS( T*C*S ).LE.TOL ) {

            // Deflation is possible.

            Z( J ) = TAU
            Z( JLAM ) = ZERO

            // Record the appropriate Givens rotation

            GIVPTR = GIVPTR + 1
            GIVCOL( 1, GIVPTR ) = INDXQ( INDX( JLAM ) )
            GIVCOL( 2, GIVPTR ) = INDXQ( INDX( J ) )
            GIVNUM( 1, GIVPTR ) = C
            GIVNUM( 2, GIVPTR ) = S
            csrot(QSIZ, Q( 1, INDXQ( INDX( JLAM ) ) ), 1, Q( 1, INDXQ( INDX( J ) ) ), 1, C, S );
            T = D( JLAM )*C*C + D( J )*S*S
            D( J ) = D( JLAM )*S*S + D( J )*C*C
            D( JLAM ) = T
            K2 = K2 - 1
            I = 1
            } // 80
            if ( K2+I.LE.N ) {
               if ( D( JLAM ).LT.D( INDXP( K2+I ) ) ) {
                  INDXP( K2+I-1 ) = INDXP( K2+I )
                  INDXP( K2+I ) = JLAM
                  I = I + 1
                  GO TO 80
               } else {
                  INDXP( K2+I-1 ) = JLAM
               }
            } else {
               INDXP( K2+I-1 ) = JLAM
            }
            JLAM = J
         } else {
            K = K + 1
            W( K ) = Z( JLAM )
            DLAMBDA( K ) = D( JLAM )
            INDXP( K ) = JLAM
            JLAM = J
         }
      }
      GO TO 70
      } // 90

      // Record the last eigenvalue.

      K = K + 1
      W( K ) = Z( JLAM )
      DLAMBDA( K ) = D( JLAM )
      INDXP( K ) = JLAM

      } // 100

      // Sort the eigenvalues and corresponding eigenvectors into DLAMBDA
      // and Q2 respectively.  The eigenvalues/vectors which were not
      // deflated go into the first K slots of DLAMBDA and Q2 respectively,
      // while those which were deflated go into the last N - K slots.

      for (J = 1; J <= N; J++) { // 110
         JP = INDXP( J )
         DLAMBDA( J ) = D( JP )
         PERM( J ) = INDXQ( INDX( JP ) )
         ccopy(QSIZ, Q( 1, PERM( J ) ), 1, Q2( 1, J ), 1 );
      } // 110

      // The deflated eigenvalues and their corresponding vectors go back
      // into the last N - K slots of D and Q respectively.

      if ( K.LT.N ) {
         scopy(N-K, DLAMBDA( K+1 ), 1, D( K+1 ), 1 );
         clacpy('A', QSIZ, N-K, Q2( 1, K+1 ), LDQ2, Q( 1, K+1 ), LDQ );
      }

      RETURN

      // End of CLAED8

      }
