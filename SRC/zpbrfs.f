      SUBROUTINE ZPBRFS( UPLO, N, KD, NRHS, AB, LDAB, AFB, LDAFB, B, LDB, X, LDX, FERR, BERR, WORK, RWORK, INFO )

*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             UPLO;
      int                INFO, KD, LDAB, LDAFB, LDB, LDX, N, NRHS;
      // ..
      // .. Array Arguments ..
      double             BERR( * ), FERR( * ), RWORK( * );
      COMPLEX*16         AB( LDAB, * ), AFB( LDAFB, * ), B( LDB, * ), WORK( * ), X( LDX, * )
      // ..

*  =====================================================================

      // .. Parameters ..
      int                ITMAX;
      const              ITMAX = 5 ;
      double             ZERO;
      const              ZERO = 0.0D+0 ;
      COMPLEX*16         ONE
      const              ONE = ( 1.0D+0, 0.0D+0 ) ;
      double             TWO;
      const              TWO = 2.0D+0 ;
      double             THREE;
      const              THREE = 3.0D+0 ;
      // ..
      // .. Local Scalars ..
      bool               UPPER;
      int                COUNT, I, J, K, KASE, L, NZ;
      double             EPS, LSTRES, S, SAFE1, SAFE2, SAFMIN, XK;
      COMPLEX*16         ZDUM
      // ..
      // .. Local Arrays ..
      int                ISAVE( 3 );
      // ..
      // .. External Subroutines ..
      // EXTERNAL XERBLA, ZAXPY, ZCOPY, ZHBMV, ZLACN2, ZPBTRS
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, DBLE, DIMAG, MAX, MIN
      // ..
      // .. External Functions ..
      bool               LSAME;
      double             DLAMCH;
      // EXTERNAL LSAME, DLAMCH
      // ..
      // .. Statement Functions ..
      double             CABS1;
      // ..
      // .. Statement Function definitions ..
      CABS1( ZDUM ) = ABS( DBLE( ZDUM ) ) + ABS( DIMAG( ZDUM ) )
      // ..
      // .. Executable Statements ..

      // Test the input parameters.

      INFO = 0
      UPPER = LSAME( UPLO, 'U' )
      if ( .NOT.UPPER .AND. .NOT.LSAME( UPLO, 'L' ) ) {
         INFO = -1
      } else if ( N.LT.0 ) {
         INFO = -2
      } else if ( KD.LT.0 ) {
         INFO = -3
      } else if ( NRHS.LT.0 ) {
         INFO = -4
      } else if ( LDAB.LT.KD+1 ) {
         INFO = -6
      } else if ( LDAFB.LT.KD+1 ) {
         INFO = -8
      } else if ( LDB.LT.MAX( 1, N ) ) {
         INFO = -10
      } else if ( LDX.LT.MAX( 1, N ) ) {
         INFO = -12
      }
      if ( INFO.NE.0 ) {
         xerbla('ZPBRFS', -INFO );
         RETURN
      }

      // Quick return if possible

      if ( N.EQ.0 .OR. NRHS.EQ.0 ) {
         DO 10 J = 1, NRHS
            FERR( J ) = ZERO
            BERR( J ) = ZERO
   10    CONTINUE
         RETURN
      }

      // NZ = maximum number of nonzero elements in each row of A, plus 1

      NZ = MIN( N+1, 2*KD+2 )
      EPS = DLAMCH( 'Epsilon' )
      SAFMIN = DLAMCH( 'Safe minimum' )
      SAFE1 = NZ*SAFMIN
      SAFE2 = SAFE1 / EPS

      // Do for each right hand side

      DO 140 J = 1, NRHS

         COUNT = 1
         LSTRES = THREE
   20    CONTINUE

         // Loop until stopping criterion is satisfied.

         // Compute residual R = B - A * X

         zcopy(N, B( 1, J ), 1, WORK, 1 );
         zhbmv(UPLO, N, KD, -ONE, AB, LDAB, X( 1, J ), 1, ONE, WORK, 1 );

         // Compute componentwise relative backward error from formula

         // max(i) ( abs(R(i)) / ( abs(A)*abs(X) + abs(B) )(i) )

         // where abs(Z) is the componentwise absolute value of the matrix
         // or vector Z.  If the i-th component of the denominator is less
         // than SAFE2, then SAFE1 is added to the i-th components of the
         // numerator and denominator before dividing.

         DO 30 I = 1, N
            RWORK( I ) = CABS1( B( I, J ) )
   30    CONTINUE

         // Compute abs(A)*abs(X) + abs(B).

         if ( UPPER ) {
            DO 50 K = 1, N
               S = ZERO
               XK = CABS1( X( K, J ) )
               L = KD + 1 - K
               DO 40 I = MAX( 1, K-KD ), K - 1
                  RWORK( I ) = RWORK( I ) + CABS1( AB( L+I, K ) )*XK
                  S = S + CABS1( AB( L+I, K ) )*CABS1( X( I, J ) )
   40          CONTINUE
               RWORK( K ) = RWORK( K ) + ABS( DBLE( AB( KD+1, K ) ) )* XK + S
   50       CONTINUE
         } else {
            DO 70 K = 1, N
               S = ZERO
               XK = CABS1( X( K, J ) )
               RWORK( K ) = RWORK( K ) + ABS( DBLE( AB( 1, K ) ) )*XK
               L = 1 - K
               DO 60 I = K + 1, MIN( N, K+KD )
                  RWORK( I ) = RWORK( I ) + CABS1( AB( L+I, K ) )*XK
                  S = S + CABS1( AB( L+I, K ) )*CABS1( X( I, J ) )
   60          CONTINUE
               RWORK( K ) = RWORK( K ) + S
   70       CONTINUE
         }
         S = ZERO
         DO 80 I = 1, N
            if ( RWORK( I ).GT.SAFE2 ) {
               S = MAX( S, CABS1( WORK( I ) ) / RWORK( I ) )
            } else {
               S = MAX( S, ( CABS1( WORK( I ) )+SAFE1 ) / ( RWORK( I )+SAFE1 ) )
            }
   80    CONTINUE
         BERR( J ) = S

         // Test stopping criterion. Continue iterating if
            // 1) The residual BERR(J) is larger than machine epsilon, and
            // 2) BERR(J) decreased by at least a factor of 2 during the
               // last iteration, and
            // 3) At most ITMAX iterations tried.

         if ( BERR( J ).GT.EPS .AND. TWO*BERR( J ).LE.LSTRES .AND. COUNT.LE.ITMAX ) {

            // Update solution and try again.

            zpbtrs(UPLO, N, KD, 1, AFB, LDAFB, WORK, N, INFO );
            zaxpy(N, ONE, WORK, 1, X( 1, J ), 1 );
            LSTRES = BERR( J )
            COUNT = COUNT + 1
            GO TO 20
         }

         // Bound error from formula

         // norm(X - XTRUE) / norm(X) .le. FERR =
         // norm( abs(inv(A))*
            // ( abs(R) + NZ*EPS*( abs(A)*abs(X)+abs(B) ))) / norm(X)

         // where
           // norm(Z) is the magnitude of the largest component of Z
           // inv(A) is the inverse of A
           // abs(Z) is the componentwise absolute value of the matrix or
              // vector Z
           // NZ is the maximum number of nonzeros in any row of A, plus 1
           // EPS is machine epsilon

         // The i-th component of abs(R)+NZ*EPS*(abs(A)*abs(X)+abs(B))
         // is incremented by SAFE1 if the i-th component of
         // abs(A)*abs(X) + abs(B) is less than SAFE2.

         // Use ZLACN2 to estimate the infinity-norm of the matrix
            // inv(A) * diag(W),
         // where W = abs(R) + NZ*EPS*( abs(A)*abs(X)+abs(B) )))

         DO 90 I = 1, N
            if ( RWORK( I ).GT.SAFE2 ) {
               RWORK( I ) = CABS1( WORK( I ) ) + NZ*EPS*RWORK( I )
            } else {
               RWORK( I ) = CABS1( WORK( I ) ) + NZ*EPS*RWORK( I ) + SAFE1
            }
   90    CONTINUE

         KASE = 0
  100    CONTINUE
         zlacn2(N, WORK( N+1 ), WORK, FERR( J ), KASE, ISAVE );
         if ( KASE.NE.0 ) {
            if ( KASE.EQ.1 ) {

               // Multiply by diag(W)*inv(A**H).

               zpbtrs(UPLO, N, KD, 1, AFB, LDAFB, WORK, N, INFO );
               DO 110 I = 1, N
                  WORK( I ) = RWORK( I )*WORK( I )
  110          CONTINUE
            } else if ( KASE.EQ.2 ) {

               // Multiply by inv(A)*diag(W).

               DO 120 I = 1, N
                  WORK( I ) = RWORK( I )*WORK( I )
  120          CONTINUE
               zpbtrs(UPLO, N, KD, 1, AFB, LDAFB, WORK, N, INFO );
            }
            GO TO 100
         }

         // Normalize error.

         LSTRES = ZERO
         DO 130 I = 1, N
            LSTRES = MAX( LSTRES, CABS1( X( I, J ) ) )
  130    CONTINUE
         IF( LSTRES.NE.ZERO ) FERR( J ) = FERR( J ) / LSTRES

  140 CONTINUE

      RETURN

      // End of ZPBRFS

      }
