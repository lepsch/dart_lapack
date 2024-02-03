      SUBROUTINE DSPRFS( UPLO, N, NRHS, AP, AFP, IPIV, B, LDB, X, LDX, FERR, BERR, WORK, IWORK, INFO )

*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             UPLO;
      int                INFO, LDB, LDX, N, NRHS;
      // ..
      // .. Array Arguments ..
      int                IPIV( * ), IWORK( * );
      double             AFP( * ), AP( * ), B( LDB, * ), BERR( * ), FERR( * ), WORK( * ), X( LDX, * );
      // ..

*  =====================================================================

      // .. Parameters ..
      int                ITMAX;
      const              ITMAX = 5 ;
      double             ZERO;
      const              ZERO = 0.0D+0 ;
      double             ONE;
      const              ONE = 1.0D+0 ;
      double             TWO;
      const              TWO = 2.0D+0 ;
      double             THREE;
      const              THREE = 3.0D+0 ;
      // ..
      // .. Local Scalars ..
      bool               UPPER;
      int                COUNT, I, IK, J, K, KASE, KK, NZ;
      double             EPS, LSTRES, S, SAFE1, SAFE2, SAFMIN, XK;
      // ..
      // .. Local Arrays ..
      int                ISAVE( 3 );
      // ..
      // .. External Subroutines ..
      // EXTERNAL DAXPY, DCOPY, DLACN2, DSPMV, DSPTRS, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, MAX
      // ..
      // .. External Functions ..
      bool               LSAME;
      double             DLAMCH;
      // EXTERNAL LSAME, DLAMCH
      // ..
      // .. Executable Statements ..

      // Test the input parameters.

      INFO = 0
      UPPER = LSAME( UPLO, 'U' )
      if ( .NOT.UPPER .AND. .NOT.LSAME( UPLO, 'L' ) ) {
         INFO = -1
      } else if ( N.LT.0 ) {
         INFO = -2
      } else if ( NRHS.LT.0 ) {
         INFO = -3
      } else if ( LDB.LT.MAX( 1, N ) ) {
         INFO = -8
      } else if ( LDX.LT.MAX( 1, N ) ) {
         INFO = -10
      }
      if ( INFO.NE.0 ) {
         CALL XERBLA( 'DSPRFS', -INFO )
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

      NZ = N + 1
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

         CALL DCOPY( N, B( 1, J ), 1, WORK( N+1 ), 1 )
         CALL DSPMV( UPLO, N, -ONE, AP, X( 1, J ), 1, ONE, WORK( N+1 ), 1 )

         // Compute componentwise relative backward error from formula

         // max(i) ( abs(R(i)) / ( abs(A)*abs(X) + abs(B) )(i) )

         // where abs(Z) is the componentwise absolute value of the matrix
         // or vector Z.  If the i-th component of the denominator is less
        t // han SAFE2, then SAFE1 is added to the i-th components of the
         // numerator and denominator before dividing.

         DO 30 I = 1, N
            WORK( I ) = ABS( B( I, J ) )
   30    CONTINUE

         // Compute abs(A)*abs(X) + abs(B).

         KK = 1
         if ( UPPER ) {
            DO 50 K = 1, N
               S = ZERO
               XK = ABS( X( K, J ) )
               IK = KK
               DO 40 I = 1, K - 1
                  WORK( I ) = WORK( I ) + ABS( AP( IK ) )*XK
                  S = S + ABS( AP( IK ) )*ABS( X( I, J ) )
                  IK = IK + 1
   40          CONTINUE
               WORK( K ) = WORK( K ) + ABS( AP( KK+K-1 ) )*XK + S
               KK = KK + K
   50       CONTINUE
         } else {
            DO 70 K = 1, N
               S = ZERO
               XK = ABS( X( K, J ) )
               WORK( K ) = WORK( K ) + ABS( AP( KK ) )*XK
               IK = KK + 1
               DO 60 I = K + 1, N
                  WORK( I ) = WORK( I ) + ABS( AP( IK ) )*XK
                  S = S + ABS( AP( IK ) )*ABS( X( I, J ) )
                  IK = IK + 1
   60          CONTINUE
               WORK( K ) = WORK( K ) + S
               KK = KK + ( N-K+1 )
   70       CONTINUE
         }
         S = ZERO
         DO 80 I = 1, N
            if ( WORK( I ).GT.SAFE2 ) {
               S = MAX( S, ABS( WORK( N+I ) ) / WORK( I ) )
            } else {
               S = MAX( S, ( ABS( WORK( N+I ) )+SAFE1 ) / ( WORK( I )+SAFE1 ) )
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

            CALL DSPTRS( UPLO, N, 1, AFP, IPIV, WORK( N+1 ), N, INFO )
            CALL DAXPY( N, ONE, WORK( N+1 ), 1, X( 1, J ), 1 )
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

         // Use DLACN2 to estimate the infinity-norm of the matrix
            // inv(A) * diag(W),
         // where W = abs(R) + NZ*EPS*( abs(A)*abs(X)+abs(B) )))

         DO 90 I = 1, N
            if ( WORK( I ).GT.SAFE2 ) {
               WORK( I ) = ABS( WORK( N+I ) ) + NZ*EPS*WORK( I )
            } else {
               WORK( I ) = ABS( WORK( N+I ) ) + NZ*EPS*WORK( I ) + SAFE1
            }
   90    CONTINUE

         KASE = 0
  100    CONTINUE
         CALL DLACN2( N, WORK( 2*N+1 ), WORK( N+1 ), IWORK, FERR( J ), KASE, ISAVE )
         if ( KASE.NE.0 ) {
            if ( KASE.EQ.1 ) {

               // Multiply by diag(W)*inv(A**T).

               CALL DSPTRS( UPLO, N, 1, AFP, IPIV, WORK( N+1 ), N, INFO )
               DO 110 I = 1, N
                  WORK( N+I ) = WORK( I )*WORK( N+I )
  110          CONTINUE
            } else if ( KASE.EQ.2 ) {

               // Multiply by inv(A)*diag(W).

               DO 120 I = 1, N
                  WORK( N+I ) = WORK( I )*WORK( N+I )
  120          CONTINUE
               CALL DSPTRS( UPLO, N, 1, AFP, IPIV, WORK( N+1 ), N, INFO )
            }
            GO TO 100
         }

         // Normalize error.

         LSTRES = ZERO
         DO 130 I = 1, N
            LSTRES = MAX( LSTRES, ABS( X( I, J ) ) )
  130    CONTINUE
         IF( LSTRES.NE.ZERO ) FERR( J ) = FERR( J ) / LSTRES

  140 CONTINUE

      RETURN

      // End of DSPRFS

      }
