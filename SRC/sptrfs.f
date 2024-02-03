      SUBROUTINE SPTRFS( N, NRHS, D, E, DF, EF, B, LDB, X, LDX, FERR, BERR, WORK, INFO )

*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                INFO, LDB, LDX, N, NRHS;
      // ..
      // .. Array Arguments ..
      REAL               B( LDB, * ), BERR( * ), D( * ), DF( * ), E( * ), EF( * ), FERR( * ), WORK( * ), X( LDX, * )
      // ..

*  =====================================================================

      // .. Parameters ..
      int                ITMAX;
      const              ITMAX = 5 ;
      REAL               ZERO
      const              ZERO = 0.0E+0 ;
      REAL               ONE
      const              ONE = 1.0E+0 ;
      REAL               TWO
      const              TWO = 2.0E+0 ;
      REAL               THREE
      const              THREE = 3.0E+0 ;
      // ..
      // .. Local Scalars ..
      int                COUNT, I, IX, J, NZ;
      REAL               BI, CX, DX, EPS, EX, LSTRES, S, SAFE1, SAFE2, SAFMIN
      // ..
      // .. External Subroutines ..
      // EXTERNAL SAXPY, SPTTRS, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, MAX
      // ..
      // .. External Functions ..
      int                ISAMAX;
      REAL               SLAMCH
      // EXTERNAL ISAMAX, SLAMCH
      // ..
      // .. Executable Statements ..

      // Test the input parameters.

      INFO = 0
      IF( N.LT.0 ) THEN
         INFO = -1
      ELSE IF( NRHS.LT.0 ) THEN
         INFO = -2
      ELSE IF( LDB.LT.MAX( 1, N ) ) THEN
         INFO = -8
      ELSE IF( LDX.LT.MAX( 1, N ) ) THEN
         INFO = -10
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'SPTRFS', -INFO )
         RETURN
      END IF

      // Quick return if possible

      IF( N.EQ.0 .OR. NRHS.EQ.0 ) THEN
         DO 10 J = 1, NRHS
            FERR( J ) = ZERO
            BERR( J ) = ZERO
   10    CONTINUE
         RETURN
      END IF

      // NZ = maximum number of nonzero elements in each row of A, plus 1

      NZ = 4
      EPS = SLAMCH( 'Epsilon' )
      SAFMIN = SLAMCH( 'Safe minimum' )
      SAFE1 = NZ*SAFMIN
      SAFE2 = SAFE1 / EPS

      // Do for each right hand side

      DO 90 J = 1, NRHS

         COUNT = 1
         LSTRES = THREE
   20    CONTINUE

         // Loop until stopping criterion is satisfied.

         // Compute residual R = B - A * X.  Also compute
         // abs(A)*abs(x) + abs(b) for use in the backward error bound.

         IF( N.EQ.1 ) THEN
            BI = B( 1, J )
            DX = D( 1 )*X( 1, J )
            WORK( N+1 ) = BI - DX
            WORK( 1 ) = ABS( BI ) + ABS( DX )
         } else {
            BI = B( 1, J )
            DX = D( 1 )*X( 1, J )
            EX = E( 1 )*X( 2, J )
            WORK( N+1 ) = BI - DX - EX
            WORK( 1 ) = ABS( BI ) + ABS( DX ) + ABS( EX )
            DO 30 I = 2, N - 1
               BI = B( I, J )
               CX = E( I-1 )*X( I-1, J )
               DX = D( I )*X( I, J )
               EX = E( I )*X( I+1, J )
               WORK( N+I ) = BI - CX - DX - EX
               WORK( I ) = ABS( BI ) + ABS( CX ) + ABS( DX ) + ABS( EX )
   30       CONTINUE
            BI = B( N, J )
            CX = E( N-1 )*X( N-1, J )
            DX = D( N )*X( N, J )
            WORK( N+N ) = BI - CX - DX
            WORK( N ) = ABS( BI ) + ABS( CX ) + ABS( DX )
         END IF

         // Compute componentwise relative backward error from formula

         // max(i) ( abs(R(i)) / ( abs(A)*abs(X) + abs(B) )(i) )

         // where abs(Z) is the componentwise absolute value of the matrix
         // or vector Z.  If the i-th component of the denominator is less
        t // han SAFE2, then SAFE1 is added to the i-th components of the
         // numerator and denominator before dividing.

         S = ZERO
         DO 40 I = 1, N
            IF( WORK( I ).GT.SAFE2 ) THEN
               S = MAX( S, ABS( WORK( N+I ) ) / WORK( I ) )
            } else {
               S = MAX( S, ( ABS( WORK( N+I ) )+SAFE1 ) / ( WORK( I )+SAFE1 ) )
            END IF
   40    CONTINUE
         BERR( J ) = S

         // Test stopping criterion. Continue iterating if
            // 1) The residual BERR(J) is larger than machine epsilon, and
            // 2) BERR(J) decreased by at least a factor of 2 during the
               // last iteration, and
            // 3) At most ITMAX iterations tried.

         IF( BERR( J ).GT.EPS .AND. TWO*BERR( J ).LE.LSTRES .AND. COUNT.LE.ITMAX ) THEN

            // Update solution and try again.

            CALL SPTTRS( N, 1, DF, EF, WORK( N+1 ), N, INFO )
            CALL SAXPY( N, ONE, WORK( N+1 ), 1, X( 1, J ), 1 )
            LSTRES = BERR( J )
            COUNT = COUNT + 1
            GO TO 20
         END IF

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

         DO 50 I = 1, N
            IF( WORK( I ).GT.SAFE2 ) THEN
               WORK( I ) = ABS( WORK( N+I ) ) + NZ*EPS*WORK( I )
            } else {
               WORK( I ) = ABS( WORK( N+I ) ) + NZ*EPS*WORK( I ) + SAFE1
            END IF
   50    CONTINUE
         IX = ISAMAX( N, WORK, 1 )
         FERR( J ) = WORK( IX )

         // Estimate the norm of inv(A).

         // Solve M(A) * x = e, where M(A) = (m(i,j)) is given by

            // m(i,j) =  abs(A(i,j)), i = j,
            // m(i,j) = -abs(A(i,j)), i .ne. j,

         // and e = [ 1, 1, ..., 1 ]**T.  Note M(A) = M(L)*D*M(L)**T.

         // Solve M(L) * x = e.

         WORK( 1 ) = ONE
         DO 60 I = 2, N
            WORK( I ) = ONE + WORK( I-1 )*ABS( EF( I-1 ) )
   60    CONTINUE

         // Solve D * M(L)**T * x = b.

         WORK( N ) = WORK( N ) / DF( N )
         DO 70 I = N - 1, 1, -1
            WORK( I ) = WORK( I ) / DF( I ) + WORK( I+1 )*ABS( EF( I ) )
   70    CONTINUE

         // Compute norm(inv(A)) = max(x(i)), 1<=i<=n.

         IX = ISAMAX( N, WORK, 1 )
         FERR( J ) = FERR( J )*ABS( WORK( IX ) )

         // Normalize error.

         LSTRES = ZERO
         DO 80 I = 1, N
            LSTRES = MAX( LSTRES, ABS( X( I, J ) ) )
   80    CONTINUE
         IF( LSTRES.NE.ZERO ) FERR( J ) = FERR( J ) / LSTRES

   90 CONTINUE

      RETURN

      // End of SPTRFS

      }
