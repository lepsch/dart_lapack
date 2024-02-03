      SUBROUTINE ZPTRFS( UPLO, N, NRHS, D, E, DF, EF, B, LDB, X, LDX, FERR, BERR, WORK, RWORK, INFO )

*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             UPLO;
      int                INFO, LDB, LDX, N, NRHS;
      // ..
      // .. Array Arguments ..
      double             BERR( * ), D( * ), DF( * ), FERR( * ), RWORK( * );
      COMPLEX*16         B( LDB, * ), E( * ), EF( * ), WORK( * ), X( LDX, * );
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
      int                COUNT, I, IX, J, NZ;
      double             EPS, LSTRES, S, SAFE1, SAFE2, SAFMIN;
      COMPLEX*16         BI, CX, DX, EX, ZDUM
      // ..
      // .. External Functions ..
      bool               LSAME;
      int                IDAMAX;
      double             DLAMCH;
      // EXTERNAL LSAME, IDAMAX, DLAMCH
      // ..
      // .. External Subroutines ..
      // EXTERNAL XERBLA, ZAXPY, ZPTTRS
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, DBLE, DCMPLX, DCONJG, DIMAG, MAX
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
      IF( .NOT.UPPER .AND. .NOT.LSAME( UPLO, 'L' ) ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( NRHS.LT.0 ) THEN
         INFO = -3
      ELSE IF( LDB.LT.MAX( 1, N ) ) THEN
         INFO = -9
      ELSE IF( LDX.LT.MAX( 1, N ) ) THEN
         INFO = -11
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'ZPTRFS', -INFO )
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
      EPS = DLAMCH( 'Epsilon' )
      SAFMIN = DLAMCH( 'Safe minimum' )
      SAFE1 = NZ*SAFMIN
      SAFE2 = SAFE1 / EPS

      // Do for each right hand side

      DO 100 J = 1, NRHS

         COUNT = 1
         LSTRES = THREE
   20    CONTINUE

         // Loop until stopping criterion is satisfied.

         // Compute residual R = B - A * X.  Also compute
         // abs(A)*abs(x) + abs(b) for use in the backward error bound.

         IF( UPPER ) THEN
            IF( N.EQ.1 ) THEN
               BI = B( 1, J )
               DX = D( 1 )*X( 1, J )
               WORK( 1 ) = BI - DX
               RWORK( 1 ) = CABS1( BI ) + CABS1( DX )
            ELSE
               BI = B( 1, J )
               DX = D( 1 )*X( 1, J )
               EX = E( 1 )*X( 2, J )
               WORK( 1 ) = BI - DX - EX
               RWORK( 1 ) = CABS1( BI ) + CABS1( DX ) + CABS1( E( 1 ) )*CABS1( X( 2, J ) )
               DO 30 I = 2, N - 1
                  BI = B( I, J )
                  CX = DCONJG( E( I-1 ) )*X( I-1, J )
                  DX = D( I )*X( I, J )
                  EX = E( I )*X( I+1, J )
                  WORK( I ) = BI - CX - DX - EX
                  RWORK( I ) = CABS1( BI ) + CABS1( E( I-1 ) )*CABS1( X( I-1, J ) ) + CABS1( DX ) + CABS1( E( I ) )* CABS1( X( I+1, J ) )
   30          CONTINUE
               BI = B( N, J )
               CX = DCONJG( E( N-1 ) )*X( N-1, J )
               DX = D( N )*X( N, J )
               WORK( N ) = BI - CX - DX
               RWORK( N ) = CABS1( BI ) + CABS1( E( N-1 ) )* CABS1( X( N-1, J ) ) + CABS1( DX )
            END IF
         ELSE
            IF( N.EQ.1 ) THEN
               BI = B( 1, J )
               DX = D( 1 )*X( 1, J )
               WORK( 1 ) = BI - DX
               RWORK( 1 ) = CABS1( BI ) + CABS1( DX )
            ELSE
               BI = B( 1, J )
               DX = D( 1 )*X( 1, J )
               EX = DCONJG( E( 1 ) )*X( 2, J )
               WORK( 1 ) = BI - DX - EX
               RWORK( 1 ) = CABS1( BI ) + CABS1( DX ) + CABS1( E( 1 ) )*CABS1( X( 2, J ) )
               DO 40 I = 2, N - 1
                  BI = B( I, J )
                  CX = E( I-1 )*X( I-1, J )
                  DX = D( I )*X( I, J )
                  EX = DCONJG( E( I ) )*X( I+1, J )
                  WORK( I ) = BI - CX - DX - EX
                  RWORK( I ) = CABS1( BI ) + CABS1( E( I-1 ) )*CABS1( X( I-1, J ) ) + CABS1( DX ) + CABS1( E( I ) )* CABS1( X( I+1, J ) )
   40          CONTINUE
               BI = B( N, J )
               CX = E( N-1 )*X( N-1, J )
               DX = D( N )*X( N, J )
               WORK( N ) = BI - CX - DX
               RWORK( N ) = CABS1( BI ) + CABS1( E( N-1 ) )* CABS1( X( N-1, J ) ) + CABS1( DX )
            END IF
         END IF

         // Compute componentwise relative backward error from formula

         // max(i) ( abs(R(i)) / ( abs(A)*abs(X) + abs(B) )(i) )

         // where abs(Z) is the componentwise absolute value of the matrix
         // or vector Z.  If the i-th component of the denominator is less
        t // han SAFE2, then SAFE1 is added to the i-th components of the
         // numerator and denominator before dividing.

         S = ZERO
         DO 50 I = 1, N
            IF( RWORK( I ).GT.SAFE2 ) THEN
               S = MAX( S, CABS1( WORK( I ) ) / RWORK( I ) )
            ELSE
               S = MAX( S, ( CABS1( WORK( I ) )+SAFE1 ) / ( RWORK( I )+SAFE1 ) )
            END IF
   50    CONTINUE
         BERR( J ) = S

         // Test stopping criterion. Continue iterating if
            // 1) The residual BERR(J) is larger than machine epsilon, and
            // 2) BERR(J) decreased by at least a factor of 2 during the
               // last iteration, and
            // 3) At most ITMAX iterations tried.

         IF( BERR( J ).GT.EPS .AND. TWO*BERR( J ).LE.LSTRES .AND. COUNT.LE.ITMAX ) THEN

            // Update solution and try again.

            CALL ZPTTRS( UPLO, N, 1, DF, EF, WORK, N, INFO )
            CALL ZAXPY( N, DCMPLX( ONE ), WORK, 1, X( 1, J ), 1 )
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

         DO 60 I = 1, N
            IF( RWORK( I ).GT.SAFE2 ) THEN
               RWORK( I ) = CABS1( WORK( I ) ) + NZ*EPS*RWORK( I )
            ELSE
               RWORK( I ) = CABS1( WORK( I ) ) + NZ*EPS*RWORK( I ) + SAFE1
            END IF
   60    CONTINUE
         IX = IDAMAX( N, RWORK, 1 )
         FERR( J ) = RWORK( IX )

         // Estimate the norm of inv(A).

         // Solve M(A) * x = e, where M(A) = (m(i,j)) is given by

            // m(i,j) =  abs(A(i,j)), i = j,
            // m(i,j) = -abs(A(i,j)), i .ne. j,

         // and e = [ 1, 1, ..., 1 ]**T.  Note M(A) = M(L)*D*M(L)**H.

         // Solve M(L) * x = e.

         RWORK( 1 ) = ONE
         DO 70 I = 2, N
            RWORK( I ) = ONE + RWORK( I-1 )*ABS( EF( I-1 ) )
   70    CONTINUE

         // Solve D * M(L)**H * x = b.

         RWORK( N ) = RWORK( N ) / DF( N )
         DO 80 I = N - 1, 1, -1
            RWORK( I ) = RWORK( I ) / DF( I ) + RWORK( I+1 )*ABS( EF( I ) )
   80    CONTINUE

         // Compute norm(inv(A)) = max(x(i)), 1<=i<=n.

         IX = IDAMAX( N, RWORK, 1 )
         FERR( J ) = FERR( J )*ABS( RWORK( IX ) )

         // Normalize error.

         LSTRES = ZERO
         DO 90 I = 1, N
            LSTRES = MAX( LSTRES, ABS( X( I, J ) ) )
   90    CONTINUE
         IF( LSTRES.NE.ZERO ) FERR( J ) = FERR( J ) / LSTRES

  100 CONTINUE

      RETURN

      // End of ZPTRFS

      }
