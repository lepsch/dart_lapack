      SUBROUTINE ZGTRFS( TRANS, N, NRHS, DL, D, DU, DLF, DF, DUF, DU2, IPIV, B, LDB, X, LDX, FERR, BERR, WORK, RWORK, INFO )
*
*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      String             TRANS;
      int                INFO, LDB, LDX, N, NRHS;
*     ..
*     .. Array Arguments ..
      int                IPIV( * );
      double             BERR( * ), FERR( * ), RWORK( * );
      COMPLEX*16         B( LDB, * ), D( * ), DF( * ), DL( * ), DLF( * ), DU( * ), DU2( * ), DUF( * ), WORK( * ), X( LDX, * )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      int                ITMAX;
      PARAMETER          ( ITMAX = 5 )
      double             ZERO, ONE;
      PARAMETER          ( ZERO = 0.0D+0, ONE = 1.0D+0 )
      double             TWO;
      PARAMETER          ( TWO = 2.0D+0 )
      double             THREE;
      PARAMETER          ( THREE = 3.0D+0 )
*     ..
*     .. Local Scalars ..
      bool               NOTRAN;
      String             TRANSN, TRANST;
      int                COUNT, I, J, KASE, NZ;
      double             EPS, LSTRES, S, SAFE1, SAFE2, SAFMIN;
      COMPLEX*16         ZDUM
*     ..
*     .. Local Arrays ..
      int                ISAVE( 3 );
*     ..
*     .. External Subroutines ..
      // EXTERNAL XERBLA, ZAXPY, ZCOPY, ZGTTRS, ZLACN2, ZLAGTM
*     ..
*     .. Intrinsic Functions ..
      // INTRINSIC ABS, DBLE, DCMPLX, DIMAG, MAX
*     ..
*     .. External Functions ..
      bool               LSAME;
      double             DLAMCH;
      // EXTERNAL LSAME, DLAMCH
*     ..
*     .. Statement Functions ..
      double             CABS1;
*     ..
*     .. Statement Function definitions ..
      CABS1( ZDUM ) = ABS( DBLE( ZDUM ) ) + ABS( DIMAG( ZDUM ) )
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters.
*
      INFO = 0
      NOTRAN = LSAME( TRANS, 'N' )
      IF( .NOT.NOTRAN .AND. .NOT.LSAME( TRANS, 'T' ) .AND. .NOT. LSAME( TRANS, 'C' ) ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( NRHS.LT.0 ) THEN
         INFO = -3
      ELSE IF( LDB.LT.MAX( 1, N ) ) THEN
         INFO = -13
      ELSE IF( LDX.LT.MAX( 1, N ) ) THEN
         INFO = -15
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'ZGTRFS', -INFO )
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( N.EQ.0 .OR. NRHS.EQ.0 ) THEN
         DO 10 J = 1, NRHS
            FERR( J ) = ZERO
            BERR( J ) = ZERO
   10    CONTINUE
         RETURN
      END IF
*
      IF( NOTRAN ) THEN
         TRANSN = 'N'
         TRANST = 'C'
      ELSE
         TRANSN = 'C'
         TRANST = 'N'
      END IF
*
*     NZ = maximum number of nonzero elements in each row of A, plus 1
*
      NZ = 4
      EPS = DLAMCH( 'Epsilon' )
      SAFMIN = DLAMCH( 'Safe minimum' )
      SAFE1 = NZ*SAFMIN
      SAFE2 = SAFE1 / EPS
*
*     Do for each right hand side
*
      DO 110 J = 1, NRHS
*
         COUNT = 1
         LSTRES = THREE
   20    CONTINUE
*
*        Loop until stopping criterion is satisfied.
*
*        Compute residual R = B - op(A) * X,
*        where op(A) = A, A**T, or A**H, depending on TRANS.
*
         CALL ZCOPY( N, B( 1, J ), 1, WORK, 1 )
         CALL ZLAGTM( TRANS, N, 1, -ONE, DL, D, DU, X( 1, J ), LDX, ONE, WORK, N )
*
*        Compute abs(op(A))*abs(x) + abs(b) for use in the backward
*        error bound.
*
         IF( NOTRAN ) THEN
            IF( N.EQ.1 ) THEN
               RWORK( 1 ) = CABS1( B( 1, J ) ) + CABS1( D( 1 ) )*CABS1( X( 1, J ) )
            ELSE
               RWORK( 1 ) = CABS1( B( 1, J ) ) + CABS1( D( 1 ) )*CABS1( X( 1, J ) ) + CABS1( DU( 1 ) )*CABS1( X( 2, J ) )
               DO 30 I = 2, N - 1
                  RWORK( I ) = CABS1( B( I, J ) ) + CABS1( DL( I-1 ) )*CABS1( X( I-1, J ) ) + CABS1( D( I ) )*CABS1( X( I, J ) ) + CABS1( DU( I ) )*CABS1( X( I+1, J ) )
   30          CONTINUE
               RWORK( N ) = CABS1( B( N, J ) ) + CABS1( DL( N-1 ) )*CABS1( X( N-1, J ) ) + CABS1( D( N ) )*CABS1( X( N, J ) )
            END IF
         ELSE
            IF( N.EQ.1 ) THEN
               RWORK( 1 ) = CABS1( B( 1, J ) ) + CABS1( D( 1 ) )*CABS1( X( 1, J ) )
            ELSE
               RWORK( 1 ) = CABS1( B( 1, J ) ) + CABS1( D( 1 ) )*CABS1( X( 1, J ) ) + CABS1( DL( 1 ) )*CABS1( X( 2, J ) )
               DO 40 I = 2, N - 1
                  RWORK( I ) = CABS1( B( I, J ) ) + CABS1( DU( I-1 ) )*CABS1( X( I-1, J ) ) + CABS1( D( I ) )*CABS1( X( I, J ) ) + CABS1( DL( I ) )*CABS1( X( I+1, J ) )
   40          CONTINUE
               RWORK( N ) = CABS1( B( N, J ) ) + CABS1( DU( N-1 ) )*CABS1( X( N-1, J ) ) + CABS1( D( N ) )*CABS1( X( N, J ) )
            END IF
         END IF
*
*        Compute componentwise relative backward error from formula
*
*        max(i) ( abs(R(i)) / ( abs(op(A))*abs(X) + abs(B) )(i) )
*
*        where abs(Z) is the componentwise absolute value of the matrix
*        or vector Z.  If the i-th component of the denominator is less
*        than SAFE2, then SAFE1 is added to the i-th components of the
*        numerator and denominator before dividing.
*
         S = ZERO
         DO 50 I = 1, N
            IF( RWORK( I ).GT.SAFE2 ) THEN
               S = MAX( S, CABS1( WORK( I ) ) / RWORK( I ) )
            ELSE
               S = MAX( S, ( CABS1( WORK( I ) )+SAFE1 ) / ( RWORK( I )+SAFE1 ) )
            END IF
   50    CONTINUE
         BERR( J ) = S
*
*        Test stopping criterion. Continue iterating if
*           1) The residual BERR(J) is larger than machine epsilon, and
*           2) BERR(J) decreased by at least a factor of 2 during the
*              last iteration, and
*           3) At most ITMAX iterations tried.
*
         IF( BERR( J ).GT.EPS .AND. TWO*BERR( J ).LE.LSTRES .AND. COUNT.LE.ITMAX ) THEN
*
*           Update solution and try again.
*
            CALL ZGTTRS( TRANS, N, 1, DLF, DF, DUF, DU2, IPIV, WORK, N, INFO )
            CALL ZAXPY( N, DCMPLX( ONE ), WORK, 1, X( 1, J ), 1 )
            LSTRES = BERR( J )
            COUNT = COUNT + 1
            GO TO 20
         END IF
*
*        Bound error from formula
*
*        norm(X - XTRUE) / norm(X) .le. FERR =
*        norm( abs(inv(op(A)))*
*           ( abs(R) + NZ*EPS*( abs(op(A))*abs(X)+abs(B) ))) / norm(X)
*
*        where
*          norm(Z) is the magnitude of the largest component of Z
*          inv(op(A)) is the inverse of op(A)
*          abs(Z) is the componentwise absolute value of the matrix or
*             vector Z
*          NZ is the maximum number of nonzeros in any row of A, plus 1
*          EPS is machine epsilon
*
*        The i-th component of abs(R)+NZ*EPS*(abs(op(A))*abs(X)+abs(B))
*        is incremented by SAFE1 if the i-th component of
*        abs(op(A))*abs(X) + abs(B) is less than SAFE2.
*
*        Use ZLACN2 to estimate the infinity-norm of the matrix
*           inv(op(A)) * diag(W),
*        where W = abs(R) + NZ*EPS*( abs(op(A))*abs(X)+abs(B) )))
*
         DO 60 I = 1, N
            IF( RWORK( I ).GT.SAFE2 ) THEN
               RWORK( I ) = CABS1( WORK( I ) ) + NZ*EPS*RWORK( I )
            ELSE
               RWORK( I ) = CABS1( WORK( I ) ) + NZ*EPS*RWORK( I ) + SAFE1
            END IF
   60    CONTINUE
*
         KASE = 0
   70    CONTINUE
         CALL ZLACN2( N, WORK( N+1 ), WORK, FERR( J ), KASE, ISAVE )
         IF( KASE.NE.0 ) THEN
            IF( KASE.EQ.1 ) THEN
*
*              Multiply by diag(W)*inv(op(A)**H).
*
               CALL ZGTTRS( TRANST, N, 1, DLF, DF, DUF, DU2, IPIV, WORK, N, INFO )
               DO 80 I = 1, N
                  WORK( I ) = RWORK( I )*WORK( I )
   80          CONTINUE
            ELSE
*
*              Multiply by inv(op(A))*diag(W).
*
               DO 90 I = 1, N
                  WORK( I ) = RWORK( I )*WORK( I )
   90          CONTINUE
               CALL ZGTTRS( TRANSN, N, 1, DLF, DF, DUF, DU2, IPIV, WORK, N, INFO )
            END IF
            GO TO 70
         END IF
*
*        Normalize error.
*
         LSTRES = ZERO
         DO 100 I = 1, N
            LSTRES = MAX( LSTRES, CABS1( X( I, J ) ) )
  100    CONTINUE
         IF( LSTRES.NE.ZERO ) FERR( J ) = FERR( J ) / LSTRES
*
  110 CONTINUE
*
      RETURN
*
*     End of ZGTRFS
*
      END
