      SUBROUTINE SLAED8( ICOMPQ, K, N, QSIZ, D, Q, LDQ, INDXQ, RHO,
     $                   CUTPNT, Z, DLAMBDA, Q2, LDQ2, W, PERM, GIVPTR,
     $                   GIVCOL, GIVNUM, INDXP, INDX, INFO )
*
*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      INTEGER            CUTPNT, GIVPTR, ICOMPQ, INFO, K, LDQ, LDQ2, N,
     $                   QSIZ
      REAL               RHO
*     ..
*     .. Array Arguments ..
      INTEGER            GIVCOL( 2, * ), INDX( * ), INDXP( * ),
     $                   INDXQ( * ), PERM( * )
      REAL               D( * ), DLAMBDA( * ), GIVNUM( 2, * ),
     $                   Q( LDQ, * ), Q2( LDQ2, * ), W( * ), Z( * )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      REAL               MONE, ZERO, ONE, TWO, EIGHT
      PARAMETER          ( MONE = -1.0E0, ZERO = 0.0E0, ONE = 1.0E0,
     $                   TWO = 2.0E0, EIGHT = 8.0E0 )
*     ..
*     .. Local Scalars ..
*
      INTEGER            I, IMAX, J, JLAM, JMAX, JP, K2, N1, N1P1, N2
      REAL               C, EPS, S, T, TAU, TOL
*     ..
*     .. External Functions ..
      INTEGER            ISAMAX
      REAL               SLAMCH, SLAPY2
      EXTERNAL           ISAMAX, SLAMCH, SLAPY2
*     ..
*     .. External Subroutines ..
      EXTERNAL           SCOPY, SLACPY, SLAMRG, SROT, SSCAL, XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, MAX, MIN, SQRT
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters.
*
      INFO = 0
*
      IF( ICOMPQ.LT.0 .OR. ICOMPQ.GT.1 ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -3
      ELSE IF( ICOMPQ.EQ.1 .AND. QSIZ.LT.N ) THEN
         INFO = -4
      ELSE IF( LDQ.LT.MAX( 1, N ) ) THEN
         INFO = -7
      ELSE IF( CUTPNT.LT.MIN( 1, N ) .OR. CUTPNT.GT.N ) THEN
         INFO = -10
      ELSE IF( LDQ2.LT.MAX( 1, N ) ) THEN
         INFO = -14
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'SLAED8', -INFO )
         RETURN
      END IF
*
*     Need to initialize GIVPTR to O here in case of quick exit
*     to prevent an unspecified code behavior (usually sigfault)
*     when IWORK array on entry to *stedc is not zeroed
*     (or at least some IWORK entries which used in *laed7 for GIVPTR).
*
      GIVPTR = 0
*
*     Quick return if possible
*
      IF( N.EQ.0 )
     $   RETURN
*
      N1 = CUTPNT
      N2 = N - N1
      N1P1 = N1 + 1
*
      IF( RHO.LT.ZERO ) THEN
         CALL SSCAL( N2, MONE, Z( N1P1 ), 1 )
      END IF
*
*     Normalize z so that norm(z) = 1
*
      T = ONE / SQRT( TWO )
      DO 10 J = 1, N
         INDX( J ) = J
   10 CONTINUE
      CALL SSCAL( N, T, Z, 1 )
      RHO = ABS( TWO*RHO )
*
*     Sort the eigenvalues into increasing order
*
      DO 20 I = CUTPNT + 1, N
         INDXQ( I ) = INDXQ( I ) + CUTPNT
   20 CONTINUE
      DO 30 I = 1, N
         DLAMBDA( I ) = D( INDXQ( I ) )
         W( I ) = Z( INDXQ( I ) )
   30 CONTINUE
      I = 1
      J = CUTPNT + 1
      CALL SLAMRG( N1, N2, DLAMBDA, 1, 1, INDX )
      DO 40 I = 1, N
         D( I ) = DLAMBDA( INDX( I ) )
         Z( I ) = W( INDX( I ) )
   40 CONTINUE
*
*     Calculate the allowable deflation tolerance
*
      IMAX = ISAMAX( N, Z, 1 )
      JMAX = ISAMAX( N, D, 1 )
      EPS = SLAMCH( 'Epsilon' )
      TOL = EIGHT*EPS*ABS( D( JMAX ) )
*
*     If the rank-1 modifier is small enough, no more needs to be done
*     except to reorganize Q so that its columns correspond with the
*     elements in D.
*
      IF( RHO*ABS( Z( IMAX ) ).LE.TOL ) THEN
         K = 0
         IF( ICOMPQ.EQ.0 ) THEN
            DO 50 J = 1, N
               PERM( J ) = INDXQ( INDX( J ) )
   50       CONTINUE
         ELSE
            DO 60 J = 1, N
               PERM( J ) = INDXQ( INDX( J ) )
               CALL SCOPY( QSIZ, Q( 1, PERM( J ) ), 1, Q2( 1, J ), 1 )
   60       CONTINUE
            CALL SLACPY( 'A', QSIZ, N, Q2( 1, 1 ), LDQ2, Q( 1, 1 ),
     $                   LDQ )
         END IF
         RETURN
      END IF
*
*     If there are multiple eigenvalues then the problem deflates.  Here
*     the number of equal eigenvalues are found.  As each equal
*     eigenvalue is found, an elementary reflector is computed to rotate
*     the corresponding eigensubspace so that the corresponding
*     components of Z are zero in this new basis.
*
      K = 0
      K2 = N + 1
      DO 70 J = 1, N
         IF( RHO*ABS( Z( J ) ).LE.TOL ) THEN
*
*           Deflate due to small z component.
*
            K2 = K2 - 1
            INDXP( K2 ) = J
            IF( J.EQ.N )
     $         GO TO 110
         ELSE
            JLAM = J
            GO TO 80
         END IF
   70 CONTINUE
   80 CONTINUE
      J = J + 1
      IF( J.GT.N )
     $   GO TO 100
      IF( RHO*ABS( Z( J ) ).LE.TOL ) THEN
*
*        Deflate due to small z component.
*
         K2 = K2 - 1
         INDXP( K2 ) = J
      ELSE
*
*        Check if eigenvalues are close enough to allow deflation.
*
         S = Z( JLAM )
         C = Z( J )
*
*        Find sqrt(a**2+b**2) without overflow or
*        destructive underflow.
*
         TAU = SLAPY2( C, S )
         T = D( J ) - D( JLAM )
         C = C / TAU
         S = -S / TAU
         IF( ABS( T*C*S ).LE.TOL ) THEN
*
*           Deflation is possible.
*
            Z( J ) = TAU
            Z( JLAM ) = ZERO
*
*           Record the appropriate Givens rotation
*
            GIVPTR = GIVPTR + 1
            GIVCOL( 1, GIVPTR ) = INDXQ( INDX( JLAM ) )
            GIVCOL( 2, GIVPTR ) = INDXQ( INDX( J ) )
            GIVNUM( 1, GIVPTR ) = C
            GIVNUM( 2, GIVPTR ) = S
            IF( ICOMPQ.EQ.1 ) THEN
               CALL SROT( QSIZ, Q( 1, INDXQ( INDX( JLAM ) ) ), 1,
     $                    Q( 1, INDXQ( INDX( J ) ) ), 1, C, S )
            END IF
            T = D( JLAM )*C*C + D( J )*S*S
            D( J ) = D( JLAM )*S*S + D( J )*C*C
            D( JLAM ) = T
            K2 = K2 - 1
            I = 1
   90       CONTINUE
            IF( K2+I.LE.N ) THEN
               IF( D( JLAM ).LT.D( INDXP( K2+I ) ) ) THEN
                  INDXP( K2+I-1 ) = INDXP( K2+I )
                  INDXP( K2+I ) = JLAM
                  I = I + 1
                  GO TO 90
               ELSE
                  INDXP( K2+I-1 ) = JLAM
               END IF
            ELSE
               INDXP( K2+I-1 ) = JLAM
            END IF
            JLAM = J
         ELSE
            K = K + 1
            W( K ) = Z( JLAM )
            DLAMBDA( K ) = D( JLAM )
            INDXP( K ) = JLAM
            JLAM = J
         END IF
      END IF
      GO TO 80
  100 CONTINUE
*
*     Record the last eigenvalue.
*
      K = K + 1
      W( K ) = Z( JLAM )
      DLAMBDA( K ) = D( JLAM )
      INDXP( K ) = JLAM
*
  110 CONTINUE
*
*     Sort the eigenvalues and corresponding eigenvectors into DLAMBDA
*     and Q2 respectively.  The eigenvalues/vectors which were not
*     deflated go into the first K slots of DLAMBDA and Q2 respectively,
*     while those which were deflated go into the last N - K slots.
*
      IF( ICOMPQ.EQ.0 ) THEN
         DO 120 J = 1, N
            JP = INDXP( J )
            DLAMBDA( J ) = D( JP )
            PERM( J ) = INDXQ( INDX( JP ) )
  120    CONTINUE
      ELSE
         DO 130 J = 1, N
            JP = INDXP( J )
            DLAMBDA( J ) = D( JP )
            PERM( J ) = INDXQ( INDX( JP ) )
            CALL SCOPY( QSIZ, Q( 1, PERM( J ) ), 1, Q2( 1, J ), 1 )
  130    CONTINUE
      END IF
*
*     The deflated eigenvalues and their corresponding vectors go back
*     into the last N - K slots of D and Q respectively.
*
      IF( K.LT.N ) THEN
         IF( ICOMPQ.EQ.0 ) THEN
            CALL SCOPY( N-K, DLAMBDA( K+1 ), 1, D( K+1 ), 1 )
         ELSE
            CALL SCOPY( N-K, DLAMBDA( K+1 ), 1, D( K+1 ), 1 )
            CALL SLACPY( 'A', QSIZ, N-K, Q2( 1, K+1 ), LDQ2,
     $                   Q( 1, K+1 ), LDQ )
         END IF
      END IF
*
      RETURN
*
*     End of SLAED8
*
      END