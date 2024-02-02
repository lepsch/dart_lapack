      SUBROUTINE CTGEX2( WANTQ, WANTZ, N, A, LDA, B, LDB, Q, LDQ, Z,
     $                   LDZ, J1, INFO )
*
*  -- LAPACK auxiliary routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      LOGICAL            WANTQ, WANTZ
      INTEGER            INFO, J1, LDA, LDB, LDQ, LDZ, N
*     ..
*     .. Array Arguments ..
      COMPLEX            A( LDA, * ), B( LDB, * ), Q( LDQ, * ),
     $                   Z( LDZ, * )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      COMPLEX            CZERO, CONE
      PARAMETER          ( CZERO = ( 0.0E+0, 0.0E+0 ),
     $                   CONE = ( 1.0E+0, 0.0E+0 ) )
      REAL               TWENTY
      PARAMETER          ( TWENTY = 2.0E+1 )
      INTEGER            LDST
      PARAMETER          ( LDST = 2 )
      LOGICAL            WANDS
      PARAMETER          ( WANDS = .TRUE. )
*     ..
*     .. Local Scalars ..
      LOGICAL            STRONG, WEAK
      INTEGER            I, M
      REAL               CQ, CZ, EPS, SA, SB, SCALE, SMLNUM, SUM,
     $                   THRESHA, THRESHB
      COMPLEX            CDUM, F, G, SQ, SZ
*     ..
*     .. Local Arrays ..
      COMPLEX            S( LDST, LDST ), T( LDST, LDST ), WORK( 8 )
*     ..
*     .. External Functions ..
      REAL               SLAMCH
      EXTERNAL           SLAMCH
*     ..
*     .. External Subroutines ..
      EXTERNAL           CLACPY, CLARTG, CLASSQ, CROT
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, CONJG, MAX, REAL, SQRT
*     ..
*     .. Executable Statements ..
*
      INFO = 0
*
*     Quick return if possible
*
      IF( N.LE.1 )
     $   RETURN
*
      M = LDST
      WEAK = .FALSE.
      STRONG = .FALSE.
*
*     Make a local copy of selected block in (A, B)
*
      CALL CLACPY( 'Full', M, M, A( J1, J1 ), LDA, S, LDST )
      CALL CLACPY( 'Full', M, M, B( J1, J1 ), LDB, T, LDST )
*
*     Compute the threshold for testing the acceptance of swapping.
*
      EPS = SLAMCH( 'P' )
      SMLNUM = SLAMCH( 'S' ) / EPS
      SCALE = REAL( CZERO )
      SUM = REAL( CONE )
      CALL CLACPY( 'Full', M, M, S, LDST, WORK, M )
      CALL CLACPY( 'Full', M, M, T, LDST, WORK( M*M+1 ), M )
      CALL CLASSQ( M*M, WORK, 1, SCALE, SUM )
      SA = SCALE*SQRT( SUM )
      SCALE = DBLE( CZERO )
      SUM = DBLE( CONE )
      CALL CLASSQ( M*M, WORK(M*M+1), 1, SCALE, SUM )
      SB = SCALE*SQRT( SUM )
*
*     THRES has been changed from
*        THRESH = MAX( TEN*EPS*SA, SMLNUM )
*     to
*        THRESH = MAX( TWENTY*EPS*SA, SMLNUM )
*     on 04/01/10.
*     "Bug" reported by Ondra Kamenik, confirmed by Julie Langou, fixed by
*     Jim Demmel and Guillaume Revy. See forum post 1783.
*
      THRESHA = MAX( TWENTY*EPS*SA, SMLNUM )
      THRESHB = MAX( TWENTY*EPS*SB, SMLNUM )
*
*     Compute unitary QL and RQ that swap 1-by-1 and 1-by-1 blocks
*     using Givens rotations and perform the swap tentatively.
*
      F = S( 2, 2 )*T( 1, 1 ) - T( 2, 2 )*S( 1, 1 )
      G = S( 2, 2 )*T( 1, 2 ) - T( 2, 2 )*S( 1, 2 )
      SA = ABS( S( 2, 2 ) ) * ABS( T( 1, 1 ) )
      SB = ABS( S( 1, 1 ) ) * ABS( T( 2, 2 ) )
      CALL CLARTG( G, F, CZ, SZ, CDUM )
      SZ = -SZ
      CALL CROT( 2, S( 1, 1 ), 1, S( 1, 2 ), 1, CZ, CONJG( SZ ) )
      CALL CROT( 2, T( 1, 1 ), 1, T( 1, 2 ), 1, CZ, CONJG( SZ ) )
      IF( SA.GE.SB ) THEN
         CALL CLARTG( S( 1, 1 ), S( 2, 1 ), CQ, SQ, CDUM )
      ELSE
         CALL CLARTG( T( 1, 1 ), T( 2, 1 ), CQ, SQ, CDUM )
      END IF
      CALL CROT( 2, S( 1, 1 ), LDST, S( 2, 1 ), LDST, CQ, SQ )
      CALL CROT( 2, T( 1, 1 ), LDST, T( 2, 1 ), LDST, CQ, SQ )
*
*     Weak stability test: |S21| <= O(EPS F-norm((A)))
*                          and  |T21| <= O(EPS F-norm((B)))
*
      WEAK = ABS( S( 2, 1 ) ).LE.THRESHA .AND. 
     $ ABS( T( 2, 1 ) ).LE.THRESHB
      IF( .NOT.WEAK )
     $   GO TO 20
*
      IF( WANDS ) THEN
*
*        Strong stability test:
*           F-norm((A-QL**H*S*QR, B-QL**H*T*QR)) <= O(EPS*F-norm((A, B)))
*
         CALL CLACPY( 'Full', M, M, S, LDST, WORK, M )
         CALL CLACPY( 'Full', M, M, T, LDST, WORK( M*M+1 ), M )
         CALL CROT( 2, WORK, 1, WORK( 3 ), 1, CZ, -CONJG( SZ ) )
         CALL CROT( 2, WORK( 5 ), 1, WORK( 7 ), 1, CZ, -CONJG( SZ ) )
         CALL CROT( 2, WORK, 2, WORK( 2 ), 2, CQ, -SQ )
         CALL CROT( 2, WORK( 5 ), 2, WORK( 6 ), 2, CQ, -SQ )
         DO 10 I = 1, 2
            WORK( I ) = WORK( I ) - A( J1+I-1, J1 )
            WORK( I+2 ) = WORK( I+2 ) - A( J1+I-1, J1+1 )
            WORK( I+4 ) = WORK( I+4 ) - B( J1+I-1, J1 )
            WORK( I+6 ) = WORK( I+6 ) - B( J1+I-1, J1+1 )
   10    CONTINUE
         SCALE = DBLE( CZERO )
         SUM = DBLE( CONE )
         CALL CLASSQ( M*M, WORK, 1, SCALE, SUM )
         SA = SCALE*SQRT( SUM )
         SCALE = DBLE( CZERO )
         SUM = DBLE( CONE )
         CALL CLASSQ( M*M, WORK(M*M+1), 1, SCALE, SUM )
         SB = SCALE*SQRT( SUM )
         STRONG = SA.LE.THRESHA .AND. SB.LE.THRESHB
         IF( .NOT.STRONG )
     $      GO TO 20
      END IF
*
*     If the swap is accepted ("weakly" and "strongly"), apply the
*     equivalence transformations to the original matrix pair (A,B)
*
      CALL CROT( J1+1, A( 1, J1 ), 1, A( 1, J1+1 ), 1, CZ, CONJG( SZ ) )
      CALL CROT( J1+1, B( 1, J1 ), 1, B( 1, J1+1 ), 1, CZ, CONJG( SZ ) )
      CALL CROT( N-J1+1, A( J1, J1 ), LDA, A( J1+1, J1 ), LDA, CQ, SQ )
      CALL CROT( N-J1+1, B( J1, J1 ), LDB, B( J1+1, J1 ), LDB, CQ, SQ )
*
*     Set  N1 by N2 (2,1) blocks to 0
*
      A( J1+1, J1 ) = CZERO
      B( J1+1, J1 ) = CZERO
*
*     Accumulate transformations into Q and Z if requested.
*
      IF( WANTZ )
     $   CALL CROT( N, Z( 1, J1 ), 1, Z( 1, J1+1 ), 1, CZ, CONJG( SZ ) )
      IF( WANTQ )
     $   CALL CROT( N, Q( 1, J1 ), 1, Q( 1, J1+1 ), 1, CQ, CONJG( SQ ) )
*
*     Exit with INFO = 0 if swap was successfully performed.
*
      RETURN
*
*     Exit with INFO = 1 if swap was rejected.
*
   20 CONTINUE
      INFO = 1
      RETURN
*
*     End of CTGEX2
*
      END