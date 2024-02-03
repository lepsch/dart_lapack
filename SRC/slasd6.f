      SUBROUTINE SLASD6( ICOMPQ, NL, NR, SQRE, D, VF, VL, ALPHA, BETA, IDXQ, PERM, GIVPTR, GIVCOL, LDGCOL, GIVNUM, LDGNUM, POLES, DIFL, DIFR, Z, K, C, S, WORK, IWORK, INFO )
*
*  -- LAPACK auxiliary routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      INTEGER            GIVPTR, ICOMPQ, INFO, K, LDGCOL, LDGNUM, NL, NR, SQRE
      REAL               ALPHA, BETA, C, S
*     ..
*     .. Array Arguments ..
      INTEGER            GIVCOL( LDGCOL, * ), IDXQ( * ), IWORK( * ), PERM( * )       REAL               D( * ), DIFL( * ), DIFR( * ), GIVNUM( LDGNUM, * ), POLES( LDGNUM, * ), VF( * ), VL( * ), WORK( * ), Z( * )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      REAL               ONE, ZERO
      PARAMETER          ( ONE = 1.0E+0, ZERO = 0.0E+0 )
*     ..
*     .. Local Scalars ..
      INTEGER            I, IDX, IDXC, IDXP, ISIGMA, IVFW, IVLW, IW, M, N, N1, N2
      REAL               ORGNRM
*     ..
*     .. External Subroutines ..
      EXTERNAL           SCOPY, SLAMRG, SLASCL, SLASD7, SLASD8, XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, MAX
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters.
*
      INFO = 0
      N = NL + NR + 1
      M = N + SQRE
*
      IF( ( ICOMPQ.LT.0 ) .OR. ( ICOMPQ.GT.1 ) ) THEN
         INFO = -1
      ELSE IF( NL.LT.1 ) THEN
         INFO = -2
      ELSE IF( NR.LT.1 ) THEN
         INFO = -3
      ELSE IF( ( SQRE.LT.0 ) .OR. ( SQRE.GT.1 ) ) THEN
         INFO = -4
      ELSE IF( LDGCOL.LT.N ) THEN
         INFO = -14
      ELSE IF( LDGNUM.LT.N ) THEN
         INFO = -16
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'SLASD6', -INFO )
         RETURN
      END IF
*
*     The following values are for bookkeeping purposes only.  They are
*     integer pointers which indicate the portion of the workspace
*     used by a particular array in SLASD7 and SLASD8.
*
      ISIGMA = 1
      IW = ISIGMA + N
      IVFW = IW + M
      IVLW = IVFW + M
*
      IDX = 1
      IDXC = IDX + N
      IDXP = IDXC + N
*
*     Scale.
*
      ORGNRM = MAX( ABS( ALPHA ), ABS( BETA ) )
      D( NL+1 ) = ZERO
      DO 10 I = 1, N
         IF( ABS( D( I ) ).GT.ORGNRM ) THEN
            ORGNRM = ABS( D( I ) )
         END IF
   10 CONTINUE
      CALL SLASCL( 'G', 0, 0, ORGNRM, ONE, N, 1, D, N, INFO )
      ALPHA = ALPHA / ORGNRM
      BETA = BETA / ORGNRM
*
*     Sort and Deflate singular values.
*
      CALL SLASD7( ICOMPQ, NL, NR, SQRE, K, D, Z, WORK( IW ), VF, WORK( IVFW ), VL, WORK( IVLW ), ALPHA, BETA, WORK( ISIGMA ), IWORK( IDX ), IWORK( IDXP ), IDXQ, PERM, GIVPTR, GIVCOL, LDGCOL, GIVNUM, LDGNUM, C, S, INFO )
*
*     Solve Secular Equation, compute DIFL, DIFR, and update VF, VL.
*
      CALL SLASD8( ICOMPQ, K, D, Z, VF, VL, DIFL, DIFR, LDGNUM, WORK( ISIGMA ), WORK( IW ), INFO )
*
*     Report the possible convergence failure.
*
      IF( INFO.NE.0 ) THEN
         RETURN
      END IF
*
*     Save the poles if ICOMPQ = 1.
*
      IF( ICOMPQ.EQ.1 ) THEN
         CALL SCOPY( K, D, 1, POLES( 1, 1 ), 1 )
         CALL SCOPY( K, WORK( ISIGMA ), 1, POLES( 1, 2 ), 1 )
      END IF
*
*     Unscale.
*
      CALL SLASCL( 'G', 0, 0, ONE, ORGNRM, N, 1, D, N, INFO )
*
*     Prepare the IDXQ sorting permutation.
*
      N1 = K
      N2 = N - K
      CALL SLAMRG( N1, N2, D, 1, -1, IDXQ )
*
      RETURN
*
*     End of SLASD6
*
      END
