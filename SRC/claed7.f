      SUBROUTINE CLAED7( N, CUTPNT, QSIZ, TLVLS, CURLVL, CURPBM, D, Q, LDQ, RHO, INDXQ, QSTORE, QPTR, PRMPTR, PERM, GIVPTR, GIVCOL, GIVNUM, WORK, RWORK, IWORK, INFO )
*
*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      int                CURLVL, CURPBM, CUTPNT, INFO, LDQ, N, QSIZ, TLVLS;
      REAL               RHO
*     ..
*     .. Array Arguments ..
      int                GIVCOL( 2, * ), GIVPTR( * ), INDXQ( * ), IWORK( * ), PERM( * ), PRMPTR( * ), QPTR( * );
      REAL               D( * ), GIVNUM( 2, * ), QSTORE( * ), RWORK( * )
      COMPLEX            Q( LDQ, * ), WORK( * )
*     ..
*
*  =====================================================================
*
*     .. Local Scalars ..
      int                COLTYP, CURR, I, IDLMDA, INDX, INDXC, INDXP, IQ, IW, IZ, K, N1, N2, PTR;
*     ..
*     .. External Subroutines ..
      EXTERNAL           CLACRM, CLAED8, SLAED9, SLAEDA, SLAMRG, XERBLA
*     ..
*     .. Intrinsic Functions ..
      // INTRINSIC MAX, MIN
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters.
*
      INFO = 0
*
*     IF( ICOMPQ.LT.0 .OR. ICOMPQ.GT.1 ) THEN
*        INFO = -1
*     ELSE IF( N.LT.0 ) THEN
      IF( N.LT.0 ) THEN
         INFO = -1
      ELSE IF( MIN( 1, N ).GT.CUTPNT .OR. N.LT.CUTPNT ) THEN
         INFO = -2
      ELSE IF( QSIZ.LT.N ) THEN
         INFO = -3
      ELSE IF( LDQ.LT.MAX( 1, N ) ) THEN
         INFO = -9
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'CLAED7', -INFO )
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( N.EQ.0 ) RETURN
*
*     The following values are for bookkeeping purposes only.  They are
*     integer pointers which indicate the portion of the workspace
*     used by a particular array in SLAED2 and SLAED3.
*
      IZ = 1
      IDLMDA = IZ + N
      IW = IDLMDA + N
      IQ = IW + N
*
      INDX = 1
      INDXC = INDX + N
      COLTYP = INDXC + N
      INDXP = COLTYP + N
*
*     Form the z-vector which consists of the last row of Q_1 and the
*     first row of Q_2.
*
      PTR = 1 + 2**TLVLS
      DO 10 I = 1, CURLVL - 1
         PTR = PTR + 2**( TLVLS-I )
   10 CONTINUE
      CURR = PTR + CURPBM
      CALL SLAEDA( N, TLVLS, CURLVL, CURPBM, PRMPTR, PERM, GIVPTR, GIVCOL, GIVNUM, QSTORE, QPTR, RWORK( IZ ), RWORK( IZ+N ), INFO )
*
*     When solving the final problem, we no longer need the stored data,
*     so we will overwrite the data from this level onto the previously
*     used storage space.
*
      IF( CURLVL.EQ.TLVLS ) THEN
         QPTR( CURR ) = 1
         PRMPTR( CURR ) = 1
         GIVPTR( CURR ) = 1
      END IF
*
*     Sort and Deflate eigenvalues.
*
      CALL CLAED8( K, N, QSIZ, Q, LDQ, D, RHO, CUTPNT, RWORK( IZ ), RWORK( IDLMDA ), WORK, QSIZ, RWORK( IW ), IWORK( INDXP ), IWORK( INDX ), INDXQ, PERM( PRMPTR( CURR ) ), GIVPTR( CURR+1 ), GIVCOL( 1, GIVPTR( CURR ) ), GIVNUM( 1, GIVPTR( CURR ) ), INFO )
      PRMPTR( CURR+1 ) = PRMPTR( CURR ) + N
      GIVPTR( CURR+1 ) = GIVPTR( CURR+1 ) + GIVPTR( CURR )
*
*     Solve Secular Equation.
*
      IF( K.NE.0 ) THEN
         CALL SLAED9( K, 1, K, N, D, RWORK( IQ ), K, RHO, RWORK( IDLMDA ), RWORK( IW ), QSTORE( QPTR( CURR ) ), K, INFO )
         CALL CLACRM( QSIZ, K, WORK, QSIZ, QSTORE( QPTR( CURR ) ), K, Q, LDQ, RWORK( IQ ) )
         QPTR( CURR+1 ) = QPTR( CURR ) + K**2
         IF( INFO.NE.0 ) THEN
            RETURN
         END IF
*
*     Prepare the INDXQ sorting permutation.
*
         N1 = K
         N2 = N - K
         CALL SLAMRG( N1, N2, D, 1, -1, INDXQ )
      ELSE
         QPTR( CURR+1 ) = QPTR( CURR )
         DO 20 I = 1, N
            INDXQ( I ) = I
   20    CONTINUE
      END IF
*
      RETURN
*
*     End of CLAED7
*
      END
