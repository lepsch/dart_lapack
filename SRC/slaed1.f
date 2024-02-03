      SUBROUTINE SLAED1( N, D, Q, LDQ, INDXQ, RHO, CUTPNT, WORK, IWORK, INFO )
*
*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
      // .. Scalar Arguments ..
      int                CUTPNT, INFO, LDQ, N;
      REAL               RHO
      // ..
      // .. Array Arguments ..
      int                INDXQ( * ), IWORK( * );
      REAL               D( * ), Q( LDQ, * ), WORK( * )
      // ..
*
*  =====================================================================
*
      // .. Local Scalars ..
      int                COLTYP, CPP1, I, IDLMDA, INDX, INDXC, INDXP, IQ2, IS, IW, IZ, K, N1, N2;
      // ..
      // .. External Subroutines ..
      // EXTERNAL SCOPY, SLAED2, SLAED3, SLAMRG, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX, MIN
      // ..
      // .. Executable Statements ..
*
      // Test the input parameters.
*
      INFO = 0
*
      IF( N.LT.0 ) THEN
         INFO = -1
      ELSE IF( LDQ.LT.MAX( 1, N ) ) THEN
         INFO = -4
      ELSE IF( MIN( 1, N / 2 ).GT.CUTPNT .OR. ( N / 2 ).LT.CUTPNT ) THEN
         INFO = -7
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'SLAED1', -INFO )
         RETURN
      END IF
*
      // Quick return if possible
*
      IF( N.EQ.0 ) RETURN
*
      // The following values are integer pointers which indicate
     t // he portion of the workspace
      // used by a particular array in SLAED2 and SLAED3.
*
      IZ = 1
      IDLMDA = IZ + N
      IW = IDLMDA + N
      IQ2 = IW + N
*
      INDX = 1
      INDXC = INDX + N
      COLTYP = INDXC + N
      INDXP = COLTYP + N
*
*
      // Form the z-vector which consists of the last row of Q_1 and the
      // first row of Q_2.
*
      CALL SCOPY( CUTPNT, Q( CUTPNT, 1 ), LDQ, WORK( IZ ), 1 )
      CPP1 = CUTPNT + 1
      CALL SCOPY( N-CUTPNT, Q( CPP1, CPP1 ), LDQ, WORK( IZ+CUTPNT ), 1 )
*
      // Deflate eigenvalues.
*
      CALL SLAED2( K, N, CUTPNT, D, Q, LDQ, INDXQ, RHO, WORK( IZ ), WORK( IDLMDA ), WORK( IW ), WORK( IQ2 ), IWORK( INDX ), IWORK( INDXC ), IWORK( INDXP ), IWORK( COLTYP ), INFO )
*
      IF( INFO.NE.0 ) GO TO 20
*
      // Solve Secular Equation.
*
      IF( K.NE.0 ) THEN
         IS = ( IWORK( COLTYP )+IWORK( COLTYP+1 ) )*CUTPNT + ( IWORK( COLTYP+1 )+IWORK( COLTYP+2 ) )*( N-CUTPNT ) + IQ2          CALL SLAED3( K, N, CUTPNT, D, Q, LDQ, RHO, WORK( IDLMDA ), WORK( IQ2 ), IWORK( INDXC ), IWORK( COLTYP ), WORK( IW ), WORK( IS ), INFO )
         IF( INFO.NE.0 ) GO TO 20
*
      // Prepare the INDXQ sorting permutation.
*
         N1 = K
         N2 = N - K
         CALL SLAMRG( N1, N2, D, 1, -1, INDXQ )
      ELSE
         DO 10 I = 1, N
            INDXQ( I ) = I
   10    CONTINUE
      END IF
*
   20 CONTINUE
      RETURN
*
      // End of SLAED1
*
      END
