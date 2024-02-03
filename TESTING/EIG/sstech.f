      SUBROUTINE SSTECH( N, A, B, EIG, TOL, WORK, INFO )
*
*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      int                INFO, N;
      REAL               TOL
*     ..
*     .. Array Arguments ..
      REAL               A( * ), B( * ), EIG( * ), WORK( * )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      REAL               ZERO
      PARAMETER          ( ZERO = 0.0E0 )
*     ..
*     .. Local Scalars ..
      int                BPNT, COUNT, I, ISUB, J, NUML, NUMU, TPNT;
      REAL               EMIN, EPS, LOWER, MX, TUPPR, UNFLEP, UPPER
*     ..
*     .. External Functions ..
      REAL               SLAMCH
      // EXTERNAL SLAMCH
*     ..
*     .. External Subroutines ..
      // EXTERNAL SSTECT
*     ..
*     .. Intrinsic Functions ..
      // INTRINSIC ABS, MAX
*     ..
*     .. Executable Statements ..
*
*     Check input parameters
*
      INFO = 0
      IF( N.EQ.0 ) RETURN
      IF( N.LT.0 ) THEN
         INFO = -1
         RETURN
      END IF
      IF( TOL.LT.ZERO ) THEN
         INFO = -5
         RETURN
      END IF
*
*     Get machine constants
*
      EPS = SLAMCH( 'Epsilon' )*SLAMCH( 'Base' )
      UNFLEP = SLAMCH( 'Safe minimum' ) / EPS
      EPS = TOL*EPS
*
*     Compute maximum absolute eigenvalue, error tolerance
*
      MX = ABS( EIG( 1 ) )
      DO 10 I = 2, N
         MX = MAX( MX, ABS( EIG( I ) ) )
   10 CONTINUE
      EPS = MAX( EPS*MX, UNFLEP )
*
*     Sort eigenvalues from EIG into WORK
*
      DO 20 I = 1, N
         WORK( I ) = EIG( I )
   20 CONTINUE
      DO 40 I = 1, N - 1
         ISUB = 1
         EMIN = WORK( 1 )
         DO 30 J = 2, N + 1 - I
            IF( WORK( J ).LT.EMIN ) THEN
               ISUB = J
               EMIN = WORK( J )
            END IF
   30    CONTINUE
         IF( ISUB.NE.N+1-I ) THEN
            WORK( ISUB ) = WORK( N+1-I )
            WORK( N+1-I ) = EMIN
         END IF
   40 CONTINUE
*
*     TPNT points to singular value at right endpoint of interval
*     BPNT points to singular value at left  endpoint of interval
*
      TPNT = 1
      BPNT = 1
*
*     Begin loop over all intervals
*
   50 CONTINUE
      UPPER = WORK( TPNT ) + EPS
      LOWER = WORK( BPNT ) - EPS
*
*     Begin loop merging overlapping intervals
*
   60 CONTINUE
      IF( BPNT.EQ.N ) GO TO 70
      TUPPR = WORK( BPNT+1 ) + EPS
      IF( TUPPR.LT.LOWER ) GO TO 70
*
*     Merge
*
      BPNT = BPNT + 1
      LOWER = WORK( BPNT ) - EPS
      GO TO 60
   70 CONTINUE
*
*     Count singular values in interval [ LOWER, UPPER ]
*
      CALL SSTECT( N, A, B, LOWER, NUML )
      CALL SSTECT( N, A, B, UPPER, NUMU )
      COUNT = NUMU - NUML
      IF( COUNT.NE.BPNT-TPNT+1 ) THEN
*
*        Wrong number of singular values in interval
*
         INFO = TPNT
         GO TO 80
      END IF
      TPNT = BPNT + 1
      BPNT = TPNT
      IF( TPNT.LE.N ) GO TO 50
   80 CONTINUE
      RETURN
*
*     End of SSTECH
*
      END
