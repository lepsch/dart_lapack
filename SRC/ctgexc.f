      SUBROUTINE CTGEXC( WANTQ, WANTZ, N, A, LDA, B, LDB, Q, LDQ, Z, LDZ, IFST, ILST, INFO )
*
*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      bool               WANTQ, WANTZ;
      int                IFST, ILST, INFO, LDA, LDB, LDQ, LDZ, N;
*     ..
*     .. Array Arguments ..
      COMPLEX            A( LDA, * ), B( LDB, * ), Q( LDQ, * ), Z( LDZ, * )
*     ..
*
*  =====================================================================
*
*     .. Local Scalars ..
      int                HERE;
*     ..
*     .. External Subroutines ..
      EXTERNAL           CTGEX2, XERBLA
*     ..
*     .. Intrinsic Functions ..
      // INTRINSIC MAX
*     ..
*     .. Executable Statements ..
*
*     Decode and test input arguments.
      INFO = 0
      IF( N.LT.0 ) THEN
         INFO = -3
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -5
      ELSE IF( LDB.LT.MAX( 1, N ) ) THEN
         INFO = -7
      ELSE IF( LDQ.LT.1 .OR. WANTQ .AND. ( LDQ.LT.MAX( 1, N ) ) ) THEN
         INFO = -9
      ELSE IF( LDZ.LT.1 .OR. WANTZ .AND. ( LDZ.LT.MAX( 1, N ) ) ) THEN
         INFO = -11
      ELSE IF( IFST.LT.1 .OR. IFST.GT.N ) THEN
         INFO = -12
      ELSE IF( ILST.LT.1 .OR. ILST.GT.N ) THEN
         INFO = -13
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'CTGEXC', -INFO )
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( N.LE.1 ) RETURN       IF( IFST.EQ.ILST ) RETURN
*
      IF( IFST.LT.ILST ) THEN
*
         HERE = IFST
*
   10    CONTINUE
*
*        Swap with next one below
*
         CALL CTGEX2( WANTQ, WANTZ, N, A, LDA, B, LDB, Q, LDQ, Z, LDZ, HERE, INFO )
         IF( INFO.NE.0 ) THEN
            ILST = HERE
            RETURN
         END IF
         HERE = HERE + 1
         IF( HERE.LT.ILST ) GO TO 10
         HERE = HERE - 1
      ELSE
         HERE = IFST - 1
*
   20    CONTINUE
*
*        Swap with next one above
*
         CALL CTGEX2( WANTQ, WANTZ, N, A, LDA, B, LDB, Q, LDQ, Z, LDZ, HERE, INFO )
         IF( INFO.NE.0 ) THEN
            ILST = HERE
            RETURN
         END IF
         HERE = HERE - 1
         IF( HERE.GE.ILST ) GO TO 20
         HERE = HERE + 1
      END IF
      ILST = HERE
      RETURN
*
*     End of CTGEXC
*
      END
