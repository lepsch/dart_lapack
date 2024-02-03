      SUBROUTINE CTBTRS( UPLO, TRANS, DIAG, N, KD, NRHS, AB, LDAB, B, LDB, INFO )

*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             DIAG, TRANS, UPLO;
      int                INFO, KD, LDAB, LDB, N, NRHS;
      // ..
      // .. Array Arguments ..
      COMPLEX            AB( LDAB, * ), B( LDB, * )
      // ..

*  =====================================================================

      // .. Parameters ..
      COMPLEX            ZERO
      const              ZERO = ( 0.0E+0, 0.0E+0 ) ;
      // ..
      // .. Local Scalars ..
      bool               NOUNIT, UPPER;
      int                J;
      // ..
      // .. External Functions ..
      bool               LSAME;
      // EXTERNAL LSAME
      // ..
      // .. External Subroutines ..
      // EXTERNAL CTBSV, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX
      // ..
      // .. Executable Statements ..

      // Test the input parameters.

      INFO = 0
      NOUNIT = LSAME( DIAG, 'N' )
      UPPER = LSAME( UPLO, 'U' )
      IF( .NOT.UPPER .AND. .NOT.LSAME( UPLO, 'L' ) ) THEN
         INFO = -1
      ELSE IF( .NOT.LSAME( TRANS, 'N' ) .AND. .NOT. LSAME( TRANS, 'T' ) .AND. .NOT.LSAME( TRANS, 'C' ) ) THEN
         INFO = -2
      ELSE IF( .NOT.NOUNIT .AND. .NOT.LSAME( DIAG, 'U' ) ) THEN
         INFO = -3
      ELSE IF( N.LT.0 ) THEN
         INFO = -4
      ELSE IF( KD.LT.0 ) THEN
         INFO = -5
      ELSE IF( NRHS.LT.0 ) THEN
         INFO = -6
      ELSE IF( LDAB.LT.KD+1 ) THEN
         INFO = -8
      ELSE IF( LDB.LT.MAX( 1, N ) ) THEN
         INFO = -10
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'CTBTRS', -INFO )
         RETURN
      END IF

      // Quick return if possible

      IF( N.EQ.0 ) RETURN

      // Check for singularity.

      IF( NOUNIT ) THEN
         IF( UPPER ) THEN
            DO 10 INFO = 1, N
               IF( AB( KD+1, INFO ).EQ.ZERO ) RETURN
   10       CONTINUE
         ELSE
            DO 20 INFO = 1, N
               IF( AB( 1, INFO ).EQ.ZERO ) RETURN
   20       CONTINUE
         END IF
      END IF
      INFO = 0

      // Solve A * X = B,  A**T * X = B,  or  A**H * X = B.

      DO 30 J = 1, NRHS
         CALL CTBSV( UPLO, TRANS, DIAG, N, KD, AB, LDAB, B( 1, J ), 1 )
   30 CONTINUE

      RETURN

      // End of CTBTRS

      }
