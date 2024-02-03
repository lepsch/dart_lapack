      REAL FUNCTION CLA_SYRCOND_X( UPLO, N, A, LDA, AF, LDAF, IPIV, X, INFO, WORK, RWORK )

*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             UPLO;
      int                N, LDA, LDAF, INFO;
      // ..
      // .. Array Arguments ..
      int                IPIV( * );
      COMPLEX            A( LDA, * ), AF( LDAF, * ), WORK( * ), X( * )
      REAL               RWORK( * )
      // ..

*  =====================================================================

      // .. Local Scalars ..
      int                KASE;
      REAL               AINVNM, ANORM, TMP
      int                I, J;
      bool               UP, UPPER;
      COMPLEX            ZDUM
      // ..
      // .. Local Arrays ..
      int                ISAVE( 3 );
      // ..
      // .. External Functions ..
      bool               LSAME;
      // EXTERNAL LSAME
      // ..
      // .. External Subroutines ..
      // EXTERNAL CLACN2, CSYTRS, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, MAX
      // ..
      // .. Statement Functions ..
      REAL               CABS1
      // ..
      // .. Statement Function Definitions ..
      CABS1( ZDUM ) = ABS( REAL( ZDUM ) ) + ABS( AIMAG( ZDUM ) )
      // ..
      // .. Executable Statements ..

      CLA_SYRCOND_X = 0.0E+0

      INFO = 0
      UPPER = LSAME( UPLO, 'U' )
      if ( .NOT.UPPER .AND. .NOT.LSAME( UPLO, 'L' ) ) {
         INFO = -1
      } else if ( N.LT.0 ) {
         INFO = -2
      } else if ( LDA.LT.MAX( 1, N ) ) {
         INFO = -4
      } else if ( LDAF.LT.MAX( 1, N ) ) {
         INFO = -6
      }
      if ( INFO.NE.0 ) {
         xerbla('CLA_SYRCOND_X', -INFO );
         RETURN
      }
      UP = .FALSE.
      IF ( LSAME( UPLO, 'U' ) ) UP = .TRUE.

      // Compute norm of op(A)*op2(C).

      ANORM = 0.0
      if ( UP ) {
         for (I = 1; I <= N; I++) {
            TMP = 0.0E+0
            for (J = 1; J <= I; J++) {
               TMP = TMP + CABS1( A( J, I ) * X( J ) )
            END DO
            DO J = I+1, N
               TMP = TMP + CABS1( A( I, J ) * X( J ) )
            END DO
            RWORK( I ) = TMP
            ANORM = MAX( ANORM, TMP )
         END DO
      } else {
         for (I = 1; I <= N; I++) {
            TMP = 0.0E+0
            for (J = 1; J <= I; J++) {
               TMP = TMP + CABS1( A( I, J ) * X( J ) )
            END DO
            DO J = I+1, N
               TMP = TMP + CABS1( A( J, I ) * X( J ) )
            END DO
            RWORK( I ) = TMP
            ANORM = MAX( ANORM, TMP )
         END DO
      }

      // Quick return if possible.

      if ( N.EQ.0 ) {
         CLA_SYRCOND_X = 1.0E+0
         RETURN
      } else if ( ANORM .EQ. 0.0E+0 ) {
         RETURN
      }

      // Estimate the norm of inv(op(A)).

      AINVNM = 0.0E+0

      KASE = 0
   10 CONTINUE
      clacn2(N, WORK( N+1 ), WORK, AINVNM, KASE, ISAVE );
      if ( KASE.NE.0 ) {
         if ( KASE.EQ.2 ) {

            // Multiply by R.

            for (I = 1; I <= N; I++) {
               WORK( I ) = WORK( I ) * RWORK( I )
            END DO

            if ( UP ) {
               csytrs('U', N, 1, AF, LDAF, IPIV, WORK, N, INFO );
            } else {
               csytrs('L', N, 1, AF, LDAF, IPIV, WORK, N, INFO );
            ENDIF

            // Multiply by inv(X).

            for (I = 1; I <= N; I++) {
               WORK( I ) = WORK( I ) / X( I )
            END DO
         } else {

            // Multiply by inv(X**T).

            for (I = 1; I <= N; I++) {
               WORK( I ) = WORK( I ) / X( I )
            END DO

            if ( UP ) {
               csytrs('U', N, 1, AF, LDAF, IPIV, WORK, N, INFO );
            } else {
               csytrs('L', N, 1, AF, LDAF, IPIV, WORK, N, INFO );
            }

            // Multiply by R.

            for (I = 1; I <= N; I++) {
               WORK( I ) = WORK( I ) * RWORK( I )
            END DO
         }
         GO TO 10
      }

      // Compute the estimate of the reciprocal condition number.

      IF( AINVNM .NE. 0.0E+0 ) CLA_SYRCOND_X = 1.0E+0 / AINVNM

      RETURN

      // End of CLA_SYRCOND_X

      }
