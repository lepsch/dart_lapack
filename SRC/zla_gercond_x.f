      double           FUNCTION ZLA_GERCOND_X( TRANS, N, A, LDA, AF, LDAF, IPIV, X, INFO, WORK, RWORK );

*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             TRANS;
      int                N, LDA, LDAF, INFO;
      // ..
      // .. Array Arguments ..
      int                IPIV( * );
      COMPLEX*16         A( LDA, * ), AF( LDAF, * ), WORK( * ), X( * )
      double             RWORK( * );
      // ..

*  =====================================================================

      // .. Local Scalars ..
      bool               NOTRANS;
      int                KASE;
      double             AINVNM, ANORM, TMP;
      int                I, J;
      COMPLEX*16         ZDUM
      // ..
      // .. Local Arrays ..
      int                ISAVE( 3 );
      // ..
      // .. External Functions ..
      bool               LSAME;
      // EXTERNAL LSAME
      // ..
      // .. External Subroutines ..
      // EXTERNAL ZLACN2, ZGETRS, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, MAX, REAL, DIMAG
      // ..
      // .. Statement Functions ..
      double             CABS1;
      // ..
      // .. Statement Function Definitions ..
      CABS1( ZDUM ) = ABS( DBLE( ZDUM ) ) + ABS( DIMAG( ZDUM ) )
      // ..
      // .. Executable Statements ..

      ZLA_GERCOND_X = 0.0D+0

      INFO = 0
      NOTRANS = LSAME( TRANS, 'N' )
      if ( .NOT. NOTRANS .AND. .NOT. LSAME( TRANS, 'T' ) .AND. .NOT. LSAME( TRANS, 'C' ) ) {
         INFO = -1
      } else if ( N.LT.0 ) {
         INFO = -2
      } else if ( LDA.LT.MAX( 1, N ) ) {
         INFO = -4
      } else if ( LDAF.LT.MAX( 1, N ) ) {
         INFO = -6
      }
      if ( INFO.NE.0 ) {
         xerbla('ZLA_GERCOND_X', -INFO );
         RETURN
      }

      // Compute norm of op(A)*op2(C).

      ANORM = 0.0D+0
      if ( NOTRANS ) {
         for (I = 1; I <= N; I++) {
            TMP = 0.0D+0
            for (J = 1; J <= N; J++) {
               TMP = TMP + CABS1( A( I, J ) * X( J ) )
            END DO
            RWORK( I ) = TMP
            ANORM = MAX( ANORM, TMP )
         END DO
      } else {
         for (I = 1; I <= N; I++) {
            TMP = 0.0D+0
            for (J = 1; J <= N; J++) {
               TMP = TMP + CABS1( A( J, I ) * X( J ) )
            END DO
            RWORK( I ) = TMP
            ANORM = MAX( ANORM, TMP )
         END DO
      }

      // Quick return if possible.

      if ( N.EQ.0 ) {
         ZLA_GERCOND_X = 1.0D+0
         RETURN
      } else if ( ANORM .EQ. 0.0D+0 ) {
         RETURN
      }

      // Estimate the norm of inv(op(A)).

      AINVNM = 0.0D+0

      KASE = 0
      } // 10
      zlacn2(N, WORK( N+1 ), WORK, AINVNM, KASE, ISAVE );
      if ( KASE.NE.0 ) {
         if ( KASE.EQ.2 ) {
            // Multiply by R.
            for (I = 1; I <= N; I++) {
               WORK( I ) = WORK( I ) * RWORK( I )
            END DO

            if ( NOTRANS ) {
               zgetrs('No transpose', N, 1, AF, LDAF, IPIV, WORK, N, INFO );
            } else {
               zgetrs('Conjugate transpose', N, 1, AF, LDAF, IPIV, WORK, N, INFO );
            }

            // Multiply by inv(X).

            for (I = 1; I <= N; I++) {
               WORK( I ) = WORK( I ) / X( I )
            END DO
         } else {

            // Multiply by inv(X**H).

            for (I = 1; I <= N; I++) {
               WORK( I ) = WORK( I ) / X( I )
            END DO

            if ( NOTRANS ) {
               zgetrs('Conjugate transpose', N, 1, AF, LDAF, IPIV, WORK, N, INFO );
            } else {
               zgetrs('No transpose', N, 1, AF, LDAF, IPIV, WORK, N, INFO );
            }

            // Multiply by R.

            for (I = 1; I <= N; I++) {
               WORK( I ) = WORK( I ) * RWORK( I )
            END DO
         }
         GO TO 10
      }

      // Compute the estimate of the reciprocal condition number.

      IF( AINVNM .NE. 0.0D+0 ) ZLA_GERCOND_X = 1.0D+0 / AINVNM

      RETURN

      // End of ZLA_GERCOND_X

      }
