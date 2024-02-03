      REAL FUNCTION CLA_GERCOND_X( TRANS, N, A, LDA, AF, LDAF, IPIV, X, INFO, WORK, RWORK )

*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             TRANS;
      int                N, LDA, LDAF, INFO;
      // ..
      // .. Array Arguments ..
      int                IPIV( * );
      COMPLEX            A( LDA, * ), AF( LDAF, * ), WORK( * ), X( * )
      REAL               RWORK( * )
      // ..

*  =====================================================================

      // .. Local Scalars ..
      bool               NOTRANS;
      int                KASE;
      REAL               AINVNM, ANORM, TMP
      int                I, J;
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
      // EXTERNAL CLACN2, CGETRS, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, MAX, REAL, AIMAG
      // ..
      // .. Statement Functions ..
      REAL               CABS1
      // ..
      // .. Statement Function Definitions ..
      CABS1( ZDUM ) = ABS( REAL( ZDUM ) ) + ABS( AIMAG( ZDUM ) )
      // ..
      // .. Executable Statements ..

      CLA_GERCOND_X = 0.0E+0

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
         CALL XERBLA( 'CLA_GERCOND_X', -INFO )
         RETURN
      }

      // Compute norm of op(A)*op2(C).

      ANORM = 0.0
      if ( NOTRANS ) {
         DO I = 1, N
            TMP = 0.0E+0
            DO J = 1, N
               TMP = TMP + CABS1( A( I, J ) * X( J ) )
            END DO
            RWORK( I ) = TMP
            ANORM = MAX( ANORM, TMP )
         END DO
      } else {
         DO I = 1, N
            TMP = 0.0E+0
            DO J = 1, N
               TMP = TMP + CABS1( A( J, I ) * X( J ) )
            END DO
            RWORK( I ) = TMP
            ANORM = MAX( ANORM, TMP )
         END DO
      }

      // Quick return if possible.

      if ( N.EQ.0 ) {
         CLA_GERCOND_X = 1.0E+0
         RETURN
      } else if ( ANORM .EQ. 0.0E+0 ) {
         RETURN
      }

      // Estimate the norm of inv(op(A)).

      AINVNM = 0.0E+0

      KASE = 0
   10 CONTINUE
      CALL CLACN2( N, WORK( N+1 ), WORK, AINVNM, KASE, ISAVE )
      if ( KASE.NE.0 ) {
         if ( KASE.EQ.2 ) {
            // Multiply by R.
            DO I = 1, N
               WORK( I ) = WORK( I ) * RWORK( I )
            END DO

            if ( NOTRANS ) {
               CALL CGETRS( 'No transpose', N, 1, AF, LDAF, IPIV, WORK, N, INFO )
            } else {
               CALL CGETRS( 'Conjugate transpose', N, 1, AF, LDAF, IPIV, WORK, N, INFO )
            ENDIF

            // Multiply by inv(X).

            DO I = 1, N
               WORK( I ) = WORK( I ) / X( I )
            END DO
         } else {

            // Multiply by inv(X**H).

            DO I = 1, N
               WORK( I ) = WORK( I ) / X( I )
            END DO

            if ( NOTRANS ) {
               CALL CGETRS( 'Conjugate transpose', N, 1, AF, LDAF, IPIV, WORK, N, INFO )
            } else {
               CALL CGETRS( 'No transpose', N, 1, AF, LDAF, IPIV, WORK, N, INFO )
            }

            // Multiply by R.

            DO I = 1, N
               WORK( I ) = WORK( I ) * RWORK( I )
            END DO
         }
         GO TO 10
      }

      // Compute the estimate of the reciprocal condition number.

      IF( AINVNM .NE. 0.0E+0 ) CLA_GERCOND_X = 1.0E+0 / AINVNM

      RETURN

      // End of CLA_GERCOND_X

      }
